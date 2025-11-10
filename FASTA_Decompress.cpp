#include "FASTA_Decompress.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <future>
#include <mutex>

#include "./zstd/lib/zstd.h"      // presumes zstd library is installed

#include "utils/helper.h"

#if !defined(__arm__) && !defined(__aarch64__) && !defined(__ARM_ARCH) && !defined(NO_AVX)
#include <immintrin.h>
const static __m256i DNA_AVX2_SHUFFLE_MASK = _mm256_set_epi32(0x07070707, 0x06060606, 0x05050505, 0x04040404, 0x03030303, 0x02020202, 0x01010101, 0x00000000);
const static __m256i DNA_AVX2_LO_MASK = _mm256_set1_epi16(0b0000110000000011);
constexpr static int32_t DNA_AVX2_LUT_I32 = (int32_t)('A') | ((int32_t)('C') << 8) | ((int32_t)('T') << 16) | ((int32_t)('G') << 24);
const static __m256i DNA_AVX2_LUT = _mm256_set_epi32((int32_t)('G'), (int32_t)('T'), (int32_t)('C'), DNA_AVX2_LUT_I32, (int32_t)('G'), (int32_t)('T'), (int32_t)('C'), DNA_AVX2_LUT_I32);
#endif

const static bool DISABLE_PARALLEL_OUTPUT = true;

std::mutex mtx;

// Only for version with packing into one byte
uint8_t FASTA_Decompress::pack_dna_no_lut(int32_t chunk) {
    return ((((chunk >> 1) & 0x03030303) * 0x40100401) >> 24) & 0xFF;
}

void FASTA_Decompress::init_LUT() {
    const uint8_t symbols[] = "ACGT";
    lut = (uint32_t *)malloc(256*sizeof(uint32_t));
    uint8_t chunk[4];
    for (int64_t i = 0; i < 4; i++) {
        for (int64_t j = 0; j < 4; j++) {
            for (int64_t k = 0; k < 4; k++) {
                for (int64_t l = 0; l < 4; l++) {
                    chunk[0] = symbols[i];
                    chunk[1] = symbols[j];
                    chunk[2] = symbols[k];
                    chunk[3] = symbols[l];
                    uint8_t packed = pack_dna_no_lut(*((uint32_t *)chunk));
                    chunk[0] = symbols[l];
                    chunk[1] = symbols[k];
                    chunk[2] = symbols[j];
                    chunk[3] = symbols[i];
                    memcpy(lut + packed, chunk, 4);
                }
            }
        }
    }
}

uint8_t *FASTA_Decompress::open_stream(char *decompressed_dest, uint8_t *buff, int64_t buff_size, int64_t dest_size) {
    if (buff[0] == NO_CODER)
        return buff + 1;
    size_t tmp = ZSTD_decompress(decompressed_dest, dest_size, buff + 1, buff_size - 1);
    auto err = ZSTD_isError(tmp);
    auto name = ZSTD_getErrorName(tmp);
    return (uint8_t*) decompressed_dest;
}

#ifndef NO_AVX
/*
 * Efficient technique taken from
 * https://github.com/Daniel-Liu-c0deb0t/cute-nucleotides/blob/master/src/n_to_bits.rs
 * Function bits_to_n_shuffle. Requires AVX2.
 * We translated the original Rust code to C++.
 */
inline void FASTA_Decompress::unpack_dna_avx2(char *dest, uint8_t *src) {
    __m256i v = _mm256_set1_epi64x(*((int64_t*) (src)));

    // duplicate each byte four times
    __m256i v1 = _mm256_shuffle_epi8(v, DNA_AVX2_SHUFFLE_MASK);

    // separately right shift each 16-bit chunk by 0 or 4 bits
    __m256i v2 = _mm256_srli_epi16(v1, 4);

    // merge together shifted chunks
    __m256i v_combined = _mm256_blend_epi16(v1, v2, 0b10101010);

    // only keep two bits in each byte (either 0b0011 or 0b1100)
    __m256i v_filtered = _mm256_and_si256(v_combined, DNA_AVX2_LO_MASK);

    // use lookup table to convert nucleotide bits to bytes
    __m256i v_final = _mm256_shuffle_epi8(DNA_AVX2_LUT, v_filtered);

    // store the result in the output vector
    _mm256_storeu_si256((__m256i*)(dest), v_final);
}
#endif

void FASTA_Decompress::decompress_block(
    int32_t b,
    size_t block_start
) {
    int t_id = get_thread_id(b);
    block_meta_t &block_meta = thread_blocks_meta[t_id];
    if (thread_buf[t_id] == nullptr) {
        thread_buf[t_id] = new char[8 + BUFFER_SIZE + 8]; // bytes for additional EOL and case routines
        thread_block[t_id] = new char[ffc_header.max_block_size + 8];
        int64_t number_of_case_chunks = (ffc_header.max_block_size + (CASE_CHUNK_SIZE - 1)) / CASE_CHUNK_SIZE;
        thread_case_mask[t_id] = new char[STREAM_HEADER_BYTES + number_of_case_chunks * 8];
        thread_dna[t_id] = new char[STREAM_HEADER_BYTES + ffc_header.max_block_size / 4];
        thread_raw[t_id] = new char[STREAM_HEADER_BYTES + ffc_header.max_block_size];
        thread_mix[t_id] = new char[STREAM_HEADER_BYTES + ffc_header.max_block_size];
        thread_subblocks_meta[t_id] = new uint32_t[STREAM_HEADER_BYTES + ffc_header.max_block_size / 2];
    }
    char* buf = thread_buf[t_id] + 8;
    char* block = thread_block[t_id];
    int chunk_size = ffc_header.chunk_size;

    if (!sequential_input) {
        if (thread_all_streams[t_id] == nullptr)
            thread_all_streams[t_id] = new char[STREAM_HEADER_BYTES * 5 + ffc_header.max_block_size +
                (ffc_header.max_block_size >> COMPRESSED_BLOCK_MARGIN_FRACTION_ORDER) +
                (ffc_header.case_compression_flag == 0 ? ffc_header.max_block_size / 8 : 0)];
        ifstream input_file(input_filename, std::ios::binary | std::ios::in);
        input_file.seekg(block_start, std::ios::beg);
        input_file.read(thread_all_streams[t_id], thread_blocks_meta[t_id].block_compressed_length);
    }
    uint8_t *source_stream = (uint8_t *) thread_all_streams[t_id];

    const int64_t CASE_CHUNK_SIZE = chunk_size * 8;
    int64_t number_of_chunks = (block_meta.block_length + (CASE_CHUNK_SIZE - 1)) / CASE_CHUNK_SIZE;
    uint8_t* case_mask_buff = source_stream;
    uint8_t *case_mask_stream = (uint8_t *)open_stream(
        thread_case_mask[t_id],
        case_mask_buff,
        block_meta.case_mask_compressed_size,
        number_of_chunks * 8
    );

    uint8_t* raw_buff = case_mask_buff + block_meta.case_mask_compressed_size;
    uint8_t* raw_stream = open_stream(
        thread_raw[t_id],
        raw_buff,
        block_meta.raw_stream_compressed_size,
        block_meta.raw_stream_size
    );

    uint8_t* dna_buff = raw_buff + block_meta.raw_stream_compressed_size;
    uint8_t* dna_stream = open_stream(
        thread_dna[t_id],
        dna_buff,
        block_meta.dna_stream_compressed_size,
        block_meta.dna_stream_size
    );

    uint8_t* mix_buff = dna_buff + block_meta.dna_stream_compressed_size;
    uint8_t* mix_stream = open_stream(
        thread_mix[t_id],
        mix_buff,
        block_meta.mix_stream_compressed_size,
        block_meta.mix_stream_size
    );

    uint8_t* subblocks_meta_buff = mix_buff + block_meta.mix_stream_compressed_size;
    uint32_t *subblocks_meta = (uint32_t*) open_stream(
        (char*) thread_subblocks_meta[t_id],
        subblocks_meta_buff,
        block_meta.subblocks_meta_compressed_size,
        block_meta.subblocks_count * sizeof(uint32_t)
    );

    int32_t buf_offset = 0;
    int32_t block_offset = 0;
    int32_t case_offset = 0;
    int32_t case_stream_offset = 0;
    int32_t raw_stream_offset = 0;
    int32_t dna_stream_offset = 0;
    int32_t mix_stream_offset = 0;

    size_t seq_line_length = block_meta.seq_line_length ? block_meta.seq_line_length : block_meta.block_length;
    int64_t predicted_first_EOL_offset = block_meta.first_EOL_offset;
    if (block_meta.first_EOL_offset < 0)
        predicted_first_EOL_offset = INT64_MAX;
    for (int64_t i = 0; i < block_meta.subblocks_count; i++) {
        uint32_t subblock_type = subblocks_meta[i] & SUBBLOCK_TYPE_MASK;
        uint32_t subblock_size = subblocks_meta[i] & SUBBLOCK_SIZE_MASK;

        bool split_subblock = buf_offset + subblock_size > BUFFER_SIZE;
        if (split_subblock) {
            uint32_t split_pos = (BUFFER_SIZE - buf_offset) / 8 * 8;
            subblocks_meta[i--] = subblock_type | (subblock_size - split_pos);
            subblock_size = split_pos;
        }

        uint32_t* dest_stream_32b;
        int64_t j = 0;
        switch (subblock_type) {
            case DNA_SUBBLOCK_TYPE:
#ifndef NO_AVX
                for (j = 0; j < subblock_size / 32; j++) {
                     unpack_dna_avx2(buf + buf_offset + (j << 5), dna_stream + dna_stream_offset);
                     dna_stream_offset += 8;
                }
#endif
                dest_stream_32b = (uint32_t*) (buf + buf_offset + (j << 5));
                for (j = j * 4; j < subblock_size / 8; j++) {
                    *(dest_stream_32b++) = lut[dna_stream[dna_stream_offset++]]; // unpacking uint8 -> uint32
                    *(dest_stream_32b++) = lut[dna_stream[dna_stream_offset++]];
                }
            break;
            case MIX_SUBBLOCK_TYPE:
                memcpy(buf + buf_offset, mix_stream + mix_stream_offset, subblock_size);
                mix_stream_offset += subblock_size;
            break;
            case NNN_SUBBLOCK_TYPE:
                memset(buf + buf_offset, 'N', subblock_size);
            break;
        }
        if (subblock_type != RAW_SUBBLOCK_TYPE)
            buf_offset += subblock_size;

        if (subblock_type == RAW_SUBBLOCK_TYPE || buf_offset >= BUFFER_SIZE - 8 || i == block_meta.subblocks_count - 1) {

            int32_t buf_pos = 0;
            while (buf_offset - buf_pos > predicted_first_EOL_offset) {
                memcpy(block + block_offset, buf + buf_pos, predicted_first_EOL_offset);
                buf_pos += predicted_first_EOL_offset;
                block_offset += predicted_first_EOL_offset;
                block[block_offset++] = '\n';
                predicted_first_EOL_offset = seq_line_length;
            }
            int len = buf_offset - buf_pos;
            predicted_first_EOL_offset -= len;
            memcpy(block + block_offset, buf + buf_pos, len);
            block_offset += len;

            if (subblock_type == RAW_SUBBLOCK_TYPE) {
                char* rawPtr = (char*) raw_stream + raw_stream_offset;
                memcpy(block + block_offset, rawPtr, subblock_size);
                raw_stream_offset += subblock_size;
                block_offset += subblock_size;
                if (!split_subblock)
                    block[block_offset++] = '\n';
                predicted_first_EOL_offset = seq_line_length;
            }
            buf_offset = 0;

            int32_t length = block_offset - case_offset;
            length = (length / 64) * 64;
            PgHelpers::convert_to_proper_case(case_mask_stream + case_stream_offset, block + case_offset, length);
            case_offset += length;
            case_stream_offset += length / 8;
        }
    }

    int64_t ceil_quantized_length = ((block_meta.block_length - case_offset + CASE_CHUNK_SIZE - 1) / CASE_CHUNK_SIZE) *
        CASE_CHUNK_SIZE;
    PgHelpers::convert_to_proper_case(case_mask_stream + case_stream_offset, block + case_offset, ceil_quantized_length);

    if (block_offset < block_meta.block_length)
        block[block_offset] = '\n';
    if (!sequential_output) {
        if (DISABLE_PARALLEL_OUTPUT) {
            std::unique_lock<std::mutex> lock(mtx);
            outStream->seekp(block_meta.block_start);
            outStream->write(block, block_meta.block_length);
            lock.unlock();
        } else {
            fstream os(output_filename, ios::out | std::ios::binary | std::ios::in);
            os.seekp(block_meta.block_start);
            os.write(block, block_meta.block_length);
            os.close();
        }

    }
}

void FASTA_Decompress::decompress_parallel(const string input_filename, const string output_filename) {
    this->input_filename = input_filename;
    this->output_filename = output_filename;

    istream* inStream;
    if (input_filename == STANDARD_IO_POSIX_ALIAS) {
#ifdef __MINGW32__
        if (_setmode(_fileno(stdin), _O_BINARY) == -1) {
            fprintf(stderr, "ERROR: switching cin to binary mode (errCode: %d)\n", strerror(errno));
            exit(EXIT_FAILURE);
        }
#endif
        inStream = &cin;
        sequential_input = true;
    } else {
        inStream = new ifstream(input_filename, ios_base::in | ios::binary);
        if (!*inStream) {
            cerr << "Error: unable to open input file: " << input_filename << endl;
            exit(1);
        }
    }

    int64_t ffc_size = inStream->tellg();
    ffc_stats_t ffc_stats;
    if (!sequential_input) {
        inStream->seekg( 0, std::ios::end );
        ffc_size = inStream->tellg() - ffc_size;
        inStream->seekg( ffc_size - sizeof(ffc_stats_t), std::ios::beg );
        inStream->read((char*) &ffc_stats, sizeof(ffc_stats_t));
        inStream->seekg( 0, std::ios::beg );

        if (input_filename != STANDARD_IO_POSIX_ALIAS){
            uint64_t space_in_bytes_required = ffc_stats.original_file_size;
            PgHelpers::createFolders(output_filename);
            if (!PgHelpers::has_enough_space(output_filename, space_in_bytes_required)) {
                *PgHelpers::appout << "Not enough space for output file." << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    if (output_filename == STANDARD_IO_POSIX_ALIAS) {
#ifdef __MINGW32__
        if (_setmode(_fileno(stdout), _O_BINARY) == -1)
            fprintf(stderr, "WARNING: switching cout to binary mode failed (errCode: %d)\n", strerror(errno));
#endif
        outStream = &cout;
        sequential_output = true;
    } else {
        bool output_exists = output_filename == "/dev/null" ? false : std::filesystem::exists(output_filename);
        if (output_exists) {
            size_t output_size = std::filesystem::file_size(output_filename);
            if (output_size > ffc_stats.original_file_size)
                std::filesystem::resize_file(output_filename, ffc_stats.original_file_size);
        } else
            PgHelpers::createFolders(output_filename);
        outStream = new ofstream(output_filename, ios::out | ios::binary |
            (output_exists ? std::ios::in : std::ios::trunc));
        if (!*outStream) {
            cerr << "Error: unable to open output file: " << output_filename << endl;
            cerr.flush();
            exit(1);
        }
        if (!output_exists && !sequential_output) {
            outStream->seekp(ffc_stats.original_file_size - 1);
            outStream->write("", 1);
            outStream->seekp(0);
        }
    }

    read_ffc_header(inStream, &ffc_header);
    if (ffc_header.ffc_header != 0x6366662e) {
        cerr << "ERROR: incorrect header: " << ffc_header.ffc_header << endl;
        cerr.flush();
        exit(EXIT_FAILURE);
    }

    string orig_filename;
    orig_filename.resize(ffc_header.orig_filename_length);
    inStream->read(orig_filename.data(), ffc_header.orig_filename_length);

    if (!sequential_input && PgHelpers::numberOfThreads > ffc_stats.blocks_count)
        PgHelpers::numberOfThreads = ffc_stats.blocks_count ? ffc_stats.blocks_count : 1;

    thread_blocks_meta.resize(PgHelpers::numberOfThreads);
    thread_buf.resize(PgHelpers::numberOfThreads, nullptr);
    thread_block.resize(PgHelpers::numberOfThreads, nullptr);
    thread_all_streams.resize(PgHelpers::numberOfThreads, nullptr);
    thread_case_mask.resize(PgHelpers::numberOfThreads, nullptr);
    thread_raw.resize(PgHelpers::numberOfThreads, nullptr);
    thread_dna.resize(PgHelpers::numberOfThreads, nullptr);
    thread_mix.resize(PgHelpers::numberOfThreads, nullptr);
    thread_subblocks_meta.resize(PgHelpers::numberOfThreads, nullptr);

    std::vector<std::future<void>> block_threads;
    int64_t blocks_count = 0;
    int64_t b = 0;
    size_t next_block_meta_pos = FFC_HEADER_SIZE + ffc_header.orig_filename_length;
    do {
        if (b <= blocks_count) {
            int t_id = b % PgHelpers::numberOfThreads;
            if (!sequential_input)
                inStream->seekg(next_block_meta_pos, std::ios::beg);
            read_block_meta_data(inStream, &thread_blocks_meta[t_id]);
            if (thread_blocks_meta[t_id].block_length > 0) {
                blocks_count++;
                if (sequential_input) {
                    if (thread_all_streams[t_id] == nullptr)
                        thread_all_streams[t_id] = new char[STREAM_HEADER_BYTES * 5 + ffc_header.max_block_size +
                         (ffc_header.max_block_size >> COMPRESSED_BLOCK_MARGIN_FRACTION_ORDER)];
                    inStream->read(thread_all_streams[t_id], thread_blocks_meta[t_id].block_compressed_length);
                }
                block_threads.emplace_back(
                        std::async(
                        launch::async,
                        &FASTA_Decompress::decompress_block,
                        this,
                        b,
                        next_block_meta_pos + FFC_BLOCK_META_SIZE
                    )
                );
            }
            next_block_meta_pos += FFC_BLOCK_META_SIZE + thread_blocks_meta[t_id].block_compressed_length;
        }

        int64_t i = b - (PgHelpers::numberOfThreads - 1);
        if (i >= 0 && i < blocks_count) {
            if (block_threads[i].valid()) {
                block_threads[i].wait();
                if (sequential_output) {
                    int ti_id = get_thread_id(i);
                    outStream->write(thread_block[ti_id], thread_blocks_meta[ti_id].block_length);
                }
            } else {
                throw "ERROR: Thread is not valid";
            }
        }
    } while (++b < blocks_count + PgHelpers::numberOfThreads);

    if (input_filename != STANDARD_IO_POSIX_ALIAS)
        delete inStream;
    if (output_filename != STANDARD_IO_POSIX_ALIAS)
        delete outStream;
    for (int t = 0; t < PgHelpers::numberOfThreads; t++) {
        delete [] thread_buf[t];
        delete [] thread_block[t];
        delete [] thread_all_streams[t];
        delete [] thread_case_mask[t];
        delete [] thread_raw[t];
        delete [] thread_dna[t];
        delete [] thread_mix[t];
        delete [] thread_subblocks_meta[t];
    }
}

void FASTA_Decompress::read_ffc_header(istream* inStream, const ffc_header_t* ffc_header)
{
    inStream->read((char*)&ffc_header->ffc_header, sizeof(ffc_header->ffc_header));
    inStream->read((char*)&ffc_header->version, sizeof(ffc_header->version));
    inStream->read((char*)&ffc_header->chunk_size, sizeof(ffc_header->chunk_size));
    inStream->read((char*)&ffc_header->max_block_size, sizeof(ffc_header->max_block_size));
    inStream->read((char*)&ffc_header->case_compression_flag, sizeof(ffc_header->case_compression_flag));
    inStream->read((char*)&ffc_header->raw_compression_flag, sizeof(ffc_header->raw_compression_flag));
    inStream->read((char*)&ffc_header->dna_compression_flag, sizeof(ffc_header->dna_compression_flag));
    inStream->read((char*)&ffc_header->mix_compression_flag, sizeof(ffc_header->mix_compression_flag));
    inStream->read((char*)&ffc_header->subblocks_meta_compression_flag, sizeof(ffc_header->subblocks_meta_compression_flag));
    inStream->read((char*)&ffc_header->crc, sizeof(ffc_header->crc));
    inStream->read((char*)&ffc_header->orig_file_timestamp, sizeof(ffc_header->orig_file_timestamp));
    inStream->read((char*)&ffc_header->orig_filename_length, sizeof(ffc_header->orig_filename_length));
}

void FASTA_Decompress::read_block_meta_data(istream* inStream, const block_meta_t* block_meta)
{
    inStream->read((char*) &block_meta->block_start, sizeof(block_meta->block_start));
    inStream->read((char*) &block_meta->block_length, sizeof(block_meta->block_length));
    inStream->read((char*) &block_meta->block_compressed_length, sizeof(block_meta->block_compressed_length));
    inStream->read((char*) &block_meta->case_mask_compressed_size, sizeof(block_meta->case_mask_compressed_size));
    inStream->read((char*) &block_meta->raw_stream_size, sizeof(block_meta->raw_stream_size));
    inStream->read((char*) &block_meta->raw_stream_compressed_size, sizeof(block_meta->raw_stream_compressed_size));
    inStream->read((char*) &block_meta->dna_stream_size, sizeof(block_meta->dna_stream_size));
    inStream->read((char*) &block_meta->dna_stream_compressed_size, sizeof(block_meta->dna_stream_compressed_size));
    inStream->read((char*) &block_meta->mix_stream_size, sizeof(block_meta->mix_stream_size));
    inStream->read((char*) &block_meta->mix_stream_compressed_size, sizeof(block_meta->mix_stream_compressed_size));
    inStream->read((char*) &block_meta->subblocks_count, sizeof(block_meta->subblocks_count));
    inStream->read((char*) &block_meta->subblocks_meta_compressed_size, sizeof(block_meta->subblocks_meta_compressed_size));
    inStream->read((char*) &block_meta->first_EOL_offset, sizeof(block_meta->first_EOL_offset));
    inStream->read((char*) &block_meta->seq_line_length, sizeof(block_meta->seq_line_length));
    inStream->read((char*) &block_meta->headers_count, sizeof(block_meta->headers_count));
}

int decompress(
    string input_filename, 
    string output_filename
) {
    if (output_filename == STANDARD_IO_POSIX_ALIAS) {
        PgHelpers::appout = &null_stream;
    }

    if (input_filename != STANDARD_IO_POSIX_ALIAS && !std::filesystem::exists(input_filename)) {
        *PgHelpers::appout << "Input file does not exist: " << input_filename << endl << flush;
        return EXIT_FAILURE;
    }

    *PgHelpers::appout << "Input file: " << input_filename  << endl << flush;
    *PgHelpers::appout << "Output file: " << output_filename << endl << flush;
    *PgHelpers::appout << "Threads: " << PgHelpers::numberOfThreads << endl << flush;

    auto stopWatch = PgHelpers::time_checkpoint();

    FASTA_Decompress decompressor;
    decompressor.decompress_parallel(input_filename, output_filename);

    *PgHelpers::appout << "Decompression finished in: " << PgHelpers::time_millis(stopWatch) << " msec." << endl << flush;

    return EXIT_SUCCESS;
}
