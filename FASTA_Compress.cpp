#include "FASTA_Compress.h"
#include "version.hpp"

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <future>
#include <filesystem>

#include "./zstd/lib/zstd.h"      // presumes zstd library is installed

block_meta_t FASTA_Compress::init_block(int32_t b, char* prev_block) {
    block_meta_t block_meta;
    int prev_t_id = get_thread_id(b - 1);
    block_meta.block_start = b ? thread_blocks_meta[prev_t_id].block_start + thread_blocks_meta[prev_t_id].block_length : 0;

    int32_t prev_block_starts_flag = b ? thread_blocks_meta[prev_t_id].first_EOL_offset : BLOCK_STARTS_JUST_AFTER_EOL;
    char* prev_end = prev_block + (b ? thread_blocks_meta[prev_t_id].block_length : 0);
    char* ptr = prev_end;
    char* guard = prev_block;
    while (ptr-- != guard) {
        if (*ptr == '\n') {
            ++ptr;
            if (ptr == prev_end)
                block_meta.first_EOL_offset = BLOCK_STARTS_JUST_AFTER_EOL;
            else if (*ptr == header_symbol)
                block_meta.first_EOL_offset = BLOCK_STARTS_INSIDE_HEADER;
            else
                block_meta.first_EOL_offset = BLOCK_STARTS_INSIDE_SEQUENCE;
            break;
        }
    }
    if (++ptr == guard)
        block_meta.first_EOL_offset = b && (*ptr == header_symbol && prev_block_starts_flag == BLOCK_STARTS_JUST_AFTER_EOL) ?
            BLOCK_STARTS_INSIDE_HEADER : prev_block_starts_flag;
    return block_meta;
}

int32_t FASTA_Compress::remove_EOLs_and_find_tail(char *seq, int32_t& seq_len, int8_t chunk_size,
    int32_t& seq_line_length, int32_t& max_seq_without_eols_length, bool ignore_first_EOL_period, int32_t declared_first_eol_pos) {
    if (seq_len == 0)
        return 0;
    int32_t last_eol_pos = -1;
    int32_t count = 0;
    int32_t i = 0;
    int32_t last_move_pos = 0;
    int32_t tail_len = 0;

    while (i + 7 < seq_len) {
        uint64_t tmp;
        memcpy(&tmp, seq + i, 8);

        uint64_t top_bits = ~(tmp ^ EIGHT_EOLS) & 0x8080808080808080;
        uint64_t lower_bits = ((~(tmp ^ EIGHT_EOLS) & 0x7F7F7F7F7F7F7F7F) + 0x0101010101010101) & 0x8080808080808080;
        if (top_bits & lower_bits) { // scan those 8 chars (from x) one by one, as there's at least one '\n' in them
            memmove(seq + last_move_pos - count, seq + last_move_pos, i - last_move_pos);
            for (int64_t i_end = i + 8; i < i_end; i++) {
                if (seq[i] == '\n') {
                    if (!ignore_first_EOL_period || last_eol_pos != -1) {
                        int32_t current_line_length = i - last_eol_pos - 1;
                        bool eol_before_header = i == seq_len - 1;
                        if (seq_line_length == 0) {
                            if (current_line_length < max_seq_without_eols_length) {
                                tail_len = seq_len - i;
                                seq_len = i;
                                break;
                            }
                            seq_line_length = eol_before_header ? 0 : current_line_length;
                        } else if (seq_line_length != current_line_length && (!eol_before_header || seq_line_length < current_line_length)) {
                            while (current_line_length-- > seq_line_length) {
                                i--;
                                seq[i] = seq[i - count];
                            }
                            tail_len = seq_len - i;
                            seq_len = i;
                            break;
                        }
                    }
                    last_eol_pos = i < seq_len - 1 ? i : last_eol_pos;
                    count++;
                }
                else
                    seq[i - count] = seq[i];
            }
            last_move_pos = i;
        }
        else {
            i += 8;
        }
    }
    memmove(seq + last_move_pos - count, seq + last_move_pos, i - last_move_pos);

    for (; !tail_len && i < seq_len; i++) {
        if (seq[i] == '\n') {
            if (!ignore_first_EOL_period || last_eol_pos != -1) {
                int32_t current_line_length = i - last_eol_pos - 1;
                bool eol_before_header = i == seq_len - 1;
                if (seq_line_length == 0) {
                    if (current_line_length < max_seq_without_eols_length) {
                        tail_len = seq_len - i;
                        seq_len = i;
                        break;
                    }
                    seq_line_length = eol_before_header ? 0 : current_line_length;
                } else if (seq_line_length != current_line_length && (!eol_before_header || seq_line_length < current_line_length)) {
                    while (current_line_length-- > seq_line_length) {
                        i--;
                        seq[i] = seq[i - count];
                    }
                    tail_len = seq_len - i;
                    seq_len = i;
                    break;
                }
            }
            last_eol_pos = i < seq_len - 1 ? i : last_eol_pos;
            count++;
        } else
            seq[i - count] = seq[i];
    }
    if (!tail_len && seq[i - 1] != '\n') {
        int32_t current_line_length = i - last_eol_pos - 1;
        if (seq_line_length < current_line_length) {
            while (current_line_length-- > seq_line_length) {
                i--;
                seq[i] = seq[i - count];
            }
            tail_len = seq_len - i;
            seq_len = i;
        }
    }
    seq_len -= count;
    int32_t seq_end = (seq_len / chunk_size) * chunk_size;
    int32_t tail_end = i;
    if (i && seq[i - 1] == '\n')
        i--;

    while (seq_len > seq_end) {
        while (declared_first_eol_pos < i-- && i == last_eol_pos) {
            seq[i] = '\n';
            last_eol_pos -= seq_line_length == 0 ? last_eol_pos : seq_line_length + 1;
        }
        seq[i] = seq[--seq_len];
    }
    if (i && i - 1 == last_eol_pos) {
        seq[--i] = '\n';
    }
    tail_len += tail_end - i;

    if (seq_line_length == 0 && max_seq_without_eols_length <= seq_len)
        max_seq_without_eols_length = seq_len - (ignore_first_EOL_period ? declared_first_eol_pos : 0);

    return tail_len;
}

template<bool binary>
bool FASTA_Compress::all_ACGT(uint64_t chunk) {
    uint64_t result = (
        ((((~(chunk ^ 0x4141414141414141)) & 0x7F7F7F7F7F7F7F7F) + 0x0101010101010101) & 0x8080808080808080) |
        ((((~(chunk ^ 0x4343434343434343)) & 0x7F7F7F7F7F7F7F7F) + 0x0101010101010101) & 0x8080808080808080) |
        ((((~(chunk ^ 0x4747474747474747)) & 0x7F7F7F7F7F7F7F7F) + 0x0101010101010101) & 0x8080808080808080) |
        ((((~(chunk ^ 0x5454545454545454)) & 0x7F7F7F7F7F7F7F7F) + 0x0101010101010101) & 0x8080808080808080)
    );
    return binary ? ((result & ~(chunk & 0x8080808080808080)) == 0x8080808080808080) : (result == 0x8080808080808080);
}

inline uint8_t FASTA_Compress::pack_dna_naive(char *chunk) {
    uint8_t encoded = ((chunk[0] & 0x6) >> 1) | ((chunk[1] & 0x6) << 1) | ((chunk[2] & 0x6) << 3) | ((chunk[3] & 0x6) << 5);
    return encoded;
}

template<typename E>
int64_t FASTA_Compress::compress_stream(int t_id, char* dest, vector<E>& stream, int32_t compression_flag, int32_t window_size_log) {
    size_t stream_header_size_in_E = (STREAM_HEADER_BYTES + sizeof(E) - 1) / sizeof(E);
    size_t compressed_size = STREAM_HEADER_BYTES + (stream.size() - stream_header_size_in_E) * sizeof(E);
    size_t uncompressed_bytes = (stream.size() - 1) * sizeof(E);
    if (compression_flag) {
        size_t compressed_size_bound = ZSTD_compressBound(uncompressed_bytes);
        ZSTD_CCtx_setParameter(thread_zsdt_cctx[t_id], ZSTD_c_enableLongDistanceMatching, window_size_log != 0);
        if (window_size_log != 0) {
            ZSTD_CCtx_setParameter(thread_zsdt_cctx[t_id], ZSTD_c_windowLog, window_size_log);
        }
        ZSTD_CCtx_setParameter(thread_zsdt_cctx[t_id], ZSTD_c_compressionLevel, compression_flag);
        compressed_size = ZSTD_compress2(thread_zsdt_cctx[t_id], dest + STREAM_HEADER_BYTES, compressed_size_bound,
            stream.data() + stream_header_size_in_E, uncompressed_bytes, compression_flag);
        compressed_size++;
    } else {
        memcpy(dest + STREAM_HEADER_BYTES, stream.data() + stream_header_size_in_E, uncompressed_bytes);
    }
    dest[0] = compression_flag ? ZSTD_CODER : NO_CODER;
    return compressed_size;
}

int64_t FASTA_Compress::compress_stream(int t_id, char* dest, char* src, size_t srcSize, int32_t compression_flag, int32_t window_size_log) {
    size_t compressed_size = srcSize;
    size_t uncompressed_bytes = srcSize - 1;
    if (compression_flag) {
        size_t compressed_size_bound = ZSTD_compressBound(uncompressed_bytes);
        ZSTD_CCtx_setParameter(thread_zsdt_cctx[t_id], ZSTD_c_enableLongDistanceMatching, window_size_log != 0);
        if (window_size_log != 0) {
            ZSTD_CCtx_setParameter(thread_zsdt_cctx[t_id], ZSTD_c_windowLog, window_size_log);
        }
        ZSTD_CCtx_setParameter(thread_zsdt_cctx[t_id], ZSTD_c_compressionLevel, compression_flag);
        compressed_size = ZSTD_compress2(thread_zsdt_cctx[t_id], dest + STREAM_HEADER_BYTES, compressed_size_bound,
            src + STREAM_HEADER_BYTES, uncompressed_bytes);
        compressed_size++;
    } else {
        memcpy(dest + STREAM_HEADER_BYTES, src + STREAM_HEADER_BYTES, uncompressed_bytes);
    }
    dest[0] = compression_flag ? ZSTD_CODER : NO_CODER;
    return compressed_size;
}

template<bool irregular_EOLs_mode>
void FASTA_Compress::compress_block(int32_t b) {
    int t_id = get_thread_id(b);
    block_meta_t& block_meta = thread_blocks_meta[t_id];
    if (!irregular_EOLs_mode && header_symbol != '>') {
        block_meta.seq_line_length = 0;
        compress_block<true>(b);
        return;
    }
    block_meta.headers_count = 0;
    size_t max_block_compressed_length = STREAM_HEADER_BYTES * 5 + ffc_header.max_block_size +
         (ffc_header.max_block_size >> COMPRESSED_BLOCK_MARGIN_FRACTION_ORDER) +
         (ffc_header.case_compression_flag == 0 ? ffc_header.max_block_size / 8 : 0);
    if (thread_dna[t_id] == nullptr) {
        thread_all_streams[t_id] = new char[max_block_compressed_length];
        int64_t number_of_case_chunks = (ffc_header.max_block_size + (CASE_CHUNK_SIZE - 1)) / CASE_CHUNK_SIZE;
        thread_case_mask[t_id] = new char[STREAM_HEADER_BYTES + number_of_case_chunks * 8];
        thread_dna[t_id] = new char[STREAM_HEADER_BYTES + ffc_header.max_block_size / 4];
        thread_raw[t_id] = new char[STREAM_HEADER_BYTES + ffc_header.max_block_size];
        thread_mix[t_id] = new char[STREAM_HEADER_BYTES + ffc_header.max_block_size];
        thread_subblocks_meta[t_id] = new uint32_t[STREAM_HEADER_BYTES + ffc_header.max_block_size / 2];
        thread_zsdt_cctx[t_id] = ZSTD_createCCtx();
    }
    char* input_data = thread_block[t_id];
    char* case_mask = thread_case_mask[t_id];
    char* raw = thread_raw[t_id];
    char* dna = thread_dna[t_id];
    char* mix = thread_mix[t_id];
    uint32_t* subblocks_meta = thread_subblocks_meta[t_id];
    size_t raw_pos = 0;
    size_t dna_pos = 0;
    size_t mix_pos = 0;
    size_t subblocks_pos = 0;
    for (int i = 0 ; i < STREAM_HEADER_BYTES; i++) {
        raw[raw_pos++] = 0;
        dna[dna_pos++] = 0;
        mix[mix_pos++] = 0;
        subblocks_meta[subblocks_pos++] = 0;
    }

    char* ptr = input_data;
    int64_t block_length = block_meta.block_length;

    int32_t first_EOL_pos = find(ptr, ptr + block_length, '\n') - ptr;
    int32_t max_seq_without_eols_length = MIN_SEQ_LINE_LENGTH;
    char* guard = ptr + block_length;

    uint64_t *casePtr = (uint64_t*) (case_mask + 1);
    char* caseChunkPtr = ptr;

    int32_t block_start_status = block_meta.first_EOL_offset;
    if (block_start_status == BLOCK_STARTS_INSIDE_HEADER) {
        uint32_t header_suffix_length = first_EOL_pos;
        memcpy(raw + raw_pos, ptr, header_suffix_length);
        raw_pos += header_suffix_length;
        ptr += first_EOL_pos + 1; // skipping the EOL (if exists)
        subblocks_meta[subblocks_pos++] = RAW_SUBBLOCK_TYPE | header_suffix_length;
    }

    block_meta.first_EOL_offset = irregular_EOLs_mode ? INCONSISTENT_EOLS : first_EOL_pos;

    while (ptr < guard) {
        while (caseChunkPtr < ptr) {
            memset(casePtr++, 0, 8);
            caseChunkPtr += CASE_CHUNK_SIZE;
        }

        char* header_pos = find(ptr, guard, header_symbol);
        while (header_pos < guard &&
            (header_pos == input_data ? block_start_status != BLOCK_STARTS_JUST_AFTER_EOL : header_pos[-1] != '\n' && header_pos[-1] != '\0'))
            header_pos = find(header_pos + 1, guard, header_symbol);

        // Stage 1 - Uppercase sequence except for a short tail

        int64_t to_uppercase_length = header_pos - caseChunkPtr;
        int64_t number_of_case_chunks = to_uppercase_length / CASE_CHUNK_SIZE;
        for (int64_t i = 0; i < number_of_case_chunks; i++) {
            uint64_t flags = PgHelpers::pack_case_flags_naive_8(caseChunkPtr);
            memcpy(casePtr++, &flags, 8);
            for (size_t j = 0; j < 8; ++j) {
                uint64_t one_chunk;
                memcpy(&one_chunk, caseChunkPtr + j * 8, 8);
                one_chunk &= 0xDFDFDFDFDFDFDFDF;
                memcpy(caseChunkPtr + j * 8, &one_chunk, 8);
            }
            caseChunkPtr += CASE_CHUNK_SIZE;
        }

        // Stage 2 - EOL removal and storing line length + 1

        int32_t sequence_without_tail_length = header_pos - ptr;
        int32_t tail_with_eols_length = 0;
        if (!irregular_EOLs_mode) {
            bool ignore_first_EOL_period = block_meta.first_EOL_offset >= ptr - input_data;
            int32_t declared_first_EOL_offset = block_meta.first_EOL_offset >= ptr - input_data ?
                block_meta.first_EOL_offset - (ptr - input_data) : -1;
            tail_with_eols_length = remove_EOLs_and_find_tail(ptr, sequence_without_tail_length, CHUNK_SIZE,
                block_meta.seq_line_length, max_seq_without_eols_length, ignore_first_EOL_period, declared_first_EOL_offset);
        } else {
            tail_with_eols_length = sequence_without_tail_length % CHUNK_SIZE;
        }

        int64_t number_of_chunks = sequence_without_tail_length / CHUNK_SIZE;
        uint64_t* chunkPtr = (uint64_t*) ptr;
        uint64_t* guardPtr = chunkPtr + number_of_chunks;
        uint64_t* subPtr = chunkPtr;
        uint64_t* encPtr = subPtr;

        // Stage 3 - Packing sequence

        while (subPtr < guardPtr) {
            --encPtr;
            while (++encPtr < guardPtr && all_ACGT<true>(*(encPtr))) {
                dna[dna_pos++] = pack_dna_naive((char*) encPtr);
                dna[dna_pos++] = pack_dna_naive(((char*) (encPtr)) + 4);
            }
            uint32_t dna_length = encPtr - subPtr;
            if (dna_length) {
                subblocks_meta[subblocks_pos++] = DNA_SUBBLOCK_TYPE | (dna_length * CHUNK_SIZE);
                subPtr = encPtr;
            }
            --encPtr;
            while (++encPtr < guardPtr && *encPtr == 0x4E4E4E4E4E4E4E4E)
                ;
            uint32_t nnn_length = encPtr - subPtr;
            if (nnn_length) {
                subblocks_meta[subblocks_pos++] = NNN_SUBBLOCK_TYPE | (nnn_length * 8);
                subPtr = encPtr;
            }
            --encPtr;
            while (++encPtr < guardPtr && !all_ACGT<true>(*(encPtr)) && *encPtr != 0x4E4E4E4E4E4E4E4E)
                ;
            uint32_t mix_length = encPtr - subPtr;
            memcpy(mix + mix_pos, subPtr, CHUNK_SIZE * mix_length);
            mix_pos += CHUNK_SIZE * mix_length;

            if (mix_length) {
                subblocks_meta[subblocks_pos++] = MIX_SUBBLOCK_TYPE | (mix_length * CHUNK_SIZE);
                subPtr = encPtr;
            }

            if (dna_length + mix_length + nnn_length == 0) {
                *PgHelpers::verboseout << "Not supported symbol in chunk at position " << (subPtr - chunkPtr) * CHUNK_SIZE << endl;
                exit(EXIT_FAILURE);
            }
        }

        uint32_t raw_length = tail_with_eols_length;
        memcpy(raw + raw_pos, header_pos - tail_with_eols_length, tail_with_eols_length);
        raw_pos += tail_with_eols_length;
        ptr = header_pos;
        if (ptr < guard) {
            block_meta.headers_count++;
            char* eol_pos = find(ptr + 1, guard, '\n');
            uint32_t header_length = eol_pos - ptr;
            raw_length += header_length;
            memcpy(raw + raw_pos, ptr, header_length);
            raw_pos += header_length;
            ptr = eol_pos + 1; // skipping the EOL (if exists)
        }
        if (raw_length)
            subblocks_meta[subblocks_pos++] = RAW_SUBBLOCK_TYPE | raw_length;
    }

    int64_t total_number_of_case_chunks = (block_length + (CASE_CHUNK_SIZE - 1)) / CASE_CHUNK_SIZE;
    while (caseChunkPtr < input_data + total_number_of_case_chunks * CASE_CHUNK_SIZE) {
        memset(casePtr++, 0, 8);
        caseChunkPtr += CASE_CHUNK_SIZE;
    }

    block_meta.raw_stream_size = raw_pos - STREAM_HEADER_BYTES;
    block_meta.dna_stream_size = dna_pos - STREAM_HEADER_BYTES;
    block_meta.mix_stream_size = mix_pos - STREAM_HEADER_BYTES;
    block_meta.subblocks_count = subblocks_pos - STREAM_HEADER_BYTES;

    char* streamsPtr = thread_all_streams[t_id];

    block_meta.case_mask_compressed_size = compress_stream(t_id, streamsPtr, case_mask, total_number_of_case_chunks * 8 + STREAM_HEADER_BYTES, ffc_header.case_compression_flag);
    streamsPtr += block_meta.case_mask_compressed_size;
    block_meta.raw_stream_compressed_size = compress_stream(t_id, streamsPtr, raw, raw_pos, ffc_header.raw_compression_flag);
    streamsPtr += block_meta.raw_stream_compressed_size;
    block_meta.dna_stream_compressed_size = compress_stream(t_id, streamsPtr, dna, dna_pos, ffc_header.dna_compression_flag, long_matching_window_size_log);
    streamsPtr += block_meta.dna_stream_compressed_size;
    block_meta.mix_stream_compressed_size = compress_stream(t_id, streamsPtr, mix, mix_pos, ffc_header.mix_compression_flag);
    streamsPtr += block_meta.mix_stream_compressed_size;
    block_meta.subblocks_meta_compressed_size = compress_stream(t_id, streamsPtr, ((char*) subblocks_meta) + (sizeof(uint32_t) - STREAM_HEADER_BYTES),
        (subblocks_pos - STREAM_HEADER_BYTES) * sizeof(uint32_t) + STREAM_HEADER_BYTES, ffc_header.subblocks_meta_compression_flag);;

    block_meta.block_compressed_length = (uint32_t) block_meta.case_mask_compressed_size +
        (uint32_t) block_meta.raw_stream_compressed_size + (uint32_t)  block_meta.dna_stream_compressed_size +
        (uint32_t) block_meta.mix_stream_compressed_size + (uint32_t) block_meta.subblocks_meta_compressed_size;

    if (block_meta.block_compressed_length > max_block_compressed_length) {
        cerr << "Internal error (try increasing compression level)." << endl;
        exit(EXIT_FAILURE);
    }
}

void FASTA_Compress::write_ffc_header(ostream* outStream, const ffc_header_t* ffc_header)
{
    outStream->write((char*)&ffc_header->ffc_header, sizeof(ffc_header->ffc_header));
    outStream->write((char*)&ffc_header->version, sizeof(ffc_header->version));
    outStream->write((char*)&ffc_header->chunk_size, sizeof(ffc_header->chunk_size));
    outStream->write((char*)&ffc_header->max_block_size, sizeof(ffc_header->max_block_size));
    outStream->write((char*)&ffc_header->case_compression_flag, sizeof(ffc_header->case_compression_flag));
    outStream->write((char*)&ffc_header->raw_compression_flag, sizeof(ffc_header->raw_compression_flag));
    outStream->write((char*)&ffc_header->dna_compression_flag, sizeof(ffc_header->dna_compression_flag));
    outStream->write((char*)&ffc_header->mix_compression_flag, sizeof(ffc_header->mix_compression_flag));
    outStream->write((char*)&ffc_header->subblocks_meta_compression_flag, sizeof(ffc_header->subblocks_meta_compression_flag));
    outStream->write((char*)&ffc_header->crc, sizeof(ffc_header->crc));
    outStream->write((char*)&ffc_header->orig_file_timestamp, sizeof(ffc_header->orig_file_timestamp));
    outStream->write((char*)&ffc_header->orig_filename_length, sizeof(ffc_header->orig_filename_length));
}

void FASTA_Compress::write_block_meta_data(ostream* outStream, const block_meta_t* block_meta)
{
    outStream->write((char*) &block_meta->block_start, sizeof(block_meta->block_start));
    outStream->write((char*) &block_meta->block_length, sizeof(block_meta->block_length));
    outStream->write((char*) &block_meta->block_compressed_length, sizeof(block_meta->block_compressed_length));
    outStream->write((char*) &block_meta->case_mask_compressed_size, sizeof(block_meta->case_mask_compressed_size));
    outStream->write((char*) &block_meta->raw_stream_size, sizeof(block_meta->raw_stream_size));
    outStream->write((char*) &block_meta->raw_stream_compressed_size, sizeof(block_meta->raw_stream_compressed_size));
    outStream->write((char*) &block_meta->dna_stream_size, sizeof(block_meta->dna_stream_size));
    outStream->write((char*) &block_meta->dna_stream_compressed_size, sizeof(block_meta->dna_stream_compressed_size));
    outStream->write((char*) &block_meta->mix_stream_size, sizeof(block_meta->mix_stream_size));
    outStream->write((char*) &block_meta->mix_stream_compressed_size, sizeof(block_meta->mix_stream_compressed_size));
    outStream->write((char*) &block_meta->subblocks_count, sizeof(block_meta->subblocks_count));
    outStream->write((char*) &block_meta->subblocks_meta_compressed_size, sizeof(block_meta->subblocks_meta_compressed_size));
    outStream->write((char*) &block_meta->first_EOL_offset, sizeof(block_meta->first_EOL_offset));
    outStream->write((char*) &block_meta->seq_line_length, sizeof(block_meta->seq_line_length));
    outStream->write((char*) &block_meta->headers_count, sizeof(block_meta->headers_count));
}

void FASTA_Compress::compress_parallel(const string input_filename, const string output_filename) {
    istream* inStream;
    if (input_filename == STANDARD_IO_POSIX_ALIAS) {
#ifdef __MINGW32__
        if (_setmode(_fileno(stdin), _O_BINARY) == -1) {
            fprintf(stderr, "ERROR: switching cin to binary mode (errCode: %d)\n", strerror(errno));
            exit(EXIT_FAILURE);
        }
#endif
        inStream = &cin;
    } else {
        inStream = new ifstream(input_filename, ios_base::in | ios::binary);
        if (!*inStream) {
            cerr << "Error: unable to open input file: " << input_filename << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    header_symbol = '>';
    if (ffc_header.max_block_size % CHUNK_SIZE) {
        cerr << "Aborting: incorrect block_size: " << ffc_header.max_block_size << endl;
        cerr.flush();
        exit(EXIT_FAILURE);
    }

    ffc_stats.streams_size = 0;
    ffc_stats.blocks_count = 0;
    ffc_stats.original_file_size = 0;

    thread_blocks_meta.resize(PgHelpers::numberOfThreads);
    thread_block.resize(PgHelpers::numberOfThreads, nullptr);
    thread_all_streams.resize(PgHelpers::numberOfThreads, nullptr);
    thread_case_mask.resize(PgHelpers::numberOfThreads, nullptr);
    thread_raw.resize(PgHelpers::numberOfThreads, nullptr);
    thread_dna.resize(PgHelpers::numberOfThreads, nullptr);
    thread_mix.resize(PgHelpers::numberOfThreads, nullptr);
    thread_subblocks_meta.resize(PgHelpers::numberOfThreads, nullptr);
    thread_zsdt_cctx.resize(PgHelpers::numberOfThreads, nullptr);

    int64_t headers_count = 0;

    ostream* outStream;
    if (output_filename == STANDARD_IO_POSIX_ALIAS) {
#ifdef __MINGW32__
        if (_setmode(_fileno(stdout), _O_BINARY) == -1)
            fprintf(stderr, "WARNING: switching cout to binary mode failed (errCode: %d)\n", strerror(errno));
#endif
        outStream = &cout;
    } else {
        outStream = new ofstream(output_filename, ios::out | ios::binary | ios::trunc);
        if (!*outStream) {
            cerr << "Error: unable to open output file: " << output_filename << endl;
            cerr.flush();
            exit(EXIT_FAILURE);
        }
    }

    write_ffc_header(outStream, &ffc_header);
    outStream->write(input_filename.c_str(), ffc_header.orig_filename_length);
    int64_t compressed_size = FFC_HEADER_SIZE + ffc_header.orig_filename_length;

    const bool adaptive_compression = ffc_header.dna_compression_flag == ADAPTIVE_DNA_COMPRESSION_LEVEL;
    if (adaptive_compression)
        ffc_header.dna_compression_flag = 1;
    block_meta_t next_block_meta = init_block(0, nullptr);
    std::vector<std::future<void>> block_threads;
    int64_t b = 0;
    do {
        if (inStream->peek() != EOF)
            ffc_stats.blocks_count++;
        if (b < ffc_stats.blocks_count) {
            int t_id = get_thread_id(b);
            if (thread_block[t_id] == nullptr)
                thread_block[t_id] = new char[ffc_header.max_block_size];
            inStream->read(thread_block[t_id], ffc_header.max_block_size);
            thread_blocks_meta[t_id] = next_block_meta;
            thread_blocks_meta[t_id].block_length = inStream->gcount();
            ffc_stats.original_file_size += thread_blocks_meta[t_id].block_length;
            if (inStream->peek() != EOF)
                next_block_meta = init_block(b + 1, thread_block[t_id]);
            block_threads.emplace_back(
                        std::async(
                        launch::async,
                        &FASTA_Compress::compress_block<false>,
                        this,
                        b
                    )
                );
        }

        int64_t i = b - (PgHelpers::numberOfThreads - 1);
        if (i >= 0 && i < ffc_stats.blocks_count) {
            if (block_threads[i].valid()) {
                block_threads[i].wait();
                int ti_id = get_thread_id(i % PgHelpers::numberOfThreads);
                if (i == 0 && adaptive_compression) {
                    if (((float) thread_blocks_meta[ti_id].dna_stream_size) / thread_blocks_meta[ti_id].dna_stream_compressed_size
                        < ADAPTIVE_DNA_THRESHOLD_RATIO) {
                        *PgHelpers::appout << "Disabling DNA compression." << endl;
                        ffc_header.dna_compression_flag = 0;
                    }
                }

                headers_count += thread_blocks_meta[ti_id].headers_count;
                write_block_meta_data(outStream, &thread_blocks_meta[ti_id]);
                outStream->write(thread_all_streams[ti_id], thread_blocks_meta[ti_id].block_compressed_length);
                compressed_size += FFC_BLOCK_META_SIZE + thread_blocks_meta[ti_id].block_compressed_length;
                ffc_stats.streams_size += thread_blocks_meta[ti_id].block_compressed_length;
            } else {
                throw "ERROR: Thread is not valid";
            }
        }
    } while (++b < ffc_stats.blocks_count + PgHelpers::numberOfThreads);

    if (input_filename != STANDARD_IO_POSIX_ALIAS)
        delete inStream;

    block_meta_t terminate_block_meta;
    write_block_meta_data(outStream, &terminate_block_meta);
    ffc_stats.number_of_sequences = headers_count;
    outStream->write((char*) &ffc_stats, sizeof(ffc_stats_t));
    compressed_size += FFC_BLOCK_META_SIZE + sizeof(ffc_stats_t);

    if (output_filename != STANDARD_IO_POSIX_ALIAS)
        delete outStream;

    for (int t = 0; t < PgHelpers::numberOfThreads; t++) {
        delete [] thread_block[t];
        delete [] thread_all_streams[t];
        delete [] thread_case_mask[t];
        delete [] thread_dna[t];
        delete [] thread_raw[t];
        delete [] thread_mix[t];
        delete [] thread_subblocks_meta[t];
        ZSTD_freeCCtx(thread_zsdt_cctx[t]);
    }

    *PgHelpers::appout << "Input file size: " << ffc_stats.original_file_size << " bytes" << endl;
    *PgHelpers::appout << "Output file size: " << compressed_size << " bytes" << endl;
}


int compress(
    string input_filename, 
    string output_filename,
    int level,
    int block_size_order,
    bool longMatchingFlag
) {
    ffc_header_t ffc_header;
    ffc_header.ffc_header = 0x6366662e;
    ffc_header.version = FFC_VERSION;

    ffc_header.orig_filename_length = input_filename.length();

    ffc_header.case_compression_flag = 1;
    ffc_header.raw_compression_flag = 1;
    ffc_header.dna_compression_flag = ADAPTIVE_DNA_COMPRESSION_LEVEL;
    ffc_header.mix_compression_flag = 1;
    ffc_header.subblocks_meta_compression_flag = 1;
    ffc_header.chunk_size = CHUNK_SIZE;

    if (output_filename == STANDARD_IO_POSIX_ALIAS) {
        PgHelpers::appout = &null_stream;
    }

    if (input_filename != STANDARD_IO_POSIX_ALIAS && !std::filesystem::exists(input_filename)) {
        *PgHelpers::appout << "Input file does not exist: " << input_filename << endl << flush;
        return EXIT_FAILURE;
    }

    *PgHelpers::appout << "Input file: " << input_filename  << endl;
    *PgHelpers::appout << "Output file: " << output_filename << endl;
    *PgHelpers::appout << "Threads: " << PgHelpers::numberOfThreads << endl;
    *PgHelpers::appout << "Level: " << (level == ADAPTIVE_DNA_COMPRESSION_LEVEL ? "adaptive" : to_string(level)) << endl;
    *PgHelpers::appout << "Block size order: " << block_size_order << endl;
    *PgHelpers::appout << "Long matching: " << (longMatchingFlag ? "enabled" : "disabled")<< endl;

    ffc_header.max_block_size = 1 << block_size_order;
    if (block_size_order == 30)
        ffc_header.max_block_size -= CASE_CHUNK_SIZE;
    if (level != ADAPTIVE_DNA_COMPRESSION_LEVEL) {
        ffc_header.case_compression_flag = level;
        ffc_header.raw_compression_flag = level;
        ffc_header.dna_compression_flag = level;
        ffc_header.mix_compression_flag = level;
        ffc_header.subblocks_meta_compression_flag = level;
    }

    int64_t input_size = std::filesystem::file_size(input_filename);
    int64_t min_blocks_count = input_size / ffc_header.max_block_size;
    uint64_t space_in_bytes_required = input_size + FFC_HEADER_SIZE + ffc_header.orig_filename_length +
        (min_blocks_count << 2) * FFC_BLOCK_META_SIZE;
    if (output_filename != STANDARD_IO_POSIX_ALIAS) {
        PgHelpers::createFolders(output_filename);
        if (!PgHelpers::has_enough_space(output_filename, space_in_bytes_required)) {
            *PgHelpers::appout << "Not enough space for output file." << endl;
            return EXIT_FAILURE;
        }
    }
    auto stopWatch = PgHelpers::time_checkpoint();

    FASTA_Compress compressor(ffc_header, longMatchingFlag ? block_size_order : 0);
    compressor.compress_parallel(
        input_filename,
        output_filename
    );
     *PgHelpers::appout << "Compression finished in: " << PgHelpers::time_millis(stopWatch) << " msec." << endl;

    return 0;
}

