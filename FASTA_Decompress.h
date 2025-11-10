#ifndef FASTA_DECOMPRESS_H
#define FASTA_DECOMPRESS_H

#include <fstream>
#include <vector>
#include <string>
#include "ffc_structures.hpp"
#include "utils/helper.h"

using namespace std;

class FASTA_Decompress {
private:
    constexpr static int64_t BUFFER_SIZE = 1 << 16;

    bool sequential_input = false;
    bool sequential_output = false;

    ffc_header_t ffc_header;
    int64_t ffc_size;
    string input_filename;
    string output_filename;
    std::ostream* outStream;

    uint32_t *lut;
    uint8_t pack_dna_no_lut(int32_t chunk);
    void init_LUT();
    void unpack_dna_avx2(char *dest, uint8_t *src);

    uint8_t *open_stream(char *decompressed_dest, uint8_t *buff, int64_t buff_size, int64_t dest_size);

    vector<block_meta_t> thread_blocks_meta;
    vector<char*> thread_buf;
    vector<char*> thread_block;

    vector<char*> thread_all_streams;
    vector<char*> thread_case_mask;
    vector<char*> thread_raw;
    vector<char*> thread_dna;
    vector<char*> thread_mix;
    vector<uint32_t*> thread_subblocks_meta;

    int get_thread_id(int32_t b) const {
        return b % PgHelpers::numberOfThreads;
    }

    void decompress_block(
        int32_t b,
        size_t block_start
    );

public:
    FASTA_Decompress() {
        init_LUT();
    }

    ~FASTA_Decompress() {
        free(lut);
    }

    void decompress_parallel(string input_filename, string output_file);

    static void read_ffc_header(istream* inStream, const ffc_header_t* ffc_header);
    static void read_block_meta_data(istream* inStream, const block_meta_t* block_meta);

};

int decompress(string input_filename, string output_filename);

#endif //FASTA_DECOMPRESS_H
