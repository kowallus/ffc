#ifndef FASTACOMPRESS_H
#define FASTACOMPRESS_H

#include <string>
#include <vector>
#include <cstring>
#include "ffc_structures.hpp"
#include "utils/helper.h"

using namespace std;

class FASTA_Compress {
private:

    constexpr static int8_t BLOCK_STARTS_JUST_AFTER_EOL = -1;
    constexpr static int8_t BLOCK_STARTS_INSIDE_HEADER = -2;
    constexpr static int8_t BLOCK_STARTS_INSIDE_SEQUENCE = -3;

    constexpr static int32_t INCONSISTENT_EOLS = -1;

    constexpr static uint64_t EIGHT_EOLS = 0x0A0A0A0A0A0A0A0A;

    constexpr static float ADAPTIVE_DNA_THRESHOLD_RATIO = 1.25;

    template<typename E>
    int64_t compress_stream(char* dest, vector<E>& stream, int32_t compression_flag);
    int64_t compress_stream(char* dest, char* src, size_t srcSize, int32_t compression_flag);

    // enum block type mapping
    enum class BLOCK_TYPE_MAPPING {
        DNA = 0,
        RAW = 1,
        MIX = 2,
        NNN = 3,
    };

    char header_symbol = '>';

    block_meta_t init_block(int32_t b, char* prev_block);

    inline int64_t remove_EOLs_and_find_eols_period(char *seq, int64_t& seq_len, int8_t chunk_size,
        int64_t& tail_len, bool ignore_first_EOL_period, int32_t declared_first_eol_pos);
    void restore_EOLs(char* seq, int64_t chunked_len, int64_t first_eol_pos, int64_t last_eol_pos, int64_t eols_period);
    void restore_EOLs_in_sequences(char *seq, int64_t len, int64_t first_eol_pos, int64_t eols_period,
        bool block_starts_with_header);

    template<bool binary>
    inline bool all_ACGT(uint64_t chunk);
    inline uint8_t pack_dna_naive(char *chunk);

    ffc_header_t ffc_header;
    ffc_stats_t ffc_stats;

    vector<block_meta_t> thread_blocks_meta;
    vector<char*> thread_block;
    vector<char*> thread_case_mask;
    vector<char*> thread_raw;
    vector<char*> thread_dna;
    vector<char*> thread_mix;
    vector<uint32_t*> thread_subblocks_meta;
    vector<char*> thread_all_streams;

    int get_thread_id(int32_t b) const {
        return b % PgHelpers::numberOfThreads;
    }

    template<bool irregular_EOLs_mode>
    void compress_block(
        int32_t b
    );

    static void write_ffc_header(ostream* outStream, const ffc_header_t* ffc_header);
    static void write_block_meta_data(ostream* outStream, const block_meta_t* block_meta);

public:

    FASTA_Compress(const ffc_header_t& ffc_header): ffc_header(ffc_header) {
    }

    ~FASTA_Compress() {
    }

    void compress_parallel(string input_filename, string output_filename);

};

int compress(
    string input_filename, 
    string output_filename,
    int level,
    int block_size_order
);

#endif //FASTACOMPRESS_H
