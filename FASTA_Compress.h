#ifndef FASTACOMPRESS_H
#define FASTACOMPRESS_H

#include <string>
#include <vector>
#include <cstring>
#include "ffc_structures.hpp"
#include "lib/zstd.h"
#include "utils/helper.h"

using namespace std;

class FASTA_Compress {
private:

    constexpr static int8_t BLOCK_STARTS_JUST_AFTER_EOL = -1;
    constexpr static int8_t BLOCK_STARTS_INSIDE_HEADER = -2;
    constexpr static int8_t BLOCK_STARTS_INSIDE_SEQUENCE = -3;

    constexpr static int32_t INCONSISTENT_EOLS = -1;

    constexpr static int32_t MIN_SEQ_LINE_LENGTH = 4;

    constexpr static uint64_t EIGHT_EOLS = 0x0A0A0A0A0A0A0A0A;

    constexpr static float ADAPTIVE_DNA_THRESHOLD_RATIO = 1.25;

    template<typename E>
    int64_t compress_stream(int t_id, char* dest, vector<E>& stream, int32_t compression_flag, int32_t window_size_log = 0);
    int64_t compress_stream(int t_id, char* dest, char* src, size_t srcSize, int32_t compression_flag, int32_t window_size_log = 0);

    // enum block type mapping
    enum class BLOCK_TYPE_MAPPING {
        DNA = 0,
        RAW = 1,
        MIX = 2,
        NNN = 3,
    };

    char header_symbol = '>';

    int long_matching_window_size_log;

    block_meta_t init_block(int32_t b, char* prev_block);

    inline int32_t remove_EOLs_and_find_tail(char *seq, int32_t& seq_len, int8_t chunk_size,
        int32_t& seq_line_length, int32_t& max_seq_without_eols_length, bool ignore_first_EOL_period, int32_t declared_first_eol_pos);

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
    vector<ZSTD_CCtx_s*> thread_zsdt_cctx;

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

    FASTA_Compress(const ffc_header_t& ffc_header, int long_matching_window_size_log): ffc_header(ffc_header),
        long_matching_window_size_log(long_matching_window_size_log) {
    }

    ~FASTA_Compress() {
    }

    void compress_parallel(string input_filename, string output_filename);

};

int compress(
    string input_filename, 
    string output_filename,
    int level,
    int block_size_order,
    bool longMatchingFlag
);

#endif //FASTACOMPRESS_H
