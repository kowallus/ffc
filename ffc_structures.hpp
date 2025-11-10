#ifndef FFC_STRUCTURES_HPP
#define FFC_STRUCTURES_HPP

#include <cstdint>

const static std::string FFC_EXTENSION = ".ffc";

const static uint8_t NO_CODER = 0;
const static uint8_t ZSTD_CODER = 7;

constexpr static int32_t ADAPTIVE_DNA_COMPRESSION_LEVEL = -1;
constexpr static int64_t DEFAULT_BLOCK_SIZE_ORDER = 22;
constexpr static int64_t CHUNK_SIZE = 8;
constexpr static int64_t CASE_CHUNK_SIZE = 64;
constexpr static int8_t STREAM_HEADER_BYTES = 1;

constexpr static int8_t COMPRESSED_BLOCK_MARGIN_FRACTION_ORDER = 5;

constexpr static int64_t FFC_HEADER_SIZE = 56;

struct ffc_header_t {
    int64_t ffc_header = 0;
    uint32_t version = 0; // 0xMMNNPPPP (MM - major (0-255), NN - minor (0-255), PPPP - patch (0-65535))
    int32_t chunk_size = 0;
    int32_t max_block_size = 0;
    int32_t case_compression_flag = 0;
    int32_t raw_compression_flag = 0;
    int32_t dna_compression_flag = 0;
    int32_t mix_compression_flag = 0;
    int32_t subblocks_meta_compression_flag = 0;
    uint32_t crc = 0;
    int64_t orig_file_timestamp = 0;
    int32_t orig_filename_length = 0;
};

constexpr static int64_t FFC_BLOCK_META_SIZE = 64;

struct block_meta_t {

    int64_t block_start = 0;
    int32_t block_length = 0;
    uint32_t block_compressed_length = 0;

    int32_t case_mask_compressed_size = 0;
    int32_t raw_stream_size = 0;
    int32_t raw_stream_compressed_size = 0;
    int32_t dna_stream_size = 0;
    int32_t dna_stream_compressed_size = 0;
    int32_t mix_stream_size = 0;
    int32_t mix_stream_compressed_size = 0;
    int32_t subblocks_count = 0;
    int32_t subblocks_meta_compressed_size = 0;
    int32_t first_EOL_offset = 0;
    int32_t seq_line_length = 0;

    int32_t headers_count = 0;
};

struct ffc_stats_t {
    int64_t blocks_count = 0;
    int64_t original_file_size = 0;
    int64_t number_of_sequences = 0;
    int64_t streams_size = 0;
};

constexpr uint32_t RAW_SUBBLOCK_TYPE = 0 << 30;
constexpr uint32_t DNA_SUBBLOCK_TYPE = 1 << 30;
constexpr uint32_t MIX_SUBBLOCK_TYPE = 2 << 30;
constexpr uint32_t NNN_SUBBLOCK_TYPE = 3 << 30;

constexpr uint32_t SUBBLOCK_TYPE_MASK = 1 << 30 | 1 << 31;
constexpr uint32_t SUBBLOCK_SIZE_MASK = ~SUBBLOCK_TYPE_MASK;

typedef uint32_t size_and_type_stream_meta_t; // size on 30-bits (largest block), two most signif. bits for type: 0 - DNA, 1 - RAW, 2 - MIX, 3 - NNN

#endif // FFC_STRUCTURES_HPP
