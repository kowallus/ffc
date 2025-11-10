#!/usr/bin/env python3
"""
FFC Decompressor
Compliant with FFC File Format Specification v1.0
"""

import struct
from typing import BinaryIO, Tuple, List, Optional
from dataclasses import dataclass
import os
import sys
import argparse

try:
    import zstandard as zstd
    ZSTD_AVAILABLE = True
except ImportError:
    ZSTD_AVAILABLE = False
    zstd = None

# Constants
FFC_MAGIC = 0x6366662e
NO_CODER = 0
ZSTD_CODER = 7

# Subblock types
RAW_SUBBLOCK_TYPE = 0x00000000
DNA_SUBBLOCK_TYPE = 0x40000000
MIX_SUBBLOCK_TYPE = 0x80000000
NNN_SUBBLOCK_TYPE = 0xC0000000

TYPE_MASK = 0xC0000000
SIZE_MASK = 0x3FFFFFFF

@dataclass
class FFCHeader:
    """FFC file header structure"""
    magic: int
    version: int
    chunk_size: int
    max_block_size: int
    case_compression_flag: int
    raw_compression_flag: int
    dna_compression_flag: int
    mix_compression_flag: int
    subblocks_meta_compression_flag: int
    crc: int
    orig_file_timestamp: int
    orig_filename_length: int
    orig_filename: str

@dataclass
class BlockMetadata:
    """Block metadata structure"""
    block_start: int
    block_size: int
    block_compressed_size: int
    case_mask_compressed_size: int
    raw_stream_size: int
    raw_stream_compressed_size: int
    dna_stream_size: int
    dna_stream_compressed_size: int
    mix_stream_size: int
    mix_stream_compressed_size: int
    subblocks_count: int
    subblocks_meta_compressed_size: int
    first_eol_offset: int
    seq_line_length: int
    seq_headers_count: int

@dataclass
class CompressionStats:
    """Compression statistics structure"""
    blocks_count: int
    original_file_size: int
    number_of_sequences: int
    streams_size: int

class FFCDecompressor:
    """FFC Decompressor v1.1"""

    def __init__(self, debug=False):
        if ZSTD_AVAILABLE:
            self.zstd_ctx = zstd.ZstdDecompressor()
        else:
            self.zstd_ctx = None
        self.debug = debug

    def debug_print(self, msg):
        """Print debug message if debug mode is enabled"""
        if self.debug:
            print(f"DEBUG: {msg}")

    def read_header(self, file: BinaryIO) -> FFCHeader:
        """
        Read FFC file header.

        Header structure (56 bytes fixed part):
        Offset | Size | Type   | Field
        -------|------|--------|----------------------------------
        0      | 8    | int64  | magic (0x6366662e)
        8      | 4    | uint32 | version
        12     | 4    | int32  | chunk_size
        16     | 4    | int32  | max_block_size
        20     | 4    | int32  | case_compression_flag
        24     | 4    | int32  | raw_compression_flag
        28     | 4    | int32  | dna_compression_flag
        32     | 4    | int32  | mix_compression_flag
        36     | 4    | int32  | subblocks_meta_compression_flag
        40     | 4    | uint32 | crc
        44     | 8    | int64  | orig_file_timestamp
        52     | 4    | int32  | orig_filename_length
        56     | var  | text   | orig_filename (if length > 0)
        """
        # Read fixed-size header (56 bytes)
        header_data = file.read(56)
        if len(header_data) != 56:
            raise ValueError("Invalid FFC file: incomplete header")

        # Unpack header fields according to spec
        # Format: <QIiiiiiiiIQi
        # Q = uint64 (magic), I = uint32 (version), 
        # 7x i = int32 (chunk_size, max_block_size, flags...),
        # I = uint32 (crc), Q = uint64 (timestamp), i = int32 (filename_length)
        fields = struct.unpack('<QIiiiiiiiIQi', header_data)

        header = FFCHeader(
            magic=fields[0],
            version=fields[1],
            chunk_size=fields[2],
            max_block_size=fields[3],
            case_compression_flag=fields[4],
            raw_compression_flag=fields[5],
            dna_compression_flag=fields[6],
            mix_compression_flag=fields[7],
            subblocks_meta_compression_flag=fields[8],
            crc=fields[9],
            orig_file_timestamp=fields[10],
            orig_filename_length=fields[11],
            orig_filename=""
        )

        # Verify magic number
        if header.magic != FFC_MAGIC:
            raise ValueError(f"Invalid FFC magic number: 0x{header.magic:x}, expected 0x{FFC_MAGIC:x}")

        # Read original filename if present
        if header.orig_filename_length > 0:
            filename_data = file.read(header.orig_filename_length)
            if len(filename_data) != header.orig_filename_length:
                raise ValueError("Invalid FFC file: incomplete filename")
            header.orig_filename = filename_data.decode('utf-8', errors='replace')

        self.debug_print(f"Header: filename='{header.orig_filename}', version=0x{header.version:08x}")

        return header

    def read_block_metadata(self, file: BinaryIO) -> Optional[BlockMetadata]:
        """
        Read block metadata (64 bytes).
        Returns None if terminator block (all zeros).

        Format: <qiIiiiiiiiiiiii
        q = int64 (block_start)
        i = int32 (block_size)
        I = uint32 (block_compressed_size)
        13x i = int32 (remaining fields)
        """
        metadata_data = file.read(64)
        if len(metadata_data) != 64:
            raise ValueError("Invalid FFC file: incomplete block metadata")

        fields = struct.unpack('<qiIiiiiiiiiiiii', metadata_data)

        # Check for terminator block (all zeros)
        if all(field == 0 for field in fields):
            self.debug_print("Terminator block found")
            return None

        return BlockMetadata(*fields)

    def decode_stream(self, file: BinaryIO, compressed_size: int) -> bytes:
        """
        Decode stream: read compression header (1 byte) and optionally decompress.

        Compression methods:
        - 0 (NO_CODER): uncompressed
        - 7 (ZSTD_CODER): zstandard compressed
        """
        if compressed_size == 0:
            return b''

        # Read 1-byte compression method header
        compression_header = file.read(1)
        if len(compression_header) != 1:
            raise ValueError("Invalid stream: missing compression header")

        compression_method = compression_header[0]
        data_size = compressed_size - 1

        # Read stream data
        stream_data = file.read(data_size)
        if len(stream_data) != data_size:
            raise ValueError(f"Invalid stream: incomplete data")

        # Decompress based on method
        if compression_method == NO_CODER:
            return stream_data
        elif compression_method == ZSTD_CODER:
            if not ZSTD_AVAILABLE:
                raise RuntimeError("zstandard library not available for decompression")
            return self.zstd_ctx.decompress(stream_data)
        else:
            raise ValueError(f"Unknown compression method: {compression_method}")

    def unpack_dna_stream(self, dna_data: bytes) -> str:
        """
        Unpack DNA stream from 2-bit packed format.

        Encoding: A=0, C=1, T=2, G=3
        Each byte contains 4 nucleotides (2 bits each):
        - bits 0-1: 1st nucleotide
        - bits 2-3: 2nd nucleotide
        - bits 4-5: 3rd nucleotide
        - bits 6-7: 4th nucleotide
        """
        dna_chars = []
        nucleotides = ['A', 'C', 'T', 'G']

        for byte in dna_data:
            # Extract 4 nucleotides from each byte
            for shift in range(0, 8, 2):
                nucleotide_idx = (byte >> shift) & 0x03
                dna_chars.append(nucleotides[nucleotide_idx])

        return ''.join(dna_chars)

    def read_subblocks_meta(self, meta_data: bytes) -> List[Tuple[int, int]]:
        """
        Parse subblocks metadata.

        Each entry is a 32-bit value:
        - bits 30-31: subblock type (00=RAW, 01=DNA, 10=MIX, 11=NNN)
        - bits 0-29: subblock size in bytes
        """
        subblocks = []
        for i in range(0, len(meta_data), 4):
            if i + 4 > len(meta_data):
                break

            meta_entry = struct.unpack('<I', meta_data[i:i+4])[0]
            subblock_type = meta_entry & TYPE_MASK
            subblock_size = meta_entry & SIZE_MASK

            subblocks.append((subblock_type, subblock_size))

        return subblocks

    def apply_case_flags(self, data: bytearray, case_flags: bytes):
        """
        Apply case flags to data according to spec v1.0.

        Each 8-byte chunk of case flags corresponds to 64 bytes of data.
        Uses special bit manipulation as described in spec pseudocode.
        """
        if not case_flags or all(b == 0 for b in case_flags):
            self.debug_print("Case flags: all zeros, skipping")
            return

        self.debug_print(f"Applying case flags: {len(case_flags)} bytes to {len(data)} bytes")

        # Process each 8-byte chunk of case flags
        for chunk_idx in range(0, len(case_flags), 8):
            # Get 8-byte case flag chunk (pad if needed)
            if chunk_idx + 8 > len(case_flags):
                flags_chunk = case_flags[chunk_idx:]
                flags_chunk += b'\x00' * (8 - len(flags_chunk))
            else:
                flags_chunk = case_flags[chunk_idx:chunk_idx + 8]

            # Unpack as 64-bit integer
            flags = struct.unpack('<Q', flags_chunk)[0]
            data_base = (chunk_idx // 8) * 64

            # Apply case flag transformations according to spec
            masks = [
                (0x0101010101010101, 5),  # result[0]
                (0x0202020202020202, 4),  # result[1]
                (0x0404040404040404, 3),  # result[2]
                (0x0808080808080808, 2),  # result[3]
                (0x1010101010101010, 1),  # result[4]
                (0x2020202020202020, 0),  # result[5]
                (0x4040404040404040, -1), # result[6]
                (0x8080808080808080, -2), # result[7]
            ]

            for result_idx, (mask, shift) in enumerate(masks):
                pos = data_base + result_idx * 8

                if pos + 8 > len(data):
                    break

                # Read 8-byte chunk from data
                chunk_bytes = struct.unpack('<Q', bytes(data[pos:pos+8]))[0]

                # Apply transformation
                if shift >= 0:
                    modified = chunk_bytes | ((flags & mask) << shift)
                else:
                    modified = chunk_bytes | ((flags & mask) >> (-shift))

                # Write back
                data[pos:pos+8] = struct.pack('<Q', modified)

    def reconstruct_block(self, raw_data: bytes, dna_chars: str, mix_data: bytes,
                         subblocks: List[Tuple[int, int]],
                         first_eol_offset: int, seq_line_length: int,
                         block_size: int) -> bytearray:
        """
        Reconstruct block from subblocks with EOL insertion.

        Key rules:
        1. Build data WITHOUT automatic EOLs first
        2. Automatic EOL after each RAW subblock if it fits within block_size
        3. Periodic EOLs every seq_line_length bytes, ONLY in DNA/MIX/NNN subblocks
        4. first_eol_offset is position in data without auto-EOLs
        """

        self.debug_print(f"\nReconstruct: first_eol={first_eol_offset}, seq_line_length={seq_line_length}, block_size={block_size}")

        # PHASE 1: Build data WITHOUT automatic EOLs after RAW
        data_without_auto_eol = bytearray()
        position_to_type = []  # Track subblock type for each position
        raw_pos = dna_pos = mix_pos = 0
        subblock_end_positions = []  # Track where each subblock ends
        current_pos = 0

        for i, (subblock_type, subblock_size) in enumerate(subblocks):
            if subblock_type == RAW_SUBBLOCK_TYPE:
                # RAW subblock: copy from raw stream
                if raw_pos + subblock_size <= len(raw_data):
                    chunk = raw_data[raw_pos:raw_pos + subblock_size]
                    data_without_auto_eol.extend(chunk)
                    position_to_type.extend([RAW_SUBBLOCK_TYPE] * len(chunk))
                    raw_pos += subblock_size
                    # Mark end position for automatic EOL
                    subblock_end_positions.append((current_pos + len(chunk), True))
                    self.debug_print(f"  Subblock {i}: RAW({subblock_size})")

            elif subblock_type == DNA_SUBBLOCK_TYPE:
                # DNA subblock: unpack from DNA stream
                if dna_pos + subblock_size <= len(dna_chars):
                    chunk = dna_chars[dna_pos:dna_pos + subblock_size]
                    data_without_auto_eol.extend(chunk.encode('ascii'))
                    position_to_type.extend([DNA_SUBBLOCK_TYPE] * subblock_size)
                    dna_pos += subblock_size
                    subblock_end_positions.append((current_pos + subblock_size, False))
                    self.debug_print(f"  Subblock {i}: DNA({subblock_size})")

            elif subblock_type == MIX_SUBBLOCK_TYPE:
                # MIX subblock: copy from MIX stream
                if mix_pos + subblock_size <= len(mix_data):
                    chunk = mix_data[mix_pos:mix_pos + subblock_size]
                    data_without_auto_eol.extend(chunk)
                    position_to_type.extend([MIX_SUBBLOCK_TYPE] * len(chunk))
                    mix_pos += subblock_size
                    subblock_end_positions.append((current_pos + len(chunk), False))
                    self.debug_print(f"  Subblock {i}: MIX({subblock_size})")

            elif subblock_type == NNN_SUBBLOCK_TYPE:
                # NNN subblock: generate N's
                data_without_auto_eol.extend(b'N' * subblock_size)
                position_to_type.extend([NNN_SUBBLOCK_TYPE] * subblock_size)
                subblock_end_positions.append((current_pos + subblock_size, False))
                self.debug_print(f"  Subblock {i}: NNN({subblock_size})")

            current_pos += subblock_size

        self.debug_print(f"Data without auto-EOL: {len(data_without_auto_eol)} bytes")

        # PHASE 2: Calculate ALL EOL positions
        eol_positions = set()

        # First EOL in block (only if in DNA/MIX/NNN)
        if 0 <= first_eol_offset < len(position_to_type):
            if position_to_type[first_eol_offset] in [DNA_SUBBLOCK_TYPE, MIX_SUBBLOCK_TYPE, NNN_SUBBLOCK_TYPE]:
                eol_positions.add(first_eol_offset)
                self.debug_print(f"First EOL at position {first_eol_offset}")

                # Periodic EOLs from first EOL (only in DNA/MIX/NNN)
                if seq_line_length > 0:
                    # Find next automatic EOL position (to know when to stop)
                    next_auto_eol = None
                    for end_pos, needs_eol in subblock_end_positions:
                        if needs_eol and end_pos > first_eol_offset:
                            next_auto_eol = end_pos
                            break

                    pos = first_eol_offset + seq_line_length
                    while pos < len(data_without_auto_eol):
                        # Stop at next automatic EOL
                        if next_auto_eol is not None and pos >= next_auto_eol:
                            break

                        # Only add if position is in DNA/MIX/NNN (not RAW)
                        if pos < len(position_to_type):
                            if position_to_type[pos] in [DNA_SUBBLOCK_TYPE, MIX_SUBBLOCK_TYPE, NNN_SUBBLOCK_TYPE]:
                                eol_positions.add(pos)
                                self.debug_print(f"Periodic EOL at {pos}")

                        pos += seq_line_length

        # Automatic EOLs after each RAW subblock
        for end_pos, needs_eol in subblock_end_positions:
            if needs_eol and end_pos <= len(data_without_auto_eol):
                eol_positions.add(end_pos)
                self.debug_print(f"Auto-EOL after RAW at {end_pos}")

                # Periodic EOLs from this auto-EOL (only in DNA/MIX/NNN)
                if seq_line_length > 0:
                    # Find next automatic EOL
                    next_auto_eol = None
                    for next_end_pos, next_needs_eol in subblock_end_positions:
                        if next_needs_eol and next_end_pos > end_pos:
                            next_auto_eol = next_end_pos
                            break

                    pos = end_pos + seq_line_length
                    while pos < len(data_without_auto_eol):
                        if next_auto_eol is not None and pos >= next_auto_eol:
                            break

                        # Only in DNA/MIX/NNN
                        if pos < len(position_to_type):
                            if position_to_type[pos] in [DNA_SUBBLOCK_TYPE, MIX_SUBBLOCK_TYPE, NNN_SUBBLOCK_TYPE]:
                                eol_positions.add(pos)
                                self.debug_print(f"Periodic EOL after auto at {pos}")

                        pos += seq_line_length

        # PHASE 3: Build final result with EOLs inserted
        result = bytearray()
        for i, byte in enumerate(data_without_auto_eol):
            # Insert EOL before this position if needed
            if i in eol_positions:
                result.append(ord('\n'))
            result.append(byte)

        # Insert last automatic EOL at the end of the block if needed
        if len(result) < block_size and len(data_without_auto_eol) in eol_positions:
            result.append(ord('\n'))
            
        self.debug_print(f"Final block: {len(result)} bytes (expected: {block_size})")

        return result

    def decompress_block(self, file: BinaryIO, block_meta: BlockMetadata) -> bytes:
        """Decompress a single block"""

        self.debug_print(f"\n=== Block: size={block_meta.block_size} ===")

        # Decode all streams
        case_mask_data = self.decode_stream(file, block_meta.case_mask_compressed_size)
        raw_stream_data = self.decode_stream(file, block_meta.raw_stream_compressed_size)
        dna_stream_data = self.decode_stream(file, block_meta.dna_stream_compressed_size)
        mix_stream_data = self.decode_stream(file, block_meta.mix_stream_compressed_size)
        subblocks_meta_data = self.decode_stream(file, block_meta.subblocks_meta_compressed_size)

        # Parse subblocks
        subblocks = self.read_subblocks_meta(subblocks_meta_data)

        # Unpack DNA stream
        dna_chars = self.unpack_dna_stream(dna_stream_data) if dna_stream_data else ""

        # Reconstruct block with EOL insertion
        result = self.reconstruct_block(
            raw_stream_data, dna_chars, mix_stream_data, subblocks,
            block_meta.first_eol_offset, block_meta.seq_line_length, block_meta.block_size
        )

        # Apply case flags (after EOL insertion)
        self.apply_case_flags(result, case_mask_data)

        return bytes(result)

    def read_compression_stats(self, file: BinaryIO) -> CompressionStats:
        """Read compression statistics (32 bytes at end of file)"""
        stats_data = file.read(32)
        if len(stats_data) != 32:
            raise ValueError("Invalid FFC file: incomplete compression statistics")

        fields = struct.unpack('<QQQQ', stats_data)
        return CompressionStats(*fields)

    def decompress_file(self, input_path: str, output_path: str) -> CompressionStats:
        """Decompress FFC file"""

        with open(input_path, 'rb') as infile:
            # Read header
            header = self.read_header(infile)
            print(f"Decompressing: {header.orig_filename or 'unknown'}")
            print(f"  Version: 0x{header.version:08x}")
            print(f"  Chunk size: {header.chunk_size}")
            print(f"  Max block size: {header.max_block_size}")

            # Process blocks
            decompressed_data = bytearray()
            block_count = 0

            while True:
                block_meta = self.read_block_metadata(infile)
                if block_meta is None:  # Terminator block
                    break

                block_data = self.decompress_block(infile, block_meta)
                decompressed_data.extend(block_data)
                block_count += 1

                if block_count % 10 == 0 or self.debug:
                    print(f"Block {block_count}: {len(block_data)} bytes")

            # Read compression statistics
            stats = self.read_compression_stats(infile)

            # Write output
            with open(output_path, 'wb') as outfile:
                outfile.write(decompressed_data)

            print(f"\nDecompression complete:")
            print(f"  Blocks processed: {block_count}")
            print(f"  Output size: {len(decompressed_data)} bytes")
            print(f"  Expected size: {stats.original_file_size} bytes")

            if len(decompressed_data) != stats.original_file_size:
                print(f"  WARNING: Size mismatch!")

            return stats

def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(description='FFC Decompressor v1.1')
    parser.add_argument('input', help='Input FFC file')
    parser.add_argument('output', help='Output file')
    parser.add_argument('--debug', '-d', action='store_true', help='Enable debug output')

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: File not found: {args.input}")
        sys.exit(1)

    decompressor = FFCDecompressor(debug=args.debug)

    if not ZSTD_AVAILABLE:
        print("Warning: zstandard library not available")
        print("         Only uncompressed streams will be supported")

    try:
        decompressor.decompress_file(args.input, args.output)
        print(f"\nSuccess! Output written to: {args.output}")
    except Exception as e:
        print(f"\nError: {e}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
