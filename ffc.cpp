#include <iostream>
#include <string>
#include <stdexcept>
#include <thread>
#include <tuple>
#include "version.hpp"
#include "CLI/CLI.hpp"
#include "FASTA_Compress.h"
#include "FASTA_Decompress.h"

#define OPTION_NOT_SET -1
#define DEF_NUM_COMPRESSION_THREADS 12
#define DEF_NUM_DECOMPRESSION_THREADS 4

std::tuple<uint8_t, uint8_t, uint16_t> decode_version(uint32_t version) {
    uint8_t major = (version >> 24) & 0xFF;
    uint8_t minor = (version >> 16) & 0xFF;
    uint16_t patch = version & 0xFFFF;
    return std::make_tuple(major, minor, patch);
}

void printVersion() {
    auto [major, minor, patch] = decode_version(FFC_VERSION);
    std::cout << "Fast FASTA Compressor (ffc), v"
              << static_cast<int>(major) << "." 
              << static_cast<int>(minor) << "." 
              << patch
              << ", (c) Sz. Grabowski, T. Kowalski, R. Susik, 2025" << std::endl;
}

bool outputFileOk(const std::string &outputFile, const std::string &inputFile, bool forceOverwrite) {
    if (outputFile == STANDARD_IO_POSIX_ALIAS)
        return true;
    if (std::filesystem::exists(outputFile)) {
        if (!forceOverwrite) {
            std::cerr << "Error: Output file already exists: " << outputFile << ". Use -f/--force to overwrite." << std::endl;
            return false;
        }
    }
    if (outputFile == inputFile) {
        std::cerr << "Error: Refusing to open an output file which is the same as input file." << std::endl;
        return false;
    }
    return true;
}

std::string deriveOutputFilename(const std::string& inputFile) {
    if (inputFile.ends_with(FFC_EXTENSION)) {
        return inputFile.substr(0, inputFile.length() - FFC_EXTENSION.length());
    }
    throw std::runtime_error("Unknown file extension (" + FFC_EXTENSION + 
                           " expected). Can't derive output filename. Use -o <outFilename>.");
}

int main(int argc, char *argv[]) {
    CLI::App app{"Fast FASTA Compressor (ffc)"};

    std::string inputFile, outputFile;
    bool 
        decompressFlag = false, 
        forceOverwriteFlag = false;
    int level = OPTION_NOT_SET, blockSize = OPTION_NOT_SET, thread = OPTION_NOT_SET;

    int maxNumThreads = std::thread::hardware_concurrency();
    int defNumCompressionThreads = maxNumThreads < DEF_NUM_COMPRESSION_THREADS ?
        maxNumThreads : DEF_NUM_COMPRESSION_THREADS;
    int defNumDecompressionThreads = maxNumThreads < DEF_NUM_DECOMPRESSION_THREADS ?
        maxNumThreads : DEF_NUM_DECOMPRESSION_THREADS;

    PgHelpers::numberOfThreads = maxNumThreads;

    app.add_flag("-d,--decompress", decompressFlag, "Decompress mode");
    app.add_flag("-f,--force", forceOverwriteFlag, "Overwrite the output file if exists");
    app.add_option("-i,--input", inputFile, "Input file, use hyphen symbol (-) for stdin")->type_name("FNAME");
    app.add_option("-o,--output", outputFile, "Output file, use hyphen symbol (-) for stdout")->type_name("FNAME");
    app.add_option("-l,--level", level, "Backend compr. level, default: adaptive")
        ->check(CLI::Range(0, 22));
    app.add_option("-b,--block", blockSize, "Block size order, default: " + std::to_string(DEFAULT_BLOCK_SIZE_ORDER))
        ->check(CLI::Range(20, 30));
    app.add_option("-t,--threads", thread, "Number of threads, default: " + std::to_string(defNumCompressionThreads)
        + "c / " + std::to_string(defNumDecompressionThreads) + "d")
        ->check(CLI::PositiveNumber);
    app.get_formatter()->column_width(30);

    bool showVersion = false;
    app.add_flag("-v,--version", showVersion, "Show version information");

    std::vector<std::string> positionalArgs;
    app.add_option("file", positionalArgs, "Input file, use hyphen symbol (-) for stdin")->type_name("FNAME");;


    app.allow_extras(false);

    CLI11_PARSE(app, argc, argv);

    if (showVersion) {
        printVersion();
        return EXIT_SUCCESS;
    }

    if (!inputFile.empty() && !positionalArgs.empty()) {
        std::cerr << "Error: Cannot specify both -i/--input and a positional argument for the input file." << std::endl;
        return EXIT_FAILURE;
    }

    if (positionalArgs.size() > 1) {
        std::cerr << "Error: Too many positional arguments provided." << std::endl;
        return EXIT_FAILURE;
    }

    if (inputFile.empty() && !positionalArgs.empty()) {
        inputFile = positionalArgs[0];
    }

    if (inputFile.empty()) {
        std::cerr << "Error: Input file is required" << std::endl;
        return EXIT_FAILURE;
    }

    PgHelpers::numberOfThreads = (thread == OPTION_NOT_SET) ?
        (decompressFlag ? defNumDecompressionThreads : defNumCompressionThreads) : thread;

    if (decompressFlag) {
        if (level != OPTION_NOT_SET || blockSize != OPTION_NOT_SET) {
            std::cerr << "Error: Compression options (-l, -b) are not allowed in decompression mode" << std::endl;
            return EXIT_FAILURE;
        }
        if (outputFile.empty()) {
            try {
                outputFile = deriveOutputFilename(inputFile);
            } catch (const std::runtime_error& e) {
                std::cerr << e.what() << std::endl;
                return EXIT_FAILURE;
            }
        }
        if (!outputFileOk(outputFile, inputFile, forceOverwriteFlag)) {
            return EXIT_FAILURE;
        }
        decompress(inputFile, outputFile);
    } else if (!decompressFlag) {
        if (outputFile.empty()) {
            outputFile = inputFile + FFC_EXTENSION;
        }
        if (!outputFileOk(outputFile, inputFile, forceOverwriteFlag)) {
            return EXIT_FAILURE;
        }
        if (level == OPTION_NOT_SET) {
            level = ADAPTIVE_DNA_COMPRESSION_LEVEL;
        }
        if (blockSize == OPTION_NOT_SET) {
            blockSize = DEFAULT_BLOCK_SIZE_ORDER;
        }
        compress(inputFile, outputFile, level, blockSize);
    }

    return EXIT_SUCCESS;
}