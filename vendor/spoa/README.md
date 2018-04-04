# Spoa

[![Latest GitHub release](https://img.shields.io/github/release/rvaser/spoa.svg)](https://github.com/rvaser/spoa/releases/latest)
[![Build status for c++/clang++](https://travis-ci.org/rvaser/spoa.svg?branch=master)](https://travis-ci.org/rvaser/spoa)
[![Published in Genome Research](https://img.shields.io/badge/published%20in-Genome%20Research-blue.svg)](https://doi.org/10.1101/gr.214270.116)

Spoa (SIMD POA) is a c++ implementation of the partial order alignment (POA) algorithm (as described in 10.1093/bioinformatics/18.3.452) which is used to generate consensus sequences (as described in 10.1093/bioinformatics/btg109). It supports three alignment modes: local (Smith-Waterman), global (Needleman-Wunsch) and semi-global alignment (overlap). It supports Intel SSE4.1+ and AVX2 (marginally faster due to high latency shifts) vectorization.

## Dependencies

### Linux

Application uses following software:

1. gcc 4.8+ or clang 3.4+
2. cmake 3.2+

## Installation

CmakeLists is provided in the project root folder. By running the following commands:

```bash
git clone --recursive https://github.com/rvaser/spoa spoa
cd spoa
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
a library named `libspoa.a` will appear in the `build/lib` directory. If you want the spoa executable, run the following two commands:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -Dspoa_build_executable=ON ..
make
```
which will place an executable named `spoa` in `build/bin` directory.

Optionally, you can run `sudo make install` to install spoa library (and executable) to your machine.

***Note***: if you omitted `--recursive` from `git clone`, run `git submodule init` and `git submodule update` before proceeding with compilation.

## Usage

Usage of spoa is as following:

    spoa [options ...] <sequences>

        <sequences>
            input file in FASTA/FASTQ format containing sequences

        options:
            -m, --match <int>
                default: 5
                score for matching bases
            -x, --mismatch <int>
                default: -4
                score for mismatching bases
            -g, --gap <int>
                default: -8
                gap penalty (must be negative)
            -l, --algorithm <int>
                default: 0
                alignment mode:
                    0 - local (Smith-Waterman)
                    1 - global (Needleman-Wunsch)
                    2 - semi-global
            -r, --result <int>
                default: 0
                result mode:
                    0 - consensus
                    1 - multiple sequence alignment
                    2 - 0 & 1
            --version
                prints the version number
            -h, --help
                prints the usage

### Library

Simple library usage can be seen in the following `example.cpp` file. This code shows how to get consensus and multiple sequence alignment for a set of sequences without quality values.

```cpp
#include "spoa/spoa.hpp"

int main(int argc, char** argv) {

    std::vector<std::string> sequences = {
        "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
        "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
        "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
        "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
        "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
        "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
    };

    auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(atoi(argv[1])),
        atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));

    auto graph = spoa::createGraph();

    for (const auto& it: sequences) {
        auto alignment = alignment_engine->align_sequence_with_graph(it, graph);
        graph->add_alignment(alignment, it);
    }

    std::string consensus = graph->generate_consensus();

    fprintf(stderr, "Consensus (%zu)\n", consensus.size());
    fprintf(stderr, "%s\n", consensus.c_str());

    std::vector<std::string> msa;
    graph->generate_multiple_sequence_alignment(msa);

    fprintf(stderr, "Multiple sequence alignment\n");
    for (const auto& it: msa) {
        fprintf(stderr, "%s\n", it.c_str());
    }

    return 0;
}
```

This code can be compiled from spoa root directory with:
```bash
g++ example.cpp -std=c++11 -Iinclude/ -Lbuild/lib/ -lspoa -o example
```
or with the following command if spoa was installed beforehand:
```bash
g++ example.cpp -std=c++11 -lspoa -o example
```

The executable can be run with:
```bash
./example 0 5 -4 -8
```

On the other hand, if you are using `cmake` you can add spoa to your project by adding commands `add_subdirectory(vendor/spoa EXCLUDE_FROM_ALL)` and `target_link_libraries(your_exe spoa)` to your main CMakeLists file.

## Contact information

For additional information, help and bug reports please send an email to: robert.vaser@fer.hr.

## Acknowledgement

This work has been supported in part by Croatian Science Foundation under the project UIP-11-2013-7353.
