#ifndef LINK_H
#define LINK_H

#include <utility>

struct FullProofReturn {
    char* proofptr;
    int err;
};

extern "C" {
    void hello_world();
    int initParams(int size, const char* file);
    int loadParams(bool check, const char* file);

    FullProofReturn createProofFull(int size, int input_size, int aux_size, uint64_t *values, uint64_t *input, uint64_t *aux);

    // just not going to deal with memory look from pointer to proof;
    // not that much data even leaked anyway
    std::pair<char*, int> createProofPartial(int size, uint64_t *data);
}

#endif /* link.hpp */