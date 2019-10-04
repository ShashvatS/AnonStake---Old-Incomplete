#ifndef FFT_H
#define FFT_H

#include <iostream>
#include <chrono>
#include <functional>
#include <omp.h>

#include "FR.hpp"

using std::cin;
using std::cout;
using std::ref;

constexpr size_t get_power_of_two(size_t n) {
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;

    return n;
}

/*
 Below we make use of pseudocode from [CLRS 2n Ed, pp. 864].
 Also, note that it's the caller's responsibility to multiply by 1/N.
 */
template <typename FieldT, size_t n>
void FFT(FieldT a[n], const FieldT &omega) {
    const size_t logn = log2(n);
    if (n != (1u << logn)) {

        // throw Exception("expected n == (1u << logn)");
        throw "expected n == (1u << logn)";
    }

    /* swapping in place (from Storer's book) */
    for (size_t k = 0; k < n; ++k) {
        const size_t rk = libff::bitreverse(k, logn);
        if (k < rk) std::swap(a[k], a[rk]);
    }

    size_t m = 1; // invariant: m = 2^{s-1}
    for (size_t s = 1; s <= logn; ++s) {
        // w_m is 2^s-th root of unity now
        const FieldT w_m = omega ^ (n / (2 * m));

        asm volatile("/* pre-inner */");

        for (int k = 0; k < n; k += (2 * m)) {
            FieldT w = FieldT::one();

            for (size_t j = 0; j < m; ++j) {
                const FieldT t = w * a[k + j + m];
                a[k + j + m] = a[k + j] - t;
                a[k + j] += t;
                w *= w_m;
            }
        }

        asm volatile("/* post-inner */");
        m *= 2;
    }
}

template <typename FieldT, size_t n>
void iFFT(FieldT a[n], const FieldT &omega) {
    FFT<FieldT, n>(a, omega.inverse());

    const FieldT sconst = FieldT(n).inverse();

    for (int i = 0; i < n; ++i) {
        a[i] *= sconst;
    }

    return;
}

template <typename FieldT, size_t n>
void multiply_by_coset(FieldT a[n], const FieldT &g) {
    FieldT u = g;
    for (int i = 1; i < n; ++i) {
        a[i] *= u;
        u *= g;
    }

    return;
}

template <typename FieldT, size_t n>
void cosetFFT(FieldT a[n], const FieldT &omega) {
    const FieldT &g = FieldT::multiplicative_generator;
    multiply_by_coset<FieldT, n>(a, g);

    FFT<FieldT, n>(a, omega);
}

template <typename FieldT, size_t n>
void divide_by_Z_on_coset(FieldT a[], const FieldT &divider, int m) {
    for (int i = 0; i < m; ++i) {
        a[i] *= divider;
    }

    return;
}

template <typename FieldT, size_t n>
void icosetFFT(FieldT a[n], const FieldT &omega) {
    const FieldT &g = FieldT::multiplicative_generator;
    iFFT<FieldT, n>(a, omega.inverse());

    multiply_by_coset<FieldT, n>(a, g.inverse());
    return;
}

template <typename FieldT, size_t m, size_t n>
float singleProof(FieldT polyA[n], FieldT polyB[n], FieldT polyC[n],
                  FieldT result[n], const FieldT &targetDiv,
                  const FieldT &omega) {

    auto start = std::chrono::steady_clock::now();
    double start2 = omp_get_wtime();

    iFFT<FieldT, n>(polyA, omega);
    iFFT<FieldT, n>(polyB, omega);
    iFFT<FieldT, n>(polyC, omega);

    cosetFFT<FieldT, n>(polyA, omega);
    cosetFFT<FieldT, n>(polyB, omega);
    cosetFFT<FieldT, n>(polyC, omega);

    for (int i = 0; i < n; ++i) {
        // result[i] = (polyA[i] * polyB[i]) - polyC[i];
        polyA[i] = polyA[i] * polyB[i] - polyC[i];
    }

    // divide_by_Z_on_coset<FieldT, n>(result, targetDiv, m);
    // icosetFFT<FieldT, n>(result, omega);
    divide_by_Z_on_coset<FieldT, n>(polyA, targetDiv, m);
    icosetFFT<FieldT, n>(polyA, omega);

    auto finish = std::chrono::steady_clock::now();
    double end2 = omp_get_wtime();
    double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(finish -
                                                                  start)
            .count();

    std::cout << "Proof time: " << elapsed_seconds << " or " << end2 - start2
              << '\n';
    return elapsed_seconds;
}

template <typename FieldT, size_t m>
void benchmarkProofCreation(int numProofs) {
    constexpr size_t n = get_power_of_two(m);
    FieldT **polyA, **polyB, **polyC, **result, targetDiv;
    polyA = new FieldT *[numProofs];
    polyB = new FieldT *[numProofs];
    polyC = new FieldT *[numProofs];
    result = new FieldT *[numProofs];

    for (int i = 0; i < numProofs; ++i) {
        polyA[i] = new FieldT[n];
        polyB[i] = new FieldT[n];
        polyC[i] = new FieldT[n];
        result[i] = new FieldT[n];

        for (int j = 0; j < n; ++j) {
            polyA[i][j] = FieldT::random_element();
            polyB[i][j] = FieldT::random_element();
            polyC[i][j] = FieldT::random_element();
        }
    }

    targetDiv = FieldT::random_element();
    FieldT omega = libff::get_root_of_unity<FieldT>(n);

    auto start = std::chrono::steady_clock::now();
    double start2 = omp_get_wtime();

    #pragma omp parallel for num_threads(12)
    for (int i = 0; i < numProofs; ++i) {
        singleProof<FieldT, m, n>(polyA[i], polyB[i], polyC[i], result[i],
                                  targetDiv, omega);
    }

    // do something
    auto finish = std::chrono::steady_clock::now();
    double end2 = omp_get_wtime();
    double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(finish -
                                                                  start)
            .count();
    std::cout << "\nTotal time: " << elapsed_seconds << " or " << end2 - start2 << '\n';

    for (int i = 0; i < numProofs; ++i) {
        delete[] polyA[i], polyB[i], polyC[i], result[i];
    }
    delete[] polyA, polyB, polyC, result;
}

#endif /* fft.hpp */