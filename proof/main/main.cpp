#include <gmp.h>
#include <iostream>
#include <libfqfft/evaluation_domain/get_evaluation_domain.hpp>
#include <memory>
#include <omp.h>
#include <tuple>

#include "FR.hpp"
#include "fft.hpp"
#include "link.hpp"

#ifdef TESTS
#include "tests.hpp"
#endif

using libfqfft::evaluation_domain;
using libfqfft::get_evaluation_domain;
using std::shared_ptr;

template <typename FieldT, size_t numConstraints, size_t inputSize>
std::tuple<FieldT **, FieldT **, FieldT **, FieldT **>
generateRandomData(int numProofs) {

    constexpr size_t n = get_power_of_two(numConstraints);

    FieldT **polyA, **polyB, **polyC;
    FieldT **inputs;

    polyA = new FieldT *[numProofs];
    polyB = new FieldT *[numProofs];
    polyC = new FieldT *[numProofs];
    inputs = new FieldT *[numProofs];

    for (int i = 0; i < numProofs; ++i) {
        polyA[i] = new FieldT[n];
        polyB[i] = new FieldT[n];
        polyC[i] = new FieldT[n];
        inputs[i] = new FieldT[inputSize];

        for (int j = 0; j < n; ++j) {
            polyA[i][j] = FieldT::random_element();
            polyB[i][j] = FieldT::random_element();
            polyC[i][j] = FieldT::random_element();
        }

        for (int j = 0; j < inputSize; ++j) {
            inputs[i][j] = FieldT::random_element();
        }
    }

    return std::make_tuple(polyA, polyB, polyC, inputs);
}

template <typename FieldT, size_t numConstraints>
void fftData(FieldT polyA[], FieldT polyB[], FieldT polyC[],
             const FieldT &omega, const FieldT &targetDiv) {
    constexpr size_t n = get_power_of_two(numConstraints);

    iFFT<FieldT, n>(polyA, omega);
    iFFT<FieldT, n>(polyB, omega);
    iFFT<FieldT, n>(polyC, omega);

    cosetFFT<FieldT, n>(polyA, omega);
    cosetFFT<FieldT, n>(polyB, omega);
    cosetFFT<FieldT, n>(polyC, omega);

    for (int i = 0; i < n; ++i) {
        polyA[i] = polyA[i] * polyB[i] - polyC[i];
    }

    divide_by_Z_on_coset<FieldT, n>(polyA, targetDiv, numConstraints);
    icosetFFT<FieldT, n>(polyA, omega);
}

template <typename FieldT, size_t numConstraints, size_t num_inputs,
          size_t num_aux>
void createProofs(std::tuple<FieldT **, FieldT **, FieldT **, FieldT **> data,
                  int numProofs, const FieldT &omega, const FieldT &targetDiv) {
    double start = omp_get_wtime();

    FieldT **polyA, **polyB, **polyC, **inputs;
    std::tie(polyA, polyB, polyC, inputs) = data;

#pragma omp parallel for num_threads(12)
    for (int i = 0; i < numProofs; ++i) {
        fftData<FieldT, numConstraints>(polyA[i], polyB[i], polyC[i], omega,
                                        targetDiv);
    }

    double end1 = omp_get_wtime();

#pragma omp barrier

    for (int i = 0; i < numProofs; ++i) {
        FieldT *aux = inputs[i];
        FieldT *input = inputs[i] + num_aux;
        FullProofReturn result = createProofFull(
            numConstraints, num_inputs, num_aux, (uint64_t *)polyA[i],
            (uint64_t *)input, (uint64_t *)aux);
    }

    double end2 = omp_get_wtime();

    std::cout << "Time: " << (end1 - start) << ' ' << (end2 - start) << '\n';
}

int main() {
    initFieldR();
    // initParams(1 << 17, "params.txt");
    loadParams(false, "params.txt");
    hello_world();

    const int numConstraints = 120000;
    const int totalInputs = 90000;
    const int numAux = 5;
    const int numInputs = totalInputs - numAux;

    auto domain = get_evaluation_domain<FieldR>(numConstraints);
    FieldR omega =
        libff::get_root_of_unity<FieldR>(get_power_of_two(numConstraints));
    FieldR coset = FieldR::multiplicative_generator;
    const FieldR Z_inverse_at_coset =
        domain->compute_vanishing_polynomial(coset).inverse();
    auto data = generateRandomData<FieldR, numConstraints, totalInputs>(12);

    createProofs<FieldR, numConstraints, numInputs, numAux>(data, 12, omega,
                                                            Z_inverse_at_coset);
}