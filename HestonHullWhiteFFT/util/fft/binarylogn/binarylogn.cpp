#include "HestonHullWhiteFFT/util/fft/binarylogn/binarylogn.h"

std::size_t heston_hull_white::util::internal::binary_logn(std::size_t n) {
    size_t binary_logn = 0;
    size_t k = 1;

    while (k < n)
    {
        k *= 2;
        binary_logn++;
    }


    if (k != n) {
        throw;
    }

    return binary_logn;

}