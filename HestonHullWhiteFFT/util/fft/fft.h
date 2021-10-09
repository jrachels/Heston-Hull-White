// based on gsl_fft_complex_radix2_transform

#ifndef HESTONHULLWHITEFFT_UTIL_FFT_FFT_H
#define HESTONHULLWHITEFFT_UTIL_FFT_FFT_H

#include <cmath>
#include <numbers>
#include <vector>

#include "HestonHullWhiteFFT/util/fft/binarylogn/binarylogn.h"
#include "HestonHullWhiteFFT/util/fft/bitreverseorder/bitreverseorder.h"

namespace heston_hull_white {
    namespace util {
        namespace internal {

            // assume prices is power of 2
            void fft_radix2_transform(std::vector<double>& prices);

        }
    }
}

#endif