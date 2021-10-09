#ifndef HESTONHULLWHITEFFT_UTIL_FFT_BINARYLOGN_BINARYLOGN_H
#define HESTONHULLWHITEFFT_UTIL_FFT_BINARYLOGN_BINARYLOGN_H

#include <cstddef>

namespace heston_hull_white {
    namespace util {
        namespace internal {

            std::size_t binary_logn(std::size_t n);

        }
    }
}

#endif