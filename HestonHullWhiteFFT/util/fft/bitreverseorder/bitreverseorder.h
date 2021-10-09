// based on fft_complex_bitreverse_order

#ifndef HESTONHULLWHITEFFT_UTIL_FFT_BITREVERSEORDER_BITREVERSEORDER_H
#define HESTONHULLWHITEFFT_UTIL_FFT_BITREVERSEORDER_BITREVERSEORDER_H

#include <cassert>
#include <vector>

namespace heston_hull_white {
    namespace util {
        namespace internal {


            void bitreverse_order(std::vector<double>& prices);

        }
    }
}


#endif