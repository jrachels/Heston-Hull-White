
#include <vector>

#include "HestonHullWhiteFFT/util/fft/bitreverseorder/bitreverseorder.h"

void heston_hull_white::util::internal::bitreverse_order(std::vector<double>& prices) {
    assert((prices.size() % 2) == 0); // included == 0 to be expresive
    std::size_t num_complex_entries = prices.size() / 2;
    size_t j = 0;

    for (std::size_t i = 0; i < num_complex_entries - 1; i++)
    {
        assert(num_complex_entries % 2 == 0);
        size_t k = num_complex_entries / 2;

        if (i < j)
        {
            std::swap(prices[2 * i], prices[2 * j]);
            std::swap(prices[2 * i + 1], prices[2 * j + 1]);
            //double tmp_real = prices[2 * i];
            //double tmp_imag = prices[2 * i + 1];
            //prices[2 * i] = [2 * i];
            //prices[2 * i + 1] = j;

            //const BASE tmp_real = REAL(data, stride, i);
            //const BASE tmp_imag = IMAG(data, stride, i);
            //REAL(data, stride, i) = REAL(data, stride, j);
            //IMAG(data, stride, i) = IMAG(data, stride, j);
            //REAL(data, stride, j) = tmp_real;
            //IMAG(data, stride, j) = tmp_imag;
        }

        while (k <= j)
        {
            j = j - k;
            k = k / 2;
        }

        j += k;
    }

    return;
}