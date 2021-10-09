#include "HestonHullWhiteFFT/util/fft/fft.h"


// assume prices is power of 2
void heston_hull_white::util::internal::fft_radix2_transform(std::vector<double>& prices) {
    const size_t num_prices = prices.size() / 2;
    //std::size_t logn = 0;
    //int status;

    if (num_prices == 1) /* identity operation */
    {
        return;
    }

    // throws if not power of 2
    std::size_t log_num_prices = binary_logn(num_prices);

    bitreverse_order(prices);

    std::size_t dual = 1;


    for (std::size_t bit = 0; bit < log_num_prices; bit++)
    {
        double w_real = 1.0;
        double w_imag = 0.0;

        const double theta = 2.0 * (std::numbers::pi_v<double> / (2.0 * static_cast<double>(dual)));

        const double s = std::sin(theta);
        const double t = std::sin(theta / 2.0);
        const double s2 = 2.0 * t * t;

        /* a = 0 */

        for (size_t b = 0; b < num_prices; b += 2 * dual)
        {
            const size_t i = b;
            const size_t j = b + dual;

            const double z1_real = prices[2 * j];
            const double z1_imag = prices[2 * j + 1];

            const double wd_real = z1_real;
            const double wd_imag = z1_imag;

            prices[2 * j] = prices[2 * i] - wd_real;
            prices[2 * j + 1] = prices[2 * i + 1] - wd_imag;
            prices[2 * i] += wd_real;
            prices[2 * i + 1] += wd_imag;
        }

        /* a = 1 .. (dual-1) */

        for (size_t a = 1; a < dual; a++)
        {

            /* trignometric recurrence for w-> exp(i theta) w */

            {
                const double tmp_real = w_real - s * w_imag - s2 * w_real;
                const double tmp_imag = w_imag + s * w_real - s2 * w_imag;
                w_real = tmp_real;
                w_imag = tmp_imag;
            }

            for (size_t b = 0; b < num_prices; b += 2 * dual)
            {
                const size_t i = b + a;
                const size_t j = b + a + dual;

                const double z1_real = prices[2 * j];
                const double z1_imag = prices[2 * j + 1];

                const double wd_real = w_real * z1_real - w_imag * z1_imag;
                const double wd_imag = w_real * z1_imag + w_imag * z1_real;

                prices[2 * j] = prices[2 * i] - wd_real;
                prices[2 * j + 1] = prices[2 * i + 1] - wd_imag;
                prices[2 * i] += wd_real;
                prices[2 * i + 1] += wd_imag;

                /*REAL(data, stride, j) = REAL(data, stride, i) - wd_real;
                IMAG(data, stride, j) = IMAG(data, stride, i) - wd_imag;
                REAL(data, stride, i) += wd_real;
                IMAG(data, stride, i) += wd_imag;*/
            }
        }
        dual *= 2;
    }

    return;

}