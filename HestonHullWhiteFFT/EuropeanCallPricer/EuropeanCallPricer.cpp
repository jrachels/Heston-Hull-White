#include "HestonHullWhiteFFT/EuropeanCallPricer/EuropeanCallPricer.h"

heston_hull_white::fourier::EuropeanCallPricer::EuropeanCallPricer(const HestonHullWhiteInputs& hhw_inputs, const FFTInputs& fft_inputs/*, double itsScale = 1.0 */) : \
fft_b_(std::numbers::pi_v<double> / fft_inputs.integration_epsilon_), fft_lambda_((2.0 * std::numbers::pi_v<double>) / (fft_inputs.number_grid_points_ * fft_inputs.integration_epsilon_)), \
itsScale_(hhw_inputs.itsScale_), alpha_(AlphaOptimizatizer(hhw_inputs)), prices_(2 * fft_inputs.number_grid_points_) {
	assert(fft_inputs.number_grid_points_ % 2 == 0);

	// step 1:  get optimal alpha
	//double alpha = AlphaOptimizatizer(hhw_inputs);

	// step 2: run this code
	double integration_epsilon = fft_inputs.integration_epsilon_;
	std::complex<double> evaluation_point{ 0, -(alpha_ + 1.0) };
	std::size_t num_grid_points_divided_by_2 = fft_inputs.number_grid_points_ / 2;
	// do 0 case
	std::complex<double> numerator = HestonHullWhiteCharacteristicFunction(hhw_inputs, evaluation_point);
	std::complex<double> denominator = evaluation_point * (evaluation_point + std::complex<double>{0, 1});
	std::complex<double> integrand = (numerator / denominator) * integration_epsilon;
	prices_[0] = integrand.real();	//Pricing array contains the discrete grid that will be Fast Fourier transformed
	prices_[1] = integrand.imag();
	evaluation_point += integration_epsilon;


	numerator = HestonHullWhiteCharacteristicFunction(hhw_inputs, evaluation_point);
	denominator = evaluation_point * (evaluation_point + std::complex<double>{0, 1});
	integrand = (numerator / denominator) * integration_epsilon;
	integrand *= (4.0 / 3.0);
	prices_[2] = integrand.real();	//Pricing array contains the discrete grid that will be Fast Fourier transformed
	prices_[3] = integrand.imag();
	evaluation_point += integration_epsilon;

	for (std::size_t j = 1; j < num_grid_points_divided_by_2; ++j) {
		// do even

		std::complex<double> numerator = HestonHullWhiteCharacteristicFunction(hhw_inputs, evaluation_point);
		std::complex<double> denominator = -evaluation_point * (evaluation_point + std::complex<double>{0, 1});
		std::complex<double> integrand = (numerator / denominator) * integration_epsilon;
		integrand *= (2.0 / 3.0);
		prices_[4 * j] = integrand.real();	//Pricing array contains the discrete grid that will be Fast Fourier transformed
		prices_[4 * j + 1] = integrand.imag();
		evaluation_point += integration_epsilon;

		// do odd
		numerator = HestonHullWhiteCharacteristicFunction(hhw_inputs, evaluation_point);
		denominator = evaluation_point * (evaluation_point + std::complex<double>{0, 1});
		integrand = (numerator / denominator) * integration_epsilon;
		integrand *= (4.0 / 3.0);
		prices_[4 * j + 2] = integrand.real();	//Pricing array contains the discrete grid that will be Fast Fourier transformed
		prices_[4 * j + 3] = integrand.imag();
		evaluation_point += integration_epsilon;
	}

	heston_hull_white::util::internal::fft_radix2_transform(prices_);

}

double heston_hull_white::fourier::EuropeanCallPricer::operator()(double strike_) {
	// get index of strike in grid
	// 
	// interpolate
	// 
	// return price

	double itsk = std::log(strike_ / itsScale_); //prices as units of scale
	std::size_t itsIndexOfStrike = static_cast<std::size_t>((itsk + fft_b_) / fft_lambda_); //Compare documentation, substitution below equation 3.9

	// i should cache (std::exp(-alpha_ * itsk)
	double itsTempPrice1 = itsScale_ * (std::exp(-alpha_ * itsk) / std::numbers::pi_v<double>) * prices_[2 * itsIndexOfStrike];	//Price will be interpolated
	double itsTempPrice2 = itsScale_ * (std::exp(-alpha_ * itsk) / std::numbers::pi_v<double>) * prices_[2 * itsIndexOfStrike + 2];

	double K1 = itsScale_ * std::exp(-fft_b_ + fft_lambda_ * itsIndexOfStrike); //rescaling
	double K2 = itsScale_ * std::exp(-fft_b_ + fft_lambda_ * (itsIndexOfStrike + 1));

	//itsTempPrice = ((itsTempPrice2 - itsTempPrice1) / (K2 - K1)) * (strike_ - K1) + itsTempPrice1;
	return  ((itsTempPrice2 - itsTempPrice1) / (K2 - K1)) * (strike_ - K1) + itsTempPrice1;
	//return itsTempPrice; //put call parity

}