#include "HestonHullWhiteFFT/util/interestcurve.h"

double heston_hull_white::fourier::util::InterestCurve::ZeroBondCurve(double t) const {
	assert(t >= 0);
	const std::size_t num_dates_minus_one = dates_.size() - 1;
	size_t index_of_greatest_lower_bound = IndexOfGreatestLowerBound_Date(t);

	// extrapolate if past last date
	if (index_of_greatest_lower_bound == (num_dates_minus_one)) {
		return std::exp(-ZeroBondToYield(dates_[index_of_greatest_lower_bound], discount_factor_spline_[index_of_greatest_lower_bound][0] * t));
	}
	// assume in middle. dates earlier than 0
	else {
		double delta_t = t - dates_[index_of_greatest_lower_bound];
		auto [a, b, c, d] = discount_factor_spline_[index_of_greatest_lower_bound];
		return ((d * delta_t + c) * delta_t + b) * delta_t + a;

	}
}

double heston_hull_white::fourier::util::InterestCurve::Theta(double t, double lambda, double eta) const {
	assert(t > 0);
	size_t index_of_greatest_lower_bound = IndexOfGreatestLowerBound_Date(t);
	assert(index_of_greatest_lower_bound < (dates_.size() - 1));
	double delta_t = t - dates_[index_of_greatest_lower_bound];

	auto [a, b, c, d] = discount_factor_spline_[index_of_greatest_lower_bound];


	double discount_factor = ((d * delta_t + c) * delta_t + b) * delta_t + a;
	double first_derivative_of_discount_factor = ((3.0 * d * delta_t) + (2.0 * c)) * delta_t + b;
	double second_derivative_of_discount_factor = (6.0 * d * delta_t) + (2.0 * c);

	double const_0 = 1 / lambda;
	return const_0 * (((-second_derivative_of_discount_factor * discount_factor) + std::pow(first_derivative_of_discount_factor, 2)) / (std::pow(discount_factor, 2))) + \
		(-first_derivative_of_discount_factor / discount_factor) + \
		(std::pow(const_0 * eta, 2) / (2.0) * (1.0 - std::exp(-2.0 * (lambda)*t)));

}


std::size_t heston_hull_white::fourier::util::InterestCurve::IndexOfGreatestLowerBound_Date(double date) const {
	std::vector<double>::const_iterator first = dates_.begin();
	// I try to use standard library algorithms to be optimal. no alternative to upper_bound
	std::vector<double>::const_iterator iter_to_greatest_lower_bound = std::prev(std::upper_bound(first, dates_.end(), date));
	return std::distance(first, iter_to_greatest_lower_bound);
}

double heston_hull_white::fourier::util::InterestCurve::ZeroBondToYield(double t, double ZeroBond) const //conversion function
{
	if (t < 1.0) //spot interest rates for maturities smaller than one year (e.g. LIBOR rates) are simply compounded rates.
		return (1.0 - ZeroBond) / (t * ZeroBond);
	else //annualy compounded rates
		return 1.0 / pow(ZeroBond, 1.0 / t) - 1.0;
}