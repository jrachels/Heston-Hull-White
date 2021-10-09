#include "HestonHullWhiteFFT/AlphaOptimization/alphaoptimization.h"

double heston_hull_white::fourier::internal::AlphaObjectiveFunction(const HestonHullWhiteInputs& hhw_inputs, const ap::real_1d_array& x) {
	double candidate_alpha = x(1);

	// equation 3.2
	std::complex<double> evaluation_point{ 0.0, -(candidate_alpha + 1.0) };
	std::complex<double> numerator = HestonHullWhiteCharacteristicFunction(hhw_inputs, evaluation_point);
	double denominator = candidate_alpha * (candidate_alpha + 1.0);

	return std::abs(numerator / denominator);
}


double heston_hull_white::fourier::AlphaOptimizatizer(const HestonHullWhiteInputs& hhw_inputs)
{
	const double epsG = 0.001;
	const double epsF = 0.0;
	const double epsX = 0.0;
	const int MaxIts = 100;

	ap::real_1d_array x;
	x.setbounds(1, 1);
	x(1) = 1.5;
	ap::integer_1d_array nbd;
	nbd.setbounds(1, 1);
	nbd(1) = 1;
	ap::real_1d_array lbound;
	lbound.setbounds(1, 1);
	lbound(1) = 0.01;
	ap::real_1d_array ubound;
	ubound.setbounds(1, 1);
	ubound(1) = 5.0;
	int info;

	heston_hull_white::util::internal::lbfgsbminimize(
		internal::AlphaObjectiveFunction,
		hhw_inputs,
		1,
		1,
		x,
		epsG,
		epsF,
		epsX,
		MaxIts,
		nbd,
		lbound,
		ubound,
		info);

	return x(1);
}