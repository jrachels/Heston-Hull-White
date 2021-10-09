
#include "HestonHullWhiteFFT/CharacteristicFunction/characteristicfunction.h"

std::complex<double> heston_hull_white::fourier::internal::HestonCharacteristicFunction(const HestonHullWhiteInputs& hhw_inputs, std::complex<double> evaluation_point) {
    //std::complex<double> U1 = evaluation_point * (evaluation_point + std::complex_literals::i );

    double gamma = std::pow(hhw_inputs.sigma_, 2);

    std::complex<double> alpha = -0.5 * evaluation_point * (evaluation_point + std::complex<double>{0, 1});
    //hhw_inputs.rho_s_v_
    std::complex<double> beta = hhw_inputs.kappa_ - hhw_inputs.rho_s_v_ * hhw_inputs.sigma_ * evaluation_point * std::complex<double>{0, 1};

    std::complex<double> D = std::sqrt(beta * beta - 2.0 * alpha * gamma);

    std::complex<double> G = (beta - D) / (beta + D);

    std::complex<double> B = ((beta - D) / (gamma)) * ((1.0 - std::exp(-D * hhw_inputs.maturity_)) / (1.0 - G * std::exp(-D * hhw_inputs.maturity_)));

    std::complex<double> psi = (G * std::exp(-D * hhw_inputs.maturity_) - 1.0) / (G - 1.0);

    std::complex<double> A = ((hhw_inputs.kappa_ * hhw_inputs.mean_volatility_) / (gamma)) * ((beta - D) * hhw_inputs.maturity_ - 2.0 * std::log(psi));

    return std::exp(A + B * hhw_inputs.starting_volatility_ + std::complex<double>{0, 1} *evaluation_point* (std::log(hhw_inputs.spot_price_ / hhw_inputs.itsScale_)));

}

std::complex<double> heston_hull_white::fourier::internal::HullWhiteCharacteristicFunctionOfIntegral(const HestonHullWhiteInputs& hhw_inputs, std::complex<double> evaluation_point) {

    // first compute variance of integral

    double cached_constant_0 = std::exp(-hhw_inputs.lambda_ * (hhw_inputs.maturity_));
    double cached_constant_1 = 1.0 / (hhw_inputs.lambda_);

    double variance_of_integral = (std::pow(hhw_inputs.eta_ * cached_constant_1, 2)) * (hhw_inputs.maturity_ + (2.0 * cached_constant_1 * cached_constant_0) - (0.5 * cached_constant_1 * cached_constant_0 * cached_constant_0) - (1.5 * cached_constant_1));

    // mean of integral

    //cached_constant_1* (1.0 - cached_constant_0)* (r_s - alpha(0.0)) + std::log(1 / hhw_inputs.interest_curve_->ZeroBondCurve(t)) + 0.5 * (variance_of_integral);

    double mean_of_integral = std::log(1 / hhw_inputs.interest_curve_->ZeroBondCurve(hhw_inputs.maturity_)) + 0.5 * (variance_of_integral);

    return std::exp((mean_of_integral * std::complex<double>{0, 1} - (.5) * variance_of_integral * evaluation_point)* evaluation_point);

    //return exp(getMeanOfIntegral(itsRInst, 0.0, T) * u * IU - 0.5 * getVarianceOfIntegral(0.0, T) * u * u);

}

std::complex<double> heston_hull_white::fourier::HestonHullWhiteCharacteristicFunction(const HestonHullWhiteInputs& hhw_inputs, std::complex<double> evaluation_point) {
    return internal::HullWhiteCharacteristicFunctionOfIntegral(hhw_inputs, evaluation_point + std::complex<double>{0, 1})* internal::HestonCharacteristicFunction(hhw_inputs, evaluation_point);
}