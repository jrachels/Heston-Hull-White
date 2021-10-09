#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A3_TPP
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A3_TPP

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A3.h"

template<typename GridType>
heston_hull_white::finite_difference::examples::A3<GridType>::A3(const double lambda, const double eta, const double theta_t, std::shared_ptr<const ADIDiscretization> discretization) : num_ss_(discretization->NumS()), num_vs_(discretization->NumV()), num_rs_(discretization->NumR()), condensed_a3_matrix_(num_rs_) {

	const double cached_second_partial_coefficient = std::pow(eta, 2) / 2.0;

	std::size_t num_rs_minus_1 = num_rs_ - 1;

	// r = -R_max

	auto [r_min, unused_zero_1, unused_zero_2, unused_zero_3, r_min_delta_minus_1, r_min_delta_0, r_min_delta_plus_1] = discretization->RCoefficientsByIndex(0);


	//condensed_a3_matrix_[0] = { 0 , cached_second_partial_coefficient * r_min_delta_0 - (r_min / 3.0), cached_second_partial_coefficient * (r_min_delta_minus_1 + r_min_delta_plus_1) };

	condensed_a3_matrix_[0] = { 0 , cached_second_partial_coefficient * (r_min_delta_0)-(r_min / 3.0), cached_second_partial_coefficient * (r_min_delta_minus_1 + r_min_delta_plus_1) };

	// r = R_max

	auto [r_max, unused_zero_4, unused_zero_5, unused_zero_6, r_max_delta_minus_1, r_max_delta_0, r_max_delta_plus_1] = discretization->RCoefficientsByIndex(num_rs_minus_1);

	condensed_a3_matrix_[num_rs_minus_1] = { cached_second_partial_coefficient * (r_max_delta_minus_1 + r_max_delta_plus_1) , cached_second_partial_coefficient * (r_max_delta_0)-(r_max / 3.0), 0 };


	// all other r

	const double cached_first_partial_constant = lambda * theta_t;

	for (std::size_t k = 1; k < num_rs_minus_1; ++k) {
		auto [current_r, beta_minus_1, beta_0, beta_plus_1, delta_minus_1, delta_0, delta_plus_1] = discretization->RCoefficientsByIndex(k);

		const double cached_first_partial_coefficient = cached_first_partial_constant - (lambda * current_r);


		condensed_a3_matrix_[k] = { (cached_second_partial_coefficient * delta_minus_1) + (cached_first_partial_coefficient * beta_minus_1), \
			(cached_second_partial_coefficient * delta_0) + (cached_first_partial_coefficient * beta_0) - (current_r / 3.0), \
			(cached_second_partial_coefficient * delta_plus_1) + (cached_first_partial_coefficient * beta_plus_1) };
	}

}

template<typename GridType>
GridType heston_hull_white::finite_difference::examples::A3<GridType>::TransformGrid(const GridType& vec) {
	// should assert that vec size is same as assumed A3 size;
	//const auto (num_ss, num_vs, num_rs) = vec.size();
	GridType grid(num_ss_, num_vs_, num_rs_);
	size_t num_rs_minus_1 = num_rs_ - 1;
	for (std::size_t i = 0; i < num_ss_; ++i) {
		for (std::size_t j = 0; j < num_vs_; ++j) {

			// r = R_min

			double current_r = vec.EntryByIndex(i, j, 0);
			double next_r = vec.EntryByIndex(i, j, 1);

			auto [unused_zero_1, coefficient_2, coefficient_3] = condensed_a3_matrix_[0];
			grid.EntryByIndex(i, j, 0) = current_r * coefficient_2 + next_r * coefficient_3;

			// R_min < r < R_max

			for (std::size_t k = 1; k < num_rs_minus_1; ++k) {
				double last_r = current_r;
				current_r = next_r;
				next_r = vec.EntryByIndex(i, j, k + 1);
				auto [coefficient_1, coefficient_2, coefficient_3] = condensed_a3_matrix_[k];
				grid.EntryByIndex(i, j, k) = last_r * coefficient_1 + current_r * coefficient_2 + next_r * coefficient_3;
			}

			// r = R_max
			auto [coefficient_4, coefficient_5, unused_zero_2] = condensed_a3_matrix_[num_rs_minus_1];
			grid.EntryByIndex(i, j, num_rs_minus_1) = current_r * coefficient_4 + next_r * coefficient_5;

		}
	}
	return grid;
}

template<typename GridType>
void heston_hull_white::finite_difference::examples::A3<GridType>::TransformGridInPlace(GridType& vec) {
	size_t num_rs_minus_1 = num_rs_ - 1;
	for (std::size_t i = 0; i < num_ss_; ++i) {
		for (std::size_t j = 0; j < num_vs_; ++j) {

			// r = R_min

			double current_r = vec.EntryByIndex(i, j, 0);
			double next_r = vec.EntryByIndex(i, j, 1);

			auto [unused_zero_1, coefficient_2, coefficient_3] = condensed_a3_matrix_[0];
			vec.EntryByIndex(i, j, 0) = current_r * coefficient_2 + next_r * coefficient_3;

			// R_min < r < R_max

			for (std::size_t k = 1; k < num_rs_minus_1; ++k) {
				double last_r = current_r;
				current_r = next_r;
				next_r = vec.EntryByIndex(i, j, k + 1);
				auto [coefficient_1, coefficient_2, coefficient_3] = condensed_a3_matrix_[k];
				vec.EntryByIndex(i, j, k) = last_r * coefficient_1 + current_r * coefficient_2 + next_r * coefficient_3;
			}

			// r = R_max
			auto [coefficient_4, coefficient_5, unused_zero_2] = condensed_a3_matrix_[num_rs_minus_1];
			vec.EntryByIndex(i, j, num_rs_minus_1) = current_r * coefficient_4 + next_r * coefficient_5;

		}
	}
}

#endif