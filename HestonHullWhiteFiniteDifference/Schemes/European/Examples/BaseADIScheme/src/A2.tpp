#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A2_TPP
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A2_TPP

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A2.h"

template<typename GridType>
heston_hull_white::finite_difference::examples::A2<GridType>::A2(const double mean_volatility, const double kappa, const double sigma /*vol of vol*/, std::size_t mean_vol_index /*included so I don't have to compute it twice*/, std::shared_ptr<const ADIDiscretization> discretization) : num_ss_(discretization->NumS()), num_vs_(discretization->NumV()), num_rs_(discretization->NumR()), mean_vol_index_(mean_vol_index), condensed_a2_matrix_(1, num_vs_, num_rs_) /*the S is irrelevant.*/, g_2_constants_(num_ss_) {

	assert(mean_vol_index > 1);

	const double cached_second_partial_constant = std::pow(sigma, 2) / 2.0;

	const double cached_first_partial_constant = mean_volatility * kappa;

	std::size_t num_vs_minus_1 = num_vs_ - 1;

	// v = 0

	auto [v_0, v_0_gamma_minus_1, v_0_gamma_0, v_0_gamma_plus_1, unused_zero_1, unused_zero_2, unused_zero_3, unused_zero_4, unused_zero_5, unused_zero_6] = discretization->VCoefficientsByIndex(0);


	for (std::size_t k = 0; k < num_rs_; ++k) {
		double current_r = discretization->RCoefficientsByIndex(k)[0];
		condensed_a2_matrix_.EntryByIndex(0, 0, k) = { 0, 0, (cached_first_partial_constant * v_0_gamma_minus_1) - (current_r / 3.0), (cached_first_partial_constant * v_0_gamma_0), (cached_first_partial_constant * v_0_gamma_plus_1) };
	}

	// 0 < v < mean_vol

	for (std::size_t j = 1; j < mean_vol_index_; ++j) {
		auto [current_v, unused_entry_1, unused_entry_2, unused_entry_3, beta_minus_1, beta_0, beta_plus_1, delta_minus_1, delta_0, delta_plus_1] = discretization->VCoefficientsByIndex(j);
		const double cached_first_partial_constant_coefficient = cached_first_partial_constant - (current_v * kappa);
		const double cached_second_partial_constant_coefficient = cached_second_partial_constant * current_v;
		for (std::size_t k = 0; k < num_rs_; ++k) {
			double current_r = discretization->RCoefficientsByIndex(k)[0];
			condensed_a2_matrix_.EntryByIndex(0, j, k) = { 0, (beta_minus_1 * cached_first_partial_constant_coefficient) + (delta_minus_1 * cached_second_partial_constant_coefficient), \
				(beta_0 * cached_first_partial_constant_coefficient) + (delta_0 * cached_second_partial_constant_coefficient) - (current_r / 3.0),
				(beta_plus_1 * cached_first_partial_constant_coefficient) + (delta_plus_1 * cached_second_partial_constant_coefficient), 0 };
		}

	}

	// mean_vol < v < v_{max-1}

	for (std::size_t j = mean_vol_index_; j < num_vs_minus_1; ++j) {
		auto [current_v, alpha_minus_2, alpha_minus_1, alpha_0, unused_entry_1, unused_entry_2, unused_entry_3, delta_minus_1, delta_0, delta_plus_1] = discretization->VCoefficientsByIndex(j);
		const double cached_first_partial_constant_coefficient = cached_first_partial_constant - (current_v * kappa);
		const double cached_second_partial_constant_coefficient = cached_second_partial_constant * current_v;
		for (std::size_t k = 0; k < num_rs_; ++k) {
			double current_r = discretization->RCoefficientsByIndex(k)[0];
			condensed_a2_matrix_.EntryByIndex(0, j, k) = { (alpha_minus_2 * cached_first_partial_constant_coefficient), (alpha_minus_1 * cached_first_partial_constant_coefficient) + (delta_minus_1 * cached_second_partial_constant_coefficient), \
				(alpha_0 * cached_first_partial_constant_coefficient) + (delta_0 * cached_second_partial_constant_coefficient) - (current_r / 3.0),
				(delta_plus_1 * cached_second_partial_constant_coefficient), 0 };
		}
	}

	// v = v_{max-1}
	auto [v_max_minus_1, v_max_alpha_minus_2, v_max_alpha_minus_1, v_max_alpha_0, unused_zero_7, unused_zero_8, unused_zero_9, delta_minus_1, delta_0, delta_plus_1] = discretization->VCoefficientsByIndex(num_vs_minus_1);

	const double cached_first_partial_constant_coefficient_at_v_max = cached_first_partial_constant - (v_max_minus_1 * kappa);
	const double cached_second_partial_constant_coefficient = cached_second_partial_constant * v_max_minus_1;

	for (std::size_t k = 0; k < num_rs_; ++k) {
		double current_r = discretization->RCoefficientsByIndex(k)[0];
		condensed_a2_matrix_.EntryByIndex(0, num_vs_minus_1, k) = { (v_max_alpha_minus_2 * cached_first_partial_constant_coefficient_at_v_max), (v_max_alpha_minus_1 * cached_first_partial_constant_coefficient_at_v_max) + (delta_minus_1 * cached_second_partial_constant_coefficient), \
				(v_max_alpha_0 * cached_first_partial_constant_coefficient_at_v_max) + (delta_0 * cached_second_partial_constant_coefficient) - (current_r / 3.0), 0, 0 };
	}

	// g2 = 0
	const double g_2_constant = cached_second_partial_constant_coefficient * delta_plus_1;
	for (std::size_t i = 0; i < num_ss_; ++i) {
		g_2_constants_[i] = g_2_constant * (discretization->SCoefficientsByIndex(i)[0]);
	}
}

template<typename GridType>
void heston_hull_white::finite_difference::examples::A2<GridType>::ApplyG2InPlace(GridType& grid) {
	std::size_t max_v_index = num_vs_ - 1;
	for (std::size_t i = 0; i < num_ss_; ++i) {
		for (std::size_t k = 0; k < num_rs_; ++k) {
			grid.EntryByIndex(i, max_v_index, k) += g_2_constants_[i];
		}
	}
}

template<typename GridType>
void heston_hull_white::finite_difference::examples::A2<GridType>::ApplyG2(const GridType& grid) {
	GridType new_grid = grid;

	std::size_t max_v_index = num_vs_ - 1;
	for (std::size_t i = 0; i < num_ss_; ++i) {
		for (std::size_t k = 0; k < num_rs_; ++k) {
			new_grid.EntryByIndex(i, max_v_index, k) += g_2_constants_[i];
		}
	}
	return new_grid;
}

template<typename GridType>
GridType heston_hull_white::finite_difference::examples::A2<GridType>::TransformGrid(const GridType& vec) {
	// should assert that vec size is same as assumed A2 size;
	//const auto (num_ss, num_vs, num_rs) = vec.size();
	GridType grid(num_ss_, num_vs_, num_rs_);
	size_t num_vs_minus_1 = num_vs_ - 1;
	for (std::size_t i = 0; i < num_ss_; ++i) {
		for (std::size_t k = 0; k < num_rs_; ++k) {

			// v = 0 and v = 1

			double last_v = vec.EntryByIndex(i, 0, k);
			double current_v = vec.EntryByIndex(i, 1, k);
			double next_v = vec.EntryByIndex(i, 2, k);

			auto [unused_zero_1, unused_zero_2, coefficient_1, coefficient_2, coefficient_3] = condensed_a2_matrix_.EntryByIndex(0, 0, k);
			grid.EntryByIndex(i, 0, k) = last_v * coefficient_1 + current_v * coefficient_2 + next_v * coefficient_3;

			auto [unused_zero_3, coefficient_4, coefficient_5, coefficient_6, unused_zero_4] = condensed_a2_matrix_.EntryByIndex(0, 1, k);
			grid.EntryByIndex(i, 1, k) = last_v * coefficient_4 + current_v * coefficient_5 + next_v * coefficient_6;

			// 0 < v < mean_vol

			for (std::size_t j = 2; j < mean_vol_index_; ++j) {
				last_v = current_v;
				current_v = next_v;
				next_v = vec.EntryByIndex(i, j + 1, k);

				auto [unused_zero_1, coefficient_1, coefficient_2, coefficient_3, unused_zero_2] = condensed_a2_matrix_.EntryByIndex(0, j, k);
				grid.EntryByIndex(i, j, k) = last_v * coefficient_1 + current_v * coefficient_2 + next_v * coefficient_3;
			}

			// v = meal_vol
			double before_last_v = last_v;
			last_v = current_v;
			current_v = next_v;
			next_v = vec.EntryByIndex(i, mean_vol_index_ + 1, k);

			auto [coefficient_7, coefficient_8, coefficient_9, coefficient_10, unused_zero_5] = condensed_a2_matrix_.EntryByIndex(0, mean_vol_index_, k);

			grid.EntryByIndex(i, mean_vol_index_, k) = before_last_v * coefficient_7 + last_v * coefficient_8 + current_v * coefficient_9 + next_v * coefficient_10;

			// mean_vol < v < v_max

			for (std::size_t j = mean_vol_index_ + 1; j < num_vs_minus_1; ++j) {
				before_last_v = last_v;
				last_v = current_v;
				current_v = next_v;
				next_v = vec.EntryByIndex(i, j + 1, k);

				auto [coefficient_1, coefficient_2, coefficient_3, coefficient_4, unused_zero_1] = condensed_a2_matrix_.EntryByIndex(0, j, k);

				grid.EntryByIndex(i, j, k) = before_last_v * coefficient_1 + last_v * coefficient_2 + current_v * coefficient_3 + next_v * coefficient_4;
			}

			// v = v_{max-1}

			auto [coefficient_11, coefficient_12, coefficient_13, unused_zero_6, unused_zero_7] = condensed_a2_matrix_.EntryByIndex(0, num_vs_minus_1, k);

			grid.EntryByIndex(i, num_vs_minus_1, k) = last_v * coefficient_11 + current_v * coefficient_12 + next_v * coefficient_13;

		}
	}
	return grid;
}

template<typename GridType>
void heston_hull_white::finite_difference::examples::A2<GridType>::TransformGridInPlace(GridType& vec) {

	size_t num_vs_minus_1 = num_vs_ - 1;
	for (std::size_t i = 0; i < num_ss_; ++i) {
		for (std::size_t k = 0; k < num_rs_; ++k) {

			// v = 0 and v = 1

			double last_v = vec.EntryByIndex(i, 0, k);
			double current_v = vec.EntryByIndex(i, 1, k);
			double next_v = vec.EntryByIndex(i, 2, k);

			auto [unused_zero_1, unused_zero_2, coefficient_1, coefficient_2, coefficient_3] = condensed_a2_matrix_.EntryByIndex(0, 0, k);
			vec.EntryByIndex(i, 0, k) = last_v * coefficient_1 + current_v * coefficient_2 + next_v * coefficient_3;

			auto [unused_zero_3, coefficient_4, coefficient_5, coefficient_6, unused_zero_4] = condensed_a2_matrix_.EntryByIndex(0, 1, k);
			vec.EntryByIndex(i, 1, k) = last_v * coefficient_4 + current_v * coefficient_5 + next_v * coefficient_6;

			// 0 < v < mean_vol

			for (std::size_t j = 2; j < mean_vol_index_; ++j) {
				last_v = current_v;
				current_v = next_v;
				next_v = vec.EntryByIndex(i, j + 1, k);

				auto [unused_zero_1, coefficient_1, coefficient_2, coefficient_3, unused_zero_2] = condensed_a2_matrix_.EntryByIndex(0, j, k);
				vec.EntryByIndex(i, j, k) = last_v * coefficient_1 + current_v * coefficient_2 + next_v * coefficient_3;
			}

			// v = meal_vol
			double before_last_v = last_v;
			last_v = current_v;
			current_v = next_v;
			next_v = vec.EntryByIndex(i, mean_vol_index_ + 1, k);

			auto [coefficient_7, coefficient_8, coefficient_9, coefficient_10, unused_zero_5] = condensed_a2_matrix_.EntryByIndex(0, mean_vol_index_, k);

			vec.EntryByIndex(i, mean_vol_index_, k) = before_last_v * coefficient_7 + last_v * coefficient_8 + current_v * coefficient_9 + next_v * coefficient_10;

			// mean_vol < v < v_max

			for (std::size_t j = mean_vol_index_ + 1; j < num_vs_minus_1; ++j) {
				before_last_v = last_v;
				last_v = current_v;
				current_v = next_v;
				next_v = vec.EntryByIndex(i, j + 1, k);

				auto [coefficient_1, coefficient_2, coefficient_3, coefficient_4, unused_zero_1] = condensed_a2_matrix_.EntryByIndex(0, j, k);

				vec.EntryByIndex(i, j, k) = before_last_v * coefficient_1 + last_v * coefficient_2 + current_v * coefficient_3 + next_v * coefficient_4;
			}

			// v = v_{max-1}

			auto [coefficient_11, coefficient_12, coefficient_13, unused_zero_6, unused_zero_7] = condensed_a2_matrix_.EntryByIndex(0, num_vs_minus_1, k);

			vec.EntryByIndex(i, num_vs_minus_1, k) = last_v * coefficient_11 + current_v * coefficient_12 + next_v * coefficient_13;

		}
	}
}


#endif