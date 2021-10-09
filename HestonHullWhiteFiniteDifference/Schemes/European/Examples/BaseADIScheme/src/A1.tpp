#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A1_TPP
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A1_TPP

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A1.h"

template<typename GridType>
heston_hull_white::finite_difference::examples::A1<GridType>::A1(std::shared_ptr<const ADIDiscretization> discretization) : num_ss_(discretization->NumS()), num_vs_(discretization->NumV()), num_rs_(discretization->NumR()), tridiagonal_matrix_(num_ss_, num_vs_, num_rs_), g_1_constants_(1, num_vs_, num_rs_) {
	//int num_rv_pairs = discretization->NumRVPairs();
	//int num_rv_pairs = (discretization->NumV())* (discretization->NumR());
	//int num_ss = (discretization->NumS());
	//int num_triples = num_rv_pairs * num_ss;
	//tridiagonal_matrix.reserve(3 * num_triples);
	std::size_t num_ss_minus_1 = num_ss_ - 1;

	// s = 0. This step is skipped because it is a boundary condition



	// the remainder of the ss

	for (std::size_t i = 0; i < num_ss_minus_1; ++i) {
		auto [current_s, beta_minus_1, beta_0, beta_plus_1, delta_minus_1, delta_0, delta_plus_1] = discretization->SCoefficientsByIndex(i);

		double cached_coefficient_1 = current_s * beta_minus_1;
		double cached_coefficient_2 = (current_s * beta_0) - (1.0 / 3.0); // including -(1.0/3.0) for ru term
		double cached_coefficient_3 = current_s * beta_plus_1;

		double cached_constant = (.5) * (current_s * current_s);

		for (std::size_t j = 0; j < num_vs_; ++j) {
			double current_v = discretization->VCoefficientsByIndex(j)[0];

			double cached_constant_0 = cached_constant * current_v;
			double cached_constant_1 = cached_constant_0 * delta_minus_1;
			double cached_constant_2 = cached_constant_0 * delta_0;
			double cached_constant_3 = cached_constant_0 * delta_plus_1;

			for (std::size_t k = 0; k < num_rs_; ++k) {
				double current_r = discretization->RCoefficientsByIndex(k)[0];
				tridiagonal_matrix_.EntryByIndex(i, j, k) = { current_r * cached_coefficient_1 + cached_constant_1, \
					current_r * cached_coefficient_2 + cached_constant_2, current_r * cached_coefficient_3 + cached_constant_3 };

			}
		}
	}

	// s = s_max
	auto [current_s, beta_minus_1, beta_0, beta_plus_1, delta_minus_1, delta_0, delta_plus_1] = discretization->SCoefficientsByIndex(num_ss_minus_1);
	double cached_constant_0 = (.5) * (current_s * current_s);
	for (std::size_t j = 0; j < num_vs_; ++j) {
		double current_v = discretization->VCoefficientsByIndex(j)[0];
		double cached_constant_1 = cached_constant_0 * current_v;
		for (std::size_t k = 0; k < num_rs_; ++k) {
			double current_r = discretization->RCoefficientsByIndex(k)[0];
			tridiagonal_matrix_.EntryByIndex(num_ss_minus_1, j, k) = { cached_constant_1 * (delta_minus_1 + delta_plus_1), cached_constant_1 * delta_0 - (current_r * (1.0 / 3.0)), 0 };
		}
	}

	double previous_s = discretization->SCoefficientsByIndex(num_ss_minus_1 - 1)[0];

	double cached_constant_2 = cached_constant_0 * (2.0 * delta_plus_1 * (current_s - previous_s));

	// store G_0 coefficient
	for (std::size_t j = 0; j < num_vs_; ++j) {
		double current_v = discretization->VCoefficientsByIndex(j)[0];
		double cached_constant_3 = cached_constant_2 * current_v;
		for (std::size_t k = 0; k < num_rs_; ++k) {
			double current_r = discretization->RCoefficientsByIndex(k)[0];
			g_1_constants_.EntryByIndex(0, j, k) = cached_constant_3 + current_r * current_s;
		}
	}

}



template<typename GridType>
void heston_hull_white::finite_difference::examples::A1<GridType>::ApplyG1InPlace(GridType& grid) {
	// needs access to discretization?
	std::size_t max_s_index = num_ss_ - 1;
	for (std::size_t j = 0; j < num_vs_; ++j) {
		for (std::size_t k = 0; k < num_rs_; ++k) {
			grid.EntryByIndex(max_s_index, j, k) += g_1_constants_.EntryByIndex(0, j, k);
		}
	}
}


template<typename GridType>
GridType heston_hull_white::finite_difference::examples::A1<GridType>::ApplyG1(const GridType& grid) {
	// needs access to discretization?
	GridType new_grid = grid;
	std::size_t max_s_index = num_ss_ - 1;
	for (std::size_t j = 0; j < num_vs_; ++j) {
		for (std::size_t k = 0; k < num_rs_; ++k) {
			new_grid.EntryByIndex(max_s_index, j, k) += g_1_constants_.EntryByIndex(0, j, k);
		}
	}
	return new_grid;
}


template<typename GridType>
GridType heston_hull_white::finite_difference::examples::A1<GridType>::TransformGrid(const GridType& vec) {
	// should assert that vec size is same as assumed A1 size;
	//const auto (num_ss, num_vs, num_rs) = vec.size();
	GridType grid(num_ss_, num_vs_, num_rs_);
	std::size_t num_ss_minus_1 = num_ss_ - 1;
	for (std::size_t j = 0; j < num_vs_; ++j) {
		for (std::size_t k = 0; k < num_rs_; ++k) {
			// s < s_max
			double last_s = 0;
			double current_s = vec.EntryByIndex(0, j, k);
			double next_s = vec.EntryByIndex(1, j, k);
			auto [coefficient_1, coefficient_2, coefficient_3] = tridiagonal_matrix_.EntryByIndex(0, j, k);
			grid.EntryByIndex(0, j, k) = current_s * coefficient_2 + next_s * coefficient_3;

			for (std::size_t i = 1; i < num_ss_minus_1; ++i) {
				last_s = current_s;
				current_s = next_s;
				next_s = vec.EntryByIndex(i + 1, j, k);
				auto [coefficient_1, coefficient_2, coefficient_3] = tridiagonal_matrix_.EntryByIndex(i, j, k);
				grid.EntryByIndex(i, j, k) = last_s * coefficient_1 + current_s * coefficient_2 + next_s * coefficient_3;
			}

			// s_max
			auto [coefficient_4, coefficient_5, coefficient_6] = tridiagonal_matrix_.EntryByIndex(num_ss_minus_1, j, k);
			grid.EntryByIndex(num_ss_minus_1, j, k) = current_s * coefficient_4 + next_s * coefficient_5;
		}
	}
	return grid;
}


template<typename GridType>
void heston_hull_white::finite_difference::examples::A1<GridType>::TransformGridInPlace(GridType& vec) {
	std::size_t num_ss_minus_1 = num_ss_ - 1;
	for (std::size_t j = 0; j < num_vs_; ++j) {
		for (std::size_t k = 0; k < num_rs_; ++k) {
			// s < s_max
			double last_s = 0;
			double current_s = vec.EntryByIndex(0, j, k);
			double next_s = vec.EntryByIndex(1, j, k);
			auto [coefficient_1, coefficient_2, coefficient_3] = tridiagonal_matrix_.EntryByIndex(0, j, k);
			vec.EntryByIndex(0, j, k) = current_s * coefficient_2 + next_s * coefficient_3;

			for (std::size_t i = 1; i < num_ss_minus_1; ++i) {
				last_s = current_s;
				current_s = next_s;
				next_s = vec.EntryByIndex(i + 1, j, k);
				auto [coefficient_1, coefficient_2, coefficient_3] = tridiagonal_matrix_.EntryByIndex(i, j, k);
				vec.EntryByIndex(i, j, k) = last_s * coefficient_1 + current_s * coefficient_2 + next_s * coefficient_3;
			}

			// s_max
			auto [coefficient_4, coefficient_5, coefficient_6] = tridiagonal_matrix_.EntryByIndex(num_ss_minus_1, j, k);
			vec.EntryByIndex(num_ss_minus_1, j, k) = current_s * coefficient_4 + next_s * coefficient_5;
		}
	}
}

#endif