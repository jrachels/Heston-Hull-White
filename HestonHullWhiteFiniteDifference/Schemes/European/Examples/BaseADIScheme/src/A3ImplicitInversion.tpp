#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A3IMPLICITINVERSION_TPP
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A3IMPLICITINVERSION_TPP

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A3ImplicitInversion.h"

template<typename GridType>
heston_hull_white::finite_difference::examples::A3ImplicitInversion<GridType>::A3ImplicitInversion(const A3<GridType>& a_3, double theta, double epsilon_t) : condensed_LU_decomposition_(a_3.GetMatrix()) {
	// copy a3 tridiagonal matrix
	// make modifications to a3
	double multiplicative_constant = -theta * epsilon_t;

	//add identity
	size_t mat_size = condensed_LU_decomposition_.size();
	for (int k = 0; k < mat_size; ++k) {
		auto [entry_1, entry_2, entry_3] = condensed_LU_decomposition_[k];
		condensed_LU_decomposition_[k] = { multiplicative_constant * entry_1, (multiplicative_constant * entry_2) + 1, multiplicative_constant * entry_3 };
	}

	// thomas algorithm in place
	// taken from https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	// Stores the Ws in the grid

	// TODO: this could probably be made a bit faster with reference. see structured binding
	std::array<double, 3> current_entries = condensed_LU_decomposition_[0];
	for (int k = 1; k < mat_size; ++k) {
		// auto& [entry_1, entry_2, entry_3] = condensed_LU_decomposition_[k];
		std::array<double, 3> next_entries = condensed_LU_decomposition_[k];
		double W = next_entries[0] / current_entries[1];
		next_entries[1] = next_entries[1] - W * current_entries[2];
		next_entries[0] = W;
		condensed_LU_decomposition_[k] = next_entries;
		current_entries = next_entries; // we do this 1 more time than we need to
	}


}

template<typename GridType>
GridType heston_hull_white::finite_difference::examples::A3ImplicitInversion<GridType>::SolveSystem(const GridType& vec) {

	size_t mat_size = condensed_LU_decomposition_.size();

	int num_ss = vec.NumS();
	int num_vs = vec.NumV();
	int num_rs = vec.NumR();

	// construct result
	GridType result{ num_ss, num_vs, num_rs };

	// use the Ws
	for (int i = 0; i < num_ss; ++i) {
		for (int j = 0; j < num_vs; ++j) {
			double previous = vec.EntryByIndex(i, j, 0);
			for (int k = 1; k < num_rs; ++k) {
				previous = vec.EntryByIndex(i, j, k) - (condensed_LU_decomposition_[k][0] * previous);
				result = previous;
			}
		}
	}

	size_t num_rs_minus_1 = num_rs - 1;

	// Second, reverse substitution
	for (int i = 0; i < num_ss; ++i) {
		for (int j = 0; j < num_vs; ++j) {
			double previous_result = result.EntryByIndex(i, j, num_rs_minus_1) / condensed_LU_decomposition_[num_rs_minus_1][1];
			result.EntryByIndex(i, j, num_rs_minus_1) = previous_result;
			for (int k = num_rs_minus_1 - 1; k > -1; --k) {
				auto [l_const, u_const_1, u_const_2] = condensed_LU_decomposition_[k];
				previous_result = (result.EntryByIndex(i, j, k) - u_const_2 * previous_result) / u_const_1;
				result.EntryByIndex(i, j, k) = previous_result;
			}
		}
	}

	// return result
	return result;
}

template<typename GridType>
void heston_hull_white::finite_difference::examples::A3ImplicitInversion<GridType>::SolveSystemInPlace(GridType& vec) {

	size_t mat_size = condensed_LU_decomposition_.size();

	std::size_t num_ss = vec.NumS();
	std::size_t num_vs = vec.NumV();
	std::size_t num_rs = vec.NumR();

	// use the Ws
	for (std::size_t i = 0; i < num_ss; ++i) {
		for (std::size_t j = 0; j < num_vs; ++j) {
			double previous = vec.EntryByIndex(i, j, 0);
			for (std::size_t k = 1; k < num_rs; ++k) {
				previous = vec.EntryByIndex(i, j, k) - (condensed_LU_decomposition_[k][0] * previous);
				vec.EntryByIndex(i, j, k) = previous;
			}
		}
	}

	size_t num_rs_minus_1 = num_rs - 1;

	// Second, reverse substitution
	for (std::size_t i = 0; i < num_ss; ++i) {
		for (std::size_t j = 0; j < num_vs; ++j) {
			double previous_result = vec.EntryByIndex(i, j, num_rs_minus_1) / condensed_LU_decomposition_[num_rs_minus_1][1];
			vec.EntryByIndex(i, j, num_rs_minus_1) = previous_result;
			for (std::size_t k = num_rs_minus_1; k > 0; --k) {
				auto [l_const, u_const_1, u_const_2] = condensed_LU_decomposition_[k-1];
				previous_result = (vec.EntryByIndex(i, j, k-1) - (u_const_2 * previous_result)) / u_const_1;
				vec.EntryByIndex(i, j, k-1) = previous_result;
			}
		}
	}

}

#endif