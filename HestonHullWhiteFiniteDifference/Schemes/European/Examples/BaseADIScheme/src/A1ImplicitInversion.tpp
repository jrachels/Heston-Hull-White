#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A1IMPLICITINVERSION_TPP
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A1IMPLICITINVERSION_TPP

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A1ImplicitInversion.h"


template<typename GridType>
heston_hull_white::finite_difference::examples::A1ImplicitInversion<GridType>::A1ImplicitInversion(const A1<GridType>& a_1, double theta, double epsilon_t) : LU_decomposition_(a_1.GetMatrix()) {
	// copy a1 tridiagonal matrix
	// make modifications to a1
	double multiplicative_constant = -theta * epsilon_t;

	//LU_decomposition_.MultiplyByConstant(multiplicative_constant);

	//add identity
	std::size_t grid_size = LU_decomposition_.Size();
	for (std::size_t i = 0; i < grid_size; ++i) {
		auto [entry_1, entry_2, entry_3] = LU_decomposition_.EntryByVectorIndex(i);
		LU_decomposition_.RefToEntryByVectorIndex(i) = { multiplicative_constant * entry_1, multiplicative_constant * entry_2 + 1, multiplicative_constant * entry_3 };
	}

	// thomas algorithm in place
	// taken from https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	// Stores the Ws in the grid

	std::size_t num_ss = LU_decomposition_.NumS();
	std::size_t num_vs = LU_decomposition_.NumV();
	std::size_t num_rs = LU_decomposition_.NumR();

	for (std::size_t j = 0; j < num_vs; ++j) {
		for (std::size_t k = 0; k < num_rs; ++k) {
			std::array<double, 3> current_entries = LU_decomposition_.EntryByIndex(0, j, k);
			for (std::size_t i = 1; i < num_ss; ++i) {
				std::array<double, 3> next_entries = LU_decomposition_.EntryByIndex(i, j, k);
				double W = next_entries[0] / current_entries[1];
				next_entries[1] = next_entries[1] - W * current_entries[2];
				next_entries[0] = W;
				LU_decomposition_.EntryByIndex(i, j, k) = next_entries;
				current_entries = next_entries; // we do this 1 more time than we need to
			}
		}
	}

}

// continuation of thomas's algorithm
template<typename GridType>
GridType heston_hull_white::finite_difference::examples::A1ImplicitInversion<GridType>::SolveSystem(const GridType& vec) {
	std::size_t num_ss = LU_decomposition_.NumS();
	std::size_t num_vs = LU_decomposition_.NumV();
	std::size_t num_rs = LU_decomposition_.NumR();
	// construct result
	GridType result{ num_ss, num_vs, num_rs };

	// use the Ws
	for (std::size_t j = 0; j < num_vs; ++j) {
		for (std::size_t k = 0; k < num_rs; ++k) {
			double previous = vec.EntryByIndex(0, j, k);
			for (std::size_t i = 1; i < num_ss; ++i) {
				previous = vec.EntryByIndex(i, j, k) - (LU_decomposition_.EntryByIndex(i, j, k)[0] * previous);
				result = previous;
			}
		}
	}

	// Second, reverse substitution
	for (std::size_t j = 0; j < num_vs; ++j) {
		for (std::size_t k = 0; k < num_rs; ++k) {
			double previous_result = result.EntryByIndex(num_ss - 1, j, k) / LU_decomposition_.EntryByIndex(num_ss - 1, j, k)[1];
			result.EntryByIndex(num_ss - 1, j, k) = previous_result;
			for (std::size_t i = num_ss - 2; i > -1; --i) {
				auto [l_const, u_const_1, u_const_2] = LU_decomposition_.EntryByIndex(i, j, k);
				previous_result = (result.EntryByIndex(i, j, k) - u_const_2 * previous_result) / u_const_1;
				result.EntryByIndex(i, j, k) = previous_result;
			}
		}
	}

	// return result
	return result;
}

// continuation of thomas's algorithm
template<typename GridType>
void heston_hull_white::finite_difference::examples::A1ImplicitInversion<GridType>::SolveSystemInPlace(GridType& vec) {
	std::size_t num_ss = LU_decomposition_.NumS();
	std::size_t num_vs = LU_decomposition_.NumV();
	std::size_t num_rs = LU_decomposition_.NumR();

	// use the Ws
	for (std::size_t j = 0; j < num_vs; ++j) {
		for (std::size_t k = 0; k < num_rs; ++k) {
			double previous = vec.EntryByIndex(0, j, k);
			for (std::size_t i = 1; i < num_ss; ++i) {
				previous = vec.EntryByIndex(i, j, k) - (LU_decomposition_.EntryByIndex(i, j, k)[0] * previous);
				vec.EntryByIndex(i, j, k) = previous;
			}
		}
	}

	// Second, reverse substitution
	for (std::size_t j = 0; j < num_vs; ++j) {
		for (std::size_t k = 0; k < num_rs; ++k) {
			double previous_result = vec.EntryByIndex(num_ss - 1, j, k) / LU_decomposition_.EntryByIndex(num_ss - 1, j, k)[1];
			vec.EntryByIndex(num_ss - 1, j, k) = previous_result;
			// had to shift i since it is a size_t. maybe everything should have just been ints
			for (std::size_t i = num_ss - 1; i > 0; --i) {
				auto [l_const, u_const_1, u_const_2] = LU_decomposition_.EntryByIndex(i-1, j, k);
				previous_result = (vec.EntryByIndex(i-1, j, k) - u_const_2 * previous_result) / u_const_1;
				vec.EntryByIndex(i-1, j, k) = previous_result;
			}
		}
	}

}

#endif