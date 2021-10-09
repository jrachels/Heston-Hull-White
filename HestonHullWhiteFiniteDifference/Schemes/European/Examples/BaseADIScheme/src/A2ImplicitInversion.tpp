#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A2IMPLICITINVERSION_TPP
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A2IMPLICITINVERSION_TPP

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A2ImplicitInversion.h"

template<typename GridType>
heston_hull_white::finite_difference::examples::A2ImplicitInversion<GridType>::A2ImplicitInversion(const A2<GridType>& a_2, double theta, double epsilon_t) : mean_vol_index_(a_2.MeanVolIndex()), LU_decomposition_(a_2.GetMatrix()) {

	num_vs_ = LU_decomposition_.NumV();
	num_rs_ = LU_decomposition_.NumR();

	// copy a2 matrix
	// make modifications to a2
	const double multiplicative_constant = -theta * epsilon_t;

	//LU_decomposition_.MultiplyByConstant(multiplicative_constant);

	//add identity
	std::size_t grid_size = LU_decomposition_.Size();
	for (std::size_t i = 0; i < grid_size; ++i) {
		auto [entry_1, entry_2, entry_3, entry_4, entry_5] = LU_decomposition_.EntryByVectorIndex(i);
		LU_decomposition_.RefToEntryByVectorIndex(i) = { multiplicative_constant * entry_1, multiplicative_constant * entry_2, multiplicative_constant * entry_3 + 1, \
			multiplicative_constant * entry_4, multiplicative_constant * entry_5 };
	}

	// ad hoc algorithm that is much like thomas's algorithm in place
	// Stores the Ws in the grid

	for (int k = 0; k < num_rs_; ++k) {
		// row cancellation at v = 0, 1
		std::array<double, 5> current_entries = LU_decomposition_.EntryByIndex(0, 0, k);
		std::array<double, 5> next_entries = LU_decomposition_.EntryByIndex(0, 1, k);

		double W = next_entries[1] / current_entries[2];
		next_entries[1] = W;
		next_entries[2] = next_entries[2] - W * current_entries[3];
		next_entries[3] = next_entries[3] - W * current_entries[4];
		LU_decomposition_.EntryByIndex(0, 1, k) = next_entries;

		current_entries = next_entries;

		// row cancellations at v = i, i+1 for 0<i<mean_vol_index_-2
		// at end, current_entries = LU_decomposition_.EntryByIndex(0, mean_vol_index_-2, k);
		std::size_t mean_vol_index_minus_1 = mean_vol_index_;
		for (std::size_t j = 2; j < mean_vol_index_minus_1; ++j) {
			next_entries = LU_decomposition_.EntryByIndex(0, j, k);
			double W = next_entries[1] / current_entries[2];
			next_entries[1] = W;
			next_entries[2] = next_entries[2] - W * current_entries[3];
			LU_decomposition_.EntryByIndex(0, j, k) = next_entries;
			current_entries = next_entries; // we do this 1 more time than we need to
		}

		// row cancellations at v = i, i+1, i+2 for mean_vol_index_-3 <i < num_vs_-2
		// at end, current_entries = LU_decomposition_.EntryByIndex(0, mean_vol_index_-2, k);
		std::size_t num_vs_minus_1 = num_vs_ - 1;
		for (std::size_t j = mean_vol_index_minus_1; j < num_vs_minus_1; ++j) {
			next_entries = LU_decomposition_.EntryByIndex(0, j, k);
			std::array<double, 5> next_next_entries = LU_decomposition_.EntryByIndex(0, j + 1, k);

			double W_1 = next_entries[1] / current_entries[2];
			next_entries[1] = W_1;
			next_entries[2] = next_entries[2] - W_1 * current_entries[3];
			double W_2 = next_next_entries[0] / current_entries[2];
			next_next_entries[0] = W_2;
			next_next_entries[1] = next_next_entries[1] - W_2 * current_entries[3];


			LU_decomposition_.EntryByIndex(0, j, k) = next_entries;
			current_entries = next_entries; // we do this 1 more time than we need to
			next_entries = next_next_entries;
		}

		// row cancellations at v = i, i+1 for i = num_vs_-2

		W = next_entries[1] / current_entries[2];
		next_entries[1] = W;
		next_entries[2] = next_entries[2] - W * current_entries[3];
		LU_decomposition_.EntryByIndex(0, num_vs_minus_1, k) = next_entries;

	}

}

template<typename GridType>
GridType heston_hull_white::finite_difference::examples::A2ImplicitInversion<GridType>::SolveSystem(const GridType& vec) {
	std::size_t num_ss = vec.NumS();

	// construct result
	GridType result{ num_ss, num_vs_, num_rs_ };

	// use the Ws
	for (int i = 0; i < num_ss; ++i) {
		for (int k = 0; k < num_rs_; ++k) {
			double previous = vec.EntryByIndex(i, 0, k);
			// vol  < mean_vol_
			for (int j = 1; j < mean_vol_index_; ++j) {
				previous = vec.EntryByIndex(i, j, k) - (LU_decomposition_.EntryByIndex(0, j, k)[1] * previous);
				result.EntryByIndex(i, j, k) = previous;


				/*double current = vec.EntryByIndex(i, j, k);
				result.EntryByIndex(i, j, k) = current - LU_decomposition_.EntryByIndex(0, j, k)[1] * previous;
				previous = current;*/
			}
			// vol >= mean_vol_
			double previous_previous = vec.EntryByIndex(i, mean_vol_index_ - 2, k);
			for (int j = mean_vol_index_; j < num_vs_; ++j) {
				previous = vec.EntryByIndex(i, j, k) - (LU_decomposition_.EntryByIndex(0, j, k)[0] * previous_previous + LU_decomposition_.EntryByIndex(0, j, k)[1] * previous);
				result.EntryByIndex(i, j, k) = previous;

				/*double current = vec.EntryByIndex(i, j, k);
				result.EntryByIndex(i, j, k) = current - (LU_decomposition_.EntryByIndex(0, j, k)[0] * previous_previous + LU_decomposition_.EntryByIndex(0, j, k)[1] * previous);
				previous = current;*/
			}

		}
	}

	// Second, reverse substitution

	for (int i = 0; i < num_ss; ++i) {
		for (int k = 0; k < num_rs_; ++k) {

			// v = v_max
			double previous_result = result.EntryByIndex(i, num_vs_ - 1, k) / LU_decomposition_.EntryByIndex(0, num_vs_ - 1, k)[2];
			result.EntryByIndex(i, num_vs_ - 1, k) = previous_result;

			// 0 < v < v_max
			for (int j = num_vs_ - 2; j > 0; --j) {
				auto [unused_entry_1, unused_entry_2, u_const_1, u_const_2, unused_entry_3] = LU_decomposition_.EntryByIndex(0, j, k);
				previous_result = (result.EntryByIndex(i, j, k) - u_const_2 * previous_result) / u_const_1;
				result.EntryByIndex(i, j, k) = previous_result;

			}

			// v = 0
			double previous_previous_result = result.EntryByIndex(i, 2, k);
			auto [unused_entry_1, unused_entry_2, u_const_1, u_const_2, u_const_3] = LU_decomposition_.EntryByIndex(0, 0, k);
			double value_at_v_0 = (result.EntryByIndex(i, 0, k) - (u_const_2 * previous_result + u_const_3 * previous_previous_result)) / u_const_1;
			result.EntryByIndex(i, 0, k) = value_at_v_0;

		}
	}

	// return result
	return result;
}

template<typename GridType>
void heston_hull_white::finite_difference::examples::A2ImplicitInversion<GridType>::SolveSystemInPlace(GridType& vec) {
	std::size_t num_ss = vec.NumS();

	// use the Ws
	for (int i = 0; i < num_ss; ++i) {
		for (int k = 0; k < num_rs_; ++k) {
			double previous = vec.EntryByIndex(i, 0, k);
			// vol  < mean_vol_
			for (int j = 1; j < mean_vol_index_; ++j) {
				previous = vec.EntryByIndex(i, j, k) - (LU_decomposition_.EntryByIndex(0, j, k)[1] * previous);
				vec.EntryByIndex(i, j, k) = previous;



				/*double current = vec.EntryByIndex(i, j, k);
				vec.EntryByIndex(i, j, k) = current - LU_decomposition_.EntryByIndex(0, j, k)[1] * previous;
				previous = current;*/
			}
			// vol >= mean_vol_
			double previous_previous = vec.EntryByIndex(i, mean_vol_index_ - 2, k);
			for (std::size_t j = mean_vol_index_; j < num_vs_; ++j) {
				previous = vec.EntryByIndex(i, j, k) - (LU_decomposition_.EntryByIndex(0, j, k)[0] * previous_previous + LU_decomposition_.EntryByIndex(0, j, k)[1] * previous);
				vec.EntryByIndex(i, j, k) = previous;

				/*double current = vec.EntryByIndex(i, j, k);
				vec.EntryByIndex(i, j, k) = current - (LU_decomposition_.EntryByIndex(0, j, k)[0] * previous_previous + LU_decomposition_.EntryByIndex(0, j, k)[1] * previous);
				previous = current;*/
			}

		}
	}

	// Second, reverse substitution

	for (int i = 0; i < num_ss; ++i) {
		for (int k = 0; k < num_rs_; ++k) {

			// v = v_max
			double previous_result = vec.EntryByIndex(i, num_vs_ - 1, k) / LU_decomposition_.EntryByIndex(0, num_vs_ - 1, k)[2];
			vec.EntryByIndex(i, num_vs_ - 1, k) = previous_result;

			// 0 < v < v_max
			for (std::size_t j = num_vs_ - 2; j > 0; --j) {
				auto [unused_entry_1, unused_entry_2, u_const_1, u_const_2, unused_entry_3] = LU_decomposition_.EntryByIndex(0, j, k);
				previous_result = (vec.EntryByIndex(i, j, k) - u_const_2 * previous_result) / u_const_1;
				vec.EntryByIndex(i, j, k) = previous_result;

			}

			// v = 0
			double previous_previous_result = vec.EntryByIndex(i, 2, k);
			auto [unused_entry_1, unused_entry_2, u_const_1, u_const_2, u_const_3] = LU_decomposition_.EntryByIndex(0, 0, k);
			double value_at_v_0 = (vec.EntryByIndex(i, 0, k) - (u_const_2 * previous_result + u_const_3 * previous_previous_result)) / u_const_1;
			vec.EntryByIndex(i, 0, k) = value_at_v_0;

		}
	}

}

#endif