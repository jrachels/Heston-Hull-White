#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A0_TPP
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_SRC_A0_TPP

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A0.h"

template<typename GridType>
heston_hull_white::finite_difference::examples::A0<GridType>::A0(std::shared_ptr<const ADIDiscretization> discretization, std::size_t mean_vol_index /*included so I don't have to compute it twice*/, \
	const double vol_of_vol, const double eta, const double rho_s_v, const double rho_s_r, const double rho_v_r) : \
	num_ss_(discretization->NumS()), num_vs_(discretization->NumV()), num_rs_(discretization->NumR()), mean_vol_index_(mean_vol_index), \
	sparse_A0_(num_ss_, num_vs_, num_rs_) {

	std::size_t num_ss_minus_1 = num_ss_ - 1;
	std::size_t num_rs_minus_1 = num_rs_ - 1;

	// instead of doing anything special at the boundary here, I just store whatever values I store and ignore them when I perform the transformation

	double  s_v_intermediate_partial_coefficient_0 = rho_s_v * vol_of_vol;
	double  s_r_intermediate_partial_coefficient_0 = rho_s_r * eta;
	double  v_r_intermediate_partial_coefficient_0 = rho_v_r * vol_of_vol * eta;

	for (std::size_t j = 0; j < mean_vol_index_; ++j) {
		auto [current_v, unused_entry_1, unused_entry_2, unused_entry_3, v_beta_minus_1, v_beta_0, v_beta_plus_1, unused_entry_4, unused_entry_5, unused_entry_6] = discretization->VCoefficientsByIndex(j);

		const double sqrt_v = std::sqrt(current_v);
		double  s_v_intermediate_partial_coefficient_1 = s_v_intermediate_partial_coefficient_0 * current_v;
		double  s_r_intermediate_partial_coefficient_1 = s_r_intermediate_partial_coefficient_0 * sqrt_v;

		double  v_r_partial_coefficient = v_r_intermediate_partial_coefficient_0 * sqrt_v;

		for (std::size_t i = 0; i < num_ss_minus_1; ++i) {
			auto [current_s, s_beta_minus_1, s_beta_0, s_beta_plus_1, unused_entry_7, unused_entry_8, unused_entry_9] = discretization->SCoefficientsByIndex(i);

			double  s_v_partial_coefficient = s_v_intermediate_partial_coefficient_1 * current_s;
			double  s_r_partial_coefficient = s_r_intermediate_partial_coefficient_1 * current_s;



			// r = R_min

			// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1), (-1, -1, 0), (1, -1, 0), (-1, 1, 0), (1, 1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1), (0, 1, 1)

			sparse_A0_.EntryByIndex(i, j, 0) = { v_beta_0 * s_beta_0 * s_v_partial_coefficient, /*1*/\
					(s_beta_minus_1 * v_beta_0 * s_v_partial_coefficient) /*2*/, \
					(s_beta_plus_1 * v_beta_0 * s_v_partial_coefficient) /*3*/, \
					(s_beta_0 * v_beta_minus_1 * s_v_partial_coefficient)  /*4*/, \
					(s_beta_0 * v_beta_plus_1 * s_v_partial_coefficient)  /*5*/, \
					0 /*6*/, \
					0 /*7*/, \
					(s_beta_minus_1 * v_beta_minus_1 * s_v_partial_coefficient) /*8*/, \
					(s_beta_plus_1 * v_beta_minus_1 * s_v_partial_coefficient) /*9*/, \
					(s_beta_minus_1 * v_beta_plus_1 * s_v_partial_coefficient) /*10*/, \
					(s_beta_plus_1 * v_beta_plus_1 * s_v_partial_coefficient) /*11*/, \
					0 /*12*/, \
					0/*13*/, \
					0/*14*/, \
					0/*15*/, \
					0 /*16*/, \
					0/*17*/, \
					0 /*18*/, \
					0 /*19*/ };

			// R_min < r < R_max


			for (std::size_t k = 1; k < num_rs_minus_1; ++k) {
				auto [current_r, r_beta_minus_1, r_beta_0, r_beta_plus_1, unused_entry_10, unused_entry_11, unused_entry_12] = discretization->RCoefficientsByIndex(k);

				// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1), (-1, -1, 0), (1, -1, 0), (-1, 1, 0), (1, 1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1), (0, 1, 1)
				sparse_A0_.EntryByIndex(i, j, k) = { v_beta_0 * s_beta_0 * s_v_partial_coefficient + \
					(s_beta_0 * r_beta_0 * s_r_partial_coefficient) + v_beta_0 * r_beta_0 * v_r_partial_coefficient, /*1*/\
					(s_beta_minus_1 * v_beta_0 * s_v_partial_coefficient) + (s_beta_minus_1 * r_beta_0 * s_r_partial_coefficient) /*2*/, \
					(s_beta_plus_1 * v_beta_0 * s_v_partial_coefficient) + (s_beta_plus_1 * r_beta_0 * s_r_partial_coefficient) /*3*/, \
					(s_beta_0 * v_beta_minus_1 * s_v_partial_coefficient) + (v_beta_minus_1 * r_beta_0 * v_r_partial_coefficient) /*4*/, \
					(s_beta_0 * v_beta_plus_1 * s_v_partial_coefficient) + (v_beta_plus_1 * r_beta_0 * v_r_partial_coefficient) /*5*/, \
					(s_beta_0 * r_beta_minus_1 * s_r_partial_coefficient) + (v_beta_0 * r_beta_minus_1 * v_r_partial_coefficient) /*6*/, \
					(s_beta_0 * r_beta_plus_1 * s_r_partial_coefficient) + (v_beta_0 * r_beta_plus_1 * v_r_partial_coefficient) /*7*/, \
					(s_beta_minus_1 * v_beta_minus_1 * s_v_partial_coefficient) /*8*/, \
					(s_beta_plus_1 * v_beta_minus_1 * s_v_partial_coefficient) /*9*/, \
					(s_beta_minus_1 * v_beta_plus_1 * s_v_partial_coefficient) /*10*/, \
					(s_beta_plus_1 * v_beta_plus_1 * s_v_partial_coefficient) /*11*/, \
					(s_beta_minus_1 * r_beta_minus_1 * s_r_partial_coefficient) /*12*/, \
					(s_beta_plus_1 * r_beta_minus_1 * s_r_partial_coefficient)/*13*/, \
					(s_beta_minus_1 * r_beta_plus_1 * s_r_partial_coefficient)/*14*/, \
					(s_beta_plus_1 * r_beta_plus_1 * s_r_partial_coefficient)/*15*/, \
					(v_beta_minus_1 * r_beta_minus_1 * v_r_partial_coefficient) /*16*/, \
					(v_beta_plus_1 * r_beta_minus_1 * v_r_partial_coefficient)/*17*/, \
					(v_beta_minus_1 * r_beta_plus_1 * v_r_partial_coefficient) /*18*/, \
					(v_beta_plus_1 * r_beta_plus_1 * v_r_partial_coefficient) /*19*/ };

			}


			// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1), (-1, -1, 0), (1, -1, 0), (-1, 1, 0), (1, 1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1), (0, 1, 1)

			// r = R_max

			sparse_A0_.EntryByIndex(i, j, num_rs_minus_1) = { v_beta_0 * s_beta_0 * s_v_partial_coefficient, /*1*/\
				(s_beta_minus_1 * v_beta_0 * s_v_partial_coefficient) /*2*/, \
				(s_beta_plus_1 * v_beta_0 * s_v_partial_coefficient) /*3*/, \
				(s_beta_0 * v_beta_minus_1 * s_v_partial_coefficient)  /*4*/, \
				(s_beta_0 * v_beta_plus_1 * s_v_partial_coefficient)  /*5*/, \
				0 /*6*/, \
				0 /*7*/, \
				(s_beta_minus_1 * v_beta_minus_1 * s_v_partial_coefficient) /*8*/, \
				(s_beta_plus_1 * v_beta_minus_1 * s_v_partial_coefficient) /*9*/, \
				(s_beta_minus_1 * v_beta_plus_1 * s_v_partial_coefficient) /*10*/, \
				(s_beta_plus_1 * v_beta_plus_1 * s_v_partial_coefficient) /*11*/, \
				0 /*12*/, \
				0/*13*/, \
				0/*14*/, \
				0/*15*/, \
				0 /*16*/, \
				0/*17*/, \
				0 /*18*/, \
				0 /*19*/ };

		}


		// r = R_min who cares
		sparse_A0_.EntryByIndex(num_ss_minus_1, j, 0) = { 0, /*1*/\
				0 /*2*/, \
				0/*3*/, \
				0 /*4*/, \
				0 /*5*/, \
				0 /*6*/, \
				0 /*7*/, \
				0 /*8*/, \
				0 /*9*/, \
				0 /*10*/, \
				0 /*11*/, \
				0 /*12*/, \
				0/*13*/, \
				0/*14*/, \
				0/*15*/, \
				0 /*16*/, \
				0/*17*/, \
				0 /*18*/, \
				0 /*19*/ };


		// R_min < r < R_max


		for (std::size_t k = 1; k < num_rs_minus_1; ++k) {
			auto [current_r, r_beta_minus_1, r_beta_0, r_beta_plus_1, unused_entry_10, unused_entry_11, unused_entry_12] = discretization->RCoefficientsByIndex(k);

			// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1), (-1, -1, 0), (1, -1, 0), (-1, 1, 0), (1, 1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1), (0, 1, 1)
			sparse_A0_.EntryByIndex(num_ss_minus_1, j, k) = { v_beta_0 * r_beta_0 * v_r_partial_coefficient, /*1*/\
				0/*2*/, \
				0/*3*/, \
				(v_beta_minus_1 * r_beta_0 * v_r_partial_coefficient) /*4*/, \
				(v_beta_plus_1 * r_beta_0 * v_r_partial_coefficient) /*5*/, \
				(v_beta_0 * r_beta_minus_1 * v_r_partial_coefficient) /*6*/, \
				(v_beta_0 * r_beta_plus_1 * v_r_partial_coefficient) /*7*/, \
				0 /*8*/, \
				0 /*9*/, \
				0 /*10*/, \
				0 /*11*/, \
				0 /*12*/, \
				0 /*13*/, \
				0 /*14*/, \
				0 /*15*/, \
				(v_beta_minus_1 * r_beta_minus_1 * v_r_partial_coefficient) /*16*/, \
				(v_beta_plus_1 * r_beta_minus_1 * v_r_partial_coefficient)/*17*/, \
				(v_beta_minus_1 * r_beta_plus_1 * v_r_partial_coefficient) /*18*/, \
				(v_beta_plus_1 * r_beta_plus_1 * v_r_partial_coefficient) /*19*/ };
		}

		// r = R_max who cares

		sparse_A0_.EntryByIndex(num_ss_minus_1, j, num_rs_minus_1) = { 0, /*1*/\
			0 /*2*/, \
			0/*3*/, \
			0 /*4*/, \
			0 /*5*/, \
			0 /*6*/, \
			0 /*7*/, \
			0 /*8*/, \
			0 /*9*/, \
			0 /*10*/, \
			0 /*11*/, \
			0 /*12*/, \
			0/*13*/, \
			0/*14*/, \
			0/*15*/, \
			0 /*16*/, \
			0/*17*/, \
			0 /*18*/, \
			0 /*19*/ };

	}



	for (std::size_t j = mean_vol_index_; j < num_vs_; ++j) {
		auto [current_v, v_alpha_minus_2, v_alpha_minus_1, v_alpha_0, unused_entry_1, unused_entry_2, unused_entry_3, unused_entry_4, unused_entry_5, unused_entry_6] = discretization->VCoefficientsByIndex(j);

		const double sqrt_v = std::sqrt(current_v);
		double  s_v_intermediate_partial_coefficient_1 = s_v_intermediate_partial_coefficient_0 * current_v;
		double  s_r_intermediate_partial_coefficient_1 = s_r_intermediate_partial_coefficient_0 * sqrt_v;

		double  v_r_partial_coefficient = v_r_intermediate_partial_coefficient_0 * sqrt_v;

		for (std::size_t i = 0; i < num_ss_minus_1; ++i) {
			auto [current_s, s_beta_minus_1, s_beta_0, s_beta_plus_1, unused_entry_7, unused_entry_8, unused_entry_9] = discretization->SCoefficientsByIndex(i);

			double  s_v_partial_coefficient = s_v_intermediate_partial_coefficient_1 * current_s;
			double  s_r_partial_coefficient = s_r_intermediate_partial_coefficient_1 * current_s;


			// r = R_min

			sparse_A0_.EntryByIndex(i, j, 0) = { v_alpha_0 * s_beta_0 * s_v_partial_coefficient, /*1*/\
					(s_beta_minus_1 * v_alpha_0 * s_v_partial_coefficient) /*2*/, \
					(s_beta_plus_1 * v_alpha_0 * s_v_partial_coefficient) /*3*/, \
					(s_beta_0 * v_alpha_minus_2 * s_v_partial_coefficient)  /*4*/, \
					(s_beta_0 * v_alpha_minus_1 * s_v_partial_coefficient)  /*5*/, \
					0 /*6*/, \
					0 /*7*/, \
					(s_beta_minus_1 * v_alpha_minus_2 * s_v_partial_coefficient) /*8*/, \
					(s_beta_plus_1 * v_alpha_minus_2 * s_v_partial_coefficient) /*9*/, \
					(s_beta_minus_1 * v_alpha_minus_1 * s_v_partial_coefficient) /*10*/, \
					(s_beta_plus_1 * v_alpha_minus_1 * s_v_partial_coefficient) /*11*/, \
					0 /*12*/, \
					0/*13*/, \
					0/*14*/, \
					0/*15*/, \
					0 /*16*/, \
					0/*17*/, \
					0 /*18*/, \
					0 /*19*/ };



			for (std::size_t k = 1; k < num_rs_minus_1; ++k) {
				auto [current_r, r_beta_minus_1, r_beta_0, r_beta_plus_1, unused_entry_10, unused_entry_11, unused_entry_12] = discretization->RCoefficientsByIndex(k);

				// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -2, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (-1, -2, 0), (1, -2, 0), (-1, -1, 0), (1, -1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -2, -1), (0, -1, -1), (0, -2, 1), (0, -1, 1)
				sparse_A0_.EntryByIndex(i, j, k) = { (s_beta_0 * v_alpha_0 * s_v_partial_coefficient) + \
					(s_beta_0 * r_beta_0 * s_r_partial_coefficient) + (v_alpha_0 * r_beta_0 * v_r_partial_coefficient), /*1*/\
					(s_beta_minus_1 * v_alpha_0 * s_v_partial_coefficient) + (s_beta_minus_1 * r_beta_0 * s_r_partial_coefficient) /*2*/, \
					(s_beta_plus_1 * v_alpha_0 * s_v_partial_coefficient) + (s_beta_plus_1 * r_beta_0 * s_r_partial_coefficient) /*3*/, \
					(s_beta_0 * v_alpha_minus_2 * s_v_partial_coefficient) + (v_alpha_minus_2 * r_beta_0 * v_r_partial_coefficient) /*4*/, \
					(s_beta_0 * v_alpha_minus_1 * s_v_partial_coefficient) + (v_alpha_minus_1 * r_beta_0 * v_r_partial_coefficient) /*5*/, \
					(s_beta_0 * r_beta_minus_1 * s_r_partial_coefficient) + (v_alpha_0 * r_beta_minus_1 * v_r_partial_coefficient) /*6*/, \
					(s_beta_0 * r_beta_plus_1 * s_r_partial_coefficient) + (v_alpha_0 * r_beta_plus_1 * v_r_partial_coefficient) /*7*/, \
					s_beta_minus_1 * v_alpha_minus_2 * s_v_partial_coefficient /*8*/, \
					s_beta_plus_1 * v_alpha_minus_2 * s_v_partial_coefficient /*9*/, \
					s_beta_minus_1 * v_alpha_minus_1 * s_v_partial_coefficient /*10*/, \
					s_beta_plus_1 * v_alpha_minus_1 * s_v_partial_coefficient /*11*/, \
					s_beta_minus_1 * r_beta_minus_1 * s_r_partial_coefficient /*12*/, \
					s_beta_plus_1 * r_beta_minus_1 * s_r_partial_coefficient /*13*/, \
					s_beta_minus_1 * r_beta_plus_1 * s_r_partial_coefficient /*14*/, \
					s_beta_plus_1 * r_beta_plus_1 * s_r_partial_coefficient /*15*/, \
					v_alpha_minus_2 * r_beta_minus_1 * v_r_partial_coefficient /*16*/, \
					v_alpha_minus_1 * r_beta_minus_1 * v_r_partial_coefficient /*17*/, \
					v_alpha_minus_2 * r_beta_plus_1 * v_r_partial_coefficient /*18*/, \
					v_alpha_minus_1 * r_beta_plus_1 * v_r_partial_coefficient /*19*/ };
			}

			// r = R_max

			sparse_A0_.EntryByIndex(i, j, num_rs_minus_1) = { v_alpha_0 * s_beta_0 * s_v_partial_coefficient, /*1*/\
					(s_beta_minus_1 * v_alpha_0 * s_v_partial_coefficient) /*2*/, \
					(s_beta_plus_1 * v_alpha_0 * s_v_partial_coefficient) /*3*/, \
					(s_beta_0 * v_alpha_minus_2 * s_v_partial_coefficient)  /*4*/, \
					(s_beta_0 * v_alpha_minus_1 * s_v_partial_coefficient)  /*5*/, \
					0 /*6*/, \
					0 /*7*/, \
					(s_beta_minus_1 * v_alpha_minus_2 * s_v_partial_coefficient) /*8*/, \
					(s_beta_plus_1 * v_alpha_minus_2 * s_v_partial_coefficient) /*9*/, \
					(s_beta_minus_1 * v_alpha_minus_1 * s_v_partial_coefficient) /*10*/, \
					(s_beta_plus_1 * v_alpha_minus_1 * s_v_partial_coefficient) /*11*/, \
					0 /*12*/, \
					0/*13*/, \
					0/*14*/, \
					0/*15*/, \
					0 /*16*/, \
					0/*17*/, \
					0 /*18*/, \
					0 /*19*/ };
		}

		// r = R_min who cares
		sparse_A0_.EntryByIndex(num_ss_minus_1, j, 0) = { 0, /*1*/\
				0 /*2*/, \
				0/*3*/, \
				0 /*4*/, \
				0 /*5*/, \
				0 /*6*/, \
				0 /*7*/, \
				0 /*8*/, \
				0 /*9*/, \
				0 /*10*/, \
				0 /*11*/, \
				0 /*12*/, \
				0/*13*/, \
				0/*14*/, \
				0/*15*/, \
				0 /*16*/, \
				0/*17*/, \
				0 /*18*/, \
				0 /*19*/ };


		// R_min < r < R_max


		for (std::size_t k = 1; k < num_rs_minus_1; ++k) {
			auto [current_r, r_beta_minus_1, r_beta_0, r_beta_plus_1, unused_entry_10, unused_entry_11, unused_entry_12] = discretization->RCoefficientsByIndex(k);

			// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -2, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (-1, -2, 0), (1, -2, 0), (-1, -1, 0), (1, -1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -2, -1), (0, -1, -1), (0, -2, 1), (0, -1, 1)
			sparse_A0_.EntryByIndex(num_ss_minus_1, j, k) = { v_alpha_0 * r_beta_0 * v_r_partial_coefficient, /*1*/\
				0/*2*/, \
				0/*3*/, \
				(v_alpha_minus_2 * r_beta_0 * v_r_partial_coefficient) /*4*/, \
				(v_alpha_minus_1 * r_beta_0 * v_r_partial_coefficient) /*5*/, \
				(v_alpha_0 * r_beta_minus_1 * v_r_partial_coefficient) /*6*/, \
				(v_alpha_0 * r_beta_plus_1 * v_r_partial_coefficient) /*7*/, \
				0 /*8*/, \
				0 /*9*/, \
				0 /*10*/, \
				0 /*11*/, \
				0 /*12*/, \
				0 /*13*/, \
				0 /*14*/, \
				0 /*15*/, \
				(v_alpha_minus_2 * r_beta_minus_1 * v_r_partial_coefficient) /*16*/, \
				(v_alpha_minus_1 * r_beta_minus_1 * v_r_partial_coefficient)/*17*/, \
				(v_alpha_minus_2 * r_beta_plus_1 * v_r_partial_coefficient) /*18*/, \
				(v_alpha_minus_1 * r_beta_plus_1 * v_r_partial_coefficient) /*19*/ };
		}

		// r = R_max who cares

		sparse_A0_.EntryByIndex(num_ss_minus_1, j, num_rs_minus_1) = { 0, /*1*/\
			0 /*2*/, \
			0/*3*/, \
			0 /*4*/, \
			0 /*5*/, \
			0 /*6*/, \
			0 /*7*/, \
			0 /*8*/, \
			0 /*9*/, \
			0 /*10*/, \
			0 /*11*/, \
			0 /*12*/, \
			0/*13*/, \
			0/*14*/, \
			0/*15*/, \
			0 /*16*/, \
			0/*17*/, \
			0 /*18*/, \
			0 /*19*/ };


	}

}


// To help keep track of entries, use the keys:
// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1), (-1, -1, 0), (1, -1, 0), (-1, 1, 0), (1, 1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1), (0, 1, 1)
// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -2, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (-1, -2, 0), (1, -2, 0), (-1, -1, 0), (1, -1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -2, -1), (0, -1, -1), (0, -2, 1), (0, -1, 1)
template<typename GridType>
[[nodiscard]] GridType heston_hull_white::finite_difference::examples::A0<GridType>::TransformGrid(const GridType& vec) {
	// should assert that vec size is same as assumed A0 size;
	//const auto (num_ss, num_vs, num_rs) = vec.size();
	GridType grid(num_ss_, num_vs_, num_rs_);

	std::size_t num_ss_minus_1 = num_ss_ - 1;
	std::size_t num_rs_minus_1 = num_rs_ - 1;

	// s = s_0.

	{
		// v = 0 (everything is 0 because every term has a coefficient with v in it

		for (std::size_t k = 0; k < num_rs_; ++k) {
			grid.EntryByIndex(0, 0, k) = 0;
		}

		// 0 < v < v_bar

		for (std::size_t j = 1; j < mean_vol_index_; ++j) {
			// r = R_min
			{
				auto [center_coefficient, unused_zero_1, s_plus_1_coefficient, v_minus_1_coefficient, v_plus_1_coefficient, unused_zero_2, unused_zero_3, \
					unused_zero_4, s_plus_1_v_minus_1_coefficient, unused_zero_5, s_plus_1_v_plus_1_coefficient, unused_zero_6, unused_zero_7, unused_zero_8, \
					unused_zero_9, unused_zero_10, unused_zero_11, unused_zero_12, unused_zero_13] = sparse_A0_.EntryByIndex(0, j, 0);

				grid.EntryByIndex(0, j, 0) = (center_coefficient * vec.EntryByIndex(0, j, 0)) + (s_plus_1_coefficient * vec.EntryByIndex(1, j, 0)) + \
					(v_minus_1_coefficient * vec.EntryByIndex(0, j - 1, 0)) + (v_plus_1_coefficient * vec.EntryByIndex(0, j + 1, 0)) + \
					(s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(1, j - 1, 0)) + (s_plus_1_v_plus_1_coefficient * vec.EntryByIndex(1, j + 1, 0));
			}

			// R_min < r < R_max
			// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1), (-1, -1, 0), (1, -1, 0), (-1, 1, 0), (1, 1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1), (0, 1, 1)

			for (std::size_t k = 1; k < num_rs_minus_1; ++k) {
				auto [center_coefficient, unused_zero_1, s_plus_1_coefficient, v_minus_1_coefficient, v_plus_1_coefficient, r_minus_1_coefficient, r_plus_1_coefficient, \
					unused_zero_4, s_plus_1_v_minus_1_coefficient, unused_zero_5, s_plus_1_v_plus_1_coefficient, unused_zero_6, s_plus_1_r_minus_1_coefficient, unused_zero_8, \
					s_plus_1_r_plus_1_coefficient, v_minus_1_r_minus_1_coefficient, v_plus_1_r_minus_1_coefficient, v_minus_1_r_plus_1_coefficient, v_plus_1_r_plus_1_coefficient] = sparse_A0_.EntryByIndex(0, j, k);


				grid.EntryByIndex(0, j, k) = (center_coefficient * vec.EntryByIndex(0, j, k)) + (s_plus_1_coefficient * vec.EntryByIndex(1, j, k)) + \
					(v_minus_1_coefficient * vec.EntryByIndex(0, j - 1, k)) + (v_plus_1_coefficient * vec.EntryByIndex(0, j + 1, k)) + \
					(r_minus_1_coefficient * vec.EntryByIndex(0, j, k - 1)) + (r_plus_1_coefficient * vec.EntryByIndex(0, j, k + 1)) + \
					(s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(1, j - 1, k)) + (s_plus_1_v_plus_1_coefficient * vec.EntryByIndex(1, j + 1, k)) + \
					(s_plus_1_r_minus_1_coefficient * vec.EntryByIndex(1, j, k - 1)) + (s_plus_1_r_plus_1_coefficient * vec.EntryByIndex(1, j, k + 1)) + \
					(v_minus_1_r_minus_1_coefficient * vec.EntryByIndex(0, j - 1, k - 1)) + (v_plus_1_r_minus_1_coefficient * vec.EntryByIndex(0, j + 1, k - 1)) + \
					(v_minus_1_r_plus_1_coefficient * vec.EntryByIndex(0, j - 1, k + 1)) + (v_plus_1_r_plus_1_coefficient * vec.EntryByIndex(0, j + 1, k + 1));

			}


			// r = R_max
			{
				auto [center_coefficient, unused_zero_1, s_plus_1_coefficient, v_minus_1_coefficient, v_plus_1_coefficient, unused_zero_2, unused_zero_3, \
					unused_zero_4, s_plus_1_v_minus_1_coefficient, unused_zero_5, s_plus_1_v_plus_1_coefficient, unused_zero_6, unused_zero_7, unused_zero_8, \
					unused_zero_9, unused_zero_10, unused_zero_11, unused_zero_12, unused_zero_13] = sparse_A0_.EntryByIndex(0, j, num_rs_minus_1);

				grid.EntryByIndex(0, j, num_rs_minus_1) = (center_coefficient * vec.EntryByIndex(0, j, num_rs_minus_1)) + (s_plus_1_coefficient * vec.EntryByIndex(1, j, num_rs_minus_1)) + \
					(v_minus_1_coefficient * vec.EntryByIndex(0, j - 1, num_rs_minus_1)) + (v_plus_1_coefficient * vec.EntryByIndex(0, j + 1, num_rs_minus_1)) + \
					(s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(1, j - 1, num_rs_minus_1)) + (s_plus_1_v_plus_1_coefficient * vec.EntryByIndex(1, j + 1, num_rs_minus_1));
			}

		}

		// v_bar < v 

		for (std::size_t j = mean_vol_index_; j < num_vs_; ++j) {
			// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -2, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (-1, -2, 0), (1, -2, 0), (-1, -1, 0), (1, -1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -2, -1), (0, -1, -1), (0, -2, 1), (0, -1, 1)

			// r = R_min
			{
				auto [center_coefficient, unused_zero_1, s_plus_1_coefficient, v_minus_2_coefficient, v_minus_1_coefficient, unused_zero_2, unused_zero_3, \
					unused_zero_4, s_plus_1_v_minus_2_coefficient, unused_zero_5, s_plus_1_v_minus_1_coefficient, unused_zero_6, unused_zero_7, unused_zero_8, \
					unused_zero_9, unused_zero_10, unused_zero_11, unused_zero_12, unused_zero_13] = sparse_A0_.EntryByIndex(0, j, 0);

				grid.EntryByIndex(0, j, 0) = (center_coefficient * vec.EntryByIndex(0, j, 0)) + (s_plus_1_coefficient * vec.EntryByIndex(1, j, 0)) + \
					(v_minus_2_coefficient * vec.EntryByIndex(0, j - 2, 0)) + (v_minus_1_coefficient * vec.EntryByIndex(0, j - 1, 0)) + \
					(s_plus_1_v_minus_2_coefficient * vec.EntryByIndex(1, j - 2, 0)) + (s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(1, j - 1, 0));

			}

			// R_min < r < R_max
			for (std::size_t k = 1; k < num_rs_minus_1; ++k) {
				auto [center_coefficient, unused_zero_1, s_plus_1_coefficient, v_minus_2_coefficient, v_minus_1_coefficient, r_minus_1_coefficient, r_plus_1_coefficient, \
					unused_zero_4, s_plus_1_v_minus_2_coefficient, unused_zero_5, s_plus_1_v_minus_1_coefficient, unused_zero_6, s_plus_1_r_minus_1_coefficient, unused_zero_8, \
					s_plus_1_r_plus_1_coefficient, v_minus_2_r_minus_1_coefficient, v_minus_1_r_minus_1_coefficient, v_minus_2_r_plus_1_coefficient, v_minus_1_r_plus_1_coefficient] = sparse_A0_.EntryByIndex(0, j, k);


				grid.EntryByIndex(0, j, k) = (center_coefficient * vec.EntryByIndex(0, j, k)) + (s_plus_1_coefficient * vec.EntryByIndex(1, j, k)) + \
					(v_minus_2_coefficient * vec.EntryByIndex(0, j - 2, k)) + (v_minus_1_coefficient * vec.EntryByIndex(0, j - 1, k)) + \
					(r_minus_1_coefficient * vec.EntryByIndex(0, j, k - 1)) + (r_plus_1_coefficient * vec.EntryByIndex(0, j, k + 1)) + \
					(s_plus_1_v_minus_2_coefficient * vec.EntryByIndex(1, j - 2, k)) + (s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(1, j - 1, k)) + \
					(s_plus_1_r_minus_1_coefficient * vec.EntryByIndex(1, j, k - 1)) + (s_plus_1_r_plus_1_coefficient * vec.EntryByIndex(1, j, k + 1)) + \
					(v_minus_2_r_minus_1_coefficient * vec.EntryByIndex(0, j - 2, k - 1)) + (v_minus_1_r_minus_1_coefficient * vec.EntryByIndex(0, j - 1, k - 1)) + \
					(v_minus_2_r_plus_1_coefficient * vec.EntryByIndex(0, j - 2, k + 1)) + (v_minus_1_r_plus_1_coefficient * vec.EntryByIndex(0, j - 1, k + 1));

			}

			// r = R_max
			{
				auto [center_coefficient, unused_zero_1, s_plus_1_coefficient, v_minus_2_coefficient, v_minus_1_coefficient, unused_zero_2, unused_zero_3, \
					unused_zero_4, s_plus_1_v_minus_2_coefficient, unused_zero_5, s_plus_1_v_minus_1_coefficient, unused_zero_6, unused_zero_7, unused_zero_8, \
					unused_zero_9, unused_zero_10, unused_zero_11, unused_zero_12, unused_zero_13] = sparse_A0_.EntryByIndex(0, j, num_rs_minus_1);

				grid.EntryByIndex(0, j, num_rs_minus_1) = (center_coefficient * vec.EntryByIndex(0, j, num_rs_minus_1)) + (s_plus_1_coefficient * vec.EntryByIndex(1, j, num_rs_minus_1)) + \
					(v_minus_2_coefficient * vec.EntryByIndex(0, j - 2, num_rs_minus_1)) + (v_minus_1_coefficient * vec.EntryByIndex(0, j - 1, num_rs_minus_1)) + \
					(s_plus_1_v_minus_2_coefficient * vec.EntryByIndex(1, j - 2, num_rs_minus_1)) + (s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(1, j - 1, num_rs_minus_1));
			}


		}


	} // end s = s_0


	// s_0 < s < s_max

	for (std::size_t i = 1; i < num_ss_minus_1; ++i) {

		// v = 0 (everything is 0 because every term has a coefficient with v in it)

		for (std::size_t k = 0; k < num_rs_; ++k) {
			grid.EntryByIndex(i, 0, k) = 0;
		}

		// 0 < v < v_bar
		// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1), (-1, -1, 0), (1, -1, 0), (-1, 1, 0), (1, 1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1), (0, 1, 1)
		for (std::size_t j = 1; j < mean_vol_index_; ++j) {

			// r = R_min
			{
				auto [center_coefficient, s_minus_1_coefficient, s_plus_1_coefficient, v_minus_1_coefficient, v_plus_1_coefficient, unused_zero_1, unused_zero_2, \
					s_minus_1_v_minus_1_coefficient, s_plus_1_v_minus_1_coefficient, s_minus_1_v_plus_1_coefficient, s_plus_1_v_plus_1_coefficient, \
					unused_zero_3, unused_zero_4, unused_zero_5, unused_zero_6, \
					unused_zero_7, unused_zero_8, unused_zero_9, unused_zero_10] = sparse_A0_.EntryByIndex(i, j, 0);

				grid.EntryByIndex(i, j, 0) = (center_coefficient * vec.EntryByIndex(i, j, 0)) + \
					(s_minus_1_coefficient * vec.EntryByIndex(i - 1, j, 0)) + (s_plus_1_coefficient * vec.EntryByIndex(i + 1, j, 0)) + \
					(v_minus_1_coefficient * vec.EntryByIndex(i, j - 1, 0)) + (v_plus_1_coefficient * vec.EntryByIndex(i, j + 1, 0)) + \
					(s_minus_1_v_minus_1_coefficient * vec.EntryByIndex(i - 1, j - 1, 0)) + (s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(i + 1, j - 1, 0)) + \
					(s_minus_1_v_plus_1_coefficient * vec.EntryByIndex(i - 1, j + 1, 0)) + (s_plus_1_v_plus_1_coefficient * vec.EntryByIndex(i + 1, j + 1, 0));

			}

			// R_min < r < R_max
			for (std::size_t k = 1; k < num_rs_minus_1; ++k) {
				auto [center_coefficient, s_minus_1_coefficient, s_plus_1_coefficient, v_minus_1_coefficient, v_plus_1_coefficient, r_minus_1_coefficient, r_plus_1_coefficient, \
					s_minus_1_v_minus_1_coefficient, s_plus_1_v_minus_1_coefficient, s_minus_1_v_plus_1_coefficient, s_plus_1_v_plus_1_coefficient, \
					s_minus_1_r_minus_1_coefficient, s_plus_1_r_minus_1_coefficient, s_minus_1_r_plus_1_coefficient, s_plus_1_r_plus_1_coefficient, \
					v_minus_1_r_minus_1_coefficient, v_plus_1_r_minus_1_coefficient, v_minus_1_r_plus_1_coefficient, v_plus_1_r_plus_1_coefficient] = sparse_A0_.EntryByIndex(i, j, k);

				grid.EntryByIndex(i, j, k) = (center_coefficient * vec.EntryByIndex(i, j, k)) + \
					(s_minus_1_coefficient * vec.EntryByIndex(i - 1, j, k)) + (s_plus_1_coefficient * vec.EntryByIndex(i + 1, j, k)) + \
					(v_minus_1_coefficient * vec.EntryByIndex(i, j - 1, k)) + (v_plus_1_coefficient * vec.EntryByIndex(i, j + 1, k)) + \
					(r_minus_1_coefficient * vec.EntryByIndex(i, j, k - 1)) + (r_plus_1_coefficient * vec.EntryByIndex(i, j, k + 1)) + \
					(s_minus_1_v_minus_1_coefficient * vec.EntryByIndex(i - 1, j - 1, k)) + (s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(i + 1, j - 1, k)) + \
					(s_minus_1_v_plus_1_coefficient * vec.EntryByIndex(i - 1, j + 1, k)) + (s_plus_1_v_plus_1_coefficient * vec.EntryByIndex(i + 1, j + 1, k)) + \
					(s_minus_1_r_minus_1_coefficient * vec.EntryByIndex(i - 1, j, k - 1)) + (s_plus_1_r_minus_1_coefficient * vec.EntryByIndex(i + 1, j, k - 1)) + \
					(s_minus_1_r_plus_1_coefficient * vec.EntryByIndex(i - 1, j, k + 1)) + (s_plus_1_r_plus_1_coefficient * vec.EntryByIndex(i + 1, j, k + 1)) + \
					(v_minus_1_r_minus_1_coefficient * vec.EntryByIndex(i, j - 1, k - 1)) + (v_plus_1_r_minus_1_coefficient * vec.EntryByIndex(i, j + 1, k - 1)) + \
					(v_minus_1_r_plus_1_coefficient * vec.EntryByIndex(i, j - 1, k + 1)) + (v_plus_1_r_plus_1_coefficient * vec.EntryByIndex(i, j + 1, k + 1));
			}

			// r = R_max
			{
				auto [center_coefficient, s_minus_1_coefficient, s_plus_1_coefficient, v_minus_1_coefficient, v_plus_1_coefficient, unused_zero_1, unused_zero_2, \
					s_minus_1_v_minus_1_coefficient, s_plus_1_v_minus_1_coefficient, s_minus_1_v_plus_1_coefficient, s_plus_1_v_plus_1_coefficient, \
					unused_zero_3, unused_zero_4, unused_zero_5, unused_zero_6, \
					unused_zero_7, unused_zero_8, unused_zero_9, unused_zero_10] = sparse_A0_.EntryByIndex(i, j, num_rs_minus_1);

				grid.EntryByIndex(i, j, num_rs_minus_1) = (center_coefficient * vec.EntryByIndex(i, j, num_rs_minus_1)) + \
					(s_minus_1_coefficient * vec.EntryByIndex(i - 1, j, num_rs_minus_1)) + (s_plus_1_coefficient * vec.EntryByIndex(i + 1, j, num_rs_minus_1)) + \
					(v_minus_1_coefficient * vec.EntryByIndex(i, j - 1, num_rs_minus_1)) + (v_plus_1_coefficient * vec.EntryByIndex(i, j + 1, num_rs_minus_1)) + \
					(s_minus_1_v_minus_1_coefficient * vec.EntryByIndex(i - 1, j - 1, num_rs_minus_1)) + (s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(i + 1, j - 1, num_rs_minus_1)) + \
					(s_minus_1_v_plus_1_coefficient * vec.EntryByIndex(i - 1, j + 1, num_rs_minus_1)) + (s_plus_1_v_plus_1_coefficient * vec.EntryByIndex(i + 1, j + 1, num_rs_minus_1));
			}
		}


		// v_bar < v
		// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -2, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (-1, -2, 0), (1, -2, 0), (-1, -1, 0), (1, -1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -2, -1), (0, -1, -1), (0, -2, 1), (0, -1, 1)
		for (std::size_t j = mean_vol_index_; j < num_vs_; ++j) {

			// r = R_min
			{
				auto [center_coefficient, s_minus_1_coefficient, s_plus_1_coefficient, v_minus_2_coefficient, v_minus_1_coefficient, unused_zero_1, unused_zero_2, \
					s_minus_1_v_minus_2_coefficient, s_plus_1_v_minus_2_coefficient, s_minus_1_v_minus_1_coefficient, s_plus_1_v_minus_1_coefficient, \
					unused_zero_3, unused_zero_4, unused_zero_5, unused_zero_6, \
					unused_zero_7, unused_zero_8, unused_zero_9, unused_zero_10] = sparse_A0_.EntryByIndex(i, j, 0);

				grid.EntryByIndex(i, j, 0) = (center_coefficient * vec.EntryByIndex(i, j, 0)) + \
					(s_minus_1_coefficient * vec.EntryByIndex(i - 1, j, 0)) + (s_plus_1_coefficient * vec.EntryByIndex(i + 1, j, 0)) + \
					(v_minus_2_coefficient * vec.EntryByIndex(i, j - 2, 0)) + (v_minus_1_coefficient * vec.EntryByIndex(i, j - 1, 0)) + \
					(s_minus_1_v_minus_2_coefficient * vec.EntryByIndex(i - 1, j - 2, 0)) + (s_plus_1_v_minus_2_coefficient * vec.EntryByIndex(i + 1, j - 2, 0)) + \
					(s_minus_1_v_minus_1_coefficient * vec.EntryByIndex(i - 1, j - 1, 0)) + (s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(i + 1, j - 1, 0));
			}

			// R_min < r < R_max
			for (std::size_t k = 1; k < num_rs_minus_1; ++k) {
				auto [center_coefficient, s_minus_1_coefficient, s_plus_1_coefficient, v_minus_2_coefficient, v_minus_1_coefficient, r_minus_1_coefficient, r_plus_1_coefficient, \
					s_minus_1_v_minus_2_coefficient, s_plus_1_v_minus_2_coefficient, s_minus_1_v_minus_1_coefficient, s_plus_1_v_minus_1_coefficient, \
					s_minus_1_r_minus_1_coefficient, s_plus_1_r_minus_1_coefficient, s_minus_1_r_plus_1_coefficient, s_plus_1_r_plus_1_coefficient, \
					v_minus_2_r_minus_1_coefficient, v_minus_1_r_minus_1_coefficient, v_minus_2_r_plus_1_coefficient, v_minus_1_r_plus_1_coefficient] = sparse_A0_.EntryByIndex(i, j, k);

				grid.EntryByIndex(i, j, k) = (center_coefficient * vec.EntryByIndex(i, j, k)) + \
					(s_minus_1_coefficient * vec.EntryByIndex(i - 1, j, k)) + (s_plus_1_coefficient * vec.EntryByIndex(i + 1, j, k)) + \
					(v_minus_2_coefficient * vec.EntryByIndex(i, j - 2, k)) + (v_minus_1_coefficient * vec.EntryByIndex(i, j - 1, k)) + \
					(r_minus_1_coefficient * vec.EntryByIndex(i, j, k - 1)) + (r_plus_1_coefficient * vec.EntryByIndex(i, j, k + 1)) + \
					(s_minus_1_v_minus_2_coefficient * vec.EntryByIndex(i - 1, j - 2, k)) + (s_plus_1_v_minus_2_coefficient * vec.EntryByIndex(i + 1, j - 2, k)) + \
					(s_minus_1_v_minus_1_coefficient * vec.EntryByIndex(i - 1, j - 1, k)) + (s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(i + 1, j - 1, k)) + \
					(s_minus_1_r_minus_1_coefficient * vec.EntryByIndex(i - 1, j, k - 1)) + (s_plus_1_r_minus_1_coefficient * vec.EntryByIndex(i + 1, j, k - 1)) + \
					(s_minus_1_r_plus_1_coefficient * vec.EntryByIndex(i - 1, j, k + 1)) + (s_plus_1_r_plus_1_coefficient * vec.EntryByIndex(i + 1, j, k + 1)) + \
					(v_minus_2_r_minus_1_coefficient * vec.EntryByIndex(i, j - 2, k - 1)) + (v_minus_1_r_minus_1_coefficient * vec.EntryByIndex(i, j - 1, k - 1)) + \
					(v_minus_2_r_plus_1_coefficient * vec.EntryByIndex(i, j - 2, k + 1)) + (v_minus_1_r_plus_1_coefficient * vec.EntryByIndex(i, j - 1, k + 1));

			}

			// r = R_max
			{
				auto [center_coefficient, s_minus_1_coefficient, s_plus_1_coefficient, v_minus_2_coefficient, v_minus_1_coefficient, unused_zero_1, unused_zero_2, \
					s_minus_1_v_minus_2_coefficient, s_plus_1_v_minus_2_coefficient, s_minus_1_v_minus_1_coefficient, s_plus_1_v_minus_1_coefficient, \
					unused_zero_3, unused_zero_4, unused_zero_5, unused_zero_6, \
					unused_zero_7, unused_zero_8, unused_zero_9, unused_zero_10] = sparse_A0_.EntryByIndex(i, j, num_rs_minus_1);

				grid.EntryByIndex(i, j, num_rs_minus_1) = (center_coefficient * vec.EntryByIndex(i, j, num_rs_minus_1)) + \
					(s_minus_1_coefficient * vec.EntryByIndex(i - 1, j, num_rs_minus_1)) + (s_plus_1_coefficient * vec.EntryByIndex(i + 1, j, num_rs_minus_1)) + \
					(v_minus_2_coefficient * vec.EntryByIndex(i, j - 2, num_rs_minus_1)) + (v_minus_1_coefficient * vec.EntryByIndex(i, j - 1, num_rs_minus_1)) + \
					(s_minus_1_v_minus_2_coefficient * vec.EntryByIndex(i - 1, j - 2, num_rs_minus_1)) + (s_plus_1_v_minus_2_coefficient * vec.EntryByIndex(i + 1, j - 2, num_rs_minus_1)) + \
					(s_minus_1_v_minus_1_coefficient * vec.EntryByIndex(i - 1, j - 1, num_rs_minus_1)) + (s_plus_1_v_minus_1_coefficient * vec.EntryByIndex(i + 1, j - 1, num_rs_minus_1));
			}
		}
	}

	// s = s_max

	{
		// v = 0 (everything is 0 because every term has a coefficient with v in it)

		for (std::size_t k = 0; k < num_rs_; ++k) {
			grid.EntryByIndex(num_ss_minus_1, 0, k) = 0;
		}


		// 0 < v < v_bar
		// NOTE: Here I copy the pattern I used before even though there is nothing to do at R_max and R_min.
		for (std::size_t j = 1; j < mean_vol_index_; ++j) {

			// r = R_min, everything is 0
			grid.EntryByIndex(num_ss_minus_1, j, 0) = 0;


			// R_min < r < R_max
			for (std::size_t k = 1; k < num_rs_minus_1; ++k) {

				auto [center_coefficient, unused_entry_1, unused_entry_2, v_minus_1_coefficient, v_plus_1_coefficient, r_minus_1_coefficient, r_plus_1_coefficient, \
					unused_entry_3, unused_entry_4, unused_entry_5, unused_entry_6, \
					unused_entry_7, unused_entry_8, unused_entry_9, unused_entry_10, \
					v_minus_1_r_minus_1_coefficient, v_plus_1_r_minus_1_coefficient, v_minus_1_r_plus_1_coefficient, v_plus_1_r_plus_1_coefficient] = sparse_A0_.EntryByIndex(num_ss_minus_1, j, k);

				grid.EntryByIndex(num_ss_minus_1, j, k) = (center_coefficient * vec.EntryByIndex(num_ss_minus_1, j, k)) + \
					(v_minus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j - 1, k)) + (v_plus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j + 1, k)) + \
					(r_minus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j, k - 1)) + (r_plus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j, k + 1)) + \
					(v_minus_1_r_minus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j - 1, k - 1)) + (v_plus_1_r_minus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j + 1, k - 1)) + \
					(v_minus_1_r_plus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j - 1, k + 1)) + (v_plus_1_r_plus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j + 1, k + 1));

			}

			// r = R_max
			grid.EntryByIndex(num_ss_minus_1, j, num_rs_minus_1) = 0;

		}




		// v_bar < v_max
		// NOTE: Here I copy the pattern I used before even though there is nothing to do at R_max and R_min.
		for (std::size_t j = mean_vol_index_; j < num_vs_; ++j) {

			// r = R_min, everything is 0
			grid.EntryByIndex(num_ss_minus_1, j, 0) = 0;


			// R_min < r < R_max
			for (std::size_t k = 1; k < num_rs_minus_1; ++k) {

				auto [center_coefficient, unused_entry_1, unused_entry_2, v_minus_2_coefficient, v_minus_1_coefficient, r_minus_1_coefficient, r_plus_1_coefficient, \
					unused_entry_3, unused_entry_4, unused_entry_5, unused_entry_6, \
					unused_entry_7, unused_entry_8, unused_entry_9, unused_entry_10, \
					v_minus_2_r_minus_1_coefficient, v_minus_1_r_minus_1_coefficient, v_minus_2_r_plus_1_coefficient, v_minus_1_r_plus_1_coefficient] = sparse_A0_.EntryByIndex(num_ss_minus_1, j, k);

				grid.EntryByIndex(num_ss_minus_1, j, k) = (center_coefficient * vec.EntryByIndex(num_ss_minus_1, j, k)) + \
					(v_minus_2_coefficient * vec.EntryByIndex(num_ss_minus_1, j - 2, k)) + (v_minus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j - 1, k)) + \
					(r_minus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j, k - 1)) + (r_plus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j, k + 1)) + \
					(v_minus_2_r_minus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j - 2, k - 1)) + (v_minus_1_r_minus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j - 1, k - 1)) + \
					(v_minus_2_r_plus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j - 2, k + 1)) + (v_minus_1_r_plus_1_coefficient * vec.EntryByIndex(num_ss_minus_1, j - 1, k + 1));

			}

			// r = R_max
			grid.EntryByIndex(num_ss_minus_1, j, num_rs_minus_1) = 0;

		}
	}

	return grid;
}

#endif