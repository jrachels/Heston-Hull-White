#include "HestonHullWhiteFiniteDifference/DiscretizationExamples/ADIDiscretization/ADIDiscretization.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace examples {
			void ADIDiscretization::ConstructSDiscretization(const ADIDiscretizationInputs& inputs) {
				// do ss
				std::size_t m_1 = inputs.m_1;
				double d_1 = inputs.d_1;
				double S_left = inputs.S_left;
				double S_right = inputs.S_right;

				//ss_.reserve(m_1 + 1);
				ss_.reserve(num_ss_);
				double xi_min = std::asinh(-S_left / d_1);

				double xi_int = (S_right - S_left) / (d_1);

				double xi_max = xi_int + std::asinh((inputs.S_max - S_right) / d_1);

				double delta_xi = (xi_max - xi_min) / (static_cast<double>(m_1));

				//double current_s = phi(xi_min);
				//double current_s = 0;
				 
				
				// REMOVED TO BE BOUNDARY CONDITION
				//ss_.push_back({ 0, 0, 0, 0, 0, 0, 0 });

				//double previous_s = 0; // start at 0
				xi_min += delta_xi;
				double current_s = Phi(xi_min, S_left, S_right, d_1, xi_int);
				double delta_s_current = current_s;
				//int num_betas = m_1;
				//for (int i = 1; i < num_betas; ++i) {
				for (std::size_t i = 1; i < num_ss_; ++i) { // start at 1 because we need to add s_max at end
					xi_min += delta_xi;
					double previous_s = current_s;
					current_s = Phi(xi_min, S_left, S_right, d_1, xi_int);

					double delta_s_previous = delta_s_current;
					delta_s_current = current_s - previous_s;

					ss_.push_back({ previous_s, BetaMinus1(delta_s_previous, delta_s_current), \
						Beta0(delta_s_previous, delta_s_current), BetaPlus1(delta_s_previous, delta_s_current),\
						DeltaMinus1(delta_s_previous, delta_s_current), \
						Delta0(delta_s_previous, delta_s_current), DeltaPlus1(delta_s_previous, delta_s_current) });
				}
				ss_.push_back({ current_s, 0, 0, 0, DeltaMinus1(delta_s_current, delta_s_current), Delta0(delta_s_current, delta_s_current), DeltaPlus1(delta_s_current, delta_s_current) });

			}


			void ADIDiscretization::ConstructVDiscretization(const ADIDiscretizationInputs& inputs) {
				double d_2 = inputs.d_2;
				std::size_t m_2 = inputs.m_2;
				//vs_.reserve(m_2 + 1);
				vs_.reserve(num_vs_);

				double delta_eta = (1.0 / static_cast<double>(m_2)) * (std::asinh(inputs.V_max / d_2));
				double eta = delta_eta;

				double previous_v = d_2 * std::sinh(eta);
				double delta_v_minus_1 = previous_v;

				eta += delta_eta;

				double current_v = d_2 * std::sinh(eta);

				double delta_v_current = current_v - previous_v;

				vs_.push_back({ 0, Gamma0(delta_v_minus_1, delta_v_current), \
				GammaPlus1(delta_v_minus_1, delta_v_current), GammaPlus2(delta_v_minus_1, delta_v_current),\
				0, 0, 0, 0, 0, 0 });

				vs_.push_back({ previous_v, 0, 0, 0, BetaMinus1(delta_v_minus_1, delta_v_current), \
				Beta0(delta_v_minus_1, delta_v_current), BetaPlus1(delta_v_minus_1, delta_v_current),\
				DeltaMinus1(delta_v_minus_1, delta_v_current), \
				Delta0(delta_v_minus_1, delta_v_current), DeltaPlus1(delta_v_minus_1, delta_v_current) });


				//int num_betas = m_2;
				//for (int i = 2; i < num_betas; ++i) {
				for (std::size_t i = 2; i < num_vs_; ++i) {
					eta += delta_eta;
					previous_v = current_v;
					current_v = d_2 * std::sinh(eta);

					double delta_v_minus_2 = delta_v_minus_1;

					delta_v_minus_1 = delta_v_current;

					delta_v_current = current_v - previous_v;


					vs_.push_back({ previous_v, AlphaMinus2(delta_v_minus_2, delta_v_minus_1), AlphaMinus1(delta_v_minus_2, delta_v_minus_1), \
						Alpha0(delta_v_minus_2, delta_v_minus_1), BetaMinus1(delta_v_minus_1, delta_v_current), \
						Beta0(delta_v_minus_1, delta_v_current), BetaPlus1(delta_v_minus_1, delta_v_current),\
						DeltaMinus1(delta_v_minus_1, delta_v_current), \
						Delta0(delta_v_minus_1, delta_v_current), DeltaPlus1(delta_v_minus_1, delta_v_current) });
				}

				// REMOVED TO BE BOUNDARY CONDITION

				//vs_.push_back({ current_v, AlphaMinus2(delta_v_minus_1, delta_v_current), AlphaMinus1(delta_v_minus_1, delta_v_current), Alpha0(delta_v_minus_1, delta_v_current), 0, 0, 0, 0, 0, 0 });

				
			}

			void ADIDiscretization::ConstructRDiscretization(const ADIDiscretizationInputs& inputs) {
				std::size_t m_3 = inputs.m_3;
				double d_3 = inputs.d_3;
				double c = inputs.c;
				double R_max = inputs.R_max;
				//rs_.reserve(m_3 + 1);
				rs_.reserve(num_rs_);

				double zeta = std::asinh((-R_max - c) / d_3);

				double delta_zeta = (1 / static_cast<double>(m_3)) * (std::asinh((R_max - c) / d_3) - std::asinh((-R_max - c) / d_3));
				double previous_r = -R_max;

				zeta += delta_zeta;
				double current_r = c + d_3 * std::sinh(zeta);

				double delta_r_current = current_r - previous_r;

				rs_.push_back({ previous_r, 0, 0, 0, DeltaMinus1(delta_r_current, delta_r_current), Delta0(delta_r_current, delta_r_current), DeltaPlus1(delta_r_current, delta_r_current) });

				for (std::size_t i = 1; i < m_3; ++i) {

					zeta += delta_zeta;
					double previous_r = current_r;
					current_r = c + d_3 * std::sinh(zeta);

					double delta_r_previous = delta_r_current;
					delta_r_current = current_r - previous_r;

					rs_.push_back({ previous_r, BetaMinus1(delta_r_previous, delta_r_current), Beta0(delta_r_previous, delta_r_current), \
						BetaPlus1(delta_r_previous, delta_r_current), DeltaMinus1(delta_r_previous, delta_r_current), \
						Delta0(delta_r_previous, delta_r_current), DeltaPlus1(delta_r_previous, delta_r_current) });
				}
				rs_.push_back({ current_r, 0, 0, 0, DeltaMinus1(delta_r_current, delta_r_current), Delta0(delta_r_current, delta_r_current), DeltaPlus1(delta_r_current, delta_r_current) });
			}


		}
	}
}