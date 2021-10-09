#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_CRAIGSNEYD_CRAIGSNEYD_H
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_CRAIGSNEYD_CRAIGSNEYD_H

// TODO: Add constructor where the A_i's are given as input

// TODO: store -epsilon_*theta_ instead of theta. pass that to the implicit_inversion solvers instead of theta and epsilon

// TODO: Prove that the number of copies and matrix transformations is optimal

#include <cassert>
#include <memory>

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A0.h"
#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A1.h"
#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A1ImplicitInversion.h"
#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A2.h"
#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A2ImplicitInversion.h"
#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A3.h"
#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A3Builder.h"
#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A3ImplicitInversion.h"

#include <HestonHullWhiteFiniteDifference/DiscretizationExamples/ADIDiscretization/ADIDiscretization.h>
#include "HestonHullWhiteFiniteDifference/Schemes/HestonHullWhiteADIInputs/hestonhullwhiteadiinputs.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace examples {

			template<typename GridType>
			class CraigSneyd {
			public:
				CraigSneyd(double expiry_time, double epsilon, HestonHullWhiteADIInputs inputs, std::shared_ptr<const ADIDiscretization> discretization) : CraigSneyd(expiry_time, epsilon, inputs, discretization, (discretization->GetIndexV(inputs.mean_volatility_)))
				{}

				// goes backwards, so computes prices epsilon before current time
				GridType Next(const GridType& current_prices);

			private:
				CraigSneyd(double expiry_time, double epsilon, HestonHullWhiteADIInputs inputs, std::shared_ptr<const ADIDiscretization> discretization, std::size_t mean_vol_index) : \
					current_time_(expiry_time), theta_(.5), epsilon_(epsilon), \
					a0_(discretization, mean_vol_index, inputs.sigma_, inputs.eta_, inputs.rho_s_v_, inputs.rho_s_r_, inputs.rho_v_r_), \
					a1_(discretization), a1_implicit_inversion_(a1_, theta_, epsilon_), \
					a2_(inputs.mean_volatility_, inputs.kappa_, inputs.sigma_, mean_vol_index, discretization), a2_implicit_inversion_(a2_, theta_, epsilon_), \
					a3_builder_(discretization, inputs.lambda_, inputs.eta_, inputs.interest_curve_), \
					previous_a3_(std::make_unique<A3<GridType>>(a3_builder_.GetA3(expiry_time))) {}
				
				double current_time_;
				const double theta_;
				const double epsilon_;
				A0<GridType> a0_;
				A1<GridType> a1_;
				A1ImplicitInversion<GridType> a1_implicit_inversion_;
				A2<GridType> a2_;
				A2ImplicitInversion<GridType> a2_implicit_inversion_;
				A3Builder<GridType> a3_builder_;
				std::unique_ptr<A3<GridType>> previous_a3_;
			};

		}
	}
}

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/CraigSneyd/src/craigsneyd.tpp"

#endif