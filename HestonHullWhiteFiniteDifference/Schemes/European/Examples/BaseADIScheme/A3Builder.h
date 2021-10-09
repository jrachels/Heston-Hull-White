#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A3BUILDER_H
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A3BUILDER_H

#include <cmath>
#include <memory>

#include <HestonHullWhiteFiniteDifference/DiscretizationExamples/ADIDiscretization/ADIDiscretization.h>
#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A3.h"
#include "HestonHullWhiteFiniteDifference/GridExamples/ThreeDimensionalGrid/threedimensionalgrid.h"
#include "HestonHullWhiteFiniteDifference/util/InterestCurve/interestcurve.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace examples {

			template<typename GridType>
			class A3Builder {
			public:
				A3Builder(std::shared_ptr<const ADIDiscretization> discretization, const double lambda, const double eta, std::shared_ptr<const util::InterestCurve>  interest_curve) : lambda_(lambda), eta_(eta), theta_const_(.5*(std::pow((eta/lambda), 2))), interest_curve_(std::move(interest_curve)), discretization_(std::move(discretization)) {
				}

				// complicating issue: Finite difference looks backwards, but interest curve looks forward
				A3<GridType> GetA3(double time) {
					// step 1: compute theta_t from interest_curve (will I need a bunch of theta_t's here?
					// step 2: construct A3
					return A3<GridType>(lambda_, eta_, Theta(time), discretization_);
				}

				// TO DO: should I overload this with Theta(const double time_t) to make const_0 a constexpr?
				// 
				
				double Theta(const double time_t) const {
					return interest_curve_->Theta(time_t, lambda_, eta_);
				}

			private:
				const double lambda_;
				const double eta_;
				const double theta_const_;
				std::shared_ptr<const util::InterestCurve> interest_curve_;
				std::shared_ptr<const ADIDiscretization> discretization_;

			};
			
		}
	}
}


#endif