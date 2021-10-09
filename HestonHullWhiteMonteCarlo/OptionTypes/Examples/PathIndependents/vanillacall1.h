#ifndef HESTONHULLWHITEMONTECARLO_EUROPEAN_EXAMPLES_PATHINDEPENDENTS_VANILLACALL1_H
#define HESTONHULLWHITEMONTECARLO_EUROPEAN_EXAMPLES_PATHINDEPENDENTS_VANILLACALL1_H

#include <algorithm>

namespace heston_hull_white {
	namespace monte_carlo {
			namespace examples {
				// Note: There is no reason that you can't combine calls and puts.
				template<typename StateType>
				class VanillaCall1
				{
				public:
					constexpr VanillaCall1(double strike) : strike_(strike) {}

					constexpr double Price(const StateType& state) const {
						//return state.DiscountFactor() * std::max((state.Value() - strike_), 0.0);
						return std::max((state.Value() - strike_), 0.0);
					}

				private:
					double strike_;
				};
		}
	}
}

#endif 
