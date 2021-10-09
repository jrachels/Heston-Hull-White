#ifndef HESTONHULLWHITEMONTECARLO_EUROPEAN_EXAMPLES_PATHINDEPENDENTS_VANILLACALL3_H
#define HESTONHULLWHITEMONTECARLO_EUROPEAN_EXAMPLES_PATHINDEPENDENTS_VANILLACALL3_H

#include <algorithm>
#include <array>

namespace heston_hull_white {
	namespace monte_carlo {
			namespace examples {
				// Note: There is no reason that you can't combine calls and puts.
				template<typename StateType, size_t N>
				class VanillaCall3
				{
				public:
					constexpr VanillaCall3(const std::array<double, N>& strikes) : strikes_(strikes) {}// perhaps passing by value will cause RVO for temporaries?

					constexpr std::array<double, N> Price(const StateType& state) const {
						std::array<double, N> result;
						double value = state.Value();
						double discount_factor = state.DiscountFactor();
						for (int i = 0; i < N; ++i) {
							//result[i] = (discount_factor * std::max((value - strikes_[i]), 0.0));
							result[i] = (std::max((value - strikes_[i]), 0.0));
						}
						return result;
					}

				private:
					std::array<double, N> strikes_;
				};
		}
	}
}

#endif 
