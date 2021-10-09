#ifndef HESTONHULLWHITEMONTECARLO_EUROPEAN_EXAMPLES_PATHINDEPENDENTS_VANILLACALL2_H
#define HESTONHULLWHITEMONTECARLO_EUROPEAN_EXAMPLES_PATHINDEPENDENTS_VANILLACALL2_H

#include <algorithm>
#include <vector>

namespace heston_hull_white {
	namespace monte_carlo {
			namespace examples {
				// Note: There is no reason that you can't combine calls and puts.
				template<typename StateType>
				class VanillaCall2
				{
				public:
					constexpr VanillaCall2(const std::vector<double>& strikes) : strikes_(strikes) {} // perhaps passing by value will cause RVO for temporaries?

					constexpr std::vector<double> Price(const StateType& state) const {
						std::vector<double> result;
						double value = state.Value();
						double discount_factor = state.DiscountFactor();
						std::size_t size = strikes_.size();
						result.reserve(size);
						for (int i = 0; i < size; ++i) {
							//result.push_back(discount_factor * std::max((value - strikes_[i]), 0.0));
							result.push_back(std::max((value - strikes_[i]), 0.0));
						}
						return result;
					}

				private:
					std::vector<double> strikes_;
				};
		}
	}
}
#endif 
