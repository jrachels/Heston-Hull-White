#ifndef HESTONHULLWHITEMONTECARLO_EUROPEAN_EXAMPLES_PATHDEPENDENTS_ASIAN1_H
#define HESTONHULLWHITEMONTECARLO_EUROPEAN_EXAMPLES_PATHDEPENDENTS_ASIAN1_H

#include <algorithm>

// TODO: make sure Asian1 updates correctly on the first and last tick. I'm not sure if the average includes the first and last tick.
// make sure the last tick isn't updated twice.
namespace heston_hull_white {
	namespace monte_carlo {
			namespace examples {
				template<typename StateType>
				class AsianCall1
				{
				public:
					constexpr AsianCall1(const double strike) : num_dates_(0), average_price_(0), strike_(strike) {}

					void Update(StateType& state) {
						//Welford's Method
						num_dates_ += 1;
						double current_average = average_price_;
						average_price_ = current_average + ((state.Value()) - current_average) / num_dates_;
					}

					constexpr double Price(const StateType& state) const {
						// maybe this should update first?
						// Update(state);
						/*double discount_factor = state.DiscountFactor();
						return discount_factor * std::max((average_price_ - strike_), 0.0);*/
						return std::max((average_price_ - strike_), 0.0);
					}

					/*constexpr double Price() const {

						return average_price_ - strike_;
					}*/

					void ResetAndInitialize() {
						num_dates_ = 0;
						average_price_ = 0;
					}

				private:
					int num_dates_;
					double average_price_;
					double strike_;
				};
		}
	}
}

#endif 
