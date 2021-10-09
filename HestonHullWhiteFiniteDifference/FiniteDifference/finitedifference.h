#ifndef HESTONHULLWHITEFINITEDIFFERENCE_FINITEDIFFERENCE_FINITEDIFFERENCE_H
#define HESTONHULLWHITEFINITEDIFFERENCE_FINITEDIFFERENCE_FINITEDIFFERENCE_H

// TODO: create overload with concepts that allows the intrinsic_pricer to accept time as an input variable in addition to a point of the discretization.

// TODO: Allow epsilon to vary at each step.

// TODO: Scheme should be a template parameter of the constructor? Find a work-around?

// TODO: Add a variadic template parameter to the FiniteDifference constructor that allows the addition of any other parameters that may be required required for the scheme
// and pass them to the scheme

#include <cassert>
#include <concepts>
#include <memory>
#include <stdexcept>
#include <vector>

namespace heston_hull_white {
	namespace finite_difference {

		template<typename Scheme, typename Discretization, typename Grid>
		class FiniteDifference {
		private:
			std::shared_ptr<const Discretization> discretization_;

		public:


			// in HHW, int asset_price_index, int variance_index, int short_rate_index
			template<typename... Args>
			requires requires (Grid g, Args... args) {
				{g.EntryByIndex(discretization_->GetIndex(args...))};
				//{discretization_->GetIndex(args...)};
			}
			//std::array<int, sizeof...(Args)> GetIndex(int time_index, Args... args) {
			std::array<std::size_t, sizeof...(Args)> GetIndex(Args... args) {
				return discretization_->GetIndex(args...);
			}

			// European Finite Difference (see Concept)
			template <typename SchemeInputType, typename DiscretizationInputType, typename IntrinsicValuePolicyType>
			requires std::constructible_from<Scheme, double, double, SchemeInputType, std::shared_ptr<const Discretization>>
			FiniteDifference(double expiry_time, double earliest_pricing_time, size_t num_time_steps, SchemeInputType& scheme_inputs, DiscretizationInputType& discretization_inputs, const IntrinsicValuePolicyType& intrinsic_pricer) : \
				expiry_time_(expiry_time), earliest_pricing_time_(earliest_pricing_time), epsilon_((expiry_time_ - earliest_pricing_time_) / (num_time_steps)), discretization_(std::make_shared<Discretization>(discretization_inputs)){
				if (expiry_time_ <= earliest_pricing_time_) {
					throw std::invalid_argument("The expiry time is the also the last pricing date. The earliest pricing time must be earlier than the last pricing date. The unit of time is years. 1 = 1 year.");
				}
				
				// step 2: reserve space in result
				result_.reserve(num_time_steps+1);
				// step 3: push back initial grid
				result_.push_back(discretization_->IntrinsicOptionValues<IntrinsicValuePolicyType>(intrinsic_pricer));
				// step 4: construct scheme
				Scheme scheme(expiry_time_, epsilon_, scheme_inputs, discretization_);
				// populate result_
				//double current_time = expiry_time_ - epsilon;
				for (int i = 0; i < num_time_steps; ++i) {
					//std::cout << "finished step " << i << "\n";
					result_.push_back(scheme.Next(result_[i]));//, current_time));
					//current_time -= epsilon_; // has one extra subtraction
				}
			}


			// American Finite Difference (see Concept)
			template <typename SchemeInputType, typename DiscretizationInputType, typename IntrinsicValuePolicyType>
			requires std::constructible_from<Scheme, double, double, SchemeInputType, std::shared_ptr<const Discretization>, IntrinsicValuePolicyType>
				FiniteDifference(double expiry_time, double earliest_pricing_time, size_t num_time_steps, SchemeInputType& scheme_inputs, DiscretizationInputType& discretization_inputs, const IntrinsicValuePolicyType& intrinsic_pricer) : \
				expiry_time_(expiry_time), earliest_pricing_time_(earliest_pricing_time), epsilon_((expiry_time_ - earliest_pricing_time_) / (num_time_steps)), discretization_(std::make_shared<Discretization>(discretization_inputs)) {
				if (expiry_time_ <= earliest_pricing_time_) {
					throw std::invalid_argument("The expiry time is the also the last pricing date. The earliest pricing time must be earlier than the last pricing date. The unit of time is years. 1 = 1 year.");
				}
				// step 2: reserve space in result
				result_.reserve(num_time_steps + 1);
				// step 3: push back initial grid
				result_.push_back(discretization_->IntrinsicOptionValues(intrinsic_pricer));
				// step 4: construct scheme
				Scheme scheme(expiry_time_, epsilon_, scheme_inputs, discretization_, intrinsic_pricer);
				// populate result_
				//double current_time = expiry_time_ - epsilon_;
				for (int i = 0; i < num_time_steps; ++i) {
					result_.push_back(scheme.Next(result_[i]));//, current_time));
					//current_time -= epsilon_; // has one extra subtraction
				}
			}

			// in HHW, int asset_price_index, int variance_index, int short_rate_index
			template<typename... Args>
			requires requires (Grid g, Args... args) {
				{g.EntryByIndex(args...)} -> std::convertible_to<double>;
			}
			double PriceByIndex(int time_index, Args... args) {
				return result_[time_index].EntryByIndex(args...);
			}

			double GetExpiryTime() {
				return expiry_time_;
			}

			double GetEarliestPricingTime() {
				return earliest_pricing_time_;
			}

			double GetEpsilon() {
				return epsilon_;
			}

		private:
			double expiry_time_;
			double earliest_pricing_time_;
			double epsilon_;
			std::vector<Grid> result_;
		};

	}
}

#endif