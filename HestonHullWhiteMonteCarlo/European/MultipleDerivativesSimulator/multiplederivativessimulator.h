#ifndef HESTONHULLWHITEMONTECARLO_EUROPEAN_MULTIPLEDERIVATIVESSIMULATOR_MULTIPLEDERIVATIVESSIMULATOR_H
#define HESTONHULLWHITEMONTECARLO_EUROPEAN_MULTIPLEDERIVATIVESSIMULATOR_MULTIPLEDERIVATIVESSIMULATOR_H

#include <concepts>
#include <vector>

#include "HestonHullWhiteMonteCarlo/OptionTypes/PathIndependents/pathindependents.h"
#include "HestonHullWhiteMonteCarlo/OptionTypes/PathDependents/pathdependents.h"

namespace heston_hull_white {
	namespace monte_carlo {
		namespace european {


			namespace internal {
				//template <typename... PathIndependentArgs, typename... PathDependentArgs>
				//void RequirePathIndependentPathDependentTypes(PathIndependents<PathIndependentArgs...>, PathDependents<PathDependentArgs...>) {}
				template <typename... Args>
				void AcceptsPathIndependents(PathIndependents<Args...>) {}

				template <typename... Args>
				void AcceptsPathDependents(PathDependents<Args...>) {}

			}


			// TO DO: I may choose to accept different options with different expiry dates and different epsilons for each expiry date.
			// In that case, I need to package PathIndependents_t, PathDependents_t, and epsilon_t for each expiry date t. Then the input
			// would be a collection of such packages.

			// TO DO: replace StateTypeInitializationParameters... with some kind of input parameter struct

			// TO DO: consider the possibility that I will need to call CustomAdvance if there is non-negligible "remaining time" between the last
			// standard subinterval and the maturity date

			// TO DO: Make smarter use of RVO in RunTrial
			template<typename StateType, typename PathIndependentsSpecialization, typename PathDependentsSpecialization, typename StateTypeInitializationParameters>
			requires requires (double epsilon, StateTypeInitializationParameters inp) {
				internal::AcceptsPathIndependents(std::declval<PathIndependentsSpecialization>()); 
				internal::AcceptsPathDependents(std::declval<PathDependentsSpecialization>());
				StateType(epsilon, inp);
				// line for statetype can be initialized with initialization parameters
			}
			class MultipleDerivativesSimulator {
			public:
				// I'm choosing to store by value here. I think this is significantly less expensive than constantly dereferencing. They should be small anyway
				// TODO: consider making num_subintervals an int input
				MultipleDerivativesSimulator(PathIndependentsSpecialization& path_independents, PathDependentsSpecialization& path_dependents, double maturity_date, int num_subintervals, StateTypeInitializationParameters& inputs) : \
					path_independents_{ path_independents }, path_dependents_{ path_dependents }, num_subintervals_(num_subintervals), epsilon_(maturity_date / static_cast<double>(num_subintervals)), starting_state_{ epsilon_, inputs } {
					//internal::RequirePathIndependentPathDependentTypes(path_independents, path_dependents); // verify that these are instantiations of specializations of pathindependents and dependents
				}

				template <std::uniform_random_bit_generator _Engine>
				std::vector<double> RunTrial(_Engine& rng_engine) {
					StateType state = starting_state_;
					//StateType state{input parameters, epsilon}; // it may be faster to do this once and then copy the state each trial. In this case that is unecessary.
					// alternatively I could create a Reset member function. This may mean I should create a StateTypeManager that holds the initial variables.
					// StateTypeManager should contain an instance of the StateType and have a Reset function to reset the StateType.

					//path_dependents.Reset();
					//path_dependents.Initialize(state);
					//path_dependents.ResetAndInitialize(state); // it may be faster to do this once and then copy the reset path_dependents each trial. In this case that is unecessary.
					path_dependents_.ResetAndInitialize();
					for (int i = 0; i < num_subintervals_; ++i) { // there could be issues here if epsilon doesn't exactly divide up the range. I need to thik about that
						state.StandardAdvance(rng_engine); // uses the epsilon from the constructor. Valuable since you need to cache exponentials.
						path_dependents_.Update(state);
					}
					//if (remaining_time_) {
					//	state.CustomAdvance(remaining_time_, rng_engine);
					//}
					//state.CustomAdvance(maturity_date - epsilon * num_steps, rng_engine);
					path_dependents_.Update(state);
					std::vector<double> result = path_independents_.ComputePrices(state);
					std::vector<double> remaining = path_dependents_.ComputePrices(state);
					result.insert(result.end(), remaining.begin(), remaining.end());
					
					// perhaps this step should be integrated into compute important thing here is that I don't need to recompute discount_factor.
					double discount_factor = state.DiscountFactor();
					for (int i = 0; i < result.size(); ++i) {
						double price_without_discount = result[i];
						result[i] = discount_factor * price_without_discount;
					}
					return result;
				}

				//template <std::uniform_random_bit_generator _Engine>
				//std::vector<double> RunTrials(int num_trials, input& parameters, _Engine& rng_engine) {
				//	StateType state{ input parameters, epsilon }; // it may be faster to do this once and then copy the state each trial. In this case that is unecessary.
				//	// alternatively I could create a Reset member function. This may mean I should create a StateTypeManager that holds the initial variables.
				//	// StateTypeManager should contain an instance of the StateType and have a Reset function to reset the StateType.

				//	//path_dependents.Reset();
				//	//path_dependents.Initialize(state);
				//	//path_dependents.ResetAndInitialize(state); // it may be faster to do this once and then copy the reset path_dependents each trial. In this case that is unecessary.
				//	path_dependents.ResetAndInitialize();
				//	for () { // there could be issues here if epsilon doesn't exactly divide up the range. I need to thik about that
				//		state.StandardAdvance(rng_engine); // uses the epsilon from the constructor. Valuable since you need to cache exponentials.
				//		path_dependents.Update(state);
				//	}
				//	state.CustomAdvance(maturity_date - epsilon * num_steps, rng_engine);
				//	path_dependents.Update(state);

				//	return merge(path_independents.Price(), path_dependents.Price());
				//}

			private:
				//PathIndepdents<PathIndependentOptions...> path_independents_;
				PathIndependentsSpecialization path_independents_;
				//PathDependents<PathDependentOptions...> path_dependents_;
				PathDependentsSpecialization path_dependents_;
				const int num_subintervals_;
				const double epsilon_;
				const StateType starting_state_;
				//double epsilon = 0; // this is an invalid epsilon
			};
		}
	}
}

#endif