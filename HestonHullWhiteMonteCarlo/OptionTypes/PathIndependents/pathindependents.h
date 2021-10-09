//TO DO:
// If too many template specializations are required, then there is a solution to 
// get rid of the variadic template. The solution is to use a proxy class that type 
// erases the option type. PathDependents and PathIndependents should each have one member
// std::vector<Proxy>, with each entry corresponding to what was before arguments to the
// variadic template. Then define the proxy class with something like

// TO DO: make the return type of ComputePrices an array. You have already forced compile
// time computation of the number of prices returned.


//class Proxy {
// public
//    template<class OptionType>
//    Proxy(OptionType& option) :
//        option_pointer_((void*)std::addressof(option)),
//        update_func_pointer_([](void* option_pointer, State state) {
//        ((OptionType*)option_pointer)->Update(state);
//            }),
//        get_prices_func_pointer_([](void* option_pointer) {
//                return ((OptionType*)option_pointer)->GetPrices();
//            })
//    {}
//            void Update(State state) { update_func_pointer_(option_pointer_, state); }
//            std::vector<double> GetPrices() { get_prices_func_pointer_(option_pointer_); }
// private:
//    void* option_pointer_ = 0;
//    void(*update_func_pointer_)(void*) = 0;
//    std::vector<double> (*get_prices_func_pointer_)(void*) = 0;
//};

// with similar function pointers to the constructor and destructor (if need be).
// Unfortunately, this code isn't memory safe. (In theory), this could cause issues
// if you are doing any concurrent programming. I don't see any solution that combines
// std::shared_ptr and type erasure. I suggest implementing the above, adding in the missing
// functions, and just being careful if you need to do anything with concurrency.
//
// The above Proxy class doesn't account for both double and std::vector<double>
// return types of the OptionType class.

#ifndef HESTONHULLWHITEMONTECARLO_OPTIONTYPES_PATHINDEPENDENTS_PATHINDEPENDENTS_H
#define HESTONHULLWHITEMONTECARLO_OPTIONTYPES_PATHINDEPENDENTS_PATHINDEPENDENTS_H

#include <concepts>
#include <tuple>

#include "HestonHullWhiteMonteCarlo/util/util.h"

namespace heston_hull_white {
	namespace monte_carlo {

			template<typename StateType, typename... OptionTypes>
			requires requires (StateType state, OptionTypes... options_) {
				(options_.Price(state), ...);
			}
			class PathIndependents {
			protected:
				std::tuple<OptionTypes...> options_;
			public:

				PathIndependents(OptionTypes&... options) : options_{ options... } {
				}

				std::vector<double> ComputePrices(StateType state) {
					std::vector<double> result{}; // starting size. should be expanded for each Arg that returns a vector of prices (at no additional cost)
					result.reserve(internal::ForceCompileTimeComputation<PathIndependents<StateType, OptionTypes...>::ComputeMinNumberOfPrices<OptionTypes...>()>());
					// grow the vector when return type is a vector.

					//std::apply([&](auto&... ts) {(internal::MoveOver(result, ts.Price()), ...); }, my_objects_);
					std::apply([&](auto&... ts) {(internal::MoveOver(result, ts.Price(state)), ...); }, options_);

					return result;
				}

				std::vector<double> ComputePrices(StateType state, int i) {
					std::vector<double> result{};
					result.reserve(i);
					// no need to grow the vector. If vector becomes too full, that should be an error. This case needs to be the most optimized.
					//std::apply([](auto&... ts) {(internal::MoveOver(result, ts.Price()), ...); }, my_objects_);
					std::apply([](auto&... ts) {(internal::MoveOver(result, ts.Price(state)), ...); }, options_);

					assert(result.size() == i);
					return result;
				}

			private:

				template <typename T>
				requires ((requires (T t, StateType state) { { t.Price(state)} -> std::same_as<double>; }) || (requires (T t, StateType state) { { t.Price(state) } -> std::same_as<std::vector<double>>; }))
					static constexpr std::size_t OptionReturnSize() {
					return 1;
				}

				//requires requires (T t) { { T.GetPrice() } -> internal::is_std_array<>; }
				template <typename T>
				requires requires (T t, StateType state) { { t.Price(state) } -> internal::is_std_array_concept<>; }
				static constexpr std::size_t OptionReturnSize() {
					//return std::tuple_size_v<decltype(std::declval<T&>().Price())>;
					return std::tuple_size_v<decltype(std::declval<T&>().Price(std::declval<StateType&>()))>;
				}

				template <typename... Args> requires (sizeof...(Args) == 0)
					static constexpr std::size_t ComputeMinNumberOfPrices() {
					return 0;
				}

				template<typename FirstArg, typename... RemainingArgs>
				static constexpr std::size_t ComputeMinNumberOfPrices() {
					return OptionReturnSize<FirstArg>() + ComputeMinNumberOfPrices<RemainingArgs...>(); // this type of recursion is slower, but who cares when it is compile time. The code works and I'm not going
					// to bother improving compilation speeds.
				}
			};
	}
}

#endif