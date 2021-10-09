//TO DO:
// If too many template specializations are required, then there is a solution to 
// get rid of the variadic template. The solution is to use a proxy class that type 
// erases the option type. PathDependents and PathIndependents should each have one member
// std::vector<Proxy>, with each entry corresponding to what was before arguments to the
// variadic template. Then define the proxy class with something like

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


// There could be an issue if the same option types are used with different
// return lengths for GetPrice(). I assume the return type is either a double or a fixed length
// std::array. The programmer may want to adapt the program to deal with std::vectors, but 
// then you can't determine the return length of GetPrices at runtime (which may require
// repeated resizing of a vector). Fixing this requires taking a close look at use cases
// of this algorithm.

#ifndef HESTONHULLWHITEMONTECARLO_OPTIONTYPES_PATHDEPENDENTS_PATHDEPENDENTS_H
#define HESTONHULLWHITEMONTECARLO_OPTIONTYPES_PATHDEPENDENTS_PATHDEPENDENTS_H

#include <concepts>
#include <tuple>

#include "HestonHullWhiteMonteCarlo/OptionTypes/PathIndependents/pathindependents.h"

namespace heston_hull_white {
	namespace monte_carlo {
			template<typename StateType, typename... OptionTypes>
			requires requires
				(StateType state, OptionTypes... options) {
				(options.Update(state), ...);
				(options.ResetAndInitialize(), ...);
			}
			class PathDependents : public PathIndependents<StateType, OptionTypes...> {
			public:

				PathDependents(OptionTypes&... options) : PathIndependents<StateType, OptionTypes...>{ options... } {
				}

				// idk how this works
				void Update(StateType state) {
					//std::apply([&](auto&... ts) {(ts.Update(state), ...); }, my_objects_);
					std::apply([&](auto&... ts) {(ts.Update(state), ...); }, PathIndependents<StateType, OptionTypes...>::options_);
				}

				void ResetAndInitialize() {
					//std::apply([&](auto&... ts) {(ts.ResetAndInitialize(), ...); }, my_objects_);
					std::apply([&](auto&... ts) {(ts.ResetAndInitialize(), ...); }, PathIndependents<StateType, OptionTypes...>::options_);
				}

			private:

			};
	}
}

#endif