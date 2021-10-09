#ifndef HESTONHULLWHITEMONTECARLO_UTIL_UTIL_H
#define HESTONHULLWHITEMONTECARLO_UTIL_UTIL_H

#include <array>
#include <vector>

namespace internal {
	template<class T>
	struct is_std_array : std::false_type {};
	// I think N below is deduced? If not, I need to figure out how to make it deduced
	template<class T, std::size_t N>
	struct is_std_array<std::array<T, N>> : std::true_type {};
	template<class T>
	struct is_std_array<T const> : is_std_array<T> {};
	template<class T>
	struct is_std_array<T volatile> : is_std_array<T> {};
	template<class T>
	struct is_std_array<T volatile const> : is_std_array<T> {};

	template<typename T>
	concept is_std_array_concept = is_std_array<T>::value;

	template<class T>
	struct is_std_vector_t : std::false_type {};
	template<class T>
	struct is_std_vector_t<std::vector<T>> : std::true_type {};
	template<class T>
	struct is_std_vector_t<T const> : is_std_vector_t<T> {};
	template<class T>
	struct is_std_vector_t<T volatile> : is_std_vector_t<T> {};
	template<class T>
	struct is_std_vector_t<T volatile const> : is_std_vector_t<T> {};

	template<typename T>
	concept is_std_vector_concept = is_std_vector_t<T>::value;

	template<std::size_t N>
	std::size_t ForceCompileTimeComputation() {
		return N;
	}

	//template<typename T, std::enable_if_t<std::is_same<decltype(std::declval<T&>().GetPrice()), double>::value, bool> = true>
//void MoveOver(std::vector<double>& ls, T t) {
//	ls.push_back(t);
//}
	inline void MoveOver(std::vector<double>& ls, double t) {
		ls.push_back(t);
	}

	//template <typename T, std::enable_if_t<internal::is_std_array<decltype(std::declval<T&>().GetPrice())>::value || internal::is_std_vector_t<decltype(std::declval<T&>().GetPrice())>::value, bool> = true>

	//void MoveOver(std::vector<double>& ls, T t) {
	//	ls.insert(ls.end(), t.begin(), t.end());
	//}

	template <typename T>
	requires requires (T t) {
		{t.begin() } -> std::incrementable<>;
		{t.end() } -> std::incrementable<>;
		//{t.begin() };
		//{t.end() };
	}
	inline void MoveOver(std::vector<double>& ls, T t) {
		ls.insert(ls.end(), t.begin(), t.end());
	}

}

#endif