#ifndef HESTONHULLWHITEFINITEDIFFERENCE_GRIDEXAMPLES_THREEDIMENSIONALGRID_THREEDIMENSIONALGRID_H
#define HESTONHULLWHITEFINITEDIFFERENCE_GRIDEXAMPLES_THREEDIMENSIONALGRID_THREEDIMENSIONALGRID_H

// The TheDimensionalGrid class is essentially an implementation of a three dimensional matrix.

#include <array>
#include <cassert>
#include <vector>

namespace heston_hull_white {
	namespace finite_difference {
		namespace examples {

			template<typename EntryType>
			class ThreeDimensionalGrid {
			public:
				ThreeDimensionalGrid(std::size_t num_ss, std::size_t num_vs, std::size_t num_rs) : num_ss_(num_ss), num_vs_(num_vs), num_rs_(num_rs), grid_(num_ss*num_vs*num_rs) {}

				EntryType& EntryByIndex(std::size_t s_index, std::size_t v_index, std::size_t r_index) {
					assert((0 <= s_index) && (s_index < num_ss_));
					assert((0 <= v_index) && (v_index < num_vs_));
					//assert(0 <= r_index < num_rs_);
					return grid_[num_vs_ * num_ss_ * r_index + num_vs_ * s_index + v_index];
				}

				EntryType EntryByIndex(std::array<std::size_t, 3> indices) const {
					assert((0 <= indices[0]) && (indices[0] < num_ss_));
					assert((0 <= indices[1]) && (indices[1] < num_vs_));
					//assert(0 <= r_index < num_rs_);
					return grid_[num_vs_ * (num_ss_ * indices[2] + indices[0]) + indices[1]];
				}

				EntryType EntryByIndex(std::size_t s_index, std::size_t v_index, std::size_t r_index) const {
					assert((0 <= s_index) && (s_index < num_ss_));
					assert((0 <= v_index) && (v_index < num_vs_));
					//assert(0 <= r_index < num_rs_);
					return grid_[num_vs_ * (num_ss_ * r_index + s_index) + v_index];
				}

				EntryType& RefToEntryByVectorIndex(std::size_t i) {
					return grid_[i];
				}

				EntryType EntryByVectorIndex(std::size_t i) const {
					return grid_[i];
				}

				std::size_t Size() const {
					return grid_.size();
				}

				// this is a template function because there is no multiplication defined for some Entry types
				template<typename ScalarType>
				void MultiplyByConstant(ScalarType scalar) {
					size_t grid_size = grid_.size();
					for (size_t i = 0; i < grid_size; ++i) {
						grid_[i] *= scalar;
					}
				}

				void Add(const ThreeDimensionalGrid<EntryType>& rhs) {
					std::size_t size = grid_.size();
					assert(size == rhs.Size());
					for (std::size_t i = 0; i < size; ++i) {
						grid_[i] += rhs.grid_[i];
					}
				}

				void Subtract(const ThreeDimensionalGrid<EntryType>& rhs) {
					std::size_t size = grid_.size();
					assert(size == rhs.Size());
					for (std::size_t i = 0; i < size; ++i) {
						grid_[i] -= rhs.grid_[i];
					}
				}

				std::size_t NumS() const {
					return num_ss_;
				}
				std::size_t NumV() const {
					return num_vs_;
				}
				std::size_t NumR() const {
					return num_rs_;
				}

			private:
				std::size_t num_ss_;
				std::size_t num_vs_;
				std::size_t num_rs_;
				std::vector<EntryType> grid_;
			};

		}
	}
}

#endif