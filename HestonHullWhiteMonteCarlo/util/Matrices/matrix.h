#ifndef HESTONHULLWHITEMONTECARLO_UTIL_MATRICES_MATRIX_H
#define HESTONHULLWHITEMONTECARLO_UTIL_MATRICES_MATRIX_H

#include <vector>

// TO DO: make header only?

namespace heston_hull_white {
	namespace util {
		class Matrix
		{
		public:
			explicit Matrix(int m, int n) : height_(m), width_(n), mat_(m* n) { //returns matrix of all 0s
				mat_.shrink_to_fit();
			}

			Matrix() = default;

			double& operator() (int i, int j) {
				// stored as [column1, column2, column3,..., columnn] for cache reasons.
				// remember to copy column first!
				return mat_[i + j * height_];
			}

			double operator() (int i, int j) const {
				// stored as [column1, column2, column3,..., columnn] for cache reasons.
				// remember to copy column first!
				return mat_[i + j * height_];
			}

			void ScaleBy(double scale);

			int rows() const {
				return height_;
			}

			int columns() const {
				return width_;
			}

			double MaxAbsOfEntries() const;

			Matrix Transpose() const;

		private:

			// forced const
			int height_;
			int width_;
			// changeable
			std::vector<double> mat_;
		};

	}
}
#endif
