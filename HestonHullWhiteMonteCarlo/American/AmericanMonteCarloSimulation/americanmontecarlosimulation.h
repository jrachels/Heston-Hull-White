#ifndef HESTONHULLWHITEMONTECARLO_AMERICAN_AMERICANMONTECARLOSIMULATION_AMERICANMONTECARLOSIMULATION_H
#define HESTONHULLWHITEMONTECARLO_AMERICAN_AMERICANMONTECARLOSIMULATION_AMERICANMONTECARLOSIMULATION_H

// TO DO: 
// check for integer overflow in total_num_trials_+ num_trials

// Currently the sum of squares are on the order of (num_trials)*(variance). Too many trials could overflow Sum of Squares. 
// If too many trials are performed, the algorithm should switch to using division to keep the sum of square estimate on the order of the variance.
// See the first algorithm here: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
// Downside is that division is slow.

// Create a SimulatorConcept that returns a double result rather than a vector (computes the outcome of one random experiment).
// Overload MonteCarloSimulation with that SimulatorConcept to deal with simulating a single random experiment. This is easier
// for the user and should provide a small speed boost.

// Need to decide if Sum of Square differences needs to be convered to Estimate_Variances in AddTrials function or in EstimateVariances
// function. Depends on how this is used. I chose EstimateVariances because I believe that should only be called once and the results saved.

// since the simulator returns an std::array, I am storing the size twice. I am storing it once in the std::array and a second time in the simulator.

// Either the simulator needs to be stateless (concept condition) and stored as a shared pointer in MonteCarloSimulation OR it can have state but needs to be stored as 
// a unique pointer? This is the in

// Use Control Variates. Compute price of vanilla European options with Fourier Transform method. Compute covariance with exotic options in using online covariance algorithm 
// (put in AddTrials). Use Control Variates to reduce covariance.

// Use Antithetic Variates

#include <concepts>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

//#include "HestonHullWhiteMonteCarlo/European/MonteCarloSimulation/src/combinedsimulatorconcept.h"

#include "HestonHullWhiteMonteCarlo/util/Matrices/matrix.h"
#include "HestonHullWhiteMonteCarlo/util/Matrices/matrixconcept.h"

namespace heston_hull_white {
	namespace monte_carlo {
		namespace american {

			template<typename AmericanCombSim, util::MatrixConcept MatrixType = util::Matrix>
			requires requires (AmericanCombSim a, MatrixType m, int n, std::mt19937 b) { //I can't use concepts in concepts, so my workaround is to us std::mt19937 specifically
				{a.RunTrial(m, n, b)} -> std::convertible_to<std::vector<double>>;
				{a.num_timesteps} -> std::convertible_to<int>;
			}
			class AmericanMonteCarloSimulation
			{
				//private:
					// I can't assume it is an array because the inputs may be decided at compile time
					//using MonteArray = std::array<double, std::tuple_size_v<decltype(std::declval<CombSim&>().RunTrial(std::declval<std::mt19937>())) >> ;
			public:
				template<std::uniform_random_bit_generator _Engine>//, CombinedSimulatorConcept CombSim>
				AmericanMonteCarloSimulation(_Engine& _Eng, std::shared_ptr<AmericanCombSim> CombSimulator, int num_trials) : total_num_trials_(0), simulator_{ std::move(CombSimulator) } {//, average_estimates{ simulator->NumRandomVariables() }, sum_of_squared_differences_from_mean{ simulator->NumRandomVariables() }{
					if (num_trials < 3) {
						throw std::invalid_argument(std::string("Error, attempted to run MonteCarloSimulation with less than 3 trials. Number of trials attempted: ") + std::to_string(num_trials) + \
							std::string(". Recommended number of trials is at least several thousand."));
					}
					//AddTrials(_Eng, simulator_, num_trials);
					AddTrials(_Eng, num_trials);
				}

				AmericanMonteCarloSimulation(int num_trials) : MonteCarloSimulation<std::mt19937>(std::mt19937{ std::random_device()() }, num_trials) {}

			private:
				// AddTrails should be called only once, so it is private

				//template <std::uniform_random_bit_generator _Engine, CombinedSimulatorConcept CombSim>
				template <std::uniform_random_bit_generator _Engine>
				void AddTrials(_Engine& _Eng, int num_trials) {
					// step 1: declare matrix of paths WAIT I NEED A VECTOR OF MATRIXES OF PATHS
					int num_options = CombSimulator->num_options;
					// Note: the traditional LSM approach discounts later
					std::vector<Matrix> matrices_of_discounted_option_prices = { Matrix(rows = num_paths, columns = CombSimulator->num_timesteps), CombSimulator->num_options }



					// initialize each entry as  Matrix(rows = num_paths, columns = CombSimulator->num_timesteps) 
						// STORE MATRIX AS A VECTOR OF ROWS?

					for (int trial_number = 0; trial_number < num_trials; ++trial_number) {
						simulator_->RunTrial(matrices_of_discounted_option_prices, trial_number, _Eng);
					}


					Matrix basis_vectors(row = num_trials, columns = num_basis_vectors)

					// step 2: SHOULD I DO EACH OPTION 1 AT A TIME, OR EACH TRAIL 1 AT A TIME?
					for (int step_number = num_steps; step_number > 0; --step_number) {
						// step i: fill matrix with 
						for (int option_number = 0; option_nmber < num_options, ++option_number) {
							// Learn how iterators work from that 
							laguerre()
							basis_vectors(trial_number)
							matrices_of_discounted_option_prices[option_number]
						}
					}







					//constexpr int num_random_variables = std::tuple_size_v<decltype(std::declval<CombSim&>.RunTrial()())>;
					int new_total = total_num_trials_ + num_trials;
					while (total_num_trials_ < new_total) {
						// warning: didn't check for integer overflow in total_num_trials_+ num_trials
						// iterate total_num_trials_
						total_num_trials_++;
						double total_num_trials_fl = static_cast<double>(total_num_trials_);
						// run simulator and get result
						//MonteArray simulation_result = simulator(_Eng);
						std::vector<double> simulation_result = simulator_->RunTrial(_Eng);
						// update mean and variances estimate	
						// combine everything into one loop
						const std::size_t num_random_variables = simulation_result.size();
						mean_outcomes_.resize(num_random_variables, 0);
						sum_of_squared_differences_from_mean_.resize(num_random_variables, 0);
						for (std::size_t i = 0; i < num_random_variables; ++i) {
							double result = simulation_result[i];
							double delta1 = result - mean_outcomes_[i];
							// update mean
							mean_outcomes_[i] += delta1 / total_num_trials_fl;
							double delta2 = result - mean_outcomes_[i];
							// update sd
							// Warning: too many trials could overflow Sum of Squares estimate. If too many trials are performed, the algorithm should switch
							// to using division. That is currently not implemented.

							sum_of_squared_differences_from_mean_[i] += delta1 * delta2;
						}

						//std::cout << (new_total - total_num_trials_) << " trials remaining.\n";
					}

				}


				////template <std::uniform_random_bit_generator _Engine, CombinedSimulatorConcept CombSim>
				//template <std::uniform_random_bit_generator _Engine>
				//void AddTrials(_Engine& _Eng, CombSim& simulator, int num_trials) {
				//	constexpr int num_random_variables = std::tuple_size_v<decltype(std::declval<CombSim&>()())>;

				//	while (total_num_trials_< total_num_trials_+ num_trials) {
				//		// warning: didn't check for integer overflow in total_num_trials_+ num_trials
				//		// iterate total_num_trials_
				//		total_num_trials_++;
				//		double total_num_trials_fl = static_cast<double>(total_num_trials_);
				//		// run simulator and get result
				//		//MonteArray simulation_result = simulator(_Eng);
				//		std::vector<double> simulation_result = simulator(_Eng);
				//		// update mean and variances estimate	
				//		// combine everything into one loop
				//		for (int i = 0; i < num_random_variables; ++i) {
				//			double result = simulation_result[i];
				//			double delta1 = result - mean_outcomes_[i];
				//			// update mean
				//			mean_outcomes_[i] += delta1 / total_num_trials_fl;
				//			double delta2 = result - mean_outcomes_[i];
				//			// update sd
				//			// Warning: too many trials could overflow Sum of Squares estimate. If too many trials are performed, the algorithm should switch
				//			// to using division. That is currently not implemented.

				//			sum_of_squared_differences_from_mean_[i] += delta1*delta2;
				//		}


				//	}

				//}

				//MonteArray MeanOutcomes() {
				std::vector<double> MeanOutcomes() {
					return mean_outcomes_;
				}
				//MonteArray EstimateVariances() {
				std::vector<double> EstimateVariances() {
					// This uses RVO, right?
					int size = sum_of_squared_differences_from_mean_.size();
					//MonteArray estimate_variances(size);
					std::vector<double> estimate_variances(size);
					double total_num_trials_fl = static_cast<double>(total_num_trials_);
					for (int i = 0; i < size; ++i) {
						estimate_variances[i] = sum_of_squared_differences_from_mean_[i] / total_num_trials_fl;
					}
					//return sum_of_squared_differences_from_mean_;
					return estimate_variances;
				}
				int NumTrials() {
					return total_num_trials_;
				}
				constexpr int NumRandomVariables() {
					return std::tuple_size_v<decltype(std::declval<CombSim&>()())>;
				}

				int total_num_trials_ = 0;
				// storing as shared_ptr requires simulator to be stateless
				std::shared_ptr<AmericanCombSim> simulator_;
				//MonteArray mean_outcomes;
				std::vector<double> mean_outcomes_;
				//MonteArray sum_of_squared_differences_from_mean_;
				std::vector<double> sum_of_squared_differences_from_mean_;
			};
		}
	}
}

#endif 