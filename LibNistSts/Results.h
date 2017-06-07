#pragma once

#include <vector>
#include <array>

namespace Nist
{
	// structure to hold the results
	struct Results {
		struct Frequency { double sum, sum_n, p_value; } frequency;
		struct BlockFrequency { double chi_squared, p_value; } blockfrequency;
		struct Runs { double pi, V, erfc_arg, p_value; } runs;
		struct LongestRunOfOnes { int N, M; double chi2, p_value; std::array<unsigned int, 7> nu; } longestrunofones;
		struct Rank { double p_30, p_31, p_32, F_30, F_31, F_32, N, chi_squared, p_value; } rank;
		struct Serial { double	p_value1, p_value2, psim0, psim1, psim2, del1, del2; } serial;
		struct NonOverlapping { unsigned int templates; std::vector<unsigned int> W; std::vector<double> chi2; std::vector<double> p_value; } nonoverlapping;
		struct Overlapping { std::array<unsigned int, 6> nu; double chi2, p_value; } overlapping;
		struct Universal { double p_value, sum, phi; } universal;
		struct ApproximateEntropy { std::array<double, 2> ApEn; double chi_squared, p_value; std::vector<unsigned int> P; unsigned int pp; } approximate_entropy;
		struct Cusum { int z, zrev; double sum1A, sum2A, sum1B, sum2B, p_valueA, p_valueB; } cusum;
		struct RandomExcursion { std::array<double, 8> p_value, sum; bool valid; std::array<int, 8> x, J; } random_excursion;
		struct RandomExcursionVariant { std::array<double, 18> p_value; bool valid; std::array<int, 18> x, count; } random_excursion_variant;
		struct LinearComplexity { double p_value, chi2; std::array<int, 7> nu; } linear_complexity;
		struct Dft { double	p_value, percentile, N_l, N_o, d; } dft;
	};
}