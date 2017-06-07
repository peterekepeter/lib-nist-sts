#pragma once

#include <vector>

namespace Nist
{
	// structure to hold the results
	struct Results {
		struct Frequency { double sum, sum_n, p_value; } frequency;
		struct BlockFrequency { double chi_squared, p_value; } blockfrequency;
		struct Runs { double pi, V, erfc_arg, p_value; } runs;
		struct LongestRunOfOnes { int N, M; double chi2, p_value; unsigned int	nu[7]; } longestrunofones;
		struct Rank { double p_30, p_31, p_32, F_30, F_31, F_32, N, chi_squared, p_value; } rank;
		struct Serial { double	p_value1, p_value2, psim0, psim1, psim2, del1, del2; } serial;
		struct NonOverlapping { unsigned int templates; std::vector<unsigned int> W; std::vector<double> chi2; std::vector<double> p_value; } nonoverlapping;
		struct Overlapping { unsigned int nu[6]; double chi2, p_value; } overlapping;
		struct Universal { double p_value, sum, phi; } universal;
		struct ApproximateEntropy { double ApEn[2], chi_squared, p_value; std::vector<unsigned int> P; unsigned int pp; } approximate_entropy;
		struct Cusum { int z, zrev; double sum1A, sum2A, sum1B, sum2B, p_valueA, p_valueB; } cusum;
		struct RandomExcursion { double p_value[8], sum[8]; int valid, x[8], J[8]; } random_excursion;
		struct RandomExcursionVariant { double p_value[18]; int valid, x[18], count[18]; } random_excursion_variant;
		struct LinearComplexity { double p_value, chi2; int nu[7]; } linear_complexity;
		struct Dft { double	p_value, percentile, N_l, N_o, d; } dft;
	};
}