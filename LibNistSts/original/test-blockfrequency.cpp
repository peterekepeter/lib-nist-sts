
#include "../Test.h"
#include "../common/cephes.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
B L O C K  F R E Q U E N C Y  T E S T
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
BlockFrequency(Nist::Test& test)
{
	auto& param = test.GetParameters();
	size_t M = param.blockFrequencyBlockLength;
	auto epsilon = test.GetEpsilon();
	size_t n = test.GetN();
	size_t	i, j, N, blockSum;
	double	p_value, sum, pi, v, chi_squared;

	N = n / M; 		/* # OF SUBSTRING BLOCKS      */
	sum = 0.0;

	for (i = 0; i<N; i++) {
		blockSum = 0;
		for (j = 0; j<M; j++)
			blockSum += epsilon[j + i*M];

		pi = static_cast<double>(blockSum) / static_cast<double>(M);
		v = pi - 0.5;
		sum += v*v;
		//printf("v*v: %lf\n",v*v);
		// added for printage printf("%d ",blockSum);
	}
	chi_squared = 4.0 * M * sum;
	p_value = cephes_igamc(N / 2.0, chi_squared / 2.0);

	auto& result = test.GetWriteableResults().blockfrequency;
	result.chi_squared = chi_squared;
	result.p_value = p_value;

#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.blockfrequency.chi_squared = chi_squared;
	R_.blockfrequency.p_value = p_value;
	if (BlockFrequency_v1 == BlockFrequency) R1 = R_;
	else R2 = R_;
#endif
	//printf("%lf",sum);
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t\tBLOCK FREQUENCY TEST\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(a) Chi^2           = %f\n", chi_squared);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(b) # of substrings = %d\n", N);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(c) block length    = %d\n", M);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(d) Note: %d bits were discarded.\n", n % M);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");

		fprintf(stats[TEST_BLOCK_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BLOCK_FREQUENCY]);
		fprintf(results[TEST_BLOCK_FREQUENCY], "%f\n", p_value); fflush(results[TEST_BLOCK_FREQUENCY]);
	}
#endif
#ifdef KS
	pvals.blockfrequency_pvals[pvals.seq_counter] = p_value;
#endif
}
