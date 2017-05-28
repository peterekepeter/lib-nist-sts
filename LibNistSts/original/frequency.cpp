#include <stdio.h>
#include <math.h>
#include "../Test.h"


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
F R E Q U E N C Y  T E S T
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
Frequency(Nist::Test& test)
{
	auto n = test.GetN();
	auto epsilon = test.GetEpsilon();
	double	f, s_obs, p_value, sum, sqrt2 = 1.41421356237309504880;

	sum = 0.0;
	for (size_t i = 0; i<n; i++)
		sum += 2 * static_cast<int>(epsilon[i]) - 1;
	//sum += epsilon[i];
	s_obs = fabs(sum) / sqrt(n);
	f = s_obs / sqrt2;
	p_value = erfc(f);
	//printf("Pval: %lf sum %lf \n", p_value, sum);

	//store results
	auto& frequency = test.GetWriteableResults().frequency;
	frequency.sum = sum;
	frequency.sum_n = sum / n;
	frequency.p_value = p_value;

#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.frequency.sum = sum;
	R_.frequency.sum_n = sum / n;
	R_.frequency.p_value = p_value;
	if (Frequency_v1 == Frequency) R1 = R_;
	else R2 = R_;
#endif


#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_FREQUENCY], "\t\t\t      FREQUENCY TEST\n");
		fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_FREQUENCY], "\t\t(a) The nth partial sum = %d\n", (int)sum);
		fprintf(stats[TEST_FREQUENCY], "\t\t(b) S_n/n               = %f\n", sum / n);
		fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");

		fprintf(stats[TEST_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_FREQUENCY]);
		fprintf(results[TEST_FREQUENCY], "%f\n", p_value); fflush(results[TEST_FREQUENCY]);
	}
#endif
#ifdef KS
	pvals.frequency_pvals[pvals.seq_counter] = p_value;
#endif
}
