
#include "../Test.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
R U N S  T E S T
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
Runs(Nist::Test& test)
{
	auto epsilon = test.GetEpsilon();
	auto n = test.GetN();
	auto& R_ = test.GetWriteableResults();
	size_t		S, k;
	double	pi, V, erfc_arg, p_value;

	S = 0;
	for (k = 0; k < n; k++)
		if (epsilon[k])
			S++;

	pi = (double)S / (double)n;

	if (fabs(pi - 0.5) >(2.0 / sqrt(n))) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\tPI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
		}
#endif
		p_value = 0.0;

		R_.runs.p_value = p_value;
		R_.runs.pi = pi;
		R_.runs.V = 0.0;
		R_.runs.erfc_arg = 0.0;
	}
	else {

		V = 1;
		for (k = 1; k < n; k++)
			if (epsilon[k] != epsilon[k - 1])
				V++;

		erfc_arg = fabs(V - 2.0 * n * pi * (1 - pi)) / (2.0 * pi * (1 - pi) * sqrt(2 * n));
		p_value = erfc(erfc_arg);
#ifdef SPEED
		dummy_result = p_value;
#endif
		R_.runs.p_value = p_value;
		R_.runs.pi = pi;
		R_.runs.V = V;
		R_.runs.erfc_arg = erfc_arg;

#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\tCOMPUTATIONAL INFORMATION:\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\t(a) Pi                        = %f\n", pi);
			fprintf(stats[TEST_RUNS], "\t\t(b) V_n_obs (Total # of runs) = %d\n", (int)V);
			fprintf(stats[TEST_RUNS], "\t\t(c) V_n_obs - 2 n pi (1-pi)\n");
			fprintf(stats[TEST_RUNS], "\t\t    -----------------------   = %f\n", erfc_arg);
			fprintf(stats[TEST_RUNS], "\t\t      2 sqrt(2n) pi (1-pi)\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			if (isNegative(p_value) || isGreaterThanOne(p_value))
				fprintf(stats[TEST_RUNS], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_RUNS], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_RUNS]);
		}
#endif
	}

#ifdef KS
	pvals.runs_pvals[pvals.seq_counter] = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(results[TEST_RUNS], "%f\n", p_value); fflush(results[TEST_RUNS]);
	}
#endif
}
