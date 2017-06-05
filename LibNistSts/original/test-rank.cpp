
#include <math.h>
#include "../Test.h"
#include "rank-matrix.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
R A N K  T E S T
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
Rank(Nist::Test& test)
{
	auto n = test.GetN();
	auto epsilon = test.GetEpsilon();
	auto R_ = test.GetWriteableResults();
	int			N, i, k, r;
	double		p_value, product, chi_squared, arg1, p_32, p_31, p_30, R, F_32, F_31, F_30;
	unsigned char	**matrix = create_matrix(32, 32);

	N = n / (32 * 32);
	if (N == 0) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[TEST_RANK], "\t\t\t\tRANK TEST\n");
			fprintf(stats[TEST_RANK], "\t\tError: Insuffucient # Of Bits To Define An 32x32 (%dx%d) Matrix\n", 32, 32);
		}
#endif
		p_value = 0.00;
	}
	else {
		r = 32;					/* COMPUTE PROBABILITIES */
		product = 1;
		for (i = 0; i <= r - 1; i++)
			product *= ((1.e0 - pow(2, i - 32))*(1.e0 - pow(2, i - 32))) / (1.e0 - pow(2, i - r));
		p_32 = pow(2, r*(32 + 32 - r) - 32 * 32) * product;

		r = 31;
		product = 1;
		for (i = 0; i <= r - 1; i++)
			product *= ((1.e0 - pow(2, i - 32))*(1.e0 - pow(2, i - 32))) / (1.e0 - pow(2, i - r));
		p_31 = pow(2, r*(32 + 32 - r) - 32 * 32) * product;

		p_30 = 1 - (p_32 + p_31);

		F_32 = 0;
		F_31 = 0;
		for (k = 0; k<N; k++) {			/* FOR EACH 32x32 MATRIX   */
			def_matrix(32, 32, matrix, k, epsilon);
#if (DISPLAY_MATRICES == 1)
			display_matrix(32, 32, matrix);
#endif
			R = computeRank(32, 32, matrix);
			if (R == 32)
				F_32++;			/* DETERMINE FREQUENCIES */
			if (R == 31)
				F_31++;
		}
		F_30 = (double)N - (F_32 + F_31);

		chi_squared = (pow(F_32 - N*p_32, 2) / (double)(N*p_32) +
			pow(F_31 - N*p_31, 2) / (double)(N*p_31) +
			pow(F_30 - N*p_30, 2) / (double)(N*p_30));

		arg1 = -chi_squared / 2.e0;

#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[TEST_RANK], "\t\t\t\tRANK TEST\n");
			fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
			fprintf(stats[TEST_RANK], "\t\tCOMPUTATIONAL INFORMATION:\n");
			fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
			fprintf(stats[TEST_RANK], "\t\t(a) Probability P_%d = %f\n", 32, p_32);
			fprintf(stats[TEST_RANK], "\t\t(b)             P_%d = %f\n", 31, p_31);
			fprintf(stats[TEST_RANK], "\t\t(c)             P_%d = %f\n", 30, p_30);
			fprintf(stats[TEST_RANK], "\t\t(d) Frequency   F_%d = %d\n", 32, (int)F_32);
			fprintf(stats[TEST_RANK], "\t\t(e)             F_%d = %d\n", 31, (int)F_31);
			fprintf(stats[TEST_RANK], "\t\t(f)             F_%d = %d\n", 30, (int)F_30);
			fprintf(stats[TEST_RANK], "\t\t(g) # of matrices    = %d\n", N);
			fprintf(stats[TEST_RANK], "\t\t(h) Chi^2            = %f\n", chi_squared);
			fprintf(stats[TEST_RANK], "\t\t(i) NOTE: %d BITS WERE DISCARDED.\n", n % (32 * 32));
			fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
		}
#endif

		p_value = exp(arg1);
#ifdef SPEED
		dummy_result = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			if (isNegative(p_value) || isGreaterThanOne(p_value))
				fprintf(stats[TEST_RANK], "WARNING:  P_VALUE IS OUT OF RANGE.\n");
		}
#endif

		delete_matrix(32, matrix);	/* DEALLOCATE MATRIX  */

		R_.rank.p_30 = p_30;
		R_.rank.p_31 = p_31;
		R_.rank.p_32 = p_32;
		R_.rank.F_30 = F_30;
		R_.rank.F_31 = F_31;
		R_.rank.F_32 = F_32;
		R_.rank.chi_squared = chi_squared;
		R_.rank.p_value = p_value;
		R_.rank.N = N;

	}
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_RANK], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_RANK]);
		fprintf(results[TEST_RANK], "%f\n", p_value); fflush(results[TEST_RANK]);
	}
#endif

#ifdef KS
	pvals.rank_pvals[pvals.seq_counter] = p_value;
#endif
	//printf("%lf %lf %lf \n",F_30,F_31,F_32);
}
