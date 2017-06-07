
#include "../Test.h"
#include "../common/cephes.h"

//100MB - M=1000 927sec
//		  M=5000 4500sec
void
LinearComplexity(Nist::Test& test)
{
	auto epsilon = test.GetEpsilon();
	auto n = test.GetN();
	auto& R_ = test.GetWriteableResults();
	auto M = test.GetParameters().linearComplexitySequenceLength;

	using BitSequence = unsigned char;

	int       i, ii, j, d, N, L, m, N_, parity, sign, K = 6, nu[7];
	double    p_value, T_, mean, chi2;
	double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
	BitSequence  *T = nullptr, *P = nullptr, *B_ = nullptr, *C = nullptr;

	N = (int)floor(n / M);
	if (((B_ = (BitSequence *)calloc(M, sizeof(BitSequence))) == nullptr) ||
		((C = (BitSequence *)calloc(M, sizeof(BitSequence))) == nullptr) ||
		((P = (BitSequence *)calloc(M, sizeof(BitSequence))) == nullptr) ||
		((T = (BitSequence *)calloc(M, sizeof(BitSequence))) == nullptr)) {
		printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		if (B_ != nullptr)
			free(B_);
		if (C != nullptr)
			free(C);
		if (P != nullptr)
			free(P);
		if (T != nullptr)
			free(T);
		return;
	}

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tL I N E A R  C O M P L E X I T Y\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tM (substring length)     = %d\n", M);
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tN (number of substrings) = %d\n", N);
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "        F R E Q U E N C Y                            \n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "  C0   C1   C2   C3   C4   C5   C6    CHI2    P-value\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tNote: %d bits were discarded!\n", n%M);
	}
#endif

	for (i = 0; i < K + 1; i++)
		nu[i] = 0;
	for (ii = 0; ii < N; ii++) {
		for (i = 0; i < M; i++) {
			B_[i] = 0;
			C[i] = 0;
			T[i] = 0;
			P[i] = 0;
		}
		L = 0;
		m = -1;
		d = 0;
		C[0] = 1;
		B_[0] = 1;

		/* DETERMINE LINEAR COMPLEXITY */
		N_ = 0;
		while (N_ < M) {
			d = (int)epsilon[ii*M + N_];
			for (i = 1; i <= L; i++)
				d += C[i] * epsilon[ii*M + N_ - i];


			d = d % 2;
			if (d == 1) {
				for (i = 0; i < M; i++) {
					T[i] = C[i];
					P[i] = 0;
				}
				for (j = 0; j < M; j++)
					if (B_[j] == 1)
						P[j + N_ - m] = 1;
				for (i = 0; i < M; i++)
					C[i] = (C[i] + P[i]) % 2;
				if (L <= N_ / 2) {
					L = N_ + 1 - L;
					m = N_;
					for (i = 0; i < M; i++)
						B_[i] = T[i];
				}
			}
			N_++;
		}
		//printf("L = %d \n",L);

		if ((parity = (M + 1) % 2) == 0)
			sign = -1;
		else
			sign = 1;
		mean = M / 2.0 + (9.0 + sign) / 36.0 - 1.0 / pow(2, M) * (M / 3.0 + 2.0 / 9.0);
		if ((parity = M % 2) == 0)
			sign = 1;
		else
			sign = -1;
		T_ = sign * (L - mean) + 2.0 / 9.0;

		if (T_ <= -2.5)
			nu[0]++;
		else if (T_ > -2.5 && T_ <= -1.5)
			nu[1]++;
		else if (T_ > -1.5 && T_ <= -0.5)
			nu[2]++;
		else if (T_ > -0.5 && T_ <= 0.5)
			nu[3]++;
		else if (T_ > 0.5 && T_ <= 1.5)
			nu[4]++;
		else if (T_ > 1.5 && T_ <= 2.5)
			nu[5]++;
		else
			nu[6]++;
	}
	chi2 = 0.00;
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		for (i = 0; i<K + 1; i++)
			fprintf(stats[TEST_LINEARCOMPLEXITY], "%4d ", (int)nu[i]);
	}
#endif
	for (i = 0; i<K + 1; i++)
	{
		chi2 += pow(nu[i] - N*pi[i], 2) / (N*pi[i]);
		//printf("%d ",nu[i]);
		R_.linear_complexity.nu[i] = nu[i];
	}

	p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
#ifdef SPEED
	dummy_result = p_value;
#endif

	R_.linear_complexity.chi2 = chi2;
	R_.linear_complexity.p_value = p_value;

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_LINEARCOMPLEXITY], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_LINEARCOMPLEXITY]);
		fprintf(results[TEST_LINEARCOMPLEXITY], "%f\n", p_value); fflush(results[TEST_LINEARCOMPLEXITY]);
	}
#endif

#ifdef KS
	pvals.linear_complexity_pvals[pvals.seq_counter] = p_value;
#endif
	free(B_);
	free(P);
	free(C);
	free(T);
}
