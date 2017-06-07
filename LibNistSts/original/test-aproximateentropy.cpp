
#include "../Test.h"
#include "../original/cephes.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
A P P R O X I M A T E  E N T R O P Y   T E S T
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
ApproximateEntropy(Nist::Test& test)
{
	auto epsilon = test.GetEpsilon();
	auto n = test.GetN();
	auto& R_ = test.GetWriteableResults();
	auto m = test.GetParameters().approximateEntropyBlockLength;

	int				i, j, k, r, blockSize, seqLength, powLen, index;
	unsigned int cc = 0;

	double			sum, numOfBlocks, ApEn[2], apen, chi_squared, p_value;
	unsigned int	*P;

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_APEN], "\t\t\tAPPROXIMATE ENTROPY TEST\n");
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_APEN], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_APEN], "\t\t(a) m (block length)    = %d\n", m);
	}
#endif

	seqLength = n;
	r = 0;

	R_.approximate_entropy.P.resize(2 * (((size_t)1) << (m + 1)));

	for (blockSize = m; blockSize <= m + 1; blockSize++) {
		if (blockSize == 0) {
			ApEn[0] = 0.00;
			r++;
		}
		else {
			numOfBlocks = (double)seqLength;
			powLen = (int)pow(2, blockSize + 1) - 1;
			if ((P = (unsigned int*)calloc(powLen, sizeof(unsigned int))) == NULL) {
#if defined(FILE_OUTPUT) ||  defined(KS)
				if (cmdFlags.output == 1 || cmdFlags.output == -1) {
					fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
				}
#endif
				printf("ApEn:  Insufficient memory available.\n");
				return;
			}
			for (i = 1; i<powLen - 1; i++)
				P[i] = 0;
			for (i = 0; i<numOfBlocks; i++) { /* COMPUTE FREQUENCY */
											  /*epsilon[0] = 1;
											  epsilon[1] = 1;
											  epsilon[2] = 1;*/
				k = 1;
				for (j = 0; j<blockSize; j++) {
					k <<= 1;
					if ((int)epsilon[(i + j) % seqLength] == 1)
						k++;
				}
				P[k - 1]++;
				//if(i < 100)printf(" %i ",k-1);
			}
			/* DISPLAY FREQUENCY */
			sum = 0.0;
			index = (int)pow(2, blockSize) - 1;
			for (i = 0; i<(int)pow(2, blockSize); i++) {
				if (P[index] > 0)
					sum += P[index] * log(P[index] / numOfBlocks);
				//printf("[%i: %d] ",index,P[index]);

				R_.approximate_entropy.P[cc++] = P[index];
				index++;
			}
			R_.approximate_entropy.pp = cc;

			sum /= numOfBlocks;
			ApEn[r] = sum;

			//printf("\n");
			//printf("SUM: %lf \n\n",sum);
			r++;
			free(P);
		}
	}
	apen = ApEn[0] - ApEn[1];

	chi_squared = 2.0*seqLength*(log(2) - apen);
	p_value = cephes_igamc(pow(2, m - 1), chi_squared / 2.0);

#ifdef SPEED
	dummy_result = p_value;
#endif
	R_.approximate_entropy.p_value = p_value;
	R_.approximate_entropy.chi_squared = chi_squared;
	R_.approximate_entropy.ApEn[0] = ApEn[0];
	R_.approximate_entropy.ApEn[1] = ApEn[1];

#ifdef KS
	pvals.approximate_entropy_pvals[pvals.seq_counter] = p_value;
#endif

	//printf("P-value %lf \n",p_value);
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_APEN], "\t\t(b) n (sequence length) = %d\n", seqLength);
		fprintf(stats[TEST_APEN], "\t\t(c) Chi^2               = %f\n", chi_squared);
		fprintf(stats[TEST_APEN], "\t\t(d) Phi(m)	       = %f\n", ApEn[0]);
		fprintf(stats[TEST_APEN], "\t\t(e) Phi(m+1)	       = %f\n", ApEn[1]);
		fprintf(stats[TEST_APEN], "\t\t(f) ApEn                = %f\n", apen);
		fprintf(stats[TEST_APEN], "\t\t(g) Log(2)              = %f\n", log(2.0));
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	}
#endif

	if (m > (int)(log(seqLength) / log(2) - 5)) {
		//printf("\t\tNote: The blockSize exceeds recommended value\n");
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[TEST_APEN], "\t\tNote: The blockSize = %d exceeds recommended value of %d\n", m,
				MAX(1, (int)(log(seqLength) / log(2) - 5)));
			fprintf(stats[TEST_APEN], "\t\tResults are inaccurate!\n");
			fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
		}
#endif
	}

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_APEN], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_APEN]);
		fprintf(results[TEST_APEN], "%f\n", p_value); fflush(results[TEST_APEN]);
	}
#endif
}
