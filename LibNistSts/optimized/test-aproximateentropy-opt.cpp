
#include "../Test.h"
#include "../common/cephes.h"
#include "tools.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
A P P R O X I M A T E  E N T R O P Y   T E S T
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* --------------------------------------------------------------------------

The following code is distributed under the following BSD-style license:

Copyright © 2013-2014 Marek Sys (syso@fi.muni.cz) & Zdenek Riha (zriha@fi.muni.cz).
All Rights Reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.

3. The name of the author may not be used to endorse or promote products derived
from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY AUTHORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------------- */


void
ApproximateEntropy2(Nist::Test& test)
{
	auto array = test.GetArray();
	auto n = test.GetN();
	auto& R_ = test.GetWriteableResults();
	auto m = test.GetParameters().approximateEntropyBlockLength;

	int				i, k, len;
	double			sum, numOfBlocks, ApEn[2], apen, chi_squared, p_value;
	unsigned int	*P, mask, help;

	unsigned int cc = 0;

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_APEN], "\t\t\tAPPROXIMATE ENTROPY TEST\n");
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_APEN], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_APEN], "\t\t(a) m (block length)    = %d\n", m);
	}
#endif


	numOfBlocks = n;
	m++;
	mask = (1 << m) - 1;

	len = (1 << m);

	R_.approximate_entropy.P.resize(len * 2);

	if ((P = (unsigned int*)calloc(len, sizeof(unsigned int))) == NULL) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
		}
#endif
		return;
	}
	for (i = 0; i < len; i++)
		P[i] = 0;
	for (i = 0; i < n - m + 1; i++) {		 /* COMPUTE FREQUENCY */
		++P[get_nth_block4(array, i)&mask];
		//if(i < 100)printf(" %i ",(1 << m)- 1 + Mirrored_int((get_nth_block4(array,i)&mask),m));
	}
	for (i = 1; i<m; i++) {
		k = get_nth_block4(array, n - m + i)&(mask >> i);
		//printf("%d ",k);
		k ^= (((unsigned int*)array)[0] << (m - i));
		//printf("%d ",k);
		k &= mask;
		//printf("%d ",k);
		P[k]++;
	}

	//DISPLAY FREQUENCY
	sum = 0.0;
	for (i = 0; i < len / 2; i++) {
		help = P[Mirrored_int(i, m - 1)] + P[Mirrored_int(i, m - 1) + len / 2];
		if (help > 0)
			sum += help*log(help / numOfBlocks);
		R_.approximate_entropy.P[cc++] = help;
		//printf("%i ",help);
	}


	sum /= numOfBlocks;
	ApEn[0] = sum;

	sum = 0.0;
	for (i = 0; i < len; i++) {
		if (P[Mirrored_int(i, m)] > 0)
			sum += P[Mirrored_int(i, m)] * log(P[Mirrored_int(i, m)] / numOfBlocks);
		//printf("[%d: %d] ",i,P[Mirrored_int(i,m)]);
		R_.approximate_entropy.P[cc++] = P[Mirrored_int(i, m)];
	}
	R_.approximate_entropy.pp = cc;

	sum /= numOfBlocks;
	ApEn[1] = sum;



	//printf("\n");
	//printf("SUM: %lf \n\n",sum);

	//printf("%lf %lf",ApEn[0],ApEn[1]);
	apen = ApEn[0] - ApEn[1];
	chi_squared = 2.0*n*(log(2) - apen);
	p_value = cephes_igamc(pow(2, m - 2), chi_squared / 2.0);
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
		fprintf(stats[TEST_APEN], "\t\t(b) n (sequence length) = %d\n", n);
		fprintf(stats[TEST_APEN], "\t\t(c) Chi^2               = %f\n", chi_squared);
		fprintf(stats[TEST_APEN], "\t\t(d) Phi(m)	       = %f\n", ApEn[0]);
		fprintf(stats[TEST_APEN], "\t\t(e) Phi(m+1)	       = %f\n", ApEn[1]);
		fprintf(stats[TEST_APEN], "\t\t(f) ApEn                = %f\n", apen);
		fprintf(stats[TEST_APEN], "\t\t(g) Log(2)              = %f\n", log(2.0));
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");

		if (m > (int)(log(n) / log(2) - 5)) {
			fprintf(stats[TEST_APEN], "\t\tNote: The blockSize = %d exceeds recommended value of %d\n", m,
				MAX(1, (int)(log(n) / log(2) - 5)));
			fprintf(stats[TEST_APEN], "\t\tResults are inaccurate!\n");
			fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
		}

		fprintf(stats[TEST_APEN], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_APEN]);
		fprintf(results[TEST_APEN], "%f\n", p_value); fflush(results[TEST_APEN]);
	}
#endif
	free(P);
}

void
ApproximateEntropy4(Nist::Test& test)
{
	auto array = test.GetArray();
	auto n = test.GetN();
	auto& R_ = test.GetWriteableResults();
	auto m = test.GetParameters().approximateEntropyBlockLength;

	int				i, k, len, mi;
	double			sum, numOfBlocks, ApEn[2], apen, chi_squared, p_value;
	unsigned int	mask, help;
	int				*P;

	unsigned int cc = 0;

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_APEN], "\t\t\tAPPROXIMATE ENTROPY TEST\n");
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_APEN], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_APEN], "\t\t(a) m (block length)    = %d\n", m);
	}
#endif


	numOfBlocks = n;
	m++;
	mask = (1 << m) - 1;

	len = (1 << m);

	R_.approximate_entropy.P.resize(len * 2);

	if ((P = (int*)calloc(len, sizeof(int))) == NULL) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
		}
#endif
		return;
	}
	for (i = 0; i < len; i++)
		P[i] = 0;
	Histogram(0, P, m, n, array);

	for (i = 1; i<m; i++) {
		k = get_nth_block4(array, n - m + i)&(mask >> i);
		//printf("%d ",k);
		k ^= (((unsigned int*)array)[0] << (m - i));
		//printf("%d ",k);
		k &= mask;
		//printf("%d ",k);
		P[k]++;
	}

	//DISPLAY FREQUENCY
	sum = 0.0;
	for (i = 0; i < len / 2; i++) {
		mi = Mirrored_int(i, m - 1);
		help = P[mi] + P[mi + len / 2];
		if (help > 0)
			sum += help*log(help / numOfBlocks);
		R_.approximate_entropy.P[cc++] = help;
		//printf("%i ",help);
	}


	sum /= numOfBlocks;
	ApEn[0] = sum;

	sum = 0.0;
	for (i = 0; i < len; i++) {
		mi = Mirrored_int(i, m);
		if (P[mi] > 0)
			sum += P[mi] * log(P[mi] / numOfBlocks);
		//printf("[%d: %d] ",i,P[Mirrored_int(i,m)]);
		R_.approximate_entropy.P[cc++] = P[Mirrored_int(i, m)];
	}
	R_.approximate_entropy.pp = cc;

	sum /= numOfBlocks;
	ApEn[1] = sum;



	//printf("\n");
	//printf("SUM: %lf \n\n",sum);

	//printf("%lf %lf",ApEn[0],ApEn[1]);
	apen = ApEn[0] - ApEn[1];
	chi_squared = 2.0*n*(log(2) - apen);
	p_value = cephes_igamc(pow(2, m - 2), chi_squared / 2.0);
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
		fprintf(stats[TEST_APEN], "\t\t(b) n (sequence length) = %d\n", n);
		fprintf(stats[TEST_APEN], "\t\t(c) Chi^2               = %f\n", chi_squared);
		fprintf(stats[TEST_APEN], "\t\t(d) Phi(m)	       = %f\n", ApEn[0]);
		fprintf(stats[TEST_APEN], "\t\t(e) Phi(m+1)	       = %f\n", ApEn[1]);
		fprintf(stats[TEST_APEN], "\t\t(f) ApEn                = %f\n", apen);
		fprintf(stats[TEST_APEN], "\t\t(g) Log(2)              = %f\n", log(2.0));
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");

		if (m > (int)(log(n) / log(2) - 5)) {
			fprintf(stats[TEST_APEN], "\t\tNote: The blockSize = %d exceeds recommended value of %d\n", m,
				MAX(1, (int)(log(n) / log(2) - 5)));
			fprintf(stats[TEST_APEN], "\t\tResults are inaccurate!\n");
			fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
		}

		fprintf(stats[TEST_APEN], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_APEN]);
		fprintf(results[TEST_APEN], "%f\n", p_value); fflush(results[TEST_APEN]);
	}
#endif
	free(P);
}
