
#include "../Test.h"
#include "../original/cephes.h"

static double psi2(int m, int n, const unsigned char* epsilon);

void
Serial(Nist::Test& test)
{
	auto n = test.GetN();
	auto m = test.GetParameters().serialBlockLength;
	auto epsilon = test.GetEpsilon();
	auto& R_ = test.GetWriteableResults();

	double	p_value1, p_value2, psim0, psim1, psim2, del1, del2;

	psim0 = psi2(m, n, epsilon);
	psim1 = psi2(m - 1, n, epsilon);
	psim2 = psi2(m - 2, n, epsilon);
	//printf("%lf %lf\n",psim1,psim2);
#ifdef SPEED
	dummy_result = psim1 + psim2;
#endif
	del1 = psim0 - psim1;
	del2 = psim0 - 2.0*psim1 + psim2;
	p_value1 = cephes_igamc(pow(2, m - 1) / 2, del1 / 2.0);
	p_value2 = cephes_igamc(pow(2, m - 2) / 2, del2 / 2.0);
	//printf("%lf %lf\n",p_value1,p_value2);
#ifdef SPEED
	dummy_result = p_value1 + p_value2; // +psim0 + psim1 + psim2 + del1 + del2;
#endif

	R_.serial.psim0 = psim0;
	R_.serial.psim1 = psim1;
	R_.serial.psim2 = psim2;
	R_.serial.del1 = del1;
	R_.serial.del2 = del2;
	R_.serial.p_value1 = p_value1;
	R_.serial.p_value2 = p_value2;

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_SERIAL], "\t\t\t       SERIAL TEST\n");
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_SERIAL], "\t\t COMPUTATIONAL INFORMATION:		  \n");
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_SERIAL], "\t\t(a) Block length    (m) = %d\n", m);
		fprintf(stats[TEST_SERIAL], "\t\t(b) Sequence length (n) = %d\n", n);
		fprintf(stats[TEST_SERIAL], "\t\t(c) Psi_m               = %f\n", psim0);
		fprintf(stats[TEST_SERIAL], "\t\t(d) Psi_m-1             = %f\n", psim1);
		fprintf(stats[TEST_SERIAL], "\t\t(e) Psi_m-2             = %f\n", psim2);
		fprintf(stats[TEST_SERIAL], "\t\t(f) Del_1               = %f\n", del1);
		fprintf(stats[TEST_SERIAL], "\t\t(g) Del_2               = %f\n", del2);
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");

		fprintf(stats[TEST_SERIAL], "%s\t\tp_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
		fprintf(results[TEST_SERIAL], "%f\n", p_value1);

		fprintf(stats[TEST_SERIAL], "%s\t\tp_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2); fflush(stats[TEST_SERIAL]);
		fprintf(results[TEST_SERIAL], "%f\n", p_value2); fflush(results[TEST_SERIAL]);
	}
#endif
#ifdef KS
	pvals.serial_pvals[0][pvals.seq_counter] = p_value1;
	pvals.serial_pvals[1][pvals.seq_counter] = p_value2;
#endif
}

static double psi2(int m, int n, const unsigned char* epsilon)
{
	int				i, j, k, powLen;
	double			sum, numOfBlocks;
	unsigned int	*P;

	if ((m == 0) || (m == -1))
		return 0.0;
	numOfBlocks = n;
	powLen = (int)pow(2, m + 1) - 1;
	if ((P = (unsigned int*)calloc(powLen, sizeof(unsigned int))) == NULL) {
		printf("Serial Test:  Insufficient memory available.\n");
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[TEST_SERIAL], "Serial Test:  Insufficient memory available.\n");
			fflush(stats[TEST_SERIAL]);
		}
#endif
		return 0.0;
	}
	for (i = 1; i<powLen - 1; i++)
		P[i] = 0;	  /* INITIALIZE NODES */
	for (i = 0; i<numOfBlocks; i++) {		 /* COMPUTE FREQUENCY */
		k = 1;
		for (j = 0; j<m; j++) {
			if (epsilon[(i + j) % n] == 0)
				k *= 2;
			else if (epsilon[(i + j) % n] == 1)
				k = 2 * k + 1;
		}
		P[k - 1]++;
	}
	sum = 0.0;
	for (i = (int)pow(2, m) - 1; i<(int)pow(2, m + 1) - 1; i++)
	{
		sum += pow(P[i], 2);
		//printf("%d %d ",i,P[i]);
	}

	sum = (sum * pow(2, m) / (double)n) - (double)n;
	free(P);

	//printf("%lf\n",sum);
	return sum;
}
