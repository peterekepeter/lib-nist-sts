
#include "../Test.h"
#include "../common/cephes.h"

#define MAXNUMOFTEMPLATES				148		/* APERIODIC TEMPLATES: 148=>temp_length=9 */
#define MIN(x,y)             ((x) >  (y)  ? (y)  : (x))

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
N O N O V E R L A P P I N G  T E M P L A T E  T E S T
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
NonOverlappingTemplateMatchings(Nist::Test& test)
{
	auto n = test.GetN();
	auto m = test.GetParameters().nonOverlappingTemplateBlockLength;
	auto epsilon = test.GetEpsilon();
	auto& R_ = test.GetWriteableResults();
	using BitSequence = unsigned char;

	int		numOfTemplates[100] = { 0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
		2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152 };
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must
	first be constructed, saved into files and then the corresponding
	number of nonperiodic templates for that file be stored in the m-th
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/
	unsigned int	bit, W_obs, /* nu[6],*/ *Wj = NULL;
	FILE			*fp = NULL;
	double			sum, chi2, p_value, lambda, pi[6], varWj;
	int				i, j, jj, k, match, SKIP, M, N, K = 5;
	char			directory[100];
	BitSequence		*sequence = NULL;

	N = 8;
	M = n / N;

	if ((Wj = (unsigned int*)calloc(N, sizeof(unsigned int))) == NULL) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
			fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
		}
#endif
		printf("\tNONOVERLAPPING TEMPLATES TESTS: Insufficient memory for required work space.\n");
		return;
	}
	lambda = (M - m + 1) / pow(2, m);
	varWj = M*(1.0 / pow(2.0, m) - (2.0*m - 1.0) / pow(2.0, 2.0*m));
	sprintf(directory, "templates/template%d", m);

	if (((lambda < 0.0) || (lambda == 0.0)) ||
		((fp = fopen(directory, "r")) == NULL) ||
		((sequence = (BitSequence *)calloc(m, sizeof(BitSequence))) == NULL)) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
			fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
			fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
			fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
		}
#endif
		printf("\tNONOVERLAPPING TEMPLATES TESTS ABORTED.\n");
		if (sequence != NULL)
		{
			free(sequence);
			sequence = NULL;
		}

#ifdef VERIFY_RESULTS
		R_.nonoverlapping.templates = 0;
		R_.nonoverlapping.W = NULL;
		R_.nonoverlapping.chi2 = NULL;
		R_.nonoverlapping.p_value = NULL;
#endif
	}
	else {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
			fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		}
#endif

		if (numOfTemplates[m] < MAXNUMOFTEMPLATES)
			SKIP = 1;
		else
			SKIP = (int)(numOfTemplates[m] / MAXNUMOFTEMPLATES);
		numOfTemplates[m] = (int)numOfTemplates[m] / SKIP;

#ifdef KS
		pvals.num_Nonoverlap_pvals = (int)numOfTemplates[m] / SKIP;
#endif

#ifdef VERIFY_RESULTS
		R_.nonoverlapping.templates = MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]);
		R_.nonoverlapping.W = (unsigned int*)malloc(sizeof(unsigned int)*R_.nonoverlapping.templates * 8);
		R_.nonoverlapping.chi2 = (double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		R_.nonoverlapping.p_value = (double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		if (R_.nonoverlapping.W == NULL || R_.nonoverlapping.chi2 == NULL || R_.nonoverlapping.p_value == NULL)
		{
			printf("NONOVERLAPPING TEMPLATES TEST: Cannot allocate memory"); return;
		}
#endif


		sum = 0.0;
		for (i = 0; i < 2; i++) {                      /* Compute Probabilities */
			pi[i] = exp(-lambda + i*log(lambda) - cephes_lgam(i + 1));
			sum += pi[i];
		}
		pi[0] = sum;
		for (i = 2; i <= K; i++) {                      /* Compute Probabilities */
			pi[i - 1] = exp(-lambda + i*log(lambda) - cephes_lgam(i + 1));
			sum += pi[i - 1];
		}
		pi[K] = 1 - sum;

		for (jj = 0; jj < MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]); jj++) {
			sum = 0;

			for (k = 0; k < m; k++) {
				fscanf(fp, "%d", &bit);
				sequence[k] = (unsigned char)bit;
#if defined(FILE_OUTPUT) ||  defined(KS)
				if (cmdFlags.output == 1 || cmdFlags.output == -1) {
					fprintf(stats[TEST_NONPERIODIC], "%d", sequence[k]);
				}
#endif
			}
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1) {
				fprintf(stats[TEST_NONPERIODIC], " ");
			}
#endif
			// nadbytecny kod? nu nikde nevyuzito...
			//			for ( k=0; k<=K; k++ )
			//				nu[k] = 0;
			for (i = 0; i < N; i++) {
				W_obs = 0;
				for (j = 0; j < M - m + 1; j++) {
					match = 1;
					for (k = 0; k < m; k++) {
						if ((int)sequence[k] != (int)epsilon[i*M + j + k]) {
							match = 0;
							break;
						}
					}
					if (match == 1)
					{
						W_obs++;
						j += m - 1;
					}
				}
				Wj[i] = W_obs;
			}
			sum = 0;
			chi2 = 0.0;                                   /* Compute Chi Square */
			for (i = 0; i < N; i++) {
#if defined(FILE_OUTPUT) ||  defined(KS)
				if (cmdFlags.output == 1 || cmdFlags.output == -1) {
					if (m == 10)
						fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[i]);
					else
						fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[i]);
				}
#endif
				/*if ( m == 10 )
				printf("%3d  ", Wj[i]);
				else
				printf("%4d ", Wj[i]);*/
#ifdef VERIFY_RESULTS
				R_.nonoverlapping.W[jj*N + i] = Wj[i];
#endif
				//
				chi2 += pow(((double)Wj[i] - lambda) / pow(varWj, 0.5), 2);
			}
			p_value = cephes_igamc(N / 2.0, chi2 / 2.0);
#ifdef SPEED
			dummy_result += p_value;
#endif
#ifdef VERIFY_RESULTS
			R_.nonoverlapping.chi2[jj] = chi2;
			R_.nonoverlapping.p_value[jj] = p_value;
#endif

#ifdef KS
			pvals.nonoverlapping_pvals[jj][pvals.seq_counter] = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1) {
				if (isNegative(p_value) || isGreaterThanOne(p_value))
					fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

				fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", jj);
			}
#endif
			if (SKIP > 1)
				fseek(fp, (long)(SKIP - 1) * 2 * m, SEEK_CUR);
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1) {
				fprintf(results[TEST_NONPERIODIC], "%f\n", p_value); fflush(results[TEST_NONPERIODIC]);
			}
#endif
		}
	}

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_NONPERIODIC], "\n"); fflush(stats[TEST_NONPERIODIC]);
	}
#endif
	if (sequence != NULL)
		free(sequence);

	free(Wj);
	if (fp != NULL)
		fclose(fp);

#ifdef VERIFY_RESULTS
	if (NonOverlappingTemplateMatchings_v1 == NonOverlappingTemplateMatchings) R1 = R_;
	else R2 = R_;
#endif
}
