#define _CRT_SECURE_NO_DEPRECATE
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
using namespace std;

static double a[20][3];
static double b[20][3];

int read_data(FILE *fp,	//reads each obs per line
		int *id,
		double *prob_obs_unmut,
		double *prob_obs_mut)
{
	if(3 != fscanf(fp, "%d %lf %lf\n", id, prob_obs_unmut, prob_obs_mut)) {
		return EOF;
	}
	return 1;
}

double safe_log(double num)
{
	if(num == 0)
		return -DBL_MAX;
	return log10(num);
}

double logSum(double log1, double log2)
{
	if (log2 >= log1) //same as if we compared p(x) to p(y)
	{
		double temp = log1;
		log1 = log2;
		log2 = temp;
	}
	return log1 + log10(1 + pow(10, log2 - log1));
}

double basic_logp_i(double unmut, double mut, double theta) //calculates p(obs/theta) for one line
{
	double sum1 = mut + safe_log(theta);
	double sum2 = unmut + safe_log(1 - theta);
	return logSum(sum1, sum2);
}

double calculateL(double theta)
{
	for (int i = 0; i < 20; i++)	//use another array b, since we have to calculate L multiple times
		for (int j = 0; j < 3; j++)
			b[i][j] = a[i][j];

	for (int k = 0; k < 3; k++)
	{
		int coeff = 0;
		if (k == 0)
			coeff = 1;
		else if (k == 1)
			coeff = 2;
		else coeff = 1;

		

		double bin = coeff*pow(theta, k)*pow((1 - theta), 2 - k);
		if (theta == 0)
		{
			bin = 0;
			if (k == 0)
				bin = 1;
		}
		if (theta == 1)
		{
			bin = 0;
			if (k == 2)
				bin = 1;
		}

		for (int i = 0; i < 20; i++)
		{
			b[i][k] += safe_log(bin);	//this incorporates the p(k/theta) term to the equation
		}
	}

/*	if (theta == 1)
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 3; j++)
			cout << a[i][j] << " ";
		cout << endl;
	}*/

	double L = 0;
	for (int i = 0; i < 20; i++)	//this calculates L
	{
		double logInd = b[i][0];
		for (int j = 0; j < 2; j++)
			logInd = logSum(logInd, b[i][j + 1]);
		L += logInd;
	}
	return L;
}

double calculateQ(int n) //recursive function
{
	if (n == 0)
		return -DBL_MAX;
	return logSum(calculateQ(n - 1), logSum(calculateL((n - 1) / 100), calculateL(n / 100)) + log10(1.0 / 100.0 / 2.0));
}

void compute_scores(FILE *fp) //computes the basic and adv scores
{
	int id;
	double prob_obs_unmut, prob_obs_mut;
//	double score, phi, logp;
//	score = phi = logp = 0;
	double theta = 0.2;
	double score_num = 0;
	double score_denom = 0;


	while (EOF != read_data(fp, &id, &prob_obs_unmut, &prob_obs_mut))
	{
		score_num += basic_logp_i(prob_obs_unmut, prob_obs_mut, theta);	//basic
		score_denom += basic_logp_i(prob_obs_unmut, prob_obs_mut, 0);	//basic
		for (int k = 0; k < 3; k++)	//populates a 2d array with ID on one side and kappa on the other
		{		//this array will be used later on to calculate L
			double unmut = safe_log(1 - k / 2) + prob_obs_unmut;
			double mut = safe_log(k / 2) + prob_obs_mut;
			a[id][k] += logSum(unmut, mut);
		}
	}

//	cout << logSum(prob_obs_mut + safe_log(0), prob_obs_unmut + safe_log(1)) << endl;
//	cout << 1.0 / 100 / 2.0 <<endl;
//	cout << log10(1 / 100 / 2) << endl;

//	cout << calculateL(0) << endl;
	for (int i = 0; i < 101; i++)
	cout << calculateQ(i) << " " << i << endl;
	cout << calculateL(0) << " L0" << endl;
	double adv_score = calculateQ(100) - calculateL(0);
	for (int i = 0; i <= 100; i++)
		cout << calculateQ(i) << endl;
	cout << "BASIC_SCORE=" << score_num - score_denom << endl;
	cout << "ADV_SCORE=" << adv_score << endl; 
}

int main(int argc, char *argv[])
{
	FILE *pFile;
	pFile = fopen("C:\\Users\\hanni_000\\Documents\\Visual Studio 2013\\Projects\\CS121Project1\\Debug\\seq.1.txt", "r");

	compute_scores(pFile);
	
//	print_adv_score(stdout, score);
//	print_logp(stdout, phi, logp);
}
