#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cstdio>
#include <math.h>
using namespace std;

int read_data(FILE *fp,
		int *id,
		double *prob_obs_unmut,
		double *prob_obs_mut)
{
	if(3 != fscanf(fp, "%d %lf %lf\n", id, prob_obs_unmut, prob_obs_mut)) {
		return EOF;
	}
	return 1;
}

int print_basic_score(FILE *fp, double score)
{
	if(1 != fprintf(fp, "BASIC_SCORE=%lf\n", score)) {
       // test = scoreBasic(score);
		return EOF;
	}
	return 1;
}

int print_adv_score(FILE *fp, double score)
{
	if(1 != fprintf(fp, "ADV_SCORE=%lf\n", score)) {
		return EOF;
	}
	return 1;
}

int print_logp(FILE *fp, double phi, double logp) {
	if(2 != fprintf(fp, "LOG_P_PHI %lf %lf\n", phi, logp)) {
		return EOF;
	}
	return 1;
}

double scoreBasic(double obsMutated, double obsUnmutated) {
    // score = LogP(obsMutated|theta = 0.2) - LogP(obsUnmutated|theta = 0)
    double score = 0;
    double mutatedTheta = 0.2;
    int unMutatedTheta = 0;
    
    // TODO: enter the math needed to calculate score here
    
    // log of Bayes Theorem for Mutated - log of Bayes Theorem for Unmutated
    
    return score;
    
}

double scoreAdvanced() {
    double score = 0;
    return score;
}

int main(int argc, char *argv[])
{
	int id;
    cout << "test" << endl;
    
	double prob_obs_unmut, prob_obs_mut;
	double score, phi, logp;
	score = phi = logp = 0; // These need to be computed
    
    // score = scoreBasic(prob_obs_mut, prob_obs_unmut);
    
	while(EOF != read_data(stdin, &id, &prob_obs_unmut, &prob_obs_mut)) {
		//fprintf(stdout, "%d %lf %lf\n", id, prob_obs_unmut, prob_obs_mut);
		// We should check the return values below!:
		print_basic_score(stdout, score);
		print_adv_score(stdout, score);
		print_logp(stdout, phi, logp);
	}
}
