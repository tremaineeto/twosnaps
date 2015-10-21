#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cstdio>
#include <math.h>
#include <float.h>
#include <limits>
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

double logSumP(double x, double y) {
    double result = 0;
    
    result = x + log10(1 + pow(10, (y - x)));
    
    return result;
}

double safe_log(double theta){
    
    if (theta == 0) {
        return -DBL_MAX;
    }
    else {
        return log10(theta);
    }
}

double scoreBasic(double obsUnmutated, double obsMutated, double theta) {
    
    double firstParam = obsMutated + safe_log(theta);
    double secondParam = obsUnmutated + safe_log(1-theta);
    
    if (secondParam > firstParam)
    {
        return logSumP(secondParam, firstParam);
    }
    
    return logSumP(firstParam, secondParam);
}

double scoreAdvanced() {
    double score = 0;
    return score;
}

int main(int argc, char *argv[])
{
    int id;
    
    double firstBasicScore = 0;
    double secondBasicScore = 0;
    
    double prob_obs_unmut, prob_obs_mut;
    double score = 0;
    
    while(EOF != read_data(stdin, &id, &prob_obs_unmut, &prob_obs_mut)) {
        firstBasicScore += scoreBasic(prob_obs_unmut, prob_obs_mut, 0.2);
        secondBasicScore += scoreBasic(prob_obs_unmut, prob_obs_mut, 0);
    }
    
    score = firstBasicScore - secondBasicScore;
    
    print_basic_score(stdout, score);
}
