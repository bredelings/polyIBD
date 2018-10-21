
#pragma once

#include <random>

//------------------------------------------------
// draw z probabilitiy distribution to get at effMOI
double psamp(vector<double> p);
int sampleZ(vector<double> p);

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
double runif_0_1();

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
double runif1(double a=0, double b=1.0);

//------------------------------------------------
// draw from Bernoulli(p) distribution
bool rbernoulli1(double p);

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd);

//------------------------------------------------
// draw from univariate normal distribution and reflect to interval (a,b)
double rnorm1_interval(double mean, double sd, double a, double b);

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal probability
int sample2(int a, int b);

//------------------------------------------------
// sample single value from given probability vector (that sums to pSum)
int sample1(std::vector<double> &p, double pSum=1.0);

//------------------------------------------------
// draw from gamma(shape,rate) distribution
double rgamma1(double shape, double rate);

//------------------------------------------------
// draw from beta(shape1,shape2) distribution
double rbeta1(double shape1, double shape2);

//------------------------------------------------
// probability mass of poisson distribution
double dpois1(int n, double lambda, bool returnLog=true);

//------------------------------------------------
// draw from poisson with lambda mean and variance
double rpois1(double lambda);

//------------------------------------------------
// draw from negative binomial distribution with mean lambda and variance gamma*lambda (gamma must be >1)
int rnbinom1(double lambda, double gamma);

//------------------------------------------------
// probability mass of negative binomial distribution with mean lambda and variance gamma*lambda (gamma must be >1)
double dnbinom1(int n, double lambda, double gamma, bool returnLog=true);

//------------------------------------------------
// draw from symmetric dichlet(alpha) distribution of length n
std::vector<double> rdirichlet1(double alpha, int n);
