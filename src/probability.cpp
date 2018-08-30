
#include <Rcpp.h>
#include <cmath>
#include "probability.h"
#include "misc.h"

using namespace std;

// comment this line out to use R default random functions
//#define USE_MY_RANDOM

//-- set random seed --
random_device rd;
default_random_engine generator(rd());

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
#ifdef USE_MY_RANDOM
double runif_0_1() {
    uniform_real_distribution<double> uniform_0_1(0.0,1.0);
    return(uniform_0_1(generator));
}

#else
double runif_0_1() {
    return(R::runif(0,1));
}
#endif

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
#ifdef USE_MY_RANDOM
double runif1(double a, double b) {
    uniform_real_distribution<double> uniform_a_b(a,b);
    return(uniform_a_b(generator));
}
#else
double runif1(double a, double b) {
    return(R::runif(a,b));
}
#endif

//------------------------------------------------
// draw from Bernoulli(p) distribution
bool rbernoulli1(double p) {
    bernoulli_distribution dist_bernoulli(p);
    return(dist_bernoulli(generator));
}

//------------------------------------------------
// draw from univariate normal distribution
#ifdef USE_MY_RANDOM
double rnorm1(double mean, double sd) {
    normal_distribution<double> dist_norm(mean,sd);
    return(dist_norm(generator));
}

#else
double rnorm1(double mean, double sd) {
    return(R::rnorm(mean, sd));
}
#endif

//------------------------------------------------
// draw from univariate normal distribution and reflect to interval (a,b)
double rnorm1_interval(double mean, double sd, double a, double b) {

    // draw raw value relative to a
    double ret = rnorm1(mean, sd) - a;

    // reflect off boundries at 0 and (b-a)
    if (ret<0 || ret>(b-a)) {
        // use multiple reflections to bring into range [-(b-a), 2(b-a)]
        while (ret < -(b-a)) {
            ret += 2*(b-a);
        }
        while (ret > 2*(b-a)) {
            ret -= 2*(b-a);
        }

        // use one more reflection to bring into range [0, (b-a)]
        if (ret < 0) {
            ret = -ret;
        }
        if (ret > (b-a)) {
            ret = 2*(b-a) - ret;
        }
    }

    // no longer relative to a
    ret += a;

    // don't let ret equal exactly a or b
    if (ret==a) {
        ret += UNDERFLO;
    } else if (ret==b) {
        ret -= UNDERFLO;
    }

    return(ret);
}

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal probability
int sample2(int a, int b) {
    int ret = floor(runif1(a, b+1));
    return(ret);
}

//------------------------------------------------
// sample single value from given probability vector (that sums to pSum)
int sample1(vector<double> &p, double pSum) {
    double rand = pSum*runif_0_1();
    double z = 0;
    for (int i=0; i<int(p.size()); i++) {
        z += p[i];
        if (rand<z) {
            return i+1;
        }
    }
    return(0);
}

//------------------------------------------------
// draw from gamma(shape,rate) distribution
#ifdef USE_MY_RANDOM
double rgamma1(double shape, double rate) {
    gamma_distribution<double> rgamma(shape,1.0/rate);
    double x = rgamma(generator);

    // check for zero or infinite values (catches bug present in Visual Studio 2010)
    if (x==0) {
        x = UNDERFLO;
    }
    if ((1.0/x)==0) {
        x = 1.0/UNDERFLO;
    }

    return(x);
}

#else
double rgamma1(double shape, double rate) {
    return(R::rgamma(shape, 1.0/rate));
}
#endif


//------------------------------------------------
// draw from beta(shape1,shape2) distribution
#ifdef USE_MY_RANDOM
double rbeta1(double shape1, double shape2) {
    double x1 = rgamma1(shape1,1.0);
    double x2 = rgamma1(shape2,1.0);
    return(x1/double(x1+x2));
}

#else
double rbeta1(double shape1, double shape2) {
    return(R::rbeta(shape1, shape2));
}
#endif

//------------------------------------------------
// probability mass of Poisson distribution
#ifdef USE_MY_RANDOM
double dpois1(int n, double lambda, bool returnLog) {
    double ret = n*log(lambda) - lambda - lgamma(n+1)
    if (!returnLog) {
        ret = exp(ret);
    }
    return(ret);
}

#else
double dpois1(int n, double lambda, bool returnLog) {
    return(R::dpois(n,lambda,returnLog));
}
#endif

//------------------------------------------------
// draw from negative binomial distribution with mean lambda and variance gamma*lambda (gamma must be >1)
int rnbinom1(double lambda, double gamma) {
    double ret = R::rnbinom(lambda/(gamma-1), 1/gamma);
    return(ret);
}


//------------------------------------------------
// probability mass of negative binomial distribution with mean lambda and variance gamma*lambda (gamma must be >1)
double dnbinom1(int n, double lambda, double gamma, bool returnLog) {
    return(R::dnbinom(n, lambda/(gamma-1), 1/gamma, returnLog));
}

//------------------------------------------------
// draw from symmetric dichlet(alpha) distribution of length n
vector<double> rdirichlet1(double alpha, int n) {
    vector<double> ret(n);
    double retSum = 0;
    for (int i=0; i<n; i++) {
        ret[i] = rgamma1(alpha,1.0);
        retSum += ret[i];
    }
    for (int i=0; i<n; i++) {
        ret[i] /= retSum;
    }
    return(ret);
}
