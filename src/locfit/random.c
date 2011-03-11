/*
 *   Copyright (c) 1996-2000 Lucent Technologies.
 *   See README file for details.
 */

#include "local.h"

#define PI_HALF         1.5707963267948966192313216916397514420986 /*pi/2*/
#define PI_QUARTER      0.7853981633974483096156608458198757210493 /*pi/4*/
#define EXP78           1.0129030479320018583185514777512982888868 /*e^(1/78)*/
#define PI128          40.7436654315252059568342434233636766808217 /*128/pi*/

static unsigned long cc, tv, ss=1;

double runif()
{	if (ss)
	{ WARN(("runif: No seed set."));
          return(0.0);
	}
        cc = cc * 69069;    /* congruential part */
        tv ^= tv >> 15;       /* tausworthe part */
        tv ^= tv << 17;
	return(((tv ^ cc) >> 1) / 2147483648.0);
}

void rseed(seed)
/*
  Seed should be string of at least 8 characters.
*/
char *seed;
{	ss = 0;
	tv = seed[0];
        tv = (tv<<8) | seed[1];
        tv = (tv<<8) | seed[2];
        tv = (tv<<8) | seed[3];
 	cc = seed[4];
        cc = (cc<<8) | seed[5];
        cc = (cc<<8) | seed[6];
        cc = (cc<<8) | seed[7];
        if(cc % 2 == 0)
             cc++;
}

/*
 * Gaussian random variable.
 * Reference: Kinderman & Monahan, Proceedings of
 * the ASA, Statistical Computing Section, 1975, 128-131.
 */
double rnorm(mu,s)
double mu, s;
{
	double rnormk, u, x2;

	do {
		u = runif();
		rnormk = 1.715527769 * (runif()-0.5) / u;
		x2 = rnormk * rnormk / 4;
	} while((x2>1-u) || (x2 > -log(u)));
	return(mu+s*rnormk);
}

double rexp(lb)
double lb;
{ return(-log(runif())/lb);
}

/*
 * Poisson random variable.
 * Simple algorithm for small lambda, else complex algorithm.
 * Crossover point must be at least 5 for the complex algorithm
 * to work correctly.
 * Reference: Devroye, pages 504, 511 and 516 (with corrections!)
 */
double rpois(lambda)
double lambda;
{
	static double olambda = -1, a, mu, delta, d, c1, c2, c3, c4, c5;
	double u, e, n, x, y, w, t, p, q;
	int new = lambda != olambda;

	olambda = lambda;
	if(lambda < 8) {
		if(new)
			a = exp(-lambda);
		q = 1;
		x = -1;
		do {
			q *= runif();
			x++;
		} while(q >= a);
		return(x);
	}

	if(new) {
		mu = floor(lambda);
		delta = sqrt(2 * mu * log(mu * PI128));
		delta = MAX(6.0, MIN(mu, floor(delta)));
		d = 2*mu + delta;
		c1 = sqrt(mu * PI_HALF);
		c2 = c1 + sqrt(d * PI_QUARTER) * exp(1/d);
		c3 = c2 + 1;
		c4 = c3 + EXP78;
		c5 = c4 + 2 * d * exp(-delta*(1+delta/2)/d) / delta;
	}
	while(1) {
		u = c5 * runif();
		e = -log(runif());
		if(u <= c1) {
			n = rnorm(0.0,1.0);
			x = floor(-fabs(n) * sqrt(mu));
			if(x < -mu)
				continue;
			w = n*n/2 + e + x*log(lambda/mu);
		} else if(u <= c2) {
			y = 1 + fabs(rnorm(0.0,1.0)) * sqrt(d/2);
			x = ceil(y);
			if(x > delta)
				continue;
			w = y*(y-2)/d + e + x*log(lambda/mu);
		} else if(u <= c3) {
			x = 0;
			w = e;
		} else if(u <= c4) {
			x = 1;
			w = e + log(lambda/mu);
		} else {
			y = delta - 2*d*log(runif())/delta;
			x = ceil(y);
			w = delta*(1+y/2)/d + e + x*log(lambda/mu);
		}
		w = -w;
		t = x*(x+1) / (2*mu);
		if(x >= 0 && w <= -t)
			return(x+mu);
		if(x < 0 && w > -t)
			continue;
		q = t * ((2*x+1)/(6*mu) - 1);
		if(w > q)
			continue;
		p = x+1 <= 0 ? x+1 : 0;
		p = q - t*t/(3*(mu+p));
		if(w <= p)
			return(x+mu);
		if(w <= x*log(mu) - LGAMMA(mu+x+1) + LGAMMA(mu+1))
			return(x+mu);
	}
}
