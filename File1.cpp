#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <vector>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/format.hpp>

using namespace boost::multiprecision;
#pragma hdrstop
#pragma argsused

/* clip_pow2 should be greater than number of bits of accuracy of the t_prob, but not overflow with 2^clip_prob */


//typedef double t_prob;
//#define clip_pow2 100

//typedef boost::rational<cpp_int> t_prob;
//#define clip_pow2 10000000

typedef number<cpp_bin_float<50>> t_prob;
#define clip_pow2 10000





#define n_combin_table 50

t_prob combin_table [n_combin_table] [n_combin_table];

t_prob *comb_prob_this;
t_prob *comb_prob_next;


#ifdef _WIN32
#include <tchar.h>
#else
  typedef char _TCHAR;
  #define _tmain main
#endif

void init_combin_table (void)
{	int i, j;

	combin_table [0] [0] = 1;
	for (i = 1; i < n_combin_table; i++)
	{	combin_table [i] [0] = 1;
		for (j = 1; j <= i; j++)
		{	combin_table [i] [j] = combin_table [i - 1] [j - 1] + combin_table [i - 1] [j];
		}
	}
}

t_prob two_to_n (int n) {
	t_prob r;
	t_prob two_p;

	r = 1;
	two_p = 2;
	while (n > 0) {
		if (n & 1) {
			r *= two_p;
		}
		two_p *= two_p;
        n /= 2;
	}
	return (r);
}

t_prob binom_prob (int n, int m) {
	return combin_table [m] [n] / two_to_n (m);
}

void calc_prob (void)
{   t_prob probs [n_combin_table];
	t_prob pthis;
	int iprob;
	int jprob;

	probs [0] = 0;
	probs [1] = 1;
	for (iprob = 2; iprob < n_combin_table; iprob++) {
		pthis = 0;
		for (jprob = 0; jprob < iprob; jprob++) {
			pthis += binom_prob (jprob, iprob) * probs [jprob];
		}
		pthis *= two_to_n (iprob) / (two_to_n (iprob) - 1);
		probs [iprob] = pthis;
		std::cout << " " << iprob << " " << pthis << "\n";
	}
}

void calc_prob_fast (int nmax)
{   t_prob *comb_prob_this;
	t_prob *comb_prob_next;
	t_prob *comb_prob_temp;
	t_prob *probs;
	t_prob pthis;

	int iprob;
	int jprob;

	comb_prob_this = new t_prob [nmax];
	comb_prob_next = new t_prob [nmax];
	probs = new t_prob [nmax];
	probs [0] = 0;
	probs [1] = 1;

	/* set to 0.5 verbosely so works with rational */

	comb_prob_this [0] = 1;
	comb_prob_this [0] /= 2;
	comb_prob_this [1] = comb_prob_this [0];

	for (iprob = 2; iprob < nmax; iprob++) {
		comb_prob_next [0] = comb_prob_this [0] / 2;
		for (jprob = 1; jprob < iprob; jprob++) {
			comb_prob_next [jprob] = (comb_prob_this [jprob - 1] + comb_prob_this [jprob]) / 2;
		}
		comb_prob_next [jprob] = comb_prob_this [jprob - 1] / 2;
		comb_prob_temp = comb_prob_this;
		comb_prob_this = comb_prob_next;
		comb_prob_next = comb_prob_temp;

		pthis = 0;
		for (jprob = 0; jprob < iprob; jprob++) {
			pthis += comb_prob_this [jprob] * probs [jprob];
		}
		if (iprob < clip_pow2) {
			pthis *= two_to_n (iprob) / (two_to_n (iprob) - 1);
        }
		probs [iprob] = pthis;
		std::cout << " " << iprob << " " << pthis << "\n";
	}



}

int _tmain(int argc, _TCHAR* argv[])
{   int nmax;

	sscanf (argv [1], "%d", &nmax);

    std::cout << std::setprecision (50);
//	init_combin_table ();

//	calc_prob ();
	calc_prob_fast (nmax);

	return 0;
}
