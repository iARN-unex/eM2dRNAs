/* Declaration for random number related variables and routines */

# ifndef _RAND_H_
# define _RAND_H_

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* Variable declarations for the random number generator */

class Rand
{
  private:
  	long rnd_uni_init;
  	long idum2=123456789;
  	long iy=0;
  	long iv[NTAB];

  public:
  	Rand(int seed);

	/* Function declarations for the random number generator */
	long get_rnd_uni_init();
  	double randomperc(void);
  	int rnd (int low, int high);
  	double rndreal (double low, double high);
	double rnd_uni(long int*);
};
# endif
