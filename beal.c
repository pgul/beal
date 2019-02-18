#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

#define LOG_PRECISE 0.0001
#define MAX_A		1000
#define MAX_B		1000
#define MAX_C		1000
#define MAX(a, b)	((a) > (b) ? (a) : (b))
#define MAX_BC		MAX(MAX_B, MAX_C)
#define MAX_ABC		MAX(MAX_A, MAX_BC)
#define MAX_POW		1000
#define MIN_POW		3
#define MULT_BASE   (8u*9*5*7*11*13*17*19*23)
// #define MULT_BASE	39916801
// #define HASH_BASE	999983 // it's a prime number
// #define HASH_BASE	(10*999983) // it's a prime number
// #define HASH_BASE	39916801 // it's a prime number
// #define HASH_BASE	1999993 // it's a prime number
// #define HASH_BASE	3999971 // it's a prime number
#define HASH_BASE	(9999991) // it's a prime number
#define VAL_BY_HASH	16
#define PRIMES_CNT	1000000
#define PREPARE_POWERS 1

typedef unsigned long int INT;

int n_primes;
unsigned int prime[PRIMES_CNT];
unsigned int prime_to_a_pow[PRIMES_CNT];
unsigned int prime_to_b_pow[PRIMES_CNT];
char has_mod_n[((MULT_BASE - 1) >> 3) + 1];

double logn[MAX_ABC];
unsigned int power_mod_n[MAX_BC][MAX_POW+1];
#if PREPARE_POWERS
mpz_t power[MAX_BC][MAX_POW+1];
#else
mpz_t man, mbm, mck;
#endif
mpz_t mpz_temp;
unsigned int hash_base[HASH_BASE][VAL_BY_HASH];

unsigned int a, b, c, n, m;
INT power_base_a_n, power_base_b_m, power_base_an_bm;
double an, bm, log_lower, log_upper;

unsigned int power_base_n(unsigned int a, unsigned int n)
{
	INT p;
	unsigned int s;

	if (a < MAX_BC && n <= MAX_POW)
		return power_mod_n[a][n];
	p = a % MULT_BASE;
	s = (n & 1) ? p : 1;
	while (n > 1)
	{
		p = (p * p) % MULT_BASE;
		n /= 2;
		if (n & 1)
			s = (s * p) % MULT_BASE;
	}
	return s;
}

int check_bigint(unsigned int k)
{
	int res;
#if PREPARE_POWERS
	mpz_add(mpz_temp, power[a][n], power[b][m]);
	res = mpz_cmp(mpz_temp, power[c][k]);
#else
	mpz_ui_pow_ui(man, a, n);
	mpz_ui_pow_ui(mbm, b, m);
	mpz_ui_pow_ui(mck, c, k);
	mpz_add(mpz_temp, man, mbm);
	res = mpz_cmp(mpz_temp, mck);
#endif
	return res;
}

INT check(void)
{
	double logc;
	unsigned int k;

	logc = logn[c];
	k = (unsigned int)(log_lower / logc) + 1;
	if (k < MIN_POW)
		return 0;
	if (k * logc > log_upper)
		return 0;
	if (power_base_an_bm != power_base_n(c, k))
		return 0;
	if (check_bigint(k))
	{
		printf("  - %u^%u + %u^%u = %u^%u\n", a, n, b, m, c, k);
		return 0;
	}
	printf("*** %u^%u + %u^%u = %u^%u\n", a, n, b, m, c, k);
	exit(0);
	return k;
}

double log_1_5, log_1_25, log_1_125, log_1_0625;

double upper_log(double an, double bm)
{
	double log_upper;

	if (an > bm + logn[16])
		log_upper = log_lower + log_1_0625 + LOG_PRECISE;
	else if (an > bm + logn[8])
		log_upper = log_lower + log_1_125 + LOG_PRECISE;
	else if (an > bm + logn[4])
		log_upper = log_lower + log_1_25 + LOG_PRECISE;
	else if (an > bm + logn[2])
		log_upper = log_lower + log_1_5 + LOG_PRECISE;
	else
		log_upper = log_lower + logn[2] + LOG_PRECISE;
	return log_upper;
}

void init_prime()
{
	int i, j;
	prime[0] = 2;
	n_primes = 1;
	for (i=3; i<MAX_ABC && n_primes<PRIMES_CNT; i+=2)
	{
		for (j=1; j<n_primes; j++)
		{
			if (i % prime[j] == 0)
				break;
		}
		if (j == n_primes)
			prime[n_primes++] = i;
	}
	if (n_primes == PRIMES_CNT)
		fprintf(stderr, "Increase PRIMES_CNT\n");
	printf("Max prime number %u\n", prime[n_primes-1]);
}

void init_log(void)
{
	int i;
	for (i=1; i<MAX_ABC; i++)
		logn[i] = log(i);
}

void init_power(void)
{
	int i, j, k;
	int hash;
	int has_mods = 0, many_vals = 0, collisions = 0;
	mpz_init(mpz_temp);
#if PREPARE_POWERS == 0
	mpz_init(man);
	mpz_init(mbm);
	mpz_init(mck);
#endif
	for (i=1; i<MAX_BC; i++)
	{
		power_mod_n[i][0] = 1;
#if PREPARE_POWERS
		mpz_init_set_ui(power[i][0], 1);
#endif
		for (j=1; j<=MAX_POW; j++)
		{
			unsigned int mod_n = ((INT)(power_mod_n[i][j-1]) * i) % MULT_BASE;
#if PREPARE_POWERS
			mpz_init(power[i][j]);
			mpz_mul_si(power[i][j], power[i][j-1], i);
#endif
			power_mod_n[i][j] = mod_n;
			if ((has_mod_n[mod_n >> 3] & (1 << (mod_n & 7))) == 0)
			{
				has_mods++;
				has_mod_n[mod_n >> 3] |= 1 << (mod_n & 7);
			}
			hash = mod_n % HASH_BASE;
			for (k=0; k<VAL_BY_HASH; k++) {
				if (hash_base[hash][k] == 0) {
					hash_base[hash][k] = i;
					if (k == VAL_BY_HASH - 1)
						collisions++;
					break;
				}
				if (hash_base[hash][k] == i)
					break;
			}
			if (k >= VAL_BY_HASH-1)
				many_vals++;
		}
	}
	printf("Has %u modules from %u total with %u collisions %u values out of hash\n", has_mods, MULT_BASE, collisions, many_vals);
}

int main(void)
{
	int i, j;
	INT checked_a = 0, tested_b = 0, checked_b = 0, collisions_b = 0;
	unsigned int *arr_c;

	if (MULT_BASE >= 0x100000000lu)
	{
		puts("Too large MULT_BASE");
		return 2;
	}
	log_1_5    = log(1.5);
	log_1_25   = log(1.25);
	log_1_125  = log(1.125);
	log_1_0625 = log(1.0625);
	init_prime();
	init_log();
	init_power();

	for (i=0; i<PRIMES_CNT; i++)
		prime_to_a_pow[i] = 1;
	a = 1;
	for (;;)
	{
		// get next a
		for (i=0; i<n_primes; i++)
		{
			if ((INT)a * prime[i] <= MAX_A)
			{
				a *= prime[i];
				prime_to_a_pow[i] *= prime[i];
				break;
			}
			a /= prime_to_a_pow[i];
			prime_to_a_pow[i] = 1;
		}
		if (i == n_primes)
			break;
		checked_a++;
		printf("a = %u (%lu/%u, %.1f%%)\n", a, checked_a, MAX_A, checked_a*100./MAX_A);
		// loop by b
		for (i=0; i<n_primes; i++)
			prime_to_b_pow[i] = prime_to_a_pow[i] == 1 ? 1 : 0;
		b = 1;
		for (;;)
		{
			for (i=0; i<n_primes; i++)
			{
				if (prime_to_b_pow[i] == 0)
					continue;
				if (prime[i] > a)
				{
					i = n_primes;
					break;
				}
				if ((INT)b * prime[i] <= MAX_B)
				{
					b *= prime[i];
					prime_to_b_pow[i] *= prime[i];
					break;
				}
				b /= prime_to_b_pow[i];
				prime_to_b_pow[i] = 1;
			}
			if (i == n_primes)
				break;
			if (b > a)
				continue;
			for (n = MIN_POW; n <= MAX_POW; n++)
			{
				power_base_a_n = power_base_n(a, n);
				an = n*logn[a];
				for (m = MIN_POW; m <= MAX_POW; m++)
				{
					if (a == b && m > n)
						continue;
					tested_b++;
					power_base_b_m = power_base_n(b, m);
					power_base_an_bm = power_base_a_n + power_base_b_m;
					if (power_base_an_bm >= MULT_BASE)
						power_base_an_bm -= MULT_BASE;
					if ((has_mod_n[power_base_an_bm >> 3] & (1 << (power_base_an_bm & 7))) == 0)
						continue;
					bm = m*logn[b];
					if (an > bm)
					{
						log_lower = an;
						log_upper = upper_log(an, bm);
					}
					else
					{
						log_lower = bm;
						log_upper = upper_log(bm, an);
					}
					checked_b++;
					arr_c = hash_base[power_base_an_bm % HASH_BASE];
					for (j=0; j<VAL_BY_HASH; j++)
					{
						c = arr_c[j];
						if (c == 0)
							break;
						check();
					}
					if (j < VAL_BY_HASH)
						continue;
					collisions_b++;
					for (c += 2; c<MAX_C; c+=2)
					{
						if (prime_to_b_pow[1] != 1 && c % 3 == 0)
							continue;
						check();
					}
				}
			}
		}
	}
	printf("Tested: %lu, checked %lu, collisions %lu\n", tested_b, checked_b, collisions_b);
}
