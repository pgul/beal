#include <stdio.h>
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
#define MULT_BASE   (8lu*9*5*7*11*13*17*19*23)
#define PRIMES_CNT	1000

typedef unsigned long int INT;

int n_primes;
INT prime[PRIMES_CNT];
INT prime_to_a_pow[PRIMES_CNT];
INT prime_to_b_pow[PRIMES_CNT];
INT prime_to_c_pow[PRIMES_CNT];

double logn[MAX_ABC];
INT power_mod_n[MAX_BC][MAX_POW+1];

INT a, b, c, n, m;
INT power_base_a_n, power_base_b_m, power_base_an_bm;
double an, bm, log_lower, log_upper;

INT power_base_n(INT a, INT n)
{
	INT p, s;
	if (a < MAX_BC && n <= MAX_POW)
		return power_mod_n[a][n];
	p = a % MULT_BASE;
	s = (n & 1) ? p : 1;
	while (n > 1) {
		p = (p * p) % MULT_BASE;
		n /= 2;
		if (n & 1)
			s = (s * p) % MULT_BASE;
	}
	return s;
}

int check_bigint(INT k)
{
	mpz_t man, mbm, mck, manbm;
	mpz_init(man);
	mpz_ui_pow_ui(man, a, n);
	mpz_init(mbm);
	mpz_ui_pow_ui(mbm, b, m);
	mpz_init(mck);
	mpz_ui_pow_ui(mck, c, k);
	mpz_init(manbm);
	mpz_add(manbm, man, mbm);
	return mpz_cmp(manbm, mck);
}

INT check(void)
{
	double logc;
	INT k;

	logc = logn[c];
	k = (INT)(log_lower / logc) + 1;
	if (k < MIN_POW)
		return 0;
	if (k * logc > log_upper)
		return 0;
	if (power_base_an_bm != power_base_n(c, k))
		return 0;
	if (check_bigint(k)) {
		printf("  - %lu^%lu + %lu^%lu = %lu^%lu\n", a, n, b, m, c, k);
		return 0;
	}
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
	for (i=3; i<MAX_ABC && n_primes<PRIMES_CNT; i+=2) {
		for (j=1; j<n_primes; j++) {
			if (i % prime[j] == 0)
				break;
		}
		if (j == n_primes)
			prime[n_primes++] = i;
	}
	if (n_primes == PRIMES_CNT)
		fprintf(stderr, "Increase PRIMES_CNT\n");
}

void init_log(void)
{
	int i;
	for (i=1; i<MAX_ABC; i++)
		logn[i] = log(i);
}

void init_power(void)
{
	int i, j;
	for (i=1; i<MAX_BC; i++) {
		power_mod_n[i][0] = 1;
		for (j=1; j<=MAX_POW; j++)
			power_mod_n[i][j] = (power_mod_n[i][j-1] * i) % MULT_BASE;
	}
}

int main(void)
{
	int i;

	if (MULT_BASE >= 0x100000000lu) {
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
	for (;;) {
		// get next a
		for (i=0; i<n_primes; i++) {
			if (a * prime[i] <= MAX_A) {
				a *= prime[i];
				prime_to_a_pow[i] *= prime[i];
				break;
			}
			a /= prime_to_a_pow[i];
			prime_to_a_pow[i] = 1;
		}
		if (i == n_primes)
			break;
		printf("a = %lu\n", a);
		// loop by b
		for (i=0; i<n_primes; i++)
			prime_to_b_pow[i] = prime_to_a_pow[i] == 1 ? 1 : 0;
		b = 1;
		for (;;) {
			for (i=0; i<n_primes; i++) {
				if (prime_to_b_pow[i] == 0)
					continue;
				if (b * prime[i] <= MAX_B) {
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
			for (n = MIN_POW; n <= MAX_POW; n++) {
				power_base_a_n = power_base_n(a, n);
				an = n*logn[a];
				for (m = MIN_POW; m <= MAX_POW; m++) {
					if (a == b && m > n)
						continue;
					power_base_b_m = power_base_n(b, m);
					bm = m*logn[b];
					// log(a+b) = log(a+a*k) = log(a) + log(1+k), k = b/a = exp(log(b) - log(a))
					if (an > bm) {
						log_lower = an;
						log_upper = upper_log(an, bm);
					}
					else {
						log_lower = bm;
						log_upper = upper_log(bm, an);
					}
					power_base_an_bm = (power_base_a_n + power_base_b_m) % MULT_BASE;
#if 0
					// loop by c
					for (i=0; i<n_primes; i++)
						prime_to_c_pow[i] = (prime_to_a_pow[i] == 1 && prime_to_b_pow[i] == 1) ? 1 : 0;
					c = 1;
					for (;;) {
						for (i=0; i<n_primes && prime[i] < MAX_C; i++) {
							if (prime_to_c_pow[i] == 0)
								continue;
							if (c * prime[i] <= MAX_C) {
								c *= prime[i];
								prime_to_c_pow[i] *= prime[i];
								break;
							}
							c /= prime_to_c_pow[i];
							prime_to_c_pow[i] = 1;
						}
						if (i == n_primes || prime[i] >= MAX_C)
							break;
						if (prime_to_c_pow[0] != 1) {
							// at least one of a, b, c must be even
#else
					for (c=((a+b)%2 ? 1 : 2); c<MAX_C; c+=2) {
						{
#endif
							INT k;
							if (prime_to_b_pow[1] != 1 && c % 3 == 0)
								continue;
							if ((k = check()) != 0)
								printf("%lu^%lu + %lu^%lu = %lu^%lu\n", a, n, b, m, c, k);
						}
					}
				}
			}
		}
	}
}
