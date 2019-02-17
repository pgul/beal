#include <stdio.h>
#include <math.h>
// #include <gmp.h>

#define LOG_PRECISE 0.0001
#define MAX_BASE    1000
#define MULT_BASE   (8lu*9*5*7*11*13*17*19*23)

typedef unsigned long int INT;

double logn[MAX_BASE];

INT a, b, c, n, m;
INT power_base_a_n, power_base_b_m, power_base;
double an, bm, maxlog;

INT power_base_n(INT a, INT n)
{
	INT p = a % MULT_BASE;
	INT s = (n & 1) ? p : 1;
	while (n > 1) {
		p = (p * p) % MULT_BASE;
		n /= 2;
		if (n & 1)
			s = (s * p) % MULT_BASE;
	}
	return s;
}

INT check(void)
{
	double ck;
	INT k;

	if (logn[c] == 0)
		logn[c] = log(c);
	k = (INT)(maxlog / logn[c]) + 1;
	if (k < 3)
		return 0;
	ck = k*logn[c];
	if (maxlog < ck - logn[2] - LOG_PRECISE)
		return 0;
	if (power_base != power_base_n(c, k))
		return 0;
	// TODO: check with gmp.h
	return k;
}

int main(void)
{
	if (MULT_BASE >= 0x100000000lu) {
		puts("Too large MULT_BASE");
		return 2;
	}
	logn[2] = log(2);
	for (a=1; a<100; a++) {
		if (logn[a] == 0)
			logn[a] = log(a);
		for (b=1; b<=100; b+=2) {
			if (logn[b] == 0)
				logn[b] = log(b);
			for (n=3; n<100; n++) {
				power_base_a_n = power_base_n(a, n);
				an = n*logn[a];
				for (m=3; m<100; m++) {
					power_base_b_m = power_base_n(b, m);
					bm = m*logn[b];
					maxlog = an > bm ? an : bm;
					power_base = (power_base_a_n + power_base_b_m) % MULT_BASE;
					for (c = (a%2 ? 2 : 3); c<100; c+=2)
					{
						INT k;
						if ((k = check()) != 0)
							printf("%lu^%lu + %lu^%lu = %lu^%lu\n", a, n, b, m, c, k);
					}
				}
			}
		}
	}
}
