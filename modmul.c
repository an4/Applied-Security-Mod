#include "modmul.h"

#define WINDOW_SIZE 4

// Read value of variable in hex.
void initialize_and_read(mpz_t variable) {
	mpz_init(variable);
	gmp_scanf( "%Zx" , variable);
}

int binary_to_decimal(mpz_t input, int start, int end) {
	int i;
	int result = 0;
	int g = 1;
	for(i = end; i<=start; i++) {
		int bit = mpz_tstbit(input, i);
		if(bit) {
			result += g * bit;
		}
		g *= 2;
	}
	return result;
}

void sliding_window_exponentiation(mpz_t output, mpz_t base, mpz_t exp, mpz_t modulus) {
	const int table_length = 1 << (WINDOW_SIZE - 1);
	mpz_t table[table_length];
	mpz_t square_base;
	mpz_t temp;
	int i;
	int exp_length;
	int l;
	// value in window 
	int u;

	mpz_init(temp);
	mpz_init(square_base);

	// T[0] = x
	mpz_init_set(table[0], base);
	// x^2 mod N
	mpz_mul(square_base, base, base);
	mpz_mod(square_base, square_base, modulus);
	// Compute lookup table.
	for(i = 1; i < table_length; i++) {
		mpz_init(table[i]);
		mpz_mul(table[i], table[i-1], square_base);
		mpz_mod(table[i], table[i], modulus);
	}
	// t = 1
	mpz_set_ui(output, 1);
	// |y|
	exp_length = mpz_size(exp) * mp_bits_per_limb;
	// |y| - 1
	i = exp_length - 1;
	while(i >= 0) {
		if(!mpz_tstbit(exp, i)) {
			l = i;
			u = 0;
		} else {
			l = ((i - WINDOW_SIZE + 1) > 0) ? (i - WINDOW_SIZE + 1) : 0;
			while(!mpz_tstbit(exp,l)) {
				l++;
			}
			// Set u = exp bits between i and l
			u = binary_to_decimal(exp, i, l);
		}

		// t' = t
		mpz_set(temp, output);		
		// t^p
		int p = 1 << (i - l + 1);
		if(p>0) {
			// t = t'^p (mod N)
			mpz_powm_ui(output, temp, p, modulus);
		}
		if(u != 0) {
			// t = t * T[(u-1)/2] mod N
			mpz_mul(output, output, table[(u-1)/2]);
			mpz_mod(output, output, modulus);
		}
		i = l - 1;
	}

	mpz_clear(square_base);
	mpz_clear(temp);
	for(i=0 ; i<table_length; i++) {
		mpz_clear(table[i]);
	}
}

/*
Perform stage 1:

- read each 3-tuple of N, e and m from stdin,
- compute the RSA encryption c,
- then write the ciphertext c to stdout.
*/

void stage1() {
  	// fill in this function with solution
	mpz_t N, e, m, c;

	initialize_and_read(N);

	while(!feof(stdin)) {

		initialize_and_read(e);
		initialize_and_read(m);

		mpz_init( c );

		sliding_window_exponentiation(c, m, e, N);

		gmp_printf( "%Zx\n", c );

		mpz_clear( N );
		mpz_clear( e );
		mpz_clear( m );
		mpz_clear( c );

		initialize_and_read(N);
  	}
}

/*
Perform stage 2:

- read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
- compute the RSA decryption m,
- then write the plaintext m to stdout.
*/

void stage2() {
	// fill in this function with solution
	mpz_t N, d, p, q, d_p, d_q, i_p, i_q, c, m;

	initialize_and_read(N);

	while(!feof(stdin)) {
		initialize_and_read(d);
		initialize_and_read(p);
		initialize_and_read(q);
		initialize_and_read(d_p);
		initialize_and_read(d_q);
		initialize_and_read(i_p);
		initialize_and_read(i_q);
		initialize_and_read(c);

		mpz_init(m);

		// mpz_powm(m, c, d, N);

		// CRT modular exponentiation

		mpz_t m_p;
		mpz_t m_q;
		mpz_init(m_p);
		mpz_init(m_q);

		mpz_powm(m_p, c, d_p, p);
		mpz_powm(m_q, c, d_q, q);

		mpz_mul(m_p, m_p, q);
		mpz_mod(m_p, m_p, N);
		mpz_mul(m_p, m_p, i_q);
		mpz_mod(m_p, m_p, N);

		mpz_mul(m_q, m_q, p);
		mpz_mod(m_q, m_q, N);
		mpz_mul(m_q, m_q, i_p);
		mpz_mod(m_q, m_q, N);		

		mpz_add(m, m_p, m_q);
		mpz_mod(m, m, N);

		gmp_printf( "%Zx\n", m );

		mpz_clear(N);
		mpz_clear(d);
		mpz_clear(p);
		mpz_clear(q);
		mpz_clear(d_p);
		mpz_clear(d_q);
		mpz_clear(i_p);
		mpz_clear(i_q);
		mpz_clear(c);
		mpz_clear(m_p);
		mpz_clear(m_q);

		initialize_and_read(N);
	}
}

/*
Perform stage 3:

- read each 5-tuple of p, q, g, h and m from stdin,
- compute the ElGamal encryption c = (c_1,c_2),
- then write the ciphertext c to stdout.
*/

void stage3() {

	// fill in this function with solution
	mpz_t p, q, g, h, m, c1, c2, k;

	initialize_and_read(p);

	while(!feof(stdin)) {

		initialize_and_read(q);
		initialize_and_read(g);
		initialize_and_read(h);
		initialize_and_read(m);

		mpz_init(c1);
		mpz_init(c2);
		mpz_init(k);

		int i_k = 1;
		mpz_set_ui(k, i_k);

		mpz_mod(k, k, q);
		mpz_powm(c1, g, k, p);
		mpz_pow_ui(h, h, i_k);
		mpz_mul(c2, m, h);
		mpz_mod(c2, c2, p);

		gmp_printf( "%Zx\n", c1 );
		gmp_printf( "%Zx\n", c2 );

		mpz_clear(p);
		mpz_clear(q);
		mpz_clear(g);
		mpz_clear(h);
		mpz_clear(m);

		initialize_and_read(p);

	}

}

/*
Perform stage 4:

- read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
- compute the ElGamal decryption m,
- then write the plaintext m to stdout.
*/

void stage4() {

 	// fill in this function with solution
	mpz_t p, q, g, x, c1, c2, m;

	initialize_and_read(p);

	while(!feof(stdin)) {
		
		initialize_and_read(q);
		initialize_and_read(g);
		initialize_and_read(x);
		initialize_and_read(c1);
		initialize_and_read(c2);

		mpz_init(m);

		mpz_mod(x, x, q);
		mpz_neg(x, x);
		mpz_add(x, x, q);

		mpz_powm(c1, c1, x, p);
		mpz_mul(m, c1, c2);
		mpz_mod(m, m, p);

		gmp_printf( "%Zx\n", m );

		mpz_clear(p);
		mpz_clear(q);
		mpz_clear(g);
		mpz_clear(x);
		mpz_clear(c1);
		mpz_clear(c2);
		mpz_clear(m);

		initialize_and_read(p);

	}

}

/*
The main function acts as a driver for the assignment by simply invoking
the correct function for the requested stage.
*/

int main( int argc, char* argv[] ) {
  if( argc != 2 ) {
    abort();
  }

  if     ( !strcmp( argv[ 1 ], "stage1" ) ) {
    stage1();
  }
  else if( !strcmp( argv[ 1 ], "stage2" ) ) {
    stage2();
  }
  else if( !strcmp( argv[ 1 ], "stage3" ) ) {
    stage3();
  }
  else if( !strcmp( argv[ 1 ], "stage4" ) ) {
    stage4();
  }
  else {
    abort();
  }

  return 0;
}
