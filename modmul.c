#include "modmul.h"

// Read value of variable in hex.
void initialize_and_read(mpz_t variable) {
	mpz_init(variable);
	gmp_scanf( "%Zx" , variable);
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

	// mpz_init( N );
	// gmp_scanf( "%Zx",  N );

	initialize_and_read(N);

	while(!feof(stdin)) {

		initialize_and_read(e);
		initialize_and_read(m);

		mpz_init( c );

		mpz_powm(c, m, e, N);

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

		mpz_powm(m, c, d, N);

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
