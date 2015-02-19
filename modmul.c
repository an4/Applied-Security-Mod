#include "modmul.h"
#include "assert.h"

#define WINDOW_SIZE 4
#define W 64

// Read value of variable in hex.
void initialize_and_read(mpz_t variable) {
    mpz_init(variable);
    gmp_scanf( "%Zx" , variable);
}

/** Sliding Window Exponentiation */
int binary_to_decimal(mpz_t input, int start, int end) {
    int i;
    int result = 0;
    int g = 1;
    for(i = end; i<=start; i++) {
        int bit = mpz_tstbit(input, i);
        if(bit) {
            result += g * bit;
        }
        g <<= 1;
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

/** Montgomery */

void montgomery_Mp(mpz_t output, mpz_t N, mpz_t BASE) {
    mpz_t temporary;
    mpz_init_set_ui(temporary, 1);
    while(mpz_cmp(temporary, N) <= 0) {
        mpz_mul(temporary, temporary, BASE);
    }
    mpz_init_set(output, temporary);
    mpz_clear(temporary);
}

void montgomery_Mw(mpz_t output, mpz_t N, mpz_t BASE) {
    mpz_t temporary;
    mpz_init_set_ui(temporary, 1);
    int i;
    for(i = 1; i < W-1; i++) {
        // t <- t*t*N (mod b)
        mpz_mul(temporary, temporary, temporary);
        mpz_mod(temporary, temporary, BASE);
        mpz_mul(temporary, temporary, N);
        mpz_mod(temporary, temporary, BASE);
    }
    // t = -t (mod b)
    mpz_neg(temporary, temporary);
    mpz_add(temporary, temporary, BASE);

    mpz_init_set(output, temporary);
    mpz_clear(temporary);
}

void montgomery_Mp_sq(mpz_t output, mpz_t N) {
    mpz_t temporary;
    mpz_init_set_ui(temporary, 1);
    int lN = mpz_size(N);
    int i;

    for(i = 1; i <= 2*lN*W; i++) {
        // t <- t + t (mod N)
        mpz_add(temporary, temporary, temporary);
        mpz_mod(temporary, temporary, N);
    }

    mpz_init_set(output, temporary);
    mpz_clear(temporary);
}

void montgomery_multiplication(mpz_t output,
                               mpz_t x,
                               mpz_t y,
                               mpz_t BASE,
                               mpz_t N,
                               mpz_t Mw) {
    mpz_t temporary;

    // r <- 0
    mpz_init(temporary);
    mpz_init_set_ui(temporary, 0);
    int i;
    int lN = mpz_size(N);

    mpz_t u;
    mpz_init(u);

    mpz_t t1, t2;
    mpz_init(t1);
    mpz_init(t2);

    mp_limb_t temp;
    mp_limb_t w = mpz_getlimbn(Mw, 0);

    for(i = 0; i < lN ; i++) {
        temp = mpz_getlimbn(temporary, 0) + mpz_getlimbn(y, i) * mpz_getlimbn(x, 0);
        temp = temp * w;

        mpz_set_ui(u, temp);

        mpz_mul_ui(t1, x, mpz_getlimbn(y, i));
        mpz_mul(t2, u, N);

        mpz_add(temporary, temporary, t1);
        mpz_add(temporary, temporary, t2);

        mpz_div(temporary, temporary, BASE);
    }

    if(mpz_cmp(temporary, N) > 0) {
        mpz_sub(temporary, temporary, N);
    }

    mpz_init_set(output, temporary);
    mpz_clear(temporary);
}

void montgomery(mpz_t output, mpz_t x, mpz_t y, mpz_t N) {
    // Calculate BASE value
    // BASE = 2^64
    mpz_t BASE;
    mpz_init_set_ui(BASE, 2);
    mpz_pow_ui(BASE, BASE, W);

    // Montgomery parameters
    mpz_t Mp;
    mpz_t Mw;
    mpz_t Mp_sq;

    mpz_t Mx;
    mpz_t My;

    mpz_t one;
    mpz_init_set_ui(one, 1);

    mpz_t temporary;

    // Compute Montgomery parameters.
    montgomery_Mp(Mp, N, BASE);
    montgomery_Mw(Mw, N, BASE);
    montgomery_Mp_sq(Mp_sq, N);

    montgomery_multiplication(Mx, x, Mp_sq, BASE, N, Mw);
    montgomery_multiplication(My, y, Mp_sq, BASE, N, Mw);
    montgomery_multiplication(temporary, Mx, My, BASE, N, Mw);
    montgomery_multiplication(temporary, temporary, one, BASE, N, Mw);

    mpz_init_set(output, temporary);

    mpz_clear(temporary);
    mpz_clear(Mp);
    mpz_clear(Mw);
    mpz_clear(Mp_sq);
    mpz_clear(Mx);
    mpz_clear(My);
    mpz_clear(one);
    mpz_clear(BASE);
}

/** Random */

void getSeed(gmp_randstate_t state) {
    // get value from /dev/urandom
    unsigned long long data[4];
    FILE* file;
    file = fopen("/dev/urandom", "r");
    int i;
    fread(&data, 4, sizeof(unsigned long long), file);
    fclose(file);

    // assign seed
    mpz_t seed;
    mpz_init(seed);
    mpz_import(seed, 4, 1, sizeof(data[0]), 0, 0, data);

    gmp_printf("S: %Zx\n", seed);

    // seed
    gmp_randseed(state, seed);
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

    while(!feof(stdin)) {
        initialize_and_read(N);
        initialize_and_read(e);
        initialize_and_read(m);

        if(mpz_cmp_ui(N, 0) == 0) {
            return;
        }

        mpz_init( c );

        sliding_window_exponentiation(c, m, e, N);

        gmp_printf( "%Zx\n", c );

        mpz_clear( N );
        mpz_clear( e );
        mpz_clear( m );
        mpz_clear( c );
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

    while(!feof(stdin)) {
        initialize_and_read(N);
        initialize_and_read(d);
        initialize_and_read(p);
        initialize_and_read(q);
        initialize_and_read(d_p);
        initialize_and_read(d_q);
        initialize_and_read(i_p);
        initialize_and_read(i_q);
        initialize_and_read(c);

        if(mpz_cmp_ui(N, 0) == 0) {
            return;
        }

        mpz_init(m);

        // mpz_powm(m, c, d, N);

        // CRT modular exponentiation

        mpz_t m_p;
        mpz_t m_q;
        mpz_init(m_p);
        mpz_init(m_q);

        // mpz_powm(m_p, c, d_p, p);
        sliding_window_exponentiation(m_p, c, d_p, p);
        // mpz_powm(m_q, c, d_q, q);
        sliding_window_exponentiation(m_q, c, d_q, q);

        // mpz_mul(m_p, m_p, q);
        // mpz_mod(m_p, m_p, N);
        montgomery(m_p, m_p, q, N);

        // mpz_mul(m_p, m_p, i_q);
        // mpz_mod(m_p, m_p, N);
        montgomery(m_p, m_p, i_q, N);

        // mpz_mul(m_q, m_q, p);
        // mpz_mod(m_q, m_q, N);
        montgomery(m_q, m_q, p, N);

        // mpz_mul(m_q, m_q, i_p);
        // mpz_mod(m_q, m_q, N);
        montgomery(m_q, m_q, i_p, N);

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

    // State
    gmp_randstate_t state;
    gmp_randinit_default(state);

    getSeed(state);

    while(!feof(stdin)) {
        initialize_and_read(p);
        initialize_and_read(q);
        initialize_and_read(g);
        initialize_and_read(h);
        initialize_and_read(m);

        if(mpz_cmp_ui(p, 0) == 0) {
            return;
        }

        mpz_init(c1);
        mpz_init(c2);
        mpz_init(k);
        mpz_set_ui(k, 0);

        // Random
        while(mpz_cmp_ui(k, 0) == 0){
            mpz_urandomm(k, state, q);
        }

        gmp_printf("Random k: %Zx\n", k);

        // mpz_powm(c1, g, k, p);
        sliding_window_exponentiation(c1, g, k, p);

        //mpz_powm(h, h, k, p);
        sliding_window_exponentiation(h, h, k, p);

        // mpz_mul(c2, m, h);
        // mpz_mod(c2, c2, p);
        montgomery(c2, m, h, p);

        gmp_printf( "%Zx\n", c1 );
        gmp_printf( "%Zx\n", c2 );

        mpz_clear(p);
        mpz_clear(q);
        mpz_clear(g);
        mpz_clear(h);
        mpz_clear(m);
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

    while(!feof(stdin)) {
        initialize_and_read(p);
        initialize_and_read(q);
        initialize_and_read(g);
        initialize_and_read(x);
        initialize_and_read(c1);
        initialize_and_read(c2);

        if(mpz_cmp_ui(p, 0) == 0) {
            return;
        }

        mpz_init(m);

        // -x (mod q)
        mpz_mod(x, x, q);
        mpz_neg(x, x);
        mpz_add(x, x, q);

        //mpz_powm(c1, c1, x, p);
        sliding_window_exponentiation(c1, c1, x, p);

        // mpz_mul(m, c1, c2);
        // mpz_mod(m, m, p);
        montgomery(m, c1, c2, p);

        gmp_printf( "%Zx\n", m );

        mpz_clear(p);
        mpz_clear(q);
        mpz_clear(g);
        mpz_clear(x);
        mpz_clear(c1);
        mpz_clear(c2);
        mpz_clear(m);
    }
}

/** TEST */

void Test_montgomery() {
    mpz_t Mp;
	mpz_t Mp_test;
    mpz_t N;
    mpz_t BASE;
	mpz_t Mw;
	mpz_t Mw_test;
	mpz_t Mp_sq;
	mpz_t Mp_sq_test;
	mpz_t x;
    mpz_t x_test;
    mpz_t Mx;
    mpz_t one;
    mpz_t Nx;
    mpz_t y;
    mpz_t r_test;
    mpz_t r;

	// Base 2^64
	// 18446744073709551616
	initialize_and_read(BASE);
	// 18446744073709551629
	initialize_and_read(N);
	// 0x100000000000000000000000000000000
	initialize_and_read(Mp_test);
	// 2^60
	initialize_and_read(x);
    // 2^50-1
    initialize_and_read(y);

	// Test montgomery rho
    montgomery_Mp(Mp, N, BASE);
    assert(mpz_cmp(Mp, Mp_test) == 0);

	// Generate test omega
	mpz_init(Mw_test);
	mpz_invert(Mw_test, N, BASE);
	mpz_neg(Mw_test, Mw_test);
	mpz_add(Mw_test, Mw_test, BASE);
	// 0xb13b13b13b13b13b

	// Test montgomery omega
	montgomery_Mw(Mw, N, BASE);
	assert(mpz_cmp(Mw, Mw_test) == 0);

	// Generate test rho squared
	mpz_init(Mp_sq_test);
	mpz_mul(Mp_sq_test, Mp, Mp);
	mpz_mod(Mp_sq_test, Mp_sq_test, N);
	// 6f91

	// Test montgomery rho squared
	montgomery_Mp_sq(Mp_sq, N);
	assert(mpz_cmp(Mp_sq, Mp_sq_test) == 0);

	// Generate x in montgomery form
    // x = x * p mod N
    mpz_init(x_test);
    mpz_mul(x_test, x, Mp);
    mpz_mod(x_test, x_test, N);
    // 8fffffffffffff7e

    // Test montgomery form
    montgomery_multiplication(Mx, x, Mp_sq, BASE, N, Mw);
    assert(mpz_cmp(Mx, x_test) == 0);

    // Test normal form from montgomery form
    mpz_init_set_ui(one, 1);
    montgomery_multiplication(Nx, Mx, one, BASE, N, Mw);
    assert(mpz_cmp(Nx, x) == 0);

    // Modular multiplication
    mpz_init(r_test);
    mpz_mul(r_test, x, y);
    mpz_mod(r_test, r_test, N);

    // Montgomery modular multiplication
    montgomery(r, x, y, N);
    assert(mpz_cmp(r, r_test) == 0);

    mpz_clear(r);
    mpz_clear(y);
    mpz_clear(r_test);
    mpz_clear(Nx);
    mpz_clear(one);
    mpz_clear(Mp);
    mpz_clear(Mp_test);
    mpz_clear(N);
    mpz_clear(BASE);
    mpz_clear(Mw);
    mpz_clear(Mw_test);
    mpz_clear(Mp_sq);
    mpz_clear(Mp_sq_test);
    mpz_clear(x);
    mpz_clear(x_test);
    mpz_clear(Mx);
}

void Test_exponentiation() {
    mpz_t base;
    mpz_t exp;
    mpz_t N;
    mpz_t result;
    mpz_t result_test;

    initialize_and_read(base);
    initialize_and_read(exp);
    initialize_and_read(N);

    // Modular exponentiation
    mpz_init(result_test);
    mpz_powm(result_test, base, exp, N);

    // Sliding window modular exponentiation
    sliding_window_exponentiation(result, base, exp, N);

    assert(mpz_cmp(result, result_test) == 0);

    mpz_clear(base);
    mpz_clear(exp);
    mpz_clear(N);
    mpz_clear(result);
}

void Test() {
    Test_montgomery();
    Test_exponentiation();
}

/** END TEST */

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
	else if( !strcmp( argv[ 1 ], "test" ) ) {
    Test();
    }
    else {
    abort();
    }

    return 0;
}
