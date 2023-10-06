
#ifndef CIFER_UZPIPFE_27_H
#define CIFER_UZPIPFE_27_H

#include "cifer/data/mat.h"
#include "cifer/data/vec_curve.h"

typedef struct cfe_uzpipfe_27
{
    size_t m1;
    size_t m2;
    mpz_t bound_x;
    mpz_t bound_y;
    mpz_t order;
} cfe_uzpipfe_27;

cfe_error cfe_uzpipfe_init_27(cfe_uzpipfe_27 *c, size_t m1, size_t m2, mpz_t bound_x, mpz_t bound_y);

void cfe_uzpipfe_copy_27(cfe_uzpipfe_27 *res, cfe_uzpipfe_27 *c);

void cfe_uzpipfe_free_27(cfe_uzpipfe_27 *c);

typedef struct cfe_uzpipfe_pub_key_27
{
    ECP_BN254 g1;
    ECP2_BN254 g2;
    FP12_BN254 gT;
    cfe_vec_G1 b[4];
    cfe_vec_G1 b_tilde[4];
} cfe_uzpipfe_pub_key_27;

void cfe_uzpipfe_master_pub_key_init_27(cfe_uzpipfe_pub_key_27 *master_pub_key, cfe_uzpipfe_27 *c);

void cfe_uzpipfe_master_pub_key_free_27(cfe_uzpipfe_pub_key_27 *master_pub_key);

typedef struct cfe_uzpipfe_sec_key_27
{
    cfe_vec_G2 b_star[4];
    cfe_vec_G2 b_tilde_star[4];
} cfe_uzpipfe_sec_key_27;

void cfe_uzpipfe_master_sec_key_init_27(cfe_uzpipfe_sec_key_27 *master_sec_key, cfe_uzpipfe_27 *c);

void cfe_uzpipfe_master_sec_key_free_27(cfe_uzpipfe_sec_key_27 *master_sec_key);


cfe_error cfe_uzpipfe_generate_keys_27(cfe_uzpipfe_sec_key_27 *master_sec_key, cfe_uzpipfe_pub_key_27 *master_pub_key,
                                    cfe_uzpipfe_27 *c);

typedef struct cfe_uzpipfe_ciphertext_27
{
    cfe_vec_G1 *c1;
    cfe_vec_G1 *c2;
} cfe_uzpipfe_ciphertext_27;

void cfe_uzpipfe_ciphertext_init_27(cfe_uzpipfe_ciphertext_27 *cipher, cfe_uzpipfe_27 *c);

void cfe_uzpipfe_ciphertext_free_27(cfe_uzpipfe_ciphertext_27 *cipher, cfe_uzpipfe_27 *c);

cfe_error cfe_uzpipfe_encrypt_27(cfe_uzpipfe_ciphertext_27 *cipher, cfe_vec *x, cfe_vec *w,
                              cfe_uzpipfe_sec_key_27 *master_sec_key, cfe_uzpipfe_pub_key_27 *master_pub_key,
                              cfe_uzpipfe_27 *c);

typedef struct cfe_uzpipfe_fe_key_27
{
    cfe_vec_G2 *k1;
    cfe_vec_G2 *k2;
} cfe_uzpipfe_fe_key_27;

void cfe_uzpipfe_fe_key_init_27(cfe_uzpipfe_fe_key_27 *fe_key, cfe_uzpipfe_27 *c);

void cfe_uzpipfe_fe_key_free_27(cfe_uzpipfe_fe_key_27 *fe_key,  cfe_uzpipfe_27 *c);

cfe_error cfe_uzpipfe_derive_fe_key_27(cfe_uzpipfe_fe_key_27 *fe_key, cfe_vec *y, cfe_vec *v,
                                    cfe_uzpipfe_sec_key_27 *master_sec_key, cfe_uzpipfe_pub_key_27 *master_pub_key,
                                    cfe_uzpipfe_27 *c);


cfe_error cfe_uzpipfe_decrypt_27(mpz_t res, cfe_uzpipfe_ciphertext_27 *cipher,
                              cfe_uzpipfe_fe_key_27 *fe_key, cfe_uzpipfe_pub_key_27 *master_pub_key,
                              cfe_uzpipfe_27 *c);

#endif
