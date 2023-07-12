/*
 * Copyright (c) 2018 XLAB d.o.o.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <string.h>

#include <amcl/pair_BN254.h>

#include "cifer/innerprod/fullysec/uzpipfe.h"
#include "cifer/internal/big.h"
#include "cifer/sample/uniform.h"
#include "cifer/internal/dlog.h"

cfe_error cfe_uzpipfe_init(cfe_uzpipfe *c, size_t m1, size_t m2, mpz_t bound_x, mpz_t bound_y)
{
    cfe_error err = CFE_ERR_NONE;
    mpz_t check, order;
    mpz_inits(check, order, NULL);
    mpz_from_BIG_256_56(order, (int64_t *)CURVE_Order_BN254);

    mpz_mul(check, bound_x, bound_y);
    mpz_mul_ui(check, check, m1);
    mpz_mul_ui(check, check, m2);

    if (mpz_cmp(check, order) >= 0)
    {
        err = CFE_ERR_PRECONDITION_FAILED;
        goto cleanup;
    }

    mpz_inits(c->bound_x, c->bound_y, c->order, NULL);

    c->m1 = m1;
    c->m2 = m2;
    mpz_set(c->bound_x, bound_x);
    mpz_set(c->bound_y, bound_y);
    mpz_set(c->order, order);

cleanup:
    mpz_clears(check, order, NULL);
    return err;
}

void cfe_uzpipfe_copy(cfe_uzpipfe *res, cfe_uzpipfe *c)
{
    mpz_inits(res->bound_x, res->bound_y, res->order, NULL);

    res->m1 = c->m1;
    res->m2 = c->m2;
    mpz_set(res->bound_x, c->bound_x);
    mpz_set(res->bound_y, c->bound_y);
    mpz_set(res->order, c->order);
}

void cfe_uzpipfe_free(cfe_uzpipfe *c)
{
    mpz_clears(c->bound_x, c->bound_y, c->order, NULL);
}

void cfe_uzpipfe_master_key_init(cfe_uzpipfe_sec_key *sec_key, cfe_uzpipfe *c)
{
    // nothing to do here
    ;
}

void cfe_uzpipfe_master_key_free(cfe_uzpipfe_sec_key *sec_key)
{
    // nothing to do here
    ;
}

cfe_error cfe_uzpipfe_generate_keys(cfe_uzpipfe_sec_key *sec_key, cfe_uzpipfe_pub_key *pub_key,
                                    cfe_uzpipfe *c)
{
    cfe_error err = CFE_ERR_NONE;

    mpz_t exp, det;
    BIG_256_56 exp_big;

    unsigned char key[32] = {
        0x39, 0x39, 0x35, 0x39, 0x38, 0x33, 0x37, 0x63,
        0x62, 0x66, 0x65, 0x65, 0x34, 0x36, 0x37, 0x62,
        0x32, 0x62, 0x34, 0x64, 0x62, 0x38, 0x63, 0x33,
        0x37, 0x35, 0x30, 0x66, 0x61, 0x66, 0x33, 0x61};

    unsigned char key_tilde[32] = {
        0x34, 0x34, 0x38, 0x38, 0x63, 0x65, 0x63, 0x62,
        0x35, 0x37, 0x34, 0x35, 0x65, 0x37, 0x31, 0x37,
        0x33, 0x30, 0x66, 0x33, 0x33, 0x63, 0x61, 0x65,
        0x64, 0x32, 0x63, 0x35, 0x61, 0x32, 0x36, 0x64};

    strncpy((char *)sec_key->key, (char *)key, 32);
    strncpy((char *)sec_key->key_tilde, (char *)key_tilde, 32);

    mpz_inits(exp, det, NULL);

    cfe_uniform_sample(exp, c->order);
    BIG_256_56_from_mpz(exp_big, exp);
    ECP_BN254_generator(&(pub_key->g1));
    ECP_BN254_mul(&(pub_key->g1), exp_big);

    cfe_uniform_sample(exp, c->order);
    BIG_256_56_from_mpz(exp_big, exp);
    ECP2_BN254_generator(&(pub_key->g2));
    ECP2_BN254_mul(&(pub_key->g2), exp_big);

    PAIR_BN254_ate(&(pub_key->gT), &(pub_key->g2), &(pub_key->g1));
    PAIR_BN254_fexp(&(pub_key->gT));

    return err;
}

void cfe_uzpipfe_fe_key_init(cfe_uzpipfe_fe_key *fe_key, cfe_uzpipfe *c)
{
    fe_key->k1 = malloc(c->m1 * sizeof(cfe_vec_G2));
    for (size_t i = 0; i < c->m1; i++)
    {
        cfe_vec_G2_init(&fe_key->k1[i], 4);
    }

    fe_key->k2 = malloc(c->m2 * sizeof(cfe_vec_G2));
    for (size_t i = 0; i < c->m2; i++)
    {
        cfe_vec_G2_init(&fe_key->k2[i], 4);
    }
}

void cfe_uzpipfe_fe_key_free(cfe_uzpipfe_fe_key *fe_key, cfe_uzpipfe *c)
{
    for (size_t i = 0; i < c->m1; i++)
    {
        cfe_vec_G2_free(&fe_key->k1[i]);
    }
    for (size_t i = 0; i < c->m2; i++)
    {
        cfe_vec_G2_free(&fe_key->k2[i]);
    }
}

cfe_error cfe_uzpipfe_derive_fe_key(cfe_uzpipfe_fe_key *fe_key, cfe_vec *y, cfe_vec *v,
                                    cfe_uzpipfe_sec_key *sec_key, cfe_uzpipfe_pub_key *pub_key,
                                    cfe_uzpipfe *c)
{
    cfe_error err = CFE_ERR_NONE;

    mpz_t zero, omega, gamma, gamma_tilde, det, el, sum, bound, bound_neg;
    cfe_vec gamma_v, gamma_tilde_v, k_tmp, B_y;
    cfe_mat B, B_t, B_inv, B_star, B_t_star;
    
    unsigned char key[32], key_tilde[32];

    mpz_inits(zero, omega, gamma, gamma_tilde, det, el, sum, bound, bound_neg, NULL);
    cfe_vec_inits(4, &k_tmp, &B_y, NULL);
    cfe_vec_init(&gamma_v, c->m1);
    cfe_vec_init(&gamma_tilde_v, c->m2);
    cfe_mat_inits(4, 4, &B, &B_inv, &B_star, &B_t, &B_t_star, NULL);

    strncpy((char *) key, (char *) sec_key->key, 32);
    strncpy((char *) key_tilde, (char *) sec_key->key_tilde, 32);

    mpz_set(bound, c->order);
    mpz_neg(bound_neg, bound);

    mpz_set_ui(zero, 0);
    cfe_uniform_sample(omega, c->order);
    
    cfe_uniform_sample_range_vec(&gamma_v, bound_neg, bound);
    cfe_uniform_sample_range_vec(&gamma_tilde_v, bound_neg, bound);
    cfe_vec_set(&gamma_tilde_v, zero, ((c->m2)-1));

    mpz_set(sum, zero);
    for(size_t i = 0; i < c->m1; i++){
        cfe_vec_get(el, &gamma_v, i);
        mpz_add(sum, sum, el);
    }
    for(size_t i = 0; i < c->m2; i++){
        cfe_vec_get(el, &gamma_tilde_v, i);
        mpz_add(sum, sum, el);
    }

    mpz_neg(el, sum);
    cfe_vec_set(&gamma_tilde_v, el, ((c->m2)-1));

    cfe_vec_G2 vec_g2;
    cfe_vec_G2_init(&vec_g2, 4);
    for (size_t i = 0; i < 4; i++)
    {
        ECP2_BN254_copy(&vec_g2.vec[i], &pub_key->g2);
    }

    for (size_t i = 0; i < c->m1; i++)
    {
        int first_it = 1;

        cfe_vec_set_const(&k_tmp, zero);
        cfe_vec_get(el, y, i);
        cfe_vec_set(&k_tmp, el, 0);
        cfe_vec_get(el, &gamma_v, i);
        cfe_vec_set(&k_tmp, el, 2);

        do
        {
            if (first_it)
            {
                key[31] += 1;
                first_it = 0;
            }
            else if (err != CFE_ERR_NONE)
            {
                printf("Inverse not found!\n");
                key[31] += 1;
            }

            cfe_uniform_sample_mat_det(&B, c->order, key);
            err = cfe_mat_inverse_mod_gauss(&B_inv, det, &B, c->order);
        } while (err != CFE_ERR_NONE);
        cfe_mat_transpose(&B_star, &B_inv);

        cfe_mat_mul_vec(&B_y, &B_star, &k_tmp);
        cfe_vec_mod(&B_y, &B_y, c->order);

        cfe_vec_mul_vec_G2(&fe_key->k1[i], &B_y, &vec_g2);
    }

    for (size_t i = 0; i < c->m2; i++)
    {
        int first_it = 1;

        cfe_vec_set_const(&k_tmp, zero);
        cfe_vec_get(el, v, i);
        mpz_mul(el, el, omega);
        cfe_vec_set(&k_tmp, el, 0);
        cfe_vec_get(el, &gamma_tilde_v, i);
        cfe_vec_set(&k_tmp, el, 2);

        do
        {
            if (first_it)
            {
                key_tilde[31] += 1;
                first_it = 0;
            }
            else if (err != CFE_ERR_NONE)
            {
                printf("Inverse not found!\n");
                key_tilde[31] += 1;
            }

            cfe_uniform_sample_mat_det(&B_t, c->order, key_tilde);
            err = cfe_mat_inverse_mod_gauss(&B_inv, det, &B_t, c->order);
        } while (err != CFE_ERR_NONE);
        cfe_mat_transpose(&B_t_star, &B_inv);

        cfe_mat_mul_vec(&B_y, &B_t_star, &k_tmp);
        cfe_vec_mod(&B_y, &B_y, c->order);

        cfe_vec_mul_vec_G2(&fe_key->k2[i], &B_y, &vec_g2);
    }

    mpz_clears(zero, omega, gamma, gamma_tilde, el, NULL);
    cfe_vec_frees(&k_tmp, &B_y, NULL);
    cfe_vec_G2_free(&vec_g2);
    cfe_mat_frees(&B, &B_inv, &B_star, &B_t, &B_t_star, NULL);

    return CFE_ERR_NONE;
}

void cfe_uzpipfe_ciphertext_init(cfe_uzpipfe_ciphertext *cipher, cfe_uzpipfe *c)
{
    cipher->c1 = malloc(c->m1 * sizeof(cfe_vec_G1));
    for (size_t i = 0; i < c->m1; i++)
    {
        cfe_vec_G1_init(&cipher->c1[i], 4);
    }

    cipher->c2 = malloc(c->m2 * sizeof(cfe_vec_G1));
    for (size_t i = 0; i < c->m2; i++)
    {
        cfe_vec_G1_init(&cipher->c2[i], 4);
    }
}

void cfe_uzpipfe_ciphertext_free(cfe_uzpipfe_ciphertext *cipher, cfe_uzpipfe *c)
{
    for (size_t i = 0; i < c->m1; i++)
    {
        cfe_vec_G1_free(&cipher->c1[i]);
    }
    for (size_t i = 0; i < c->m2; i++)
    {
        cfe_vec_G1_free(&cipher->c2[i]);
    }
}

cfe_error cfe_uzpipfe_encrypt(cfe_uzpipfe_ciphertext *cipher, cfe_vec *x, cfe_vec *w,
                              cfe_uzpipfe_sec_key *sec_key, cfe_uzpipfe_pub_key *pub_key,
                              cfe_uzpipfe *c)
{
    cfe_error err = CFE_ERR_NONE;

    mpz_t zero, alpha, delta, det, el;
    cfe_vec c_tmp, B_x;
    cfe_mat B, B_t, B_inv;
    
    unsigned char key[32], key_tilde[32];

    mpz_inits(zero, alpha, delta, el, det, NULL);
    cfe_vec_inits(4, &c_tmp, &B_x, NULL);
    cfe_mat_inits(4, 4, &B, &B_inv, &B_t, NULL);

    strncpy((char *) key, (char *) sec_key->key, 32);
    strncpy((char *) key_tilde, (char *) sec_key->key_tilde, 32);

    mpz_set_ui(zero, 0);
    cfe_uniform_sample(alpha, c->order);
    cfe_uniform_sample(delta, c->order);

    cfe_vec_G1 vec_g1;
    cfe_vec_G1_init(&vec_g1, 4);
    for (size_t i = 0; i < 4; i++)
    {
        ECP_BN254_copy(&vec_g1.vec[i], &pub_key->g1);
    }

    for (size_t i = 0; i < c->m1; i++)
    {
        int first_it = 1;

        cfe_vec_set_const(&c_tmp, zero);
        cfe_vec_get(el, x, i);
        cfe_vec_set(&c_tmp, el, 0);
        cfe_vec_set(&c_tmp, alpha, 2);

        do
        {
            if (first_it)
            {
                key[31] += 1;
                first_it = 0;
            }
            else if (err != CFE_ERR_NONE)
            {
                printf("Inverse not found!\n");
                key[31] += 1;
            }

            cfe_uniform_sample_mat_det(&B, c->order, key);
            err = cfe_mat_inverse_mod_gauss(&B_inv, det, &B, c->order);
        } while (err != CFE_ERR_NONE);

        cfe_mat_mul_vec(&B_x, &B, &c_tmp);
        cfe_vec_mod(&B_x, &B_x, c->order);

        cfe_vec_mul_vec_G1(&cipher->c1[i], &B_x, &vec_g1);
    }

    for (size_t i = 0; i < c->m2; i++)
    {
        int first_it = 1;

        cfe_vec_set_const(&c_tmp, zero);
        cfe_vec_get(el, w, i);
        mpz_mul(el, el, delta);
        cfe_vec_set(&c_tmp, el, 0);
        cfe_vec_set(&c_tmp, alpha, 2);

        do
        {
            if (first_it)
            {
                key_tilde[31] += 1;
                first_it = 0;
            }
            else if (err != CFE_ERR_NONE)
            {
                printf("Inverse not found!\n");
                key_tilde[31] += 1;
            }

            cfe_uniform_sample_mat_det(&B_t, c->order, key_tilde);
            err = cfe_mat_inverse_mod_gauss(&B_inv, det, &B_t, c->order);
        } while (err != CFE_ERR_NONE);

        cfe_mat_mul_vec(&B_x, &B_t, &c_tmp);
        cfe_vec_mod(&B_x, &B_x, c->order);

        cfe_vec_mul_vec_G1(&cipher->c2[i], &B_x, &vec_g1);
    }

    mpz_clears(zero, alpha, delta, el, NULL);
    cfe_vec_frees(&c_tmp, &B_x, NULL);
    cfe_vec_G1_free(&vec_g1);
    cfe_mat_frees(&B, &B_inv, &B_t, NULL);

    return CFE_ERR_NONE;
}

cfe_error cfe_uzpipfe_decrypt(mpz_t res, cfe_uzpipfe_ciphertext *cipher,
                              cfe_uzpipfe_fe_key *fe_key, cfe_uzpipfe_pub_key *pub_key,
                              cfe_uzpipfe *c)
{
    FP12_BN254 h, g;
    FP12_BN254 paired;

    FP12_BN254_copy(&g, &pub_key->gT);
    FP12_BN254_one(&h);

    for (size_t i = 0; i < c->m1; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            PAIR_BN254_ate(&paired, &fe_key->k1[i].vec[j], &cipher->c1[i].vec[j]);
            PAIR_BN254_fexp(&paired);
            FP12_BN254_mul(&h, &paired);
        }
    }

    for (size_t i = 0; i < c->m2; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            PAIR_BN254_ate(&paired, &fe_key->k2[i].vec[j], &cipher->c2[i].vec[j]);
            PAIR_BN254_fexp(&paired);
            FP12_BN254_mul(&h, &paired);
        }
    }

    mpz_t res_bound;
    mpz_init(res_bound);
    mpz_mul(res_bound, c->bound_x, c->bound_y);
    mpz_mul_ui(res_bound, res_bound, c->m1);

    cfe_error err = cfe_baby_giant_FP12_BN256_with_neg(res, &h, &g, res_bound);

    return err;
}