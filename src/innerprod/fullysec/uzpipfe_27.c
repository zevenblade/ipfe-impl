#include <string.h>

#include <amcl/pair_BN254.h>

#include "cifer/innerprod/fullysec/uzpipfe_27.h"
#include "cifer/internal/big.h"
#include "cifer/sample/uniform.h"
#include "cifer/internal/dlog.h"

cfe_error cfe_uzpipfe_init_27(cfe_uzpipfe_27 *c, size_t m1, size_t m2, mpz_t bound_x, mpz_t bound_y)
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

void cfe_uzpipfe_copy_27(cfe_uzpipfe_27 *res, cfe_uzpipfe_27 *c)
{
    mpz_inits(res->bound_x, res->bound_y, res->order, NULL);

    res->m1 = c->m1;
    res->m2 = c->m2;
    mpz_set(res->bound_x, c->bound_x);
    mpz_set(res->bound_y, c->bound_y);
    mpz_set(res->order, c->order);
}

void cfe_uzpipfe_free_27(cfe_uzpipfe_27 *c)
{
    mpz_clears(c->bound_x, c->bound_y, c->order, NULL);
}

void cfe_uzpipfe_master_pub_key_init_27(cfe_uzpipfe_pub_key_27 *master_pub_key, cfe_uzpipfe_27 *c)
{
    for (size_t i = 0; i < 4; i++)
    {
        cfe_vec_G1_init(&master_pub_key->b[i], 7);
        cfe_vec_G1_init(&master_pub_key->b_tilde[i], 7);
    }
}

void cfe_uzpipfe_master_pub_key_free_27(cfe_uzpipfe_pub_key_27 *master_pub_key)
{
    for (size_t i = 0; i < 4; i++)
    {
        cfe_vec_G1_free(&master_pub_key->b[i]);
        cfe_vec_G1_free(&master_pub_key->b_tilde[i]);
    }
}

void cfe_uzpipfe_master_sec_key_init_27(cfe_uzpipfe_sec_key_27 *master_sec_key, cfe_uzpipfe_27 *c)
{
    for (size_t i = 0; i < 4; i++)
    {
        cfe_vec_G2_init(&master_sec_key->b_star[i], 7);
        cfe_vec_G2_init(&master_sec_key->b_tilde_star[i], 7);
    }
}

void cfe_uzpipfe_master_sec_key_free_27(cfe_uzpipfe_sec_key_27 *master_sec_key)
{
    for (size_t i = 0; i < 4; i++)
    {
        cfe_vec_G2_free(&master_sec_key->b_star[i]);
        cfe_vec_G2_free(&master_sec_key->b_tilde_star[i]);
    }
}

cfe_error cfe_uzpipfe_generate_keys_27(cfe_uzpipfe_sec_key_27 *master_sec_key, cfe_uzpipfe_pub_key_27 *master_pub_key,
                                       cfe_uzpipfe_27 *c)
{
    cfe_error err = CFE_ERR_NONE;
    mpz_t exp, det;
    BIG_256_56 exp_big;
    cfe_vec row;
    cfe_vec_G1 vec_g1;
    cfe_vec_G2 vec_g2;
    cfe_mat B, B_inv, B_star, B_tilde, B_tilde_inv, B_tilde_star;

    mpz_inits(exp, det, NULL);
    cfe_vec_init(&row, 7);
    cfe_vec_G1_init(&vec_g1, 7);
    cfe_vec_G2_init(&vec_g2, 7);
    cfe_mat_inits(7, 7, &B, &B_inv, &B_star, &B_tilde, &B_tilde_inv, &B_tilde_star, NULL);

    cfe_uniform_sample(exp, c->order);
    BIG_256_56_from_mpz(exp_big, exp);
    ECP_BN254_generator(&master_pub_key->g1);
    ECP_BN254_mul(&master_pub_key->g1, exp_big);

    cfe_uniform_sample(exp, c->order);
    BIG_256_56_from_mpz(exp_big, exp);
    ECP2_BN254_generator(&master_pub_key->g2);
    ECP2_BN254_mul(&master_pub_key->g2, exp_big);

    PAIR_BN254_ate(&master_pub_key->gT, &master_pub_key->g2, &master_pub_key->g1);
    PAIR_BN254_fexp(&master_pub_key->gT);

    for (size_t i = 0; i < 7; i++)
    {
        ECP_BN254_copy(&vec_g1.vec[i], &master_pub_key->g1);
    }

    for (size_t i = 0; i < 7; i++)
    {
        ECP2_BN254_copy(&vec_g2.vec[i], &master_pub_key->g2);
    }

    do
    {
        cfe_uniform_sample_mat(&B, c->order);
        err = cfe_mat_inverse_mod_gauss(&B_inv, det, &B, c->order);
    } while (err != CFE_ERR_NONE);
    cfe_mat_transpose(&B_star, &B_inv);

    do
    {
        cfe_uniform_sample_mat(&B_tilde, c->order);
        err = cfe_mat_inverse_mod_gauss(&B_tilde_inv, det, &B_tilde, c->order);
    } while (err != CFE_ERR_NONE);
    cfe_mat_transpose(&B_tilde_star, &B_tilde_inv);

    for (size_t i = 0; i < 4; i++)
    {
        cfe_mat_get_row(&row, &B, i);
        cfe_vec_mul_vec_G1(&master_pub_key->b[i], &row, &vec_g1);

        cfe_mat_get_row(&row, &B_tilde, i);
        cfe_vec_mul_vec_G1(&master_pub_key->b_tilde[i], &row, &vec_g1);

        cfe_mat_get_row(&row, &B_star, i);
        cfe_vec_mul_vec_G2(&master_sec_key->b_star[i], &row, &vec_g2);

        cfe_mat_get_row(&row, &B_tilde_star, i);
        cfe_vec_mul_vec_G2(&master_sec_key->b_tilde_star[i], &row, &vec_g2);
    }

    mpz_clears(exp, det, NULL);
    cfe_vec_G1_free(&vec_g1);
    cfe_mat_frees(&B, &B_inv, &B_star, &B_tilde, &B_tilde_inv, &B_tilde_star, NULL);

    return err;
}

void cfe_uzpipfe_ciphertext_init_27(cfe_uzpipfe_ciphertext_27 *cipher, cfe_uzpipfe_27 *c)
{
    cipher->c1 = malloc(c->m1 * sizeof(cfe_vec_G1));
    for (size_t i = 0; i < c->m1; i++)
    {
        cfe_vec_G1_init(&cipher->c1[i], 7);
    }

    cipher->c2 = malloc(c->m2 * sizeof(cfe_vec_G1));
    for (size_t i = 0; i < c->m2; i++)
    {
        cfe_vec_G1_init(&cipher->c2[i], 7);
    }
}

void cfe_uzpipfe_ciphertext_free_27(cfe_uzpipfe_ciphertext_27 *cipher, cfe_uzpipfe_27 *c)
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

cfe_error cfe_uzpipfe_encrypt_27(cfe_uzpipfe_ciphertext_27 *cipher, cfe_vec *x, cfe_vec *w,
                                 cfe_uzpipfe_sec_key_27 *master_sec_key, cfe_uzpipfe_pub_key_27 *master_pub_key,
                                 cfe_uzpipfe_27 *c)
{
    mpz_t delta, alpha, pi, pi_tilde, el;
    BIG_256_56 el_big;
    ECP_BN254 b_1, b_2, b_3, b_4;

    mpz_inits(delta, alpha, pi, pi_tilde, el, NULL);

    cfe_uniform_sample(delta, c->order);
    // cfe_uniform_sample(alpha, c->order);
    mpz_set_ui(alpha, 1);

    for (size_t i = 0; i < c->m1; i++)
    {
        cfe_uniform_sample(pi, c->order);
        for (size_t j = 0; j < 7; j++)
        {
            mpz_set(el, pi);
            BIG_256_56_from_mpz(el_big, el);
            ECP_BN254_copy(&b_1, &master_pub_key->b[0].vec[j]);
            ECP_BN254_mul(&b_1, el_big);

            mpz_mul_ui(el, pi, i + 1);
            BIG_256_56_from_mpz(el_big, el);
            ECP_BN254_copy(&b_2, &master_pub_key->b[1].vec[j]);
            ECP_BN254_mul(&b_2, el_big);

            cfe_vec_get(el, x, i);
            mpz_mod(el, el, c->order);
            BIG_256_56_from_mpz(el_big, el);
            ECP_BN254_copy(&b_3, &master_pub_key->b[2].vec[j]);
            ECP_BN254_mul(&b_3, el_big);

            mpz_set(el, alpha);
            BIG_256_56_from_mpz(el_big, el);
            ECP_BN254_copy(&b_4, &master_pub_key->b[3].vec[j]);
            ECP_BN254_mul(&b_4, el_big);

            ECP_BN254_copy(&cipher->c1[i].vec[j], &b_1);
            ECP_BN254_add(&cipher->c1[i].vec[j], &b_2);
            ECP_BN254_add(&cipher->c1[i].vec[j], &b_3);
            ECP_BN254_add(&cipher->c1[i].vec[j], &b_4);
        }
    }

    for (size_t i = 0; i < c->m2; i++)
    {
        cfe_uniform_sample(pi_tilde, c->order);
        for (size_t j = 0; j < 7; j++)
        {
            mpz_set(el, pi_tilde);
            BIG_256_56_from_mpz(el_big, el);
            ECP_BN254_copy(&b_1, &master_pub_key->b_tilde[0].vec[j]);
            ECP_BN254_mul(&b_1, el_big);

            mpz_mul_ui(el, pi_tilde, i + 1);
            BIG_256_56_from_mpz(el_big, el);
            ECP_BN254_copy(&b_2, &master_pub_key->b_tilde[1].vec[j]);
            ECP_BN254_mul(&b_2, el_big);

            cfe_vec_get(el, w, i);
            mpz_mul(el, el, delta);
            mpz_mod(el, el, c->order);
            BIG_256_56_from_mpz(el_big, el);
            ECP_BN254_copy(&b_3, &master_pub_key->b_tilde[2].vec[j]);
            ECP_BN254_mul(&b_3, el_big);

            mpz_set(el, alpha);
            BIG_256_56_from_mpz(el_big, el);
            ECP_BN254_copy(&b_4, &master_pub_key->b_tilde[3].vec[j]);
            ECP_BN254_mul(&b_4, el_big);

            ECP_BN254_copy(&cipher->c2[i].vec[j], &b_1);
            ECP_BN254_add(&cipher->c2[i].vec[j], &b_2);
            ECP_BN254_add(&cipher->c2[i].vec[j], &b_3);
            ECP_BN254_add(&cipher->c2[i].vec[j], &b_4);
        }
    }

    return CFE_ERR_NONE;
}

void cfe_uzpipfe_fe_key_init_27(cfe_uzpipfe_fe_key_27 *fe_key, cfe_uzpipfe_27 *c)
{
    fe_key->k1 = malloc(c->m1 * sizeof(cfe_vec_G2));
    for (size_t i = 0; i < c->m1; i++)
    {
        cfe_vec_G2_init(&fe_key->k1[i], 7);
    }

    fe_key->k2 = malloc(c->m2 * sizeof(cfe_vec_G2));
    for (size_t i = 0; i < c->m2; i++)
    {
        cfe_vec_G2_init(&fe_key->k2[i], 7);
    }
}

void cfe_uzpipfe_fe_key_free_27(cfe_uzpipfe_fe_key_27 *fe_key, cfe_uzpipfe_27 *c)
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

cfe_error cfe_uzpipfe_derive_fe_key_27(cfe_uzpipfe_fe_key_27 *fe_key, cfe_vec *y, cfe_vec *v,
                                       cfe_uzpipfe_sec_key_27 *master_sec_key, cfe_uzpipfe_pub_key_27 *master_pub_key,
                                       cfe_uzpipfe_27 *c)
{
    mpz_t zero, omega, rho, rho_tilde, el, sum, bound, bound_neg, dot;
    BIG_256_56 el_big;
    ECP2_BN254 b_1, b_2, b_3, b_4;
    cfe_vec gamma_v, gamma_tilde_v, k_tmp, B_y;

    mpz_inits(zero, omega, rho, rho_tilde, el, sum, bound, bound_neg, dot, NULL);
    cfe_vec_inits(7, &k_tmp, &B_y, NULL);
    cfe_vec_init(&gamma_v, c->m1);
    cfe_vec_init(&gamma_tilde_v, c->m2);

    cfe_uniform_sample(omega, c->order);

    mpz_set(bound, c->order);
    mpz_neg(bound_neg, bound);
    mpz_set_ui(zero, 0);

    cfe_uniform_sample_range_vec(&gamma_v, bound_neg, bound);
    cfe_uniform_sample_range_vec(&gamma_tilde_v, bound_neg, bound);
    cfe_vec_set(&gamma_tilde_v, zero, ((c->m2) - 1));

    mpz_set(sum, zero);
    for (size_t i = 0; i < c->m1; i++)
    {
        cfe_vec_get(el, &gamma_v, i);
        mpz_add(sum, sum, el);
    }
    for (size_t i = 0; i < c->m2; i++)
    {
        cfe_vec_get(el, &gamma_tilde_v, i);
        mpz_add(sum, sum, el);
    }

    mpz_neg(el, sum);
    cfe_vec_set(&gamma_tilde_v, el, ((c->m2) - 1));

    for (size_t i = 0; i < c->m1; i++)
    {
        cfe_uniform_sample(rho, c->order);
        for (size_t j = 0; j < 7; j++)
        {
            mpz_mul_ui(el, rho, i + 1);
            mpz_neg(el, el);
            mpz_mod(el, el, c->order);
            BIG_256_56_from_mpz(el_big, el);
            ECP2_BN254_copy(&b_1, &master_sec_key->b_star[0].vec[j]);
            ECP2_BN254_mul(&b_1, el_big);

            mpz_set(el, rho);
            BIG_256_56_from_mpz(el_big, el);
            ECP2_BN254_copy(&b_2, &master_sec_key->b_star[1].vec[j]);
            ECP2_BN254_mul(&b_2, el_big);

            cfe_vec_get(el, y, i);
            mpz_mod(el, el, c->order);
            BIG_256_56_from_mpz(el_big, el);
            ECP2_BN254_copy(&b_3, &master_sec_key->b_star[2].vec[j]);
            ECP2_BN254_mul(&b_3, el_big);

            cfe_vec_get(el, &gamma_v, i);
            mpz_mod(el, el, c->order);
            BIG_256_56_from_mpz(el_big, el);
            ECP2_BN254_copy(&b_4, &master_sec_key->b_star[3].vec[j]);
            ECP2_BN254_mul(&b_4, el_big);

            ECP2_BN254_copy(&fe_key->k1[i].vec[j], &b_1);
            ECP2_BN254_add(&fe_key->k1[i].vec[j], &b_2);
            ECP2_BN254_add(&fe_key->k1[i].vec[j], &b_3);
            ECP2_BN254_add(&fe_key->k1[i].vec[j], &b_4);
        }
    }

    for (size_t i = 0; i < c->m2; i++)
    {
        cfe_uniform_sample(rho_tilde, c->order);
        for (size_t j = 0; j < 7; j++)
        {
            mpz_mul_ui(el, rho_tilde, i + 1);
            mpz_neg(el, el);
            mpz_mod(el, el, c->order);
            BIG_256_56_from_mpz(el_big, el);
            ECP2_BN254_copy(&b_1, &master_sec_key->b_tilde_star[0].vec[j]);
            ECP2_BN254_mul(&b_1, el_big);

            mpz_set(el, rho_tilde);
            BIG_256_56_from_mpz(el_big, el);
            ECP2_BN254_copy(&b_2, &master_sec_key->b_tilde_star[1].vec[j]);
            ECP2_BN254_mul(&b_2, el_big);

            cfe_vec_get(el, v, i);
            mpz_mul(el, el, omega);
            mpz_mod(el, el, c->order);
            BIG_256_56_from_mpz(el_big, el);
            ECP2_BN254_copy(&b_3, &master_sec_key->b_tilde_star[2].vec[j]);
            ECP2_BN254_mul(&b_3, el_big);

            cfe_vec_get(el, &gamma_tilde_v, i);
            mpz_mod(el, el, c->order);
            BIG_256_56_from_mpz(el_big, el);
            ECP2_BN254_copy(&b_4, &master_sec_key->b_tilde_star[3].vec[j]);
            ECP2_BN254_mul(&b_4, el_big);

            ECP2_BN254_copy(&fe_key->k2[i].vec[j], &b_1);
            ECP2_BN254_add(&fe_key->k2[i].vec[j], &b_2);
            ECP2_BN254_add(&fe_key->k2[i].vec[j], &b_3);
            ECP2_BN254_add(&fe_key->k2[i].vec[j], &b_4);
        }
    }

    return CFE_ERR_NONE;
}

cfe_error cfe_uzpipfe_decrypt_27(mpz_t res, cfe_uzpipfe_ciphertext_27 *cipher,
                                 cfe_uzpipfe_fe_key_27 *fe_key, cfe_uzpipfe_pub_key_27 *master_pub_key,
                                 cfe_uzpipfe_27 *c)
{
    ECP_BN254 g1;
    ECP2_BN254 g2;
    FP12_BN254 h, g;
    FP12_BN254 paired;

    ECP_BN254_copy(&g1, &master_pub_key->g1);
    ECP2_BN254_copy(&g2, &master_pub_key->g2);

    FP12_BN254_copy(&g, &master_pub_key->gT);
    FP12_BN254_one(&h);

    for (size_t i = 0; i < c->m1; i++)
    {
        for (size_t j = 0; j < 7; j++)
        {
            PAIR_BN254_ate(&paired, &fe_key->k1[i].vec[j], &cipher->c1[i].vec[j]);
            PAIR_BN254_fexp(&paired);
            FP12_BN254_mul(&h, &paired);
        }
    }

    for (size_t i = 0; i < c->m2; i++)
    {
        for (size_t j = 0; j < 7; j++)
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
    //mpz_set_ui(res_bound, 10000);

    cfe_error err = cfe_baby_giant_FP12_BN256_with_neg(res, &h, &g, res_bound);

    return err;
}