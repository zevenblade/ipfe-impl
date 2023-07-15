#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <amcl/big_256_56.h>
#include <amcl/pair_BN254.h>

#include "cifer/abe/uabipfe.h"
#include "cifer/internal/big.h"
#include "cifer/sample/uniform.h"
#include "cifer/internal/dlog.h"

#include "cifer/abe/policy.h"

#include "cifer/internal/common.h"
#include "cifer/internal/big.h"
#include "cifer/internal/hash.h"
#include "cifer/sample/uniform.h"
#include "cifer/abe/fame.h"

cfe_error cfe_uabipfe_init(cfe_uabipfe *c, size_t m1, mpz_t bound_x, mpz_t bound_y)
{
    cfe_error err = CFE_ERR_NONE;
    mpz_t check, order;
    mpz_inits(check, order, NULL);
    mpz_from_BIG_256_56(order, (int64_t *)CURVE_Order_BN254);

    mpz_mul(check, bound_x, bound_y);
    mpz_mul_ui(check, check, m1);

    if (mpz_cmp(check, order) >= 0)
    {
        err = CFE_ERR_PRECONDITION_FAILED;
        goto cleanup;
    }

    mpz_inits(c->bound_x, c->bound_y, c->order, NULL);

    c->m1 = m1;
    mpz_set(c->bound_x, bound_x);
    mpz_set(c->bound_y, bound_y);
    mpz_set(c->order, order);

cleanup:
    mpz_clears(check, order, NULL);
    return err;
}

void cfe_uabipfe_copy(cfe_uabipfe *res, cfe_uabipfe *c)
{
    mpz_inits(res->bound_x, res->bound_y, res->order, NULL);

    res->m1 = c->m1;
    mpz_set(res->bound_x, c->bound_x);
    mpz_set(res->bound_y, c->bound_y);
    mpz_set(res->order, c->order);
}

void cfe_uabipfe_free(cfe_uabipfe *c)
{
    mpz_clears(c->bound_x, c->bound_y, c->order, NULL);
}

void cfe_uabipfe_master_pub_key_init(cfe_uabipfe_master_pub_key *master_pub_key)
{
    for (size_t i = 0; i < 3; i++)
    {
        cfe_vec_G1_init(&master_pub_key->f[i], 8);
    }
    cfe_vec_G1_init(&master_pub_key->h12, 4);
    cfe_vec_G1_init(&master_pub_key->h3, 4);
}

void cfe_uabipfe_master_pub_key_free(cfe_uabipfe_master_pub_key *master_pub_key)
{
    cfe_vec_G1_free(&master_pub_key->f[0]);
    cfe_vec_G1_free(&master_pub_key->f[1]);
    cfe_vec_G1_free(&master_pub_key->f[2]);
    cfe_vec_G1_free(&master_pub_key->h12);
    cfe_vec_G1_free(&master_pub_key->h3);
}

void cfe_uabipfe_master_sec_key_init(cfe_uabipfe_master_sec_key *master_sec_key)
{
    for (size_t i = 0; i < 3; i++)
    {
        cfe_vec_init(&master_sec_key->f_star[i], 8);
    }
    for (size_t i = 0; i < 3; i++)
    {
        cfe_vec_init(&master_sec_key->h_star[i], 4);
    }

    mpz_init(master_sec_key->z);
}

void cfe_uabipfe_master_sec_key_free(cfe_uabipfe_master_sec_key *master_sec_key)
{
    for (size_t i = 0; i < 3; i++)
    {
        cfe_vec_free(&master_sec_key->f_star[i]);
    }
    for (size_t i = 0; i < 3; i++)
    {
        cfe_vec_free(&master_sec_key->h_star[i]);
    }
}

cfe_error cfe_uabipfe_generate_keys(cfe_uabipfe_master_pub_key *master_pub_key, cfe_uabipfe_master_sec_key *master_sec_key,
                                    cfe_uabipfe *c)
{
    cfe_error err = CFE_ERR_NONE;

    mpz_t det, mu;
    BIG_256_56 mu_big;
    cfe_vec col_4_1, col_4_2, col_8, sum;
    cfe_vec_G1 vec_g1_4, vec_g1_8;
    cfe_mat F, F_inv, F_star, H, H_inv, H_star;

    mpz_inits(det, mu, NULL);

    cfe_vec_inits(4, &sum, &col_4_1, &col_4_2, NULL);
    cfe_vec_inits(8, &col_8, NULL);

    cfe_vec_G1_init(&vec_g1_4, 4);
    cfe_vec_G1_init(&vec_g1_8, 8);

    cfe_mat_inits(8, 8, &F, &F_inv, &F_star, NULL);
    cfe_mat_inits(4, 4, &H, &H_inv, &H_star, NULL);

    cfe_uniform_sample(mu, c->order);
    BIG_256_56_from_mpz(mu_big, mu);

    cfe_uniform_sample(master_sec_key->z, c->order);

    ECP_BN254_generator(&master_pub_key->g1);
    ECP2_BN254_generator(&master_pub_key->g2);

    PAIR_BN254_ate(&master_pub_key->gT, &master_pub_key->g2, &master_pub_key->g1);
    PAIR_BN254_fexp(&master_pub_key->gT);

    for (size_t i = 0; i < 4; i++)
    {
        ECP_BN254_copy(&vec_g1_4.vec[i], &master_pub_key->g1);
    }

    for (size_t i = 0; i < 8; i++)
    {
        ECP_BN254_copy(&vec_g1_8.vec[i], &master_pub_key->g1);
    }

    ECP_BN254_generator(&master_sec_key->g1);
    ECP2_BN254_generator(&master_sec_key->g2);

    PAIR_BN254_ate(&master_sec_key->gT, &master_sec_key->g2, &master_sec_key->g1);
    PAIR_BN254_fexp(&master_sec_key->gT);

    do
    {
        cfe_uniform_sample_mat(&F, c->order);
        err = cfe_mat_inverse_mod_gauss(&F_inv, det, &F, c->order);
    } while (err != CFE_ERR_NONE);
    cfe_mat_transpose(&F_star, &F_inv);

    do
    {
        cfe_uniform_sample_mat(&H, c->order);
        err = cfe_mat_inverse_mod_gauss(&H_inv, det, &H, c->order);
    } while (err != CFE_ERR_NONE);
    cfe_mat_transpose(&H_star, &H_inv);

    for (size_t i = 0; i < 3; i++)
    {
        cfe_mat_get_row(&col_8, &F, i);
        cfe_vec_mul_vec_G1(&master_pub_key->f[i], &col_8, &vec_g1_8);
    }

    cfe_mat_get_row(&col_4_1, &H, 0);
    cfe_mat_get_row(&col_4_2, &H, 1);
    cfe_vec_mul_scalar(&col_4_2, &col_4_2, mu);
    cfe_vec_mod(&col_4_2, &col_4_2, c->order);
    cfe_vec_add(&sum, &col_4_1, &col_4_2);
    cfe_vec_mul_vec_G1(&master_pub_key->h12, &sum, &vec_g1_4);

    cfe_mat_get_row(&col_4_1, &H, 2);
    cfe_vec_mul_vec_G1(&master_pub_key->h3, &col_4_1, &vec_g1_4);

    ECP_BN254_copy(&master_pub_key->mu, &master_pub_key->g1);
    ECP_BN254_mul(&master_pub_key->mu, mu_big);

    for (size_t i = 0; i < 3; i++)
    {
        cfe_mat_get_row(&col_8, &F_star, i);
        cfe_vec_copy(&master_sec_key->f_star[i], &col_8);
    }

    for (size_t i = 0; i < 3; i++)
    {
        cfe_mat_get_row(&col_4_1, &H_star, i);
        cfe_vec_copy(&master_sec_key->h_star[i], &col_4_1);
    }

    mpz_clears(det, mu, NULL);
    cfe_vec_frees(&col_4_1, &col_4_2, &col_8, &sum, NULL);
    cfe_vec_G1_free(&vec_g1_4);
    cfe_vec_G1_free(&vec_g1_8);
    cfe_mat_frees(&F, &F_inv, &F_star, &H, &H_inv, &H_star, NULL);

    return err;
}

void cfe_uabipfe_ciphertext_init(cfe_uabipfe_ciphertext *cipher, cfe_msp *msp, cfe_uabipfe *c)
{
    cipher->c_j = cfe_malloc(msp->mat.rows * sizeof(cfe_vec_G1));
    for (size_t i = 0; i < msp->mat.rows; i++)
    {
        cfe_vec_G1_init(&cipher->c_j[i], 8);
    }

    cfe_vec_G1_init(&cipher->c_ipfe, 4);
    cipher->t = malloc(c->m1 * sizeof(FP12_BN254));
}

void cfe_uabipfe_ciphertext_free(cfe_uabipfe_ciphertext *cipher, cfe_msp *msp, cfe_uabipfe *c)
{
    for (size_t i = 0; i < msp->mat.rows; i++)
    {
        cfe_vec_G1_free(&cipher->c_j[i]);
    }
    cfe_vec_G1_free(&cipher->c_ipfe);
}

cfe_error cfe_uabipfe_encrypt(cfe_uabipfe_ciphertext *cipher, cfe_uabipfe_master_pub_key *master_pub_key,
                              cfe_vec *x, cfe_msp *msp, cfe_uabipfe *c)
{
    cfe_error err = CFE_ERR_NONE;

    mpz_t omega, psi, el;
    BIG_256_56 x_big, omega_big, psi_big, sigma_big, sigma_j_big;
    cfe_vec x_tmp, sigma;
    ECP_BN254 el_g1, h_t_1, h_t_2, f_t_1, f_t_2, f_t_3;
    ECP2_BN254 el_g2, u_i, s_i;
    FP12_BN254 paired;
    cfe_string str_i;

    mpz_inits(omega, psi, el, NULL);

    cfe_vec_init(&x_tmp, c->m1);
    cfe_vec_init(&sigma, msp->mat.rows);

    cfe_vec_copy(&x_tmp, x);
    cfe_vec_mod(&x_tmp, &x_tmp, c->order);

    cfe_uniform_sample_vec(&sigma, c->order);

    cfe_uniform_sample(omega, c->order);
    BIG_256_56_from_mpz(omega_big, omega);
    cfe_uniform_sample(psi, c->order);
    BIG_256_56_from_mpz(psi_big, psi);

    for (size_t j = 0; j < msp->mat.rows; j++)
    {
        for (size_t i = 0; i < 8; i++)
        {
            mpz_set(el, sigma.vec[j]);
            BIG_256_56_from_mpz(sigma_big, el);
            ECP_BN254_copy(&f_t_1, &master_pub_key->f[0].vec[i]);
            ECP_BN254_mul(&f_t_1, sigma_big);

            mpz_set(el, sigma.vec[j]);
            mpz_mul_si(el, el, j + 1);
            mpz_mul_si(el, el, -1);
            mpz_mod(el, el, c->order);
            BIG_256_56_from_mpz(sigma_j_big, el);
            ECP_BN254_copy(&f_t_2, &master_pub_key->f[1].vec[i]);
            ECP_BN254_mul(&f_t_2, sigma_j_big);

            ECP_BN254_copy(&f_t_3, &master_pub_key->f[2].vec[i]);
            ECP_BN254_mul(&f_t_3, psi_big);

            ECP_BN254_copy(&cipher->c_j[j].vec[i], &f_t_1);
            ECP_BN254_add(&cipher->c_j[j].vec[i], &f_t_2);
            ECP_BN254_add(&cipher->c_j[j].vec[i], &f_t_3);
        }
    }

    for (size_t i = 0; i < 4; i++)
    {
        ECP_BN254_copy(&h_t_1, &master_pub_key->h12.vec[i]);
        ECP_BN254_mul(&h_t_1, omega_big);

        ECP_BN254_copy(&h_t_2, &master_pub_key->h3.vec[i]);
        ECP_BN254_mul(&h_t_2, psi_big);

        ECP_BN254_copy(&cipher->c_ipfe.vec[i], &h_t_1);
        ECP_BN254_add(&cipher->c_ipfe.vec[i], &h_t_2);
    }

    for (size_t i = 0; i < c->m1; i++)
    {
        ECP_BN254_copy(&el_g1, &master_pub_key->g1);
        ECP2_BN254_copy(&el_g2, &master_pub_key->g2);

        BIG_256_56_from_mpz(x_big, x_tmp.vec[i]);
        ECP2_BN254_mul(&el_g2, x_big);

        PAIR_BN254_ate(&cipher->t[i], &el_g2, &el_g1);
        PAIR_BN254_fexp(&cipher->t[i]);

        cfe_int_to_str(&str_i, (int)i);
        cfe_hash_G2(&u_i, &str_i);
        cfe_hash_G2(&s_i, &str_i);

        ECP2_BN254_mul(&u_i, omega_big);
        PAIR_BN254_ate(&paired, &u_i, &master_pub_key->g1);
        PAIR_BN254_fexp(&paired);
        FP12_BN254_mul(&cipher->t[i], &paired);

        ECP2_BN254_mul(&s_i, omega_big);
        PAIR_BN254_ate(&paired, &s_i, &master_pub_key->mu);
        PAIR_BN254_fexp(&paired);
        FP12_BN254_mul(&cipher->t[i], &paired);
    }

    return err;
}

void cfe_uabipfe_fe_key_init(cfe_uabipfe_fe_key *fe_key, cfe_msp *msp, cfe_uabipfe *c)
{
    fe_key->k_j = cfe_malloc(msp->mat.rows * sizeof(cfe_vec_G2));
    for (size_t i = 0; i < msp->mat.rows; i++)
    {
        cfe_vec_G2_init(&fe_key->k_j[i], 8);
    }

    cfe_vec_G2_init(&fe_key->k_ipfe, 4);
}

void cfe_uabipfe_fe_key_free(cfe_uabipfe_fe_key *fe_key, cfe_msp *msp, cfe_uabipfe *c)
{
    for (size_t i = 0; i < msp->mat.rows; i++)
    {
        cfe_vec_G2_free(&fe_key->k_j[i]);
    }
    cfe_vec_G2_free(&fe_key->k_ipfe);
}

cfe_error cfe_uabipfe_derive_fe_key(cfe_uabipfe_fe_key *fe_key, cfe_vec *y, cfe_msp *msp,
                                    cfe_uabipfe_master_sec_key *sec_key, cfe_uabipfe *c)
{
    cfe_error err = CFE_ERR_NONE;

    mpz_t a_0, a_0_z, el, zero;
    ECP2_BN254 u_i, s_i, sum_u_i, sum_s_i, ss_u, ss_i, g_a_z, ss_a_z;
    cfe_vec y_tmp, k_tmp, pi_j, s_x, shares;
    cfe_vec_G2 vec_g2;
    cfe_string str_i;
    BIG_256_56 y_big, h1_big, h2_big, h3_big, a_0_z_big;

    mpz_inits(a_0, a_0_z, el, zero, NULL);

    cfe_vec_init(&y_tmp, c->m1);
    cfe_vec_init(&k_tmp, 8);
    cfe_vec_init(&pi_j, msp->mat.rows);
    cfe_vec_init(&s_x, msp->mat.cols);
    cfe_vec_init(&shares, msp->mat.rows);

    cfe_vec_G2_init(&vec_g2, 8);

    cfe_vec_copy(&y_tmp, y);
    cfe_vec_neg(&y_tmp, &y_tmp);
    cfe_vec_mod(&y_tmp, &y_tmp, c->order);

    mpz_set_ui(zero, 0);

    for (size_t i = 0; i < 8; i++)
    {
        ECP2_BN254_copy(&vec_g2.vec[i], &sec_key->g2);
    }

    cfe_uniform_sample_vec(&pi_j, c->order);

    cfe_uniform_sample(a_0, c->order);
    mpz_neg(a_0_z, a_0);
    mpz_mul(a_0_z, a_0_z, sec_key->z);
    mpz_mod(a_0_z, a_0_z, c->order);
    BIG_256_56_from_mpz(a_0_z_big, a_0_z);
    ECP2_BN254_copy(&g_a_z, &sec_key->g2);
    ECP2_BN254_mul(&g_a_z, a_0_z_big);

    cfe_uniform_sample_vec(&s_x, c->order);
    mpz_set(s_x.vec[0], a_0);
    cfe_mat_mul_vec(&shares, &msp->mat, &s_x);

    for (size_t j = 0; j < msp->mat.rows; j++)
    {
        cfe_vec_set_const(&k_tmp, zero);
        for (size_t i = 0; i < 8; i++)
        {
            mpz_set(el, pi_j.vec[j]);
            mpz_mul_ui(el, el, j + 1);
            mpz_mul(el, el, sec_key->f_star[0].vec[i]);
            mpz_mod(el, el, c->order);
            mpz_add(k_tmp.vec[i], k_tmp.vec[i], el);

            mpz_set(el, pi_j.vec[j]);
            mpz_mul(el, el, sec_key->f_star[1].vec[i]);
            mpz_mod(el, el, c->order);
            mpz_add(k_tmp.vec[i], k_tmp.vec[i], el);

            mpz_set(el, shares.vec[j]);
            mpz_mul(el, el, sec_key->z);
            mpz_mul(el, el, sec_key->f_star[2].vec[i]);
            mpz_mod(el, el, c->order);
            mpz_add(k_tmp.vec[i], k_tmp.vec[i], el);
        }
        cfe_vec_mul_vec_G2(&fe_key->k_j[j], &k_tmp, &vec_g2);
    }

    cfe_int_to_str(&str_i, (int)0);
    cfe_hash_G2(&u_i, &str_i);
    cfe_hash_G2(&s_i, &str_i);

    BIG_256_56_from_mpz(y_big, y_tmp.vec[0]);

    ECP2_BN254_copy(&sum_u_i, &u_i);
    ECP2_BN254_mul(&sum_u_i, y_big);

    ECP2_BN254_copy(&sum_s_i, &s_i);
    ECP2_BN254_mul(&sum_s_i, y_big);

    for (size_t i = 1; i < c->m1; i++)
    {
        cfe_int_to_str(&str_i, (int)i);
        cfe_hash_G2(&u_i, &str_i);
        cfe_hash_G2(&s_i, &str_i);

        BIG_256_56_from_mpz(y_big, y_tmp.vec[i]);

        ECP2_BN254_mul(&u_i, y_big);
        ECP2_BN254_add(&sum_u_i, &u_i);

        ECP2_BN254_mul(&s_i, y_big);
        ECP2_BN254_add(&sum_s_i, &s_i);
    }

    for (size_t i = 0; i < 4; i++)
    {
        ECP2_BN254_copy(&ss_u, &sum_u_i);
        ECP2_BN254_copy(&ss_i, &sum_s_i);
        ECP2_BN254_copy(&ss_a_z, &g_a_z);

        BIG_256_56_from_mpz(h1_big, sec_key->h_star[0].vec[i]);
        ECP2_BN254_mul(&ss_u, h1_big);

        BIG_256_56_from_mpz(h2_big, sec_key->h_star[1].vec[i]);
        ECP2_BN254_mul(&ss_i, h2_big);

        BIG_256_56_from_mpz(h3_big, sec_key->h_star[2].vec[i]);
        ECP2_BN254_mul(&ss_a_z, h3_big);

        ECP2_BN254_copy(&fe_key->k_ipfe.vec[i], &ss_u);
        ECP2_BN254_add(&fe_key->k_ipfe.vec[i], &ss_i);
        ECP2_BN254_add(&fe_key->k_ipfe.vec[i], &ss_a_z);
    }

    return err;
}

cfe_error cfe_uabipfe_decrypt(mpz_t res, cfe_uabipfe_ciphertext *cipher,
                              cfe_uabipfe_master_pub_key *master_pub_key,
                              cfe_uabipfe_fe_key *fe_key, cfe_vec *y,
                              cfe_msp *msp, int *owned_attrib, size_t num_attrib,
                              cfe_uabipfe *c)
{
    cfe_error err = CFE_ERR_NONE;

    int attrib_vec[msp->mat.rows];
    for (size_t i = 0; i < msp->mat.rows; i++)
    {
        attrib_vec[i] = 0;
    }

    mpz_t zero;
    ECP2_BN254 k_j_c_j;
    FP12_BN254 prod_t, t_tmp, g, prod_ipfe, paired_ipfe, prod_j, paired_j;
    BIG_256_56 y_big, alpha_big;
    cfe_vec y_tmp, row_vec, one_vec, alpha;
    cfe_mat mat_red, mat_red_t;

    mpz_inits(zero, NULL);

    FP12_BN254_copy(&g, &master_pub_key->gT);

    cfe_vec_init(&y_tmp, c->m1);
    cfe_vec_init(&one_vec, msp->mat.cols);
    cfe_vec_init(&row_vec, msp->mat.cols);

    cfe_mat_init(&mat_red, num_attrib, msp->mat.cols);
    cfe_mat_init(&mat_red_t, msp->mat.cols, num_attrib);

    cfe_vec_copy(&y_tmp, y);
    cfe_vec_mod(&y_tmp, &y_tmp, c->order);

    mpz_set_ui(zero, 0);

    for (size_t i = 0; i < num_attrib; i++)
    {
        for (size_t j = 0; j < msp->mat.rows; j++)
        {
            if (msp->row_to_attrib[j] == owned_attrib[i])
            {
                attrib_vec[j] = 1;
            }
        }
    }

    for (size_t i = 0, j = 0; i < msp->mat.rows; i++)
    {
        if (attrib_vec[i] == 1)
        {
            cfe_mat_get_row(&row_vec, &msp->mat, i);
            cfe_mat_set_vec(&mat_red, &row_vec, j);
            j++;
        }
    }
    cfe_mat_transpose(&mat_red_t, &mat_red);

    cfe_vec_set_const(&one_vec, zero);
    mpz_set_ui(one_vec.vec[0], 1);

    err = cfe_gaussian_elimination_solver(&alpha, &mat_red_t, &one_vec, c->order);
    if (err == CFE_ERR_NO_SOLUTION_EXISTS)
    {
        return CFE_ERR_INSUFFICIENT_KEYS;
    }

    FP12_BN254_one(&prod_t);
    for (size_t i = 0; i < c->m1; i++)
    {
        if (mpz_cmp_ui(y->vec[i], 0) != 0)
        {
            FP12_BN254_copy(&t_tmp, &cipher->t[i]);

            BIG_256_56_from_mpz(y_big, y_tmp.vec[i]);
            FP12_BN254_pow(&t_tmp, &t_tmp, y_big);
            FP12_BN254_mul(&prod_t, &t_tmp);
        }
    }

    FP12_BN254_one(&prod_ipfe);
    for (size_t i = 0; i < 4; i++)
    {
        PAIR_BN254_ate(&paired_ipfe, &fe_key->k_ipfe.vec[i], &cipher->c_ipfe.vec[i]);
        PAIR_BN254_fexp(&paired_ipfe);
        FP12_BN254_mul(&prod_ipfe, &paired_ipfe);
    }

    FP12_BN254_one(&prod_j);
    for (size_t j = 0, k = 0; j < msp->mat.rows; j++)
    {
        if (attrib_vec[j] == 1)
        {
            for (size_t i = 0; i < 8; i++)
            {
                ECP2_BN254_copy(&k_j_c_j, &fe_key->k_j[j].vec[i]);
                BIG_256_56_from_mpz(alpha_big, alpha.vec[k]);
                ECP2_BN254_mul(&k_j_c_j, alpha_big);
                PAIR_BN254_ate(&paired_j, &k_j_c_j, &cipher->c_j[j].vec[i]);
                PAIR_BN254_fexp(&paired_j);
                FP12_BN254_mul(&prod_j, &paired_j);
            }
            k++;
        }
    }

    FP12_BN254_mul(&prod_t, &prod_ipfe);
    FP12_BN254_mul(&prod_t, &prod_j);

    mpz_t res_bound;
    mpz_init(res_bound);
    mpz_mul(res_bound, c->bound_x, c->bound_y);
    mpz_mul_ui(res_bound, res_bound, c->m1);

    err = cfe_baby_giant_FP12_BN256_with_neg(res, &prod_t, &g, res_bound);

    return err;
}