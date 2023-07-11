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

#ifndef CIFER_UZPIPFE_H
#define CIFER_UZPIPFE_H

#include "cifer/data/mat.h"
#include "cifer/data/vec_curve.h"

/**
 * cfe_uzpipfe contains the shared choice for parameters on which
 * the functionality of the scheme depend.
 */
typedef struct cfe_uzpipfe
{
    size_t m1;
    size_t m2;
    mpz_t bound_x;
    mpz_t bound_y;
    mpz_t order;
} cfe_uzpipfe;
/**
 * Configures a new client for the fhipe scheme. It returns an error if
 * the bounds and length of vectors are too high.
 *
 * @param c A pointer to an uninitialized struct representing the scheme
 * @param m1 Length of the vector x that will be encrypted
 * @param m2 Length of the attribute vector w
 * @param bound_x Bound on the inputs of the vectors that will be encrypted
 * @param bound_y Bound on the inputs of the inner product vectors for which
 * the functional keys will be generated.
 * @return Error code
 */
cfe_error cfe_uzpipfe_init(cfe_uzpipfe *c, size_t m1, size_t m2, mpz_t bound_x, mpz_t bound_y);

/**
 * Reconstructs the scheme with the same configuration parameters from
 * an already existing fhipe scheme instance.
 *
 * @param res A pointer to an uninitialized cfe_uzpipfe struct
 * @param c A pointer to an instance of the scheme (*initialized* cfe_uzpipfe
 * struct)
 */
void cfe_uzpipfe_copy(cfe_uzpipfe *res, cfe_uzpipfe *c);

/**
 * Frees the memory occupied by the struct members. It does not free the
 * memory occupied by the struct itself.
 *
 * @param c A pointer to an instance of the scheme (*initialized* cfe_fhipe
 * struct)
 */
void cfe_uzpipfe_free(cfe_uzpipfe *c);

/**
 * cfe_fhipe_sec_key represents a master secret key in fhipe scheme.
 */
typedef struct cfe_uzpipfe_pub_key
{
    ECP_BN254 g1;
    ECP2_BN254 g2;
    FP12_BN254 gT;
} cfe_uzpipfe_pub_key;

typedef struct cfe_uzpipfe_sec_key
{
    unsigned char key[32];
    unsigned char key_tilde[32];
} cfe_uzpipfe_sec_key;

/**
 * Initializes the struct which represents the master secret key in fhipe.
 *
 * @param sec_key A pointer to an uninitialized cfe_uzpipfe_sec_key struct
 * @param c A pointer to an instance of the scheme (*initialized* cfe_uzpipfe
 * struct)
 */
void cfe_uzpipfe_master_key_init(cfe_uzpipfe_sec_key *sec_key, cfe_uzpipfe *c);

/**
 * Frees the memory occupied by the struct members. It does
 * not free the memory occupied by the struct itself.
 *
 * @param sec_key A pointer to an *initialized* cfe_fhipe_sec_key struct
 */
void cfe_uzpipfe_master_key_free(cfe_uzpipfe_sec_key *sec_key);

/**
 * Generates a master secret key and a public key for the scheme. It returns an
 * error if generating one of the parts of the secret key failed.
 *
 * @param sec_key A pointer to a cfe_uzpipfe_sec_key struct (the
 * master secret key will be stored here)
 * @param pub_key A pointer to a cfe_uzpipfe_pub_key struct (the public key will
 * be stored here)
 * @param c A pointer to an instance of the scheme (*initialized*
 * cfe_fh_multi_ipe struct)
 * @return Error code
 */
cfe_error cfe_uzpipfe_generate_keys(cfe_uzpipfe_sec_key *sec_key, cfe_uzpipfe_pub_key *pub_key,
                                    cfe_uzpipfe *c);

typedef struct cfe_uzpipfe_fe_key
{
    cfe_vec_G2 *k1;
    cfe_vec_G2 *k2;
} cfe_uzpipfe_fe_key;

/**
 * Initializes the struct which represents the functional encryption key in fhipe.
 *
 * @param fe_key A pointer to an uninitialized cfe_uzpipfe_FE_key struct
 * @param c A pointer to an instance of the scheme (*initialized* cfe_uzpipfe
 * struct)
 */
void cfe_uzpipfe_fe_key_init(cfe_uzpipfe_fe_key *fe_key, cfe_uzpipfe *c);

/**
 * Frees the memory occupied by the struct members. It does
 * not free the memory occupied by the struct itself.
 *
 * @param fe_key A pointer to an *initialized* cfe_uzpipfe_FE_key struct
 */
void cfe_uzpipfe_fe_key_free(cfe_uzpipfe_fe_key *fe_key,  cfe_uzpipfe *c);

/**
 * Takes a master secret key and input vector y, and derives the functional
 * encryption key. In case the key could not be derived, it returns an error.
 *
 * @param fe_key A pointer to a cfe_fhipe_fe_key struct (the functional
 * encryption key will be stored here)
 * @param y A pointer to the inner product vector
 * @param v A pointer to the predicate product vector
 * @param sec_key A pointer to the master secret key
 * @param pub_key A pointer to the public key
 * @param c A pointer to an instance of the scheme (*initialized* cfe_fhipe
 * struct)
 * @return Error code
 */
cfe_error cfe_uzpipfe_derive_fe_key(cfe_uzpipfe_fe_key *fe_key, cfe_vec *y, cfe_vec *v,
                                    cfe_uzpipfe_sec_key *sec_key, cfe_uzpipfe_pub_key *pub_key,
                                    cfe_uzpipfe *c);

typedef struct cfe_uzpipfe_ciphertext
{
    cfe_vec_G1 *c1;
    cfe_vec_G1 *c2;
} cfe_uzpipfe_ciphertext;

/**
 * Initializes the struct which represents the ciphertext.
 *
 * @param cipher A pointer to an uninitialized cfe_uzpipfe_ciphertext struct
 * @param c A pointer to an instance of the scheme (*initialized* cfe_uzpipfe
 * struct)
 */
void cfe_uzpipfe_ciphertext_init(cfe_uzpipfe_ciphertext *cipher, cfe_uzpipfe *c);

/**
 * Frees the memory occupied by the struct members. It does
 * not free the memory occupied by the struct itself.
 *
 * @param cipher A pointer to an *initialized* cfe_uzpipfe_ciphertext struct
 */
void cfe_uzpipfe_ciphertext_free(cfe_uzpipfe_ciphertext *cipher, cfe_uzpipfe *c);

/**
 * Encrypts input vector x with the provided master secret key. It returns a
 * ciphertext struct. If encryption failed, an error is returned.
 *
 * @param cipher A pointer to an initialized cfe_uzpipfe_ciphertext struct
 * (the resulting ciphertext will be stored here)
 * @param x A pointer to the plaintext vector
 * @param sec_key A pointer to the master secret key
 * @param pub_key A pointer to the public key
 * @param c A pointer to an instance of the scheme (*initialized* cfe_uzpipfe
 * struct)
 * @return Error code
 */
cfe_error cfe_uzpipfe_encrypt(cfe_uzpipfe_ciphertext *cipher, cfe_vec *x, cfe_vec *w,
                              cfe_uzpipfe_sec_key *sec_key, cfe_uzpipfe_pub_key *pub_key,
                              cfe_uzpipfe *c);

/**
 * Accepts the encrypted vector and functional encryption key. It returns the
 * inner product of x and y. If decryption failed, an error is returned.
 *
 * @param res The result of the decryption (the value will be stored here)
 * @param cipher A pointer to the ciphertext vector
 * @param fe_key The functional encryption key
 * @param c A pointer to an instance of the scheme (*initialized* cfe_fhipe
 * struct)
 * @return Error code
 */
cfe_error cfe_uzpipfe_decrypt(mpz_t res, cfe_uzpipfe_ciphertext *cipher,
                              cfe_uzpipfe_fe_key *fe_key, cfe_uzpipfe_pub_key *pub_key,
                              cfe_uzpipfe *c);

#endif
