
#ifndef CIFER_UABIPFE_H
#define CIFER_UABIPFE_H


#include <gmp.h>
#include <amcl/fp12_BN254.h>
#include <amcl/ecp_BN254.h>
#include <amcl/ecp2_BN254.h>

#include "cifer/data/vec.h"
#include "cifer/data/mat.h"
#include "cifer/data/vec_curve.h"

#include "cifer/abe/policy.h"
#include "cifer/internal/errors.h"

/**
 * cfe_uabipfe contains the shared choice for parameters on which
 * the functionality of the scheme depend.
 */
typedef struct cfe_uabipfe
{
    size_t m1;
    mpz_t bound_x;
    mpz_t bound_y;
    mpz_t order;
} cfe_uabipfe;

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
cfe_error cfe_uabipfe_init(cfe_uabipfe *c, size_t m1, mpz_t bound_x, mpz_t bound_y);

/**
 * Reconstructs the scheme with the same configuration parameters from
 * an already existing fhipe scheme instance.
 *
 * @param res A pointer to an uninitialized cfe_uzpipfe struct
 * @param c A pointer to an instance of the scheme (*initialized* cfe_uzpipfe
 * struct)
 */
void cfe_uabipfe_copy(cfe_uabipfe *res, cfe_uabipfe *c);

/**
 * Frees the memory occupied by the struct members. It does not free the
 * memory occupied by the struct itself.
 *
 * @param c A pointer to an instance of the scheme (*initialized* cfe_fhipe
 * struct)
 */
void cfe_uabipfe_free(cfe_uabipfe *c);

/**
 * cfe_uabipfe_master_pub_key represents a master public key in fhipe scheme.
 */
typedef struct cfe_uabipfe_master_pub_key
{
    ECP_BN254 g1;
    ECP2_BN254 g2;
    FP12_BN254 gT;
    cfe_vec_G1 f[3];
    cfe_vec_G1 h12;
    cfe_vec_G1 h3;
    ECP_BN254 mu;
} cfe_uabipfe_master_pub_key;

/**
 * Initializes the struct which represents the master public key in uabipfe.
 *
 * @param master_pub_key A pointer to an uninitialized cfe_uabipfe_master_pub_key struct
 * 
 */
void cfe_uabipfe_master_pub_key_init(cfe_uabipfe_master_pub_key *master_pub_key);

/**
 * Frees the memory occupied by the struct members. It does
 * not free the memory occupied by the struct itself.
 * 
 * @param master_pub_key A pointer to an *initialized* cfe_uabipfe_master_pub_key struct
 */
void cfe_uabipfe_master_pub_key_free(cfe_uabipfe_master_pub_key *master_pub_key);

/**
 * cfe_uabipfe_master_pub_key represents a master public key in fhipe scheme.
 */
typedef struct cfe_uabipfe_master_sec_key
{
    ECP_BN254 g1;
    ECP2_BN254 g2;
    FP12_BN254 gT;
    cfe_vec f_star[3];
    cfe_vec h_star[3];
    mpz_t z;
} cfe_uabipfe_master_sec_key;

/**
 * Initializes the struct which represents the master public key in uabipfe.
 *
 * @param master_sec_key A pointer to an uninitialized cfe_uabipfe_master_sec_key struct
 * 
 */
void cfe_uabipfe_master_sec_key_init(cfe_uabipfe_master_sec_key *master_sec_key);

/**
 * Frees the memory occupied by the struct members. It does
 * not free the memory occupied by the struct itself.
 * 
 * @param master_sec_key A pointer to an *initialized* cfe_uabipfe_master_sec_key struct
 */
void cfe_uabipfe_master_sec_key_free(cfe_uabipfe_master_sec_key *master_sec_key);

/**
 * Generates a master public key and a master secret key for the scheme. It returns an
 * error if generating one of the parts of the secret key failed.
 *
 * @param master_pub_key A pointer to a cfe_uabipfe_master_pub_key struct 
 * @param master_sec_key A pointer to a cfe_uabipfe_master_sec_key struct 
 * @param c A pointer to an instance of the scheme (*initialized*
 * cfe_uabipfe struct)
 * @return Error code
 */
cfe_error cfe_uabipfe_generate_keys(cfe_uabipfe_master_pub_key *master_pub_key, cfe_uabipfe_master_sec_key *master_sec_key,
                                cfe_uabipfe *c);


typedef struct cfe_uabipfe_ciphertext
{
    cfe_vec_G1 c_j;
    cfe_vec_G1 c_ipfe;
    FP12_BN254 *t;
} cfe_uabipfe_ciphertext;


/**
 * Initializes the struct which represents the ciphertext.
 *
 * @param cipher A pointer to an uninitialized cfe_uabipfe_ciphertext struct
 * @param c A pointer to an instance of the scheme (*initialized* cfe_uabipfe
 * struct)
 */
void cfe_uabipfe_ciphertext_init(cfe_uabipfe_ciphertext *cipher, cfe_uabipfe *c);

/**
 * Frees the memory occupied by the struct members. It does
 * not free the memory occupied by the struct itself.
 *
 * @param cipher A pointer to an *initialized* cfe_uabipfe_ciphertext struct
 * @param c A pointer to an instance of the scheme (*initialized* cfe_uabipfe
 * struct)
 */
void cfe_uabipfe_ciphertext_free(cfe_uabipfe_ciphertext *cipher, cfe_uabipfe *c);

/**
 * Encrypts input vector x with the provided master secret key. It returns a
 * ciphertext struct. If encryption failed, an error is returned.
 *
 * @param cipher A pointer to an initialized cfe_uabipfe_ciphertext struct
 * (the resulting ciphertext will be stored here)
 * @param master_pub_key A pointer to the master public key
 * @param x A pointer to the plaintext vector x
 * @param c A pointer to an instance of the scheme (*initialized* cfe_uabipfe
 * struct)
 * @return Error code
 */
cfe_error cfe_uabipfe_encrypt(cfe_uabipfe_ciphertext *cipher, cfe_uabipfe_master_pub_key *master_pub_key, 
                                cfe_vec *x, cfe_uabipfe *c);

typedef struct cfe_uabipfe_fe_key
{
    cfe_vec_G2 k_j;
    cfe_vec_G2 k_ipfe;
} cfe_uabipfe_fe_key;


/**
 * Initializes the struct which represents the functional encryption key in fhipe.
 *
 * @param fe_key A pointer to an uninitialized cfe_uabipfe_FE_key struct
 * @param c A pointer to an instance of the scheme (*initialized* cfe_uabipfe
 * struct)
 */
void cfe_uabipfe_fe_key_init(cfe_uabipfe_fe_key *fe_key, cfe_uabipfe *c);

/**
 * Frees the memory occupied by the struct members. It does
 * not free the memory occupied by the struct itself.
 *
 * @param fe_key A pointer to an *initialized* cfe_uzpipfe_FE_key struct
 * @param c A pointer to an instance of the scheme (*initialized* cfe_uabipfe
 * struct)
 */
void cfe_uabipfe_fe_key_free(cfe_uabipfe_fe_key *fe_key,  cfe_uabipfe *c);

/**
 * Takes a master secret key and input vector y, and derives the functional
 * encryption key. In case the key could not be derived, it returns an error.
 *
 * @param fe_key A pointer to a cfe_fhipe_fe_key struct (the functional
 * encryption key will be stored here)
 * @param y A pointer to the inner product vector
 * @param sec_key A pointer to the master secret key
 * @param c A pointer to an instance of the scheme (*initialized* cfe_fhipe
 * struct)
 * @return Error code
 */
cfe_error cfe_uabipfe_derive_fe_key(cfe_uabipfe_fe_key *fe_key, cfe_vec *y, 
                                    cfe_uabipfe_master_sec_key *sec_key, cfe_uabipfe *c);

/**
 * Accepts the encrypted vector and functional encryption key. It returns the
 * inner product of x and y. If decryption failed, an error is returned.
 *
 * @param res The result of the decryption (the value will be stored here)
 * @param cipher A pointer to the ciphertext vector
 * @param fe_key The functional encryption key
 * @param y A pointer to the plaintext vector y
 * @param c A pointer to an instance of the scheme (*initialized* cfe_fhipe
 * struct)
 * @return Error code
 */
cfe_error cfe_uabipfe_decrypt(mpz_t res, cfe_uabipfe_ciphertext *cipher,
                              cfe_uabipfe_master_pub_key *master_pub_key,
                              cfe_uabipfe_fe_key *fe_key, cfe_vec *y,
                              cfe_uabipfe *c);

#endif