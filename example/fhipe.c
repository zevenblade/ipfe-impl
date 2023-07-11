#include <stdio.h>
#include <gmp.h>
#include "cifer/innerprod/fullysec/fhipe.h"
#include "cifer/sample/uniform.h"

void main()
{
    size_t l = 3;
    mpz_t bound, bound_neg, xy_check, xy;
    mpz_inits(bound, bound_neg, xy_check, xy, NULL);
    mpz_set_ui(bound, 2);
    mpz_pow_ui(bound, bound, 1);
    mpz_neg(bound_neg, bound);

    cfe_fhipe fhipe;
    cfe_error err = cfe_fhipe_init(&fhipe, l, bound, bound);

    // generate master key
    cfe_fhipe_sec_key sec_key;
    cfe_fhipe_master_key_init(&sec_key, &fhipe);
    err = cfe_fhipe_generate_master_key(&sec_key, &fhipe);

    // sample a vector that will be encrypted and an inner
    // product vector
    cfe_vec x, y;
    cfe_vec_inits(l, &x, &y, NULL);
    cfe_uniform_sample_range_vec(&x, bound_neg, bound);
    cfe_uniform_sample_range_vec(&y, bound_neg, bound);
    cfe_vec_dot(xy_check, &x, &y);

    cfe_vec_print(&x);
    cfe_vec_print(&y);
    printf("\n");

    // derive a functional key for vector y
    cfe_fhipe_fe_key FE_key;
    cfe_fhipe_fe_key_init(&FE_key, &fhipe);
    err = cfe_fhipe_derive_fe_key(&FE_key, &y, &sec_key, &fhipe);

    // encrypt the vector
    cfe_fhipe_ciphertext cipher;
    cfe_fhipe_ciphertext_init(&cipher, &fhipe);
    err = cfe_fhipe_encrypt(&cipher, &x, &sec_key, &fhipe);

    // simulate a decryptor
    cfe_fhipe decryptor;
    cfe_fhipe_copy(&decryptor, &fhipe);
    // decryptor decrypts the inner-product without knowing
    // vectors x and y
    err = cfe_fhipe_decrypt(xy, &cipher, &FE_key, &decryptor);

    // check the correctness of the result
    gmp_printf("%Zd\n",  xy);
    
    if(err == CFE_ERR_NONE){
        printf("No error\n");
    } else if(err = CFE_ERR_DLOG_NOT_FOUND){
        printf("Dlog not found\n");
    } else{
        printf("Other error\n");
    }

    // clean up
    mpz_clears(bound, bound_neg, xy_check, xy, NULL);
    cfe_vec_frees(&x, &y, NULL);
    cfe_fhipe_free(&fhipe);
    cfe_fhipe_free(&decryptor);
    cfe_fhipe_master_key_free(&sec_key);
    cfe_fhipe_fe_key_free(&FE_key);
    cfe_fhipe_ciphertext_free(&cipher);
}