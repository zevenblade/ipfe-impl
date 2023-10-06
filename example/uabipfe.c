#include <stdio.h>
#include <gmp.h>
#include "cifer/abe/uabipfe.h"
#include "cifer/sample/uniform.h"

void main()
{
    cfe_error err;
    size_t m1 = 10;
    mpz_t zero, one, bound, bound_neg, xy_check, xy, wv_check;
    cfe_vec x, y;

    mpz_inits(zero, one, bound, bound_neg, xy_check, xy, wv_check, NULL);
    cfe_vec_inits(m1, &x, &y, NULL);

    mpz_set_ui(zero, 0);
    mpz_set_ui(one, 1);

    mpz_set_ui(bound, 2);
    mpz_pow_ui(bound, bound, 1);
    mpz_neg(bound_neg, bound);

    cfe_uabipfe uabipfe;
    err = cfe_uabipfe_init(&uabipfe, m1, bound, bound);

    cfe_uabipfe_master_pub_key master_pub_key;
    cfe_uabipfe_master_pub_key_init(&master_pub_key);
    cfe_uabipfe_master_sec_key master_sec_key;
    cfe_uabipfe_master_sec_key_init(&master_sec_key);
    err = cfe_uabipfe_generate_keys(&master_pub_key, &master_sec_key, &uabipfe);

    cfe_uniform_sample_range_vec(&x, bound_neg, bound);
    cfe_uniform_sample_range_vec(&y, bound_neg, bound);
    cfe_vec_dot(xy_check, &x, &y);
    
    printf("x: ");
    cfe_vec_print(&x);
    printf("\ny: ");
    cfe_vec_print(&y);
    gmp_printf("\nx * y: %Zd\n", xy_check);

    // create a msp structure out of a boolean expression representing the
    // policy specifying which attributes are needed to decrypt the ciphertext
    //char bool_exp[] = "1";
    // char bool_exp[] = "(1 OR 2) OR 3";
    // char bool_exp[] = "(1 AND 2) AND 3";
    //char bool_exp[] = "(1 AND 2) OR 3";
     char bool_exp[] = "(5 OR 3) AND ((2 OR 4) OR (1 AND 6))";
    size_t bool_exp_len = 36; // length of the boolean expression string
    cfe_msp msp;
    err = cfe_boolean_to_msp(&msp, bool_exp, bool_exp_len, false);
    // cfe_mat_print(&msp.mat);
    
    cfe_uabipfe_ciphertext cipher;
    cfe_uabipfe_ciphertext_init(&cipher, &msp, &uabipfe);
    err = cfe_uabipfe_encrypt(&cipher, &master_pub_key, &x, &msp, &uabipfe);
    
    cfe_uabipfe_fe_key fe_key;
    cfe_uabipfe_fe_key_init(&fe_key, &msp, &uabipfe);
    err = cfe_uabipfe_derive_fe_key(&fe_key, &y, &msp, &master_sec_key, &uabipfe);

    int owned_attrib[] = {3, 1, 6};

    cfe_uabipfe decryptor;
    cfe_uabipfe_copy(&decryptor, &uabipfe);
    err = cfe_uabipfe_decrypt(xy, &cipher, &master_pub_key, &fe_key, &y, &msp, owned_attrib, 3, &decryptor);
    
    
    gmp_printf("FE x * y: %Zd\n", xy);

    if (err == CFE_ERR_NONE)
    {
        printf("No error\n");
    }
    else if (err == CFE_ERR_DLOG_NOT_FOUND)
    {
        printf("Dlog not found\n");
    }
    else if (err == CFE_ERR_INSUFFICIENT_KEYS)
    {
        printf("Insufficient attributes\n");
    }
    else
    {
        printf("Other error\n");
    }
    
}