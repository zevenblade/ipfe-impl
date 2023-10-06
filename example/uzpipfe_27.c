#include <stdio.h>
#include <gmp.h>
#include "cifer/innerprod/fullysec/uzpipfe_27.h"
#include "cifer/sample/uniform.h"

void main()
{
    cfe_error err;
    size_t m1 = 100;
    size_t m2 = 30;
    mpz_t zero, one, bound, bound_neg, xy_check, xy, wv_check;
    cfe_vec x, y, w, v;

    mpz_inits(zero, one, bound, bound_neg, xy_check, xy, wv_check, NULL);
    cfe_vec_inits(m1, &x, &y, NULL);
    cfe_vec_inits(m2, &w, &v, NULL);
    
    mpz_set_ui(zero, 0);
    mpz_set_ui(one, 1);

    mpz_set_ui(bound, 2);
    mpz_pow_ui(bound, bound, 1);
    mpz_neg(bound_neg, bound);

    cfe_uzpipfe_27 uzpipfe;
    err = cfe_uzpipfe_init_27(&uzpipfe, m1, m2, bound, bound);

    cfe_uzpipfe_pub_key_27 pub_key;
    cfe_uzpipfe_sec_key_27 sec_key;
    cfe_uzpipfe_master_pub_key_init_27(&pub_key, &uzpipfe);
    cfe_uzpipfe_master_sec_key_init_27(&sec_key, &uzpipfe);
    err = cfe_uzpipfe_generate_keys_27(&sec_key, &pub_key, &uzpipfe);
    
    cfe_uniform_sample_range_vec(&x, bound_neg, bound);
    cfe_uniform_sample_range_vec(&y, bound_neg, bound);
    cfe_vec_dot(xy_check, &x, &y);

    printf("x: ");
    cfe_vec_print(&x);
    printf("\ny: ");
    cfe_vec_print(&y);
    gmp_printf("\nx * y: %Zd\n", xy_check);

    cfe_uniform_sample_range_vec(&w, bound_neg, bound);
    cfe_uniform_sample_range_vec(&v, bound_neg, bound);
    cfe_vec_set(&w, zero, (m2-1));
    cfe_vec_set(&v, zero, (m2-1));
    cfe_vec_dot(wv_check, &w, &v);
    mpz_neg(wv_check, wv_check);
    cfe_vec_set(&w, one, (m2-1));
    cfe_vec_set(&v, wv_check, (m2-1));
    cfe_vec_dot(wv_check, &w, &v);

    printf("w: ");
    cfe_vec_print(&w);
    printf("\nv: ");
    cfe_vec_print(&v);
    gmp_printf("\nw * v: %Zd\n", wv_check);

    cfe_uzpipfe_ciphertext_27 cipher;
    cfe_uzpipfe_ciphertext_init_27(&cipher, &uzpipfe);
    err = cfe_uzpipfe_encrypt_27(&cipher, &x, &w, &sec_key, &pub_key, &uzpipfe);

    cfe_uzpipfe_fe_key_27 fe_key;
    cfe_uzpipfe_fe_key_init_27(&fe_key, &uzpipfe);
    err = cfe_uzpipfe_derive_fe_key_27(&fe_key, &y, &v, &sec_key, &pub_key,
                                    &uzpipfe);

    cfe_uzpipfe_27 decryptor;
    cfe_uzpipfe_copy_27(&decryptor, &uzpipfe);
    err = cfe_uzpipfe_decrypt_27(xy, &cipher, &fe_key, &pub_key, &uzpipfe);
    
    
    gmp_printf("FE x * y: %Zd\n", xy);

    if(err == CFE_ERR_NONE){
        printf("No error\n");
    } else if(err = CFE_ERR_DLOG_NOT_FOUND){
        printf("Dlog not found\n");
    } else{
        printf("Other error\n");
    }
    

    return ;
}