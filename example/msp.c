#include <stdlib.h>
#include <stdarg.h>

#include <amcl/big_256_56.h>
#include <amcl/pair_BN254.h>

#include "cifer/abe/uabipfe.h"
#include "cifer/internal/big.h"
#include "cifer/sample/uniform.h"
#include "cifer/internal/dlog.h"

#include "cifer/internal/common.h"
#include "cifer/internal/big.h"
#include "cifer/internal/hash.h"
#include "cifer/sample/uniform.h"

#include "cifer/abe/fame.h"
#include <stdio.h>
#include <gmp.h>
#include "cifer/abe/uabipfe.h"
#include "cifer/sample/uniform.h"

void main()
{
    cfe_error err;
    
    mpz_t p;
    cfe_vec one_vec, alpha;
    cfe_mat msp_trans;

    mpz_inits(p, NULL);
    mpz_from_BIG_256_56(p, (int64_t *)CURVE_Order_BN254);

    // create a msp structure out of a boolean expression representing the
    // policy specifying which attributes are needed to decrypt the ciphertext
    //char bool_exp[] = "(1 OR 2) OR 3";
    //char bool_exp[] = "(1 AND 2) AND 3";
    char bool_exp[] = "1";
    //char bool_exp[] = "(5 OR 3) AND ((2 OR 4) OR (1 AND 6))";
    size_t bool_exp_len = 1;  // length of the boolean expression string
    cfe_msp msp;
    err = cfe_boolean_to_msp(&msp, bool_exp, bool_exp_len, false);
    cfe_mat_print(&msp.mat);

    for(size_t i = 0; i < 3; i++){
        printf("%d ", msp.row_to_attrib[i]);
    }
    printf("\n");

    cfe_vec_init(&one_vec, msp.mat.cols);
    mpz_set_ui(one_vec.vec[0], 1);
    cfe_vec_print(&one_vec);
    printf("\n");

    cfe_mat_init(&msp_trans, msp.mat.cols, msp.mat.rows);
    cfe_mat_transpose(&msp_trans, &msp.mat);


    cfe_gaussian_elimination_solver(&alpha, &msp_trans, &one_vec, p);
    cfe_vec_print(&alpha);
    printf("\n");
}