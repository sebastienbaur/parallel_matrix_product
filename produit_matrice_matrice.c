#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include "fonctions.h"

void compute_product_matrix_matrix(double A,double B,double result){
    if(A[0][1] == B[0][0] and A[0][0]==result[0][0] and result[0][1]==B[0][1]){
        int i;
        int j;
        int k;
        for(i=1;i<=result[0][0];i++){
            for(j=1;j<=result[0][1];j++){
                result[i][j] = 0
                k = 0;
                for(k=0;k<A_column_number;k++){
                    result[i][j] += A[i][k]*B[k][j];
                }
            }
        }
    }
    else{
        printf("Invalid dimensions : (%d,%d)*(%d,%d), (%d,%d) \n\n", A[0][0],A[0][1],B[0][0],B[0][1],result[0][0], result[0][1]);
    }
}
