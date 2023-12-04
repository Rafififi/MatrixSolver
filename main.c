#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "functions.h"

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        printf("Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    // <Handle the inputs here>
    const char *filename = argv[1];
    //If the file is not in the directory
    

    CSRMatrix aMatrix;
    ReadMMtoCSR(filename, &aMatrix);

    CSRMatrix *nonConsantMatrix = (CSRMatrix *)malloc(sizeof(CSRMatrix)); // allocate memory for the matrix

    // check if memory allocation was successful
    if (nonConsantMatrix == NULL)
    {
        printf("Error: memory allocation failed\n");
        exit(1);
    }

    nonConsantMatrix->num_rows = aMatrix.num_rows;
    nonConsantMatrix->num_cols = aMatrix.num_cols;
    nonConsantMatrix->num_non_zeros = aMatrix.num_non_zeros;
    // allocate memory for the matrix
    nonConsantMatrix->csr_data = (double *)malloc(nonConsantMatrix->num_non_zeros * sizeof(double));
    nonConsantMatrix->col_ind = (int *)malloc(nonConsantMatrix->num_non_zeros * sizeof(int));
    nonConsantMatrix->row_ptr = (int *)malloc((nonConsantMatrix->num_rows + 1) * sizeof(int));
    
    if (nonConsantMatrix->csr_data == NULL || nonConsantMatrix->col_ind == NULL || nonConsantMatrix->row_ptr == NULL)
    {
        printf("Error: memory allocation failed\n");
        exit(1);
    }

    memcpy(nonConsantMatrix->csr_data, aMatrix.csr_data, nonConsantMatrix->num_non_zeros * sizeof(double));
    memcpy(nonConsantMatrix->col_ind, aMatrix.col_ind, nonConsantMatrix->num_non_zeros * sizeof(int));
    memcpy(nonConsantMatrix->row_ptr, aMatrix.row_ptr, (nonConsantMatrix->num_rows + 1) * sizeof(int));
    // we need to copy as if we don't then the original matrix will be changed

    CSRTranspose(nonConsantMatrix);

    puts("Matrix read successfully");
    // Initializing all the vector b (in Ax=b)
    double *bMatrix = (double *)malloc(aMatrix.num_cols * sizeof(double));
    double *xMatrix = (double *)malloc(aMatrix.num_cols * sizeof(double));
    if (bMatrix == NULL || xMatrix == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }

    // Set all elements of b to 1
    for (int i = 0; i < aMatrix.num_cols; ++i)
    {
        bMatrix[i] = 1.0;
    }

    for (int i = 0; i < aMatrix.num_cols; ++i)
    {
        xMatrix[i] = 1.0;
    }
    triangularCheck(aMatrix);
   

    solver(aMatrix, bMatrix, xMatrix, *nonConsantMatrix);
    puts("Solver finished");

    printf("XMatix: \n");
    
    

    double *residual = (double *)malloc(aMatrix.num_cols * sizeof(double));
    if (residual == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }

    computeResidual(aMatrix, bMatrix, xMatrix, residual, *nonConsantMatrix);
    

    
    double norm = computeNorm(residual, aMatrix.num_cols);
    printf(" \nNorm: %e\n", norm);
    // <The rest of your code goes here>
    freeCSRMatrix(&aMatrix);
    freeCSRMatrix(nonConsantMatrix);
    free(nonConsantMatrix);


    free(bMatrix);
    free(xMatrix);
    free(residual);
    bMatrix = NULL;
    xMatrix = NULL;
    residual = NULL;

    printf("\nDone!\n");
    return 0;
}