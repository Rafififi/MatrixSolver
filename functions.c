#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"


void ReadMMtoCSR(const char *filename, CSRMatrix *aMatrix)
{
    FILE *fileName = fopen(filename, "r");
    fscanf(fileName, "%*[^\n]\n"); // Skip the first line
    //skip lines that start with a %
    char c = getc(fileName);
    while (c == '%')
    {
        fscanf(fileName, "%*[^\n]\n"); // Skip the line
        c = getc(fileName); // Get the next character
    }
    ungetc(c, fileName); // Goes back 1 character because the while loop goes 1 character too far

    // Read the number of rows, columns and non-zero elements
    fscanf(fileName, "%d %d %d", &aMatrix->num_rows, &aMatrix->num_cols, &aMatrix->num_non_zeros);
    
    

    int *rows = (int *)calloc(aMatrix->num_non_zeros, sizeof(int));
    int *cols = (int *)calloc(aMatrix->num_non_zeros,  sizeof(int));
    double *data = (double *)calloc(aMatrix->num_non_zeros,  sizeof(double));
    for (int i = 0; i<aMatrix->num_non_zeros; i++)
    {
        fscanf(fileName, "%d %d %lf", &rows[i], &cols[i], &data[i]);
        c = getc(fileName);
    }


    printf("Numbers have been read\n");


    aMatrix->row_ptr = (int *)calloc((aMatrix->num_rows + 1), sizeof(int));
    aMatrix->col_ind = (int *)calloc(aMatrix->num_non_zeros, sizeof(int));
    aMatrix->csr_data = (double *)calloc(aMatrix->num_non_zeros, sizeof(double));

    int value = 0;
    int col = 0;
    for (int row = 0; row <= aMatrix->num_rows; row++)
    {
        for (int i = 0; i < aMatrix->num_non_zeros; i++)
        {
            if (rows[i] == row)
            {
                value++;
                aMatrix->row_ptr[row] = value;
                aMatrix->col_ind[col] = cols[i] - 1;
                aMatrix->csr_data[col] = data[i];
                col++;
            }
        }
    }

    for (int i = 0; i <= aMatrix->num_rows; i++)
    {
        printf("%d ", aMatrix->row_ptr[i]);
    }
    printf("\n");
    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        printf("col index: %d ", aMatrix->col_ind[i]);
        printf("csr Data: %lf ", aMatrix->csr_data[i]);
        printf("\n");
    }
    
    free(rows);
    free(cols);
    free(data);
    
    fclose(fileName);

}

void solver(const CSRMatrix AMatrix, double *b, double *x)
{
    int n = AMatrix.num_rows;
    double *x_new = malloc(n * sizeof(double));
    double tol = 1e-6;
    int max_iter = 1000;

    for (int i = 0; i < n; i++)
    {
        x_new[i] = x[i];
    }

    for (int iter = 0; iter < max_iter; iter++)
    {
        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;
            for (int j = AMatrix.row_ptr[i]; j < AMatrix.row_ptr[i + 1]; j++)
            {
                if (AMatrix.col_ind[j] != i)
                {
                    sum += AMatrix.csr_data[j] * x[AMatrix.col_ind[j]];
                }
            }
            x_new[i] = (b[i] - sum) / AMatrix.csr_data[AMatrix.row_ptr[i] + i];
        }

        double error = 0.0;
        for (int i = 0; i < n; i++)
        {
            error += fabs(x_new[i] - x[i]);
        }

        if (error < tol)
        {
            break;
        }

        for (int i = 0; i < n; i++)
        {
            x[i] = x_new[i];
        }
    }

    free(x_new);
}


void spmv_csr(const CSRMatrix *AMatrix, const double *x, double *y)
{
    double *product = (double *)calloc(AMatrix->num_non_zeros, sizeof(double));
    for (int row = 0; row < AMatrix->num_rows; row++)
    {

        for (int j = AMatrix->row_ptr[row]; j < AMatrix->row_ptr[row + 1]; j++)
        {
            product[row] += AMatrix->csr_data[j] * x[AMatrix->col_ind[j]];
        }
        y[row] = product[row];
    }
    

    free(product);
    printf("spmv done\n");
}

void compute_residual(const CSRMatrix AMatrix, const double *b, const double *x, double *r)
{
    double *product = (double *)calloc(AMatrix.num_non_zeros, sizeof(double));
    spmv_csr(&AMatrix, x, product);

    for (int i = 0; i < AMatrix.num_rows; i++)
    {
        r[i] = product[i] - b[i];
    }
    
    free(product);
    printf("residual done\n");

}

double compute_norm(const double *r, int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
    {
        norm += r[i] * r[i];
    }
    norm = sqrt(norm);
    printf("norm done\n");
    return norm;

}

void free_csr_matrix(CSRMatrix *aMatrix)
{
    free(aMatrix->row_ptr);
    free(aMatrix->col_ind);
    free(aMatrix->csr_data);
    aMatrix->row_ptr = NULL;
    aMatrix->col_ind = NULL;
    aMatrix->csr_data = NULL;
    aMatrix->num_rows = 0;
    aMatrix->num_cols = 0;
    aMatrix->num_non_zeros = 0;
    printf("free done\n");
}
