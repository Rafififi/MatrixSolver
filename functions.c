#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"



void ReadMMtoCSR(const char *filename, CSRMatrix *aMatrix)
{
    FILE *fileName = fopen(filename, "r");
    if (fileName == NULL)
    {
        printf("File not found\n");
        exit(1);
    }

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


    MatrixElement *elements = (MatrixElement *)calloc(aMatrix->num_non_zeros, sizeof(MatrixElement));
    if (elements == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }
    /*
    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        printf("%d %d %lf\n", rows[i], cols[i], data[i]);
    }
    */

    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        elements[i].row = rows[i];
        elements[i].col = cols[i];
        elements[i].data = data[i];
    }

    //using qsort to sort the rows and applying the same permutation to cols and data
    qsort(elements, aMatrix->num_non_zeros, sizeof(MatrixElement), compare);
    //qsort is a built in function that sorts the array elements in ascending order
    printf("qsort done\n");


    aMatrix->row_ptr = (int *)calloc((aMatrix->num_rows + 1), sizeof(int));
    aMatrix->col_ind = (int *)calloc(aMatrix->num_non_zeros, sizeof(int));
    aMatrix->csr_data = (double *)calloc(aMatrix->num_non_zeros, sizeof(double));

    if (aMatrix->row_ptr == NULL || aMatrix->col_ind == NULL || aMatrix->csr_data == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }

    int current_row = 0;
    int value = 0;
    int col = 0;

    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        while (current_row < elements[i].row && current_row <= aMatrix->num_rows)
        {
            aMatrix->row_ptr[current_row] = value;
            current_row++;
        }
        value++;
        aMatrix->col_ind[col] = elements[i].col - 1;
        aMatrix->csr_data[col] = elements[i].data;
        col++;
    }

    // Fill the rest of row_ptr with the total number of non-zero elements
    while (current_row <= aMatrix->num_rows)
    {
        aMatrix->row_ptr[current_row] = value;
        current_row++;
    }
    if (aMatrix->num_non_zeros < 50)
    {
        printf("row_ptr: ");
        for (int i = 0; i <= aMatrix->num_rows; i++)
        {
            printf("%d ", aMatrix->row_ptr[i]);
        }
        printf("\n");
        printf("col_ind: ");
        for (int i = 0; i < aMatrix->num_non_zeros; i++)
        {
            printf("%d ", aMatrix->col_ind[i]);
        }
        printf("\n");
        printf("csr_data: ");
        for (int i = 0; i < aMatrix->num_non_zeros; i++)
        {
            printf("%lf ", aMatrix->csr_data[i]);
        }
        printf("\n");
    }
    
    free(rows);
    free(cols);
    free(data);
    free(elements);
    elements = NULL;
    fclose(fileName);

}

void solver(const CSRMatrix AMatrix, double *b, double *x)
{
    double *diagonal = (double *)calloc(AMatrix.num_rows, sizeof(double));
    int maxIterations = 1000;
    if (diagonal == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < AMatrix.num_rows; i++)
    {
        for (int j = AMatrix.row_ptr[i]; j < AMatrix.row_ptr[i + 1]; j++)
        {
            if (AMatrix.col_ind[j] == i)
            {
                diagonal[i] = AMatrix.csr_data[j];
                if (diagonal[i] == 0)
                {
                    printf("The matrix is singular\n");
                    exit(1);
                }
            }
        }
    }

    for (int i = 0; i < maxIterations; i++)
    {
        for (int j = 0; j < AMatrix.num_rows; j++)
        {
            double sum = 0;
            for (int k = AMatrix.row_ptr[j]; k < AMatrix.row_ptr[j + 1]; k++)
            {
                if (AMatrix.col_ind[k] != j)
                {
                    sum += AMatrix.csr_data[k] * x[AMatrix.col_ind[k]];
                }
            }
            x[j] = (b[j] - sum) / diagonal[j];
        }
    }

    
}


void spmv_csr(const CSRMatrix *AMatrix, const double *x, double *y)
{
    int matrixType = triangularcheck(*AMatrix);
    if (matrixType == 2)
    {
        //If the value is not on the diagonal I need to the x value that is equal to the col_index of the transpose of that value
        double *product = (double *)calloc(AMatrix->num_non_zeros, sizeof(double));
        if (product == NULL)
        {
            printf("Memory allocation failed\n");
            exit(1);
        }
        printf("Product:");
        for (int row = 0; row < AMatrix->num_rows; row++)
        {
            for (int j = AMatrix->row_ptr[row]; j < AMatrix->row_ptr[row + 1]; j++)
            {
                product[row] += AMatrix->csr_data[j] * x[AMatrix->col_ind[j]];
            }
            y[row] = product[row];
            printf(" %lf", y[row]);
        }
        free(product);
    }

    else
    {
        double *product = (double *)calloc(AMatrix->num_non_zeros, sizeof(double));
        if (product == NULL)
        {
            printf("Memory allocation failed\n");
            exit(1);
        }
        printf("Product:");
        for (int row = 0; row < AMatrix->num_rows; row++)
        {

            for (int j = AMatrix->row_ptr[row]; j < AMatrix->row_ptr[row + 1]; j++)
            {
                product[row] += AMatrix->csr_data[j] * x[AMatrix->col_ind[j]];
            }
            y[row] = product[row];
            printf(" %lf", y[row]);
        }
        free(product);
    } 
    printf("\nspmv done\n");
}

void compute_residual(const CSRMatrix AMatrix, const double *b, const double *x, double *r)
{
    double *product = (double *)calloc(AMatrix.num_non_zeros, sizeof(double));

    if (product == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }
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

int compare(const void *a, const void *b) //this is for qsort to sort the rows and applying the same permutation to cols and data
{
    MatrixElement *elementA = (MatrixElement *)a;
    MatrixElement *elementB = (MatrixElement *)b;
    return elementA->row - elementB->row;
    //It works by returning a negative, zero, or positive integer, 
    //depending on whether the first argument is less than, equal to, or greater than the second argument.
}

int triangularcheck(const CSRMatrix AMatrix)
{
    int upper = 1;
    for (int i = 0; i < AMatrix.num_rows; i++)
    {
        for (int j = AMatrix.row_ptr[i]; j < AMatrix.row_ptr[i + 1]; j++)
        {
            if (AMatrix.col_ind[j] < i)
            {
                upper = 0;
                break;
            }
        }
        if (upper == 0)
        {
            break;
        }
    }
    
    int lower = 1;
    for (int i = 0; i < AMatrix.num_rows; i++)
    {
        for (int j = AMatrix.row_ptr[i]; j < AMatrix.row_ptr[i + 1]; j++)
        {
            if (AMatrix.col_ind[j] > i)
            {
                lower = 0;
                break;
            }
        }
        if (lower == 0)
        {
            break;
        }
    }

    if (upper == 1 && lower == 1)
    {
        printf("The matrix not triangular\n");
        return 0;
    }
    else if (upper == 1)
    {
        printf("The matrix is upper triangular\n");
        return 1;
    }
    else if (lower == 1)
    {
        printf("The matrix is lower triangular\n");
        return 2;
    }
    else
    {
        printf("The matrix is not triangular\n");
        return 3;
    }
}