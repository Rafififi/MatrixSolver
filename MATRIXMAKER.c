#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

typedef struct
{
    double *csr_data;  // Array of non-zero values
    int *col_ind;      // Array of column indices
    int *row_ptr;      // Array of row pointers
    int num_non_zeros; // Number of non-zero elements
    int num_rows;      // Number of rows in matrix
    int num_cols;      // Number of columns in matrix
} CSRMatrix;

void ReadMMtoCSR(const char *filename, CSRMatrix *aMatrix)
{
    FILE *fileName = fopen(filename, "r");
    fscanf(fileName, "%*[^\n]\n"); // Skip the first line
    // skip lines that start with a %
    char c = getc(fileName);
    while (c == '%')
    {
        fscanf(fileName, "%*[^\n]\n"); // Skip the line
        c = getc(fileName);            // Get the next character
    }
    ungetc(c, fileName); // Goes back 1 character because the while loop goes 1 character too far

    // Read the number of rows, columns and non-zero elements
    fscanf(fileName, "%d %d %d", &aMatrix->num_rows, &aMatrix->num_cols, &aMatrix->num_non_zeros);

    int *rows = (int *)calloc(aMatrix->num_non_zeros, sizeof(int));
    int *cols = (int *)calloc(aMatrix->num_non_zeros, sizeof(int));
    double *data = (double *)calloc(aMatrix->num_non_zeros, sizeof(double));
    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        fscanf(fileName, "%d %d %lf", &rows[i], &cols[i], &data[i]);
        c = getc(fileName);
    }

    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        printf("Amatrix(%d, %d) = %lf\n", rows[i], cols[i], data[i]);
    }

    printf("Numbers have been read\n");
    free(rows);
    free(cols);
    free(data);
}



int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        printf("Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    // <Handle the inputs here>
    const char *filename = argv[1];

    // If the file is not in the directory
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("File not found\n");
        exit(1);
    }

    CSRMatrix aMatrix;
    ReadMMtoCSR(filename, &aMatrix);
    puts("Matrix read successfully");

    
}
