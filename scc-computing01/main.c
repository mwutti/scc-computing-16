//
//  Main.c
//  sccc01
//
//  Created by Michael Wutti on 01/05/16.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>

typedef struct YSMF {
    int m;
    int n;
    int nonZeros;
    int *a;
    int *ai;
    int *aj;
} YSMF;

/**
 * Creates a YSMF from given matrix a
 * matrix sparse matrix
 * m Number of Rows of matrix
 * n Number of Columns of matrix
 * nonZeros number of non zero elements in matrix (easier allocation of ai and aj)
 */
YSMF *initYaleMatrix(int** matrix, int m, int n, int nonZeros) {
    YSMF *yaleMatrix = (YSMF *)malloc (sizeof(YSMF ));
    if ( yaleMatrix == NULL ) {
        return NULL;
    }
    
    yaleMatrix->aj = (int *)malloc(nonZeros * sizeof(int ));
    yaleMatrix->ai = (int *)malloc((m+1) * sizeof(int ));
    yaleMatrix->a = (int *)malloc(nonZeros * sizeof(int ));
    
    if ( yaleMatrix->a == NULL ) {
        free(yaleMatrix);
        return NULL;
    }
    
    if ( yaleMatrix->aj == NULL ) {
        free(yaleMatrix->a);
        free(yaleMatrix);
        return NULL;
    }
    
    if ( yaleMatrix->ai == NULL ) {
        free(yaleMatrix->a);
        free(yaleMatrix->aj);
        free(yaleMatrix);
        return NULL;
    }
    
    int i,j;
    int indexA = 0;
    int indexAI = 1;
    int indexAJ = 0;
    
    yaleMatrix->ai[0] = 0;
    for ( i = 0; i < m; i++ ) {
        int nonZeroCount = 0;
        for (j = 0; j < n; j++) {
            if ( matrix[i][j] != 0 ) {
                yaleMatrix->a[indexA++] = matrix[i][j];
                yaleMatrix->aj[indexAJ++] = j;
                nonZeroCount++;
            }
        }
        yaleMatrix->ai[indexAI] = yaleMatrix->ai[indexAI-1] + nonZeroCount;
        indexAI++;
    }
    
    yaleMatrix->m = m;
    yaleMatrix->n = n;
    yaleMatrix->nonZeros = nonZeros;
    return yaleMatrix;
}

/**
 * Prints matrix
 * matrix matrix
 * m Number of rows of matrix
 * n Number of columns of matrix
 */
void printMatrix(int **matrix, int m, int n) {
    int i;
    int j;
    for ( i = 0 ; i < m ; i++ ) {
        for ( j = 0; j < n; j++ ) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/**
 * prints YSMF
 *
 */
void printYaleMatrix(YSMF *matrix) {
    int i;
    
    printf("A : ");
    for (i = 0; i < matrix->nonZeros; i++) {
        printf("%d ", matrix->a[i]);
    }
    printf("\n");
    
    printf("Ai: ");
    for (i = 0; i < matrix->m + 1; i++) {
        printf("%d ", matrix->ai[i]);
    }
    printf("\n");
    
    printf("Aj: ");
    for (i = 0; i < matrix->nonZeros; i++) {
        printf("%d ", matrix->aj[i]);
    }
}

/**
 * initializes a mxn sparse matrix with perc non zero values
 * m Rows of matrix
 * n Columns of matrix
 * perc percantage of non zero elements in Matrix
 */
int** initSparseMatrix(int m, int n, double perc) {
    const int nonZeroCount = m * n * perc;
    int **a = (int**)malloc(m * sizeof(int *));
    
    int i;
    for (i = 0 ; i < m; i++ ) {
        a[i] =(int *) calloc(n * sizeof(int), sizeof(int));
    }
    
    //set nonZeroCount x random value into random place in matrices
    for (i = 0; i < nonZeroCount; i++) {
        int randM, randN;
        do {
            randM = rand() % m;
            randN = rand() % n;
        } while (a[randM][randN] != 0);
        
        a[randM][randN] = rand() % 9 + 1;
    }
    return a;
}

/**
 * Converts a YSMF into a simple Matrix format
 */
int** convertFromYale(YSMF *yaleMatrix) {
    int i, j, indexAJ = 0, indexA = 0;
    int **matrix = (int**)malloc(yaleMatrix->m * sizeof(int*));
    
    for ( i = 0; i < yaleMatrix->m; i++ ) {
        matrix[i] = (int*) calloc(yaleMatrix->n * sizeof(int), sizeof(int));
    }
    
    //iterate over Ai[]
    for (i = 1; i < yaleMatrix->m + 1; i++) {
        //How Many nonZeros in row i ?
        int nonZeros = yaleMatrix->ai[i] - yaleMatrix->ai[i -1];
        
        for ( j = 0; j < nonZeros; j++ ) {
            matrix[i -1][yaleMatrix->aj[indexAJ++]] = yaleMatrix->a[indexA++];
        }
    }
    
    return matrix;
}

int** addSimple(int **a, int **b, int m, int n) {
    int **c = (int**)malloc(m * sizeof(int *));
    
    int i,j;
    for ( i = 0 ; i < m; i++ ) {
        c[i] =(int *) malloc(n * sizeof(int));
    }
    
    for ( i = 0; i < m; i++ ) {
        for ( j = 0; j < n; j++ ) {
            c[i][j] = a[i][j] + b[i][j];
        }
    }

    return c;
}

int main( int argc, const char* argv[] ) {
    //arg[0], arg[1] -> m x n
    //arg[2] perc. of non 0 values
    //arg[3] #threads
    if ( argc != 5 ) {
        printf("usage: <m> <n> <percantage of non 0 values> <#threads>\n");
        return -1;
    }
    
    const int m = atoi(argv[1]);
    const int n = atoi(argv[2]);
    const double perc = atof(argv[3]);
    const int numThreads = atoi(argv[4]);
    
    //Validation
    if ( perc < 0.0 || perc > 1.0 ) {
        printf("percentage invalid: must be a value between 0.0 and 0.99999\n");
        return -1;
    }
    
    if ( numThreads <= 0 || numThreads > 128 ) {
        printf("#Threads must be between 1 and 128");
        return -1;
    }
    
    if ( m <= 0 || n <= 0 ) {
        printf("m and n must be positive");
        return -1;
    }//Validation
    
    //Init Sparse Matrices
    srand(time(NULL));
    int **a = initSparseMatrix(m, n, perc);
    int **b = initSparseMatrix(m, n, perc);
    printMatrix(a, m, n);
    printMatrix(b, m, n);
    printf("\n");
    
    //Create YSMF
    YSMF *yaleMatrixA = initYaleMatrix(a, m, n, m * n * perc);
    YSMF *yaleMatrixB = initYaleMatrix(b, m, n, m * n * perc);
    //printYaleMatrix(yaleMatrixA);
    //int **matrix = convertFromYale(yaleMatrixA);
    //printf("\n");
    //printMatrix(matrix, m, n);
    
    //ADD matrices simple
    int **c = addSimple(a, b, m, n);
    
    printMatrix(c, m, n);
    
    return 0;
}

