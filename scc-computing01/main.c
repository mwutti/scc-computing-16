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

void printMatrix(int **pArray, int m, int n) {
    int i;
    int j;
    for ( i = 0 ; i < m ; i++ ) {
        for ( j = 0; j < n; j++ ) {
            printf("%d ", pArray[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int** initSparseMatrix(int m, int n, double perc) {
    const int nonZeroCount = m * n * perc;
    int **a = (int**)malloc(m * sizeof(int *));
    
    int i;
    for (i = 0 ; i < m; i++ ) {
        a[i] =(int *) calloc(n * sizeof(int), sizeof(int));
    }
    
    //set nonZeroCount x random value into random place in matrices
    srand(time(NULL));
    printf("nonZeroCount %d\n", nonZeroCount);
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
    int **a = initSparseMatrix(m, n, perc);
    int **b = initSparseMatrix(m, n, perc);
    printMatrix(a, m, n);
    YSMF ya;
    ya.m = m;
    ya.n = n;
    //storage has to be reallocated in case of unknown matrices
    ya.nonZeros = m * n * perc;
    //create YSMF
    ya.a = malloc(ya.nonZeros * sizeof(int));
    ya.ai = malloc((ya.m + 1) * sizeof(int));
    ya.aj = malloc(ya.n * sizeof(int));
    
    int i,j;
    int indexA = 0;
    int indexAI = 0;
    int indexAJ = 0;
    
    for ( i = 0; i < m; i++ ) {
        for (j = 0; j < n; j++) {
            if ( a[i][j] != 0 ) {
                ya.a[indexA++] = a[i][j];
                ya.aj[indexAJ++] = j;
            }
        }
    }
    
    //print oneD- Array
    for (i = 0; i < ya.nonZeros; i++) {
        printf("%d ", ya.a[i]);
    }
    return 0;
}






