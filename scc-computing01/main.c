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

//print 2dArray
void printArray(int **pArray, int m, int n) {
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

int** arrayAlloc(int m, int n) {
    int **a = (int**)malloc(m * sizeof(int *));
    int i;
    for (i = 0 ; i < m; i++ ) {
        a[i] =(int *) malloc(n * sizeof(int));
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
    //# of nonZeroValues in the matrices truncated
    const int nonZeroCount = m * n * perc;

    
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
    }
    
    //initialize Matrices in heap, stack does not hold tons of MB
    //m,n not know during compile time -> manually initialize arrays
    int **a = arrayAlloc(m,n);
    int **b = arrayAlloc(m,n);
    a = (int**)malloc(m * sizeof(int *));
    
    int i = 0;
    for (i = 0 ; i < m; i++ ) {
        a[i] =(int *) malloc(n * sizeof(int));
    }
    
    //set nonZeroCount x random value into random place in matrices
    srand(time(NULL));
    for (i = 0; i < nonZeroCount; i++) {
        int randMA = rand() % m;
        int randNA = rand() % n;
        int randMB = rand() % m;
        int randNB = rand() % n;
    
        a[randMA][randNA] = rand() % 9;
        b[randMB][randNB] = rand() % 9;
    }
    
    //printArray(a, m, n);
    //printArray(b, m, n);
    printf("finished");
    return 0;
}




