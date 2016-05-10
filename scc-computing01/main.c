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
#include <sys/timeb.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <pthread.h>
#include <stdbool.h>

typedef struct YSMF {
    int *a;
    int *ai;
    int *aj;
} YSMF;

// Struct or parallel initialization of Sparse Matrices
typedef struct INIT_STRUCT {
    int row;
    int **matrix;
} INIT_STRUCT;

typedef struct YSMF_ADD {
    int row;
    short *added;
    int **tempA;
    int *tempAI;
    int **tempAJ;
} YSMF_ADD;

int m;
int n;
double perc;
long ms;
YSMF *yaleMatrixA;
YSMF *yaleMatrixB;
YSMF *yaleMatrixC;

int numThreads;
pthread_t * threads;
pthread_mutex_t mu1;
pthread_mutex_t mu2;
pthread_mutex_t mu3;

/**
 * Creates a YSMF from given matrix a
 * matrix sparse matrix
 * m Number of Rows of matrix
 * n Number of Columns of matrix
 * nonZeros number of non zero elements in matrix (easier allocation of ai and aj)
 */
YSMF *initYaleMatrix(int** matrix, int m, int n, int nonZeros) {
    YSMF *yaleMatrix = (YSMF *)malloc(sizeof(YSMF));
    if ( yaleMatrix == NULL ) {
        return NULL;
    }
    
    yaleMatrix->aj = (int *)malloc(nonZeros * sizeof(int));
    yaleMatrix->ai = (int *)malloc((m+1) * sizeof(int));
    yaleMatrix->a = (int *)malloc(nonZeros * sizeof(int));
    
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
    for (i = 0; i < matrix->ai[m]; i++) {
        printf("%d ", matrix->a[i]);
    }
    printf("\n");
    
    printf("Ai: ");
    for (i = 0; i < m + 1; i++) {
        printf("%d ", matrix->ai[i]);
    }
    printf("\n");
    
    printf("Aj: ");
    for (i = 0; i < matrix->ai[m]; i++) {
        printf("%d ", matrix->aj[i]);
    }
    printf("\n\n");
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

void *initRow(void *init_struct) {
    INIT_STRUCT *init = (INIT_STRUCT *)init_struct;
    while ( init->row < m ) {
        pthread_mutex_lock(&mu1);
        int row = init->row++;
        if ( row < m ) {
            init->matrix[row] =(int *) calloc(n * sizeof(int), sizeof(int));
        }
        pthread_mutex_unlock(&mu1);
    }

    return (NULL);
}


int ** initSparseMatrixParallel(int m, int n, double perc) {
    const int nonZeroCount = m * n * perc;
    int **a = (int**)malloc(m * sizeof(int *));
    int i;
    
    INIT_STRUCT *init = (INIT_STRUCT *)malloc(sizeof(INIT_STRUCT *));
    init->matrix = a;
    init->row = 0;
    
    pthread_mutex_init(&mu1, NULL);
    for ( i = 0; i < numThreads ; i++ ) {
        pthread_create(&threads[i], NULL, initRow, (void *)init);
    }
    
    for ( i = 0; i < numThreads ; i++ ) {
        pthread_join(threads[i], NULL);
    }
    
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
    int **matrix = (int**)malloc(m * sizeof(int*));
    
    for ( i = 0; i < m; i++ ) {
        matrix[i] = (int*) calloc(n * sizeof(int), sizeof(int));
    }
    
    //iterate over Ai[]
    for (i = 1; i < m + 1; i++) {
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

YSMF* addYSMF(YSMF *yaleMatrixA, YSMF *yaleMatrixB){
    yaleMatrixC = (YSMF *)malloc(sizeof(YSMF));

    if ( yaleMatrixC == NULL ) {
        return NULL;
    }
    
    // optimistic allocation of a and ai with nonZerosA +nonZerosB ... reallocated when finished
    yaleMatrixC->aj = (int *)malloc((yaleMatrixA->ai[m] + yaleMatrixB->ai[m]) * sizeof(int));
    yaleMatrixC->ai = (int *)malloc((m + 1) * sizeof(int));
    yaleMatrixC->a = (int *)malloc((yaleMatrixA->ai[m] + yaleMatrixB->ai[m]) * sizeof(int));
    
    if ( yaleMatrixC->a == NULL ) {
        free(yaleMatrixC);
        return NULL;
    }
    
    if ( yaleMatrixC->aj == NULL ) {
        free(yaleMatrixC->a);
        free(yaleMatrixC);
        return NULL;
    }
    
    if ( yaleMatrixC->ai == NULL ) {
        free(yaleMatrixC->a);
        free(yaleMatrixC->aj);
        free(yaleMatrixC);
        return NULL;
    }
    
    int aiA= 1, aiB = 1, aiC = 1;      //index for rowArrays
    int ajA= 0, ajB = 0, ajC = 0;   //index for columnArrays
    int i;
    
    yaleMatrixC->ai[0] = 0;
    for ( i = 1; i < m + 1; i++ ) {
        if ( yaleMatrixA->ai[aiA] - yaleMatrixA->ai[aiA - 1] == 0 &&    // case both rows have no 0 Values
            yaleMatrixB->ai[aiB] - yaleMatrixB->ai[aiB - 1] == 0 ) {
            yaleMatrixC->ai[aiC] = yaleMatrixC->ai[aiC - 1];
            aiA++;
            aiB++;
            aiC++;
        } else if ( yaleMatrixA->ai[aiA] - yaleMatrixA->ai[aiA - 1] == 0) { // in Matrix B there is at least 1 nonZero value
            int nonZeros = yaleMatrixB->ai[aiB] - yaleMatrixB->ai[aiB - 1]; // how many nonZeros in this row?
            yaleMatrixC->ai[aiC] = yaleMatrixC->ai[aiC - 1] + nonZeros;     // number of non Zeros for this row
            aiA++;
            aiB++;
            aiC++;
            while ( nonZeros > 0 ) {                                      // for all nonZeros in this row
                yaleMatrixC->aj[ajC] = yaleMatrixB->aj[ajB];
                yaleMatrixC->a[ajC++] = yaleMatrixB->a[ajB++];
                nonZeros--;
            }
        } else if ( yaleMatrixB->ai[aiB] - yaleMatrixB->ai[aiB - 1] == 0) { // in Matrix A there is at least 1 nonZero value
            int nonZeros = yaleMatrixA->ai[aiA] - yaleMatrixA->ai[aiA - 1];
            yaleMatrixC->ai[aiC] = yaleMatrixC->ai[aiC - 1] + nonZeros;     // number of non Zeros for this row
            aiA++;
            aiB++;
            aiC++;
            while ( nonZeros > 0 ) {                                      // for all nonZeros in this row
                yaleMatrixC->aj[ajC] = yaleMatrixA->aj[ajA];
                yaleMatrixC->a[ajC++] = yaleMatrixA->a[ajA++];
                nonZeros--;
            }
        } else {        // Potential values to add
            int nonZerosA = yaleMatrixA->ai[aiA] - yaleMatrixA->ai[aiA - 1];
            int nonZerosB = yaleMatrixB->ai[aiB] - yaleMatrixB->ai[aiB - 1];
            int nonZerosForRow = 0;
            while ( nonZerosA > 0 && nonZerosB > 0 ) {    // Start iterating over nonZeros
                if ( yaleMatrixA->aj[ajA] == yaleMatrixB->aj[ajB] ) { // addition of values
                    yaleMatrixC->a[ajC] = yaleMatrixA->a[ajA++] + yaleMatrixB->a[ajB];
                    yaleMatrixC->aj[ajC++] = yaleMatrixB->aj[ajB++];
                    nonZerosForRow++;
                    nonZerosA--;
                    nonZerosB--;
                } else if ( yaleMatrixA->aj[ajA] < yaleMatrixB->aj[ajB] ) { // first add smaller value in Matrix A
                    yaleMatrixC->aj[ajC] = yaleMatrixA->aj[ajA];
                    yaleMatrixC->a[ajC++] = yaleMatrixA->a[ajA++];
                    nonZerosForRow++;
                    nonZerosA--;
                } else  { // first add smaller value in Matrix B
                    yaleMatrixC->aj[ajC] = yaleMatrixB->aj[ajB];
                    yaleMatrixC->a[ajC++] = yaleMatrixB->a[ajB++];
                    nonZerosForRow++;
                    nonZerosB--;
                }
            }
            
            while ( nonZerosA > 0 ) { //add remaining values from A
                yaleMatrixC->aj[ajC] = yaleMatrixA->aj[ajA];
                yaleMatrixC->a[ajC++] = yaleMatrixA->a[ajA++];
                nonZerosForRow++;
                nonZerosA--;
            }
            
            while ( nonZerosB > 0 ) { //add remaining values from B
                yaleMatrixC->aj[ajC] = yaleMatrixB->aj[ajB];
                yaleMatrixC->a[ajC++] = yaleMatrixB->a[ajB++];
                nonZerosForRow++;
                nonZerosB--;
            }
            yaleMatrixC->ai[aiC] = yaleMatrixC->ai[aiC - 1] + nonZerosForRow;     // number of non Zeros for this row
            aiA++;
            aiB++;
            aiC++;
        }
    }
    realloc(yaleMatrixC->a, yaleMatrixC->ai[m] * sizeof(int));
    realloc(yaleMatrixC->aj, yaleMatrixC->ai[m] * sizeof(int));
    
    return yaleMatrixC;
}

void * addRowParallel(void *ysmf_add) {
    YSMF_ADD *ysmfADD = (YSMF_ADD *)ysmf_add;
    int myRow;
    int i;
    
    while ( ysmfADD->row < m ) {
        //which row to add?
        pthread_mutex_lock(&mu1);
        myRow = ysmfADD->row++;
        pthread_mutex_unlock(&mu1);
        
        if ( myRow < m ) {
            if ( yaleMatrixA->ai[myRow + 1] - yaleMatrixA->ai[myRow] == 0 &&    // case both rows have no 0 Values
                yaleMatrixB->ai[myRow + 1] - yaleMatrixB->ai[myRow] == 0 ) {
                ysmfADD->tempAI[myRow] = 0;
                ysmfADD->added[myRow] = 1;
            } else if ( yaleMatrixA->ai[myRow + 1] - yaleMatrixA->ai[myRow] == 0) { // in Matrix B there is at least 1 nonZero value
                int nonZeros = yaleMatrixB->ai[myRow + 1] - yaleMatrixB->ai[myRow]; // how many nonZeros in this row?
                ysmfADD->tempAJ[myRow] = (int *) malloc(nonZeros * sizeof(int));
                ysmfADD->tempA[myRow] = (int *) malloc(nonZeros * sizeof(int));
                
                ysmfADD->tempAI[myRow] = nonZeros;
                
                for (i = 0; nonZeros > 0; i++) {                                      // for all nonZeros in this row
                    ysmfADD->tempAJ[myRow][i] = yaleMatrixB->aj[yaleMatrixB->ai[myRow] +i];
                    ysmfADD->tempA[myRow][i] = yaleMatrixB->a[yaleMatrixB->ai[myRow] + i];
                    nonZeros--;
                }
                ysmfADD->added[myRow] = 1;
            } else if ( yaleMatrixB->ai[myRow + 1] - yaleMatrixB->ai[myRow] == 0) { // in Matrix A there is at least 1 nonZero value
                int nonZeros = yaleMatrixA->ai[myRow + 1] - yaleMatrixA->ai[myRow];
                ysmfADD->tempAJ[myRow] = (int *)malloc(nonZeros * sizeof(int));
                ysmfADD->tempA[myRow] = (int *)malloc(nonZeros * sizeof(int));
                
                ysmfADD->tempAI[myRow] = nonZeros;
                
                for (i = 0; nonZeros > 0; i++) {                                        // for all nonZeros in this row
                    ysmfADD->tempAJ[myRow][i] = yaleMatrixA->aj[yaleMatrixA->ai[myRow] + i];
                    ysmfADD->tempA[myRow][i] = yaleMatrixA->a[yaleMatrixA->ai[myRow] + i];
                    nonZeros--;
                }
                
                ysmfADD->added[myRow] = 1;
            } else {        // Potential values to add
                int nonZerosA = yaleMatrixA->ai[myRow + 1] - yaleMatrixA->ai[myRow];
                int nonZerosB = yaleMatrixB->ai[myRow + 1] - yaleMatrixB->ai[myRow];
                int jA = 0, jB = 0;
                ysmfADD->tempAJ[myRow] = (int *)malloc((nonZerosA + nonZerosB) * sizeof(int));
                ysmfADD->tempA[myRow] = (int *)malloc((nonZerosA + nonZerosB) * sizeof(int));
                
                
                for (i = 0; (nonZerosA > 0 && nonZerosB > 0); i++ ) {    // Start iterating over nonZeros
                    if ( yaleMatrixA->aj[yaleMatrixA->ai[myRow] + jA] == yaleMatrixB->aj[yaleMatrixB->ai[myRow] + jB] ) { // addition of values
                        ysmfADD->tempA[myRow][i] = yaleMatrixA->a[yaleMatrixA->ai[myRow] + jA] + yaleMatrixB->a[yaleMatrixB->ai[myRow] + jB];
                        ysmfADD->tempAJ[myRow][i] = yaleMatrixA->aj[yaleMatrixA->ai[myRow] + jA];
                        nonZerosA--;
                        nonZerosB--;
                        jB++;
                        jA++;
                    } else if ( yaleMatrixA->aj[yaleMatrixA->ai[myRow] + jA] < yaleMatrixB->aj[yaleMatrixB->ai[myRow] + jB] ) { // first add smaller value in Matrix A
                        ysmfADD->tempAJ[myRow][i] = yaleMatrixA->aj[yaleMatrixA->ai[myRow] + jA];
                        ysmfADD->tempA[myRow][i] = yaleMatrixA->a[yaleMatrixA->ai[myRow] + jA];
                        nonZerosA--;
                        jA++;
                    } else  { // first add smaller value in Matrix B
                        ysmfADD->tempAJ[myRow][i] = yaleMatrixB->aj[yaleMatrixB->ai[myRow] + jB];
                        ysmfADD->tempA[myRow][i] = yaleMatrixB->a[yaleMatrixB->ai[myRow] + jB];
                        nonZerosB--;
                        jB++;
                    }
                }
                
                for (; nonZerosA > 0; i++ ) { //add remaining values from A
                    ysmfADD->tempAJ[myRow][i] = yaleMatrixA->aj[yaleMatrixA->ai[myRow] + jA];
                    ysmfADD->tempA[myRow][i] = yaleMatrixA->a[yaleMatrixA->ai[myRow] + jA];
                    nonZerosA--;
                    jB++;jA++;
                }
                
                for (; nonZerosB > 0; i++ ) { //add remaining values from B
                    ysmfADD->tempAJ[myRow][i] = yaleMatrixB->aj[yaleMatrixB->ai[myRow] + jB];
                    ysmfADD->tempA[myRow][i] = yaleMatrixB->a[yaleMatrixB->ai[myRow] + jB];
                    nonZerosB--;
                    jB++;
                }
                ysmfADD->tempAI[myRow] =  i;     // number of non Zeros for this row
                ysmfADD->added[myRow] = 1;
            }
        } else {
            return (NULL);
        }
    }
    return (NULL);
}

YSMF* addYSMFParallel(YSMF *yaleMatrixA, YSMF *yaleMatrixB) {
    int i;
    yaleMatrixC = (YSMF *)malloc(sizeof(YSMF));

    if ( yaleMatrixC == NULL ) {
        return NULL;
    }
    
    // optimistic allocation of a and ai with nonZerosA +nonZerosB ... reallocated when finished
    yaleMatrixC->aj = (int *)malloc((yaleMatrixA->ai[m] + yaleMatrixB->ai[m]) * sizeof(int));
    yaleMatrixC->ai = (int *)malloc((m + 1) * sizeof(int));
    yaleMatrixC->a = (int *)malloc((yaleMatrixA->ai[m] + yaleMatrixB->ai[m]) * sizeof(int));
    
    if ( yaleMatrixC->a == NULL ) {
        free(yaleMatrixC);
        return NULL;
    }
    
    if ( yaleMatrixC->aj == NULL ) {
        free(yaleMatrixC->a);
        free(yaleMatrixC);
        return NULL;
    }
    
    if ( yaleMatrixC->ai == NULL ) {
        free(yaleMatrixC->a);
        free(yaleMatrixC->aj);
        free(yaleMatrixC);
        return NULL;
    }
    yaleMatrixC->ai[0] = 0;
    //Start parallel
    
    pthread_mutex_init(&mu1, NULL);
    pthread_mutex_init(&mu2, NULL);

    YSMF_ADD *ysmfADD = (YSMF_ADD *)malloc(sizeof(YSMF_ADD));
    ysmfADD->row = 0;
    ysmfADD->tempA = (int **)malloc(m * sizeof(int*));
    ysmfADD->tempAJ = (int **)malloc(m * sizeof(int*));
    ysmfADD->tempAI = (int *)malloc(m * sizeof(int));
    ysmfADD->added = (short *)calloc(sizeof(short), m * sizeof(short));
    
    for ( i = 0; i < numThreads; i++ ) {
        pthread_create(&threads[i], NULL, addRowParallel, (void *)ysmfADD);
    }

    int j, totalAdded = 0;
    
    
    for ( i = 0; i < m; i++ ) {
        //Busy Waiting
        while ( !ysmfADD->added[i] );

        for ( j = 0; j < ysmfADD->tempAI[i]; j++ ) {
            yaleMatrixC->a[totalAdded] = ysmfADD->tempA[i][j];
            yaleMatrixC->aj[totalAdded] = ysmfADD->tempAJ[i][j];
            totalAdded++;
        }
        yaleMatrixC->ai[i+1] = totalAdded;
    }
    
    //realloc(yaleMatrixC->a, yaleMatrixC->ai[m] * sizeof(int));
    //realloc(yaleMatrixC->aj, yaleMatrixC->ai[m] * sizeof(int));
    
    return yaleMatrixC;
}

int main( int argc, const char* argv[] ) {
    //arg[0], arg[1] -> m x n
    //arg[2] perc. of non 0 values
    //arg[3] #threads
    if ( argc != 5 ) {
        printf("usage: <m> <n> <percantage of non 0 values> <#threads>\n");
        return -1;
    }
    
    m = atoi(argv[1]);
    n = atoi(argv[2]);
    perc = atof(argv[3]);
    numThreads = atoi(argv[4]);

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
    
    threads = malloc(sizeof(pthread_t) * numThreads);
    //Init Sparse Matrices
    srand(time(NULL));
    int **a = initSparseMatrix(m, n, perc);
    //Create YSMF
    printf("%lu\n", sizeof(int));
    printf("%f", (double)(8 * m * n / 1024 ));
    yaleMatrixA = initYaleMatrix(a, m, n, m * n * perc);
    free(a);
    int **b = initSparseMatrix(m, n, perc);
    yaleMatrixB = initYaleMatrix(b, m, n, m * n * perc);
    free(b);
    
    struct timeb start, end;
    ftime(&start);
    
    if ( numThreads == 1 ) {
        yaleMatrixC = addYSMF(yaleMatrixA, yaleMatrixB);
    } else {
        numThreads++;
        yaleMatrixC = addYSMFParallel(yaleMatrixA, yaleMatrixB);
    }
    ftime(&end);
    
    printf("Matrix Addition took %d ms.\n",(int) (1000.0 * (end.time - start.time) + (end.millitm - start.millitm)));
    printf("Number of Threads: %d\n", numThreads);
    printf("m x n: %d x %d\n", m, n);
    printf("Non Zero Values: %d%% \n", (int) (perc * 100));
    free(yaleMatrixA);
    free(yaleMatrixB);
    free(threads);
    return 0;
}

