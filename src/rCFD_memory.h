#ifndef RCFD_MEMORY
#define RCFD_MEMORY

/* (C) 2022
    Daniel Queteschiner
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/

// 2d int array
int **malloc_i_2d(int n1, int n2)
{
    int i, n = 0;
    int **array;
    int nbytes = (int)(sizeof(int) * n1 * n2);
    int *data = (int*)malloc(nbytes);
    nbytes = (int)(sizeof(int*) * n1);
    array = (int**)malloc(nbytes);

    for(i = 0; i < n1; ++i){
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

void free_i_2d(int **array)
{
    if (array == NULL) return;
    free(array[0]);
    free(array);
}

// 2d double array
double **malloc_r_2d(int n1, int n2)
{
    int i, n = 0;
    double **array;
    int nbytes = (int)(sizeof(double) * n1 * n2);
    double *data = (double*)malloc(nbytes);
    nbytes = (int)(sizeof(double*) * n1);
    array = (double**)malloc(nbytes);

    for(i = 0; i < n1; ++i){
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

void free_r_2d(double **array)
{
    if (array == NULL) return;
    free(array[0]);
    free(array);
}

// 3d int array
int ***malloc_i_3d(int n1, int n2, int n3)
{
    int i, j, m, n = 0;
    int **plane;
    int ***array;
    int nbytes = (int)(sizeof(int) * n1 * n2 * n3);
    int *data = (int*)malloc(nbytes);
    nbytes = (int)(sizeof(int*) * n1 * n2);
    plane = (int**)malloc(nbytes);
    nbytes = (int)(sizeof(int**) * n1);
    array = (int***)malloc(nbytes);

    for(i = 0; i < n1; ++i){
        m = i * n2;
        array[i] = &plane[m];
        for(j = 0; j < n2; ++j){
            plane[m+j] = &data[n];
            n += n3;
        }
    }
    return array;
}

void free_i_3d(int ***array)
{
    if (array == NULL) return;
    free(array[0][0]);
    free(array[0]);
    free(array);
}

// 3d double array
double ***malloc_r_3d(int n1, int n2, int n3)
{
    int i, j, m, n = 0;
    double **plane;
    double ***array;
    int nbytes = (int)(sizeof(double) * n1 * n2 * n3);
    double *data = (double*)malloc(nbytes);
    nbytes = (int)(sizeof(double*) * n1 * n2);
    plane = (double**)malloc(nbytes);
    nbytes = (int)(sizeof(double**) * n1);
    array = (double***)malloc(nbytes);

    for(i = 0; i < n1; ++i){
        m = i * n2;
        array[i] = &plane[m];
        for(j = 0; j < n2; ++j){
            plane[m+j] = &data[n];
            n += n3;
        }
    }
    return array;
}

void free_r_3d(double ***array)
{
    if (array == NULL) return;
    free(array[0][0]);
    free(array[0]);
    free(array);
}

// 4d int array
int ****malloc_i_4d(int n1, int n2, int n3, int n4)
{
    int i, j, k;
    int m1, m2, n = 0;
    int **cube;
    int ***plane;
    int ****array;
    int nbytes = (int)(sizeof(int) * n1 * n2 * n3 * n4);
    int *data = (int*)malloc(nbytes);
    nbytes = (int)(sizeof(int*)) * n1 * n2 * n3;
    cube   = (int**)malloc(nbytes);
    nbytes = (int)(sizeof(int**)) * n1 * n2;
    plane  = (int***)malloc(nbytes);
    nbytes = (int)(sizeof(int***)) * n1;
    array  = (int****)malloc(nbytes);

    for(i = 0; i < n1; ++i){
        m1 = i * n2;
        array[i] = &plane[m1];
        for(j = 0; j < n2; ++j){
            m2 = i * n2 * n3 + j * n3;
            plane[m1+j] = &cube[m2];
            for(k = 0; k < n3; ++k){
                cube[m2+k] = &data[n];
                n += n4;
            }
        }
    }
    return array;
}

void free_i_4d(int ****array)
{
    if (array == NULL) return;
    free(array[0][0][0]);
    free(array[0][0]);
    free(array[0]);
    free(array);
}

#if 0
// test with valgrind
void test_memory_allocation()
{
    int n1=2, n2=3, n3=4, n4=5;
    int i, j, k, l;
    int **iarray2d = NULL;
    int ***iarray3d = NULL;
    int ****iarray4d = NULL;
    double **rarray2d = NULL;
    double ***rarray3d = NULL;

    // test free functions with NULL pointers
    free_i_2d(iarray2d);
    free_i_3d(iarray3d);
    free_i_4d(iarray4d);
    free_r_2d(rarray2d);
    free_r_3d(rarray3d);

    // test malloc functions
    iarray2d = malloc_i_2d(n1, n2);
    iarray3d = malloc_i_3d(n1, n2, n3);
    iarray4d = malloc_i_4d(n1, n2, n3, n4);
    rarray2d = malloc_r_2d(n1, n2);
    rarray3d = malloc_r_3d(n1, n2, n3);

    // fill arrays
    for(i = 0; i < n1; ++i){
        for(j = 0; j < n2; ++j){
            iarray2d[i][j] = (i+1)*(j+1);
            rarray2d[i][j] = (i+1)*(j+1)*2.0;
            for(k = 0; k < n3; ++k){
                iarray3d[i][j][k] = (i+1)*(j+1)*(k+1);
                rarray3d[i][j][k] = (i+1)*(j+1)*(k+1)*2.0;
                for(l = 0; l< n4; ++l){
                    iarray4d[i][j][k][l] = (i+1)*(j+1)*(k+1)*(l+1);
                }
            }
        }
    }

    // print to screen
    for(i = 0; i < n1; ++i){
        for(j = 0; j < n2; ++j){
            printf("iarray2d[%d][%d] = %6d ", i, j, iarray2d[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for(i = 0; i < n1; ++i){
        for(j = 0; j < n2; ++j){
            printf("rarray2d[%d][%d] = %8.2f ", i, j, rarray2d[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for(i = 0; i < n1; ++i){
        for(j = 0; j < n2; ++j){
            for(k = 0; k < n3; ++k){
                printf("iarray3d[%d][%d][%d] = %6d ", i, j, k, iarray3d[i][j][k]);
            }
            printf("\n");
        }
    }
    printf("\n");
    for(i = 0; i < n1; ++i){
        for(j = 0; j < n2; ++j){
            for(k = 0; k < n3; ++k){
                printf("rarray3d[%d][%d][%d] = %8.2f ", i, j, k, rarray3d[i][j][k]);
            }
            printf("\n");
        }
    }
    printf("\n");
    for(i = 0; i < n1; ++i){
        for(j = 0; j < n2; ++j){
           for(k = 0; k < n3; ++k){
                for(l = 0; l< n4; ++l){
                    printf("iarray4d[%d][%d][%d][%d] = %6d ", i, j, k, l, iarray4d[i][j][k][l]);
                }
                printf("\n");
            }
        }
    }

    // test free functions with non-NULL pointers
    free_i_2d(iarray2d);
    free_i_3d(iarray3d);
    free_i_4d(iarray4d);
    free_r_2d(rarray2d);
    free_r_3d(rarray3d);
}
#endif

#endif
