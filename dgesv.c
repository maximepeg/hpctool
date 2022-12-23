#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <openblas/lapacke.h>
//#include <mkl_lapacke.h>

double *generate_matrix(int size)
{
    srand(1);
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}
double *initialize_matrix(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = 0;
    }

    return matrix;
}




int is_nearly_equal(double x, double y)
{
    const double epsilon = 1e-5 /* some small number */;
    return abs(x - y) <= epsilon * abs(x);
    // see Knuth section 4.2.2 pages 217-218
}

int check_result(double *bref, double *b, int size)
{
    int i;

    for (i = 0; i < size * size; i++)
    {
        if (!is_nearly_equal(bref[i], b[i]))
            return 0;
    }

    return 1;
}
void matmul(double *a, double *b, int n)
{
    int i,j,k;
    double sum;
    double *temp;
    temp = initialize_matrix(n);

    #pragma omp parallel for collapse(3) private(i, j, k) shared(temp, a, b)
    for(i=0;i<n;i++)
        for(k=0;k<n;k++)
            for(j=0;j<n;j++)
                temp[i*n+j] += a[i*n+k]*b[k*n+j];

    #pragma omp parallel for private(i) shared(b, temp)
    for(i=0;i<n*n;i++)
        b[i]=temp[i];

    free(temp);
}
int gjinv (double *a, int n, double *b) //source is rosettacode
{
	int i, j, k, p;
	double f, g, tol;
	if (n < 1) return -1;  /* Function Body */
	f = 0.;  /* Frobenius norm of a */
    #pragma omp parallel for collapse(2) private(i,j,g) shared(a) reduction(+: f) 
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			g = a[j+i*n];
			f += g * g;
		}
	}
	f = sqrt(f);
	tol = f * 2.2204460492503131e-016;
    #pragma omp parallel for collapse(2) private(i,j) shared(b)
	for (i = 0; i < n; ++i) {  /* Set b to identity matrix. */
		for (j = 0; j < n; ++j) {
			b[j+i*n] = (i == j) ? 1. : 0.;
		}
	}

	for (k = 0; k < n; ++k) {  /* Main loop */
		f = fabs(a[k+k*n]);  /* Find pivot. */
		p = k;

        #pragma omp parallel for schedule(static) private(i, p, g, f) shared(a, k) 
		for (i = k+1; i < n; ++i) {
			g = fabs(a[k+i*n]);
			if (g > f) {
				f = g;
				p = i;
			}
		}
		if (f < tol) return 1;  /* Matrix is singular. */
		if (p != k) {  /* Swap rows. */

            #pragma omp parallel for schedule(static) private(j, f) shared(a) 
			for (j = k; j < n; ++j) {
				f = a[j+k*n];
				a[j+k*n] = a[j+p*n];
				a[j+p*n] = f;
			}

            #pragma omp parallel for schedule(static) private(j, f) shared(b) 
			for (j = 0; j < n; ++j) {
				f = b[j+k*n];
				b[j+k*n] = b[j+p*n];
				b[j+p*n] = f;
			}
		}
		f = 1. / a[k+k*n];  /* Scale row so pivot is 1. */
        #pragma omp parallel for schedule(static) private(j) shared(a) 
		for (j = k; j < n; ++j) a[j+k*n] *= f;

        #pragma omp parallel for schedule(static) private(j) shared(b) 
		for (j = 0; j < n; ++j) b[j+k*n] *= f;

		for (i = 0; i < n; ++i) {  /* Subtract to get zeros. */
			if (i == k) continue;
			f = a[k+i*n];
            #pragma omp parallel for private(j) shared(a) 
			for (j = k; j < n; ++j) a[j+i*n] -= a[j+k*n] * f;

            #pragma omp parallel for private(j) shared(b) 
			for (j = 0; j < n; ++j) b[j+i*n] -= b[j+k*n] * f;
		}
	}
	return 0;
} /* end of gjinv */

void display(double *a, int n)
{

    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            printf("%f ", a[i + n * j]);
        printf("\n");
    }
}

double my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
  double *temp;
  double start;
  start = omp_get_wtime();
  temp = initialize_matrix(lda);
  gjinv(a, lda, temp);
  //display(temp,lda);
  //printf("\n");
  matmul(temp, b, lda);
  //display(b, lda);
  free(temp);
  printf("Wall Time elapsed %2f\n", omp_get_wtime()-start);

}


void main(int argc, char *argv[])
{

    int size = atoi(argv[1]);

    double *a, *aref;
    double *b, *bref;

    a = generate_matrix(size);
    aref = generate_matrix(size);
    b = generate_matrix(size);
    bref = generate_matrix(size);

    // Using LAPACK dgesv OpenBLAS implementation to solve the system
    int n = size, nrhs = size, lda = size, ldb = size, info;
    int *ipiv = (int *)malloc(sizeof(int) * size);

    clock_t tStart = clock();

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
    printf("Time taken by OpenBLAS LAPACK: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    int *ipiv2 = (int *)malloc(sizeof(int) * size);

    tStart = clock();
    my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb);
    printf("Time taken by my implementation: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
    /*printf("default\n");
    display(bref, size);
    printf("my implementation\n");
    display(b, size);*/
    if (check_result(bref, b, size) == 1)
        printf("Result is ok!\n");
    else
        printf("Result is wrong!\n");

    free(a);
    free(b);
    free(aref);
    free(bref);
    free(ipiv);
    free(ipiv2);
}
