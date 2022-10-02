#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <openblas/lapacke.h>
//#include <mkl_lapacke.h>

double *generate_matrix(int size)
{
  int i;
  double *matrix = (double *) malloc(sizeof(double) * size * size);


  for (i = 0; i < size * size; i++) {
    matrix[i] = rand() % 100;
  }

  return matrix;
}

double *generate_vector(int size)
{
  int i;
  double *vector = (double *) malloc(sizeof(double) * size);


  for (i = 0; i < size * size; i++) {
    vector[i] = 0;
  }

  return vector;

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

  for(i = 0; i < size*size; i++) {
    if (!is_nearly_equal(bref[i], b[i]))
      return 0;
  }

  return 1;
}

double det(double *a, int lda)
{
    double det =1.0;
    int i, j,k;
    double tmp;

    for (j=0; j<lda; j++0)
    {
      i_max = j;
      for (i=0; i<lda; i++)
        if(a[i*lda+j] > a[i_max*lda+j])
          i_max = i;
      
      if (imax!=j)
      {
        for(k=0; k<lda; k++)
        {
            tmp = a[i_max*lda+k];
            a[i_max*lda+k] = a[j*lda+k];
            a[j*lda+k] = tmp;
        }
        det *= -1;

      }
        if (abs(a[j*(lda+1)]) < 1e-12)
        {
          printf("singular matrix");
          return NAN;
        }
        for (i= j + 1; i < lda; i++)
        {
          double mult = -a[i*lda+j]/a[j*(lda+1)];
          for(k = 0; k<lda; k++)
            a[i*lda+k] += mult * a[j*lda + k];
        }
    } 

    for (i=0; i < lda; i++)
      det *= a[i*(lda+1)];
    
    return det;
}

double cramer()
{

} 

double my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{


}

void main(int argc, char *argv[])
{

  ssrand(1);

  int size = atoi(argv[1]);

  double *a, *aref;
  double *b, *bref;

  a = generate_matrix(size);
  aref = generate_matrix(size);
  b = generate_matrix(size);
  bref = generate_matrix(size);

  // Using LAPACK dgesv OpenBLAS implementation to solve the system
  int n = size, nrhs = size, lda = size, ldb = size, info;
  int *ipiv = (int *) malloc(sizeof(int) * size);

  clock_t tStart = clock();
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
  printf("Time taken by OpenBLAS LAPACK: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

  int *ipiv2 = (int *) malloc(sizeof(int) * size);

  tStart = clock();
  my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb);
  printf("Time taken by my implementation: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

  if (check_result(bref, b, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");
}
