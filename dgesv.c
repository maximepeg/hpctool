#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <openblas/lapacke.h>
//#include <mkl_lapacke.h>

double *generate_matrix(int size)
{
  srand(1);
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
    int i, j, k, i_max;
    double tmp;

    for (j=0; j<lda; j++)
    {
      i_max = j;
      for (i=0; i<lda; i++)
        if(a[i*lda+j] > a[i_max*lda+j])
          i_max = i;
      
      if (i_max!=j)
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

double cramer(double *a, int lda, double det_a, double *b, int var)
{
  double *tmp;
  int i;
  tmp = generate_matrix(lda);
  memcpy(tmp, a, sizeof(*tmp));

  for (i = 0; i < lda; i++)
    tmp[i*lda+var] = b[var*lda+i];

  double det_tmp = det(tmp, lda);  

   
  return det_tmp / det_a;
} 

double my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
  int i, j;
  double *tmp;
  double det_a;
  det_a = det(a, lda);
  for (i = 0; i < lda; i++)
    for( j = 0; j < lda; j++)
      b[i*lda + j] = cramer(a, lda, det_a, b, i);

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
