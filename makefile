dgesv_gcc: dgesv.c
	gcc dgesv.c -lopenblas -o dgesv
dgesv_icc: dgesv.c
	icc dgesv.c -mkl -o dgesv
 
