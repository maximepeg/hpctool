dgesv_gcc: dgesv.c
	gcc dgesv.c -lm -lopenblas -o dgesv
dgesv_icc: dgesv.c
	icc dgesv.c -lm -mkl -o dgesv
 
