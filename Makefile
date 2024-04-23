parallel_gaussian_quadrature:
	mpicc -O3 -o parallel_gaussian_quadrature parallel_gaussian_quadrature.c  -lm







clean:
	rm parallel_gaussian_quadrature
