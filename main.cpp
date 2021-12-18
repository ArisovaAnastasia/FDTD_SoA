#define _USE_MATH_DEFINES
#include <iostream>
#include <cstdlib>
#include "start_fdtd+pml.h"

int main() {

	//setlocale(LC_ALL, "english");

	omp_set_num_threads(16);

	double n = 4;

	double T = 1.0 / 64.0 * M_PI;
	int delta_x = 8, delta_y = 8, delta_z = 8;
	int Nx = 32, Ny = 32, Nz = 16;
	std::cout << "Time: " << T << std::endl;
	run_fdtd<double, double>(n, Nx, Ny, Nz, (double)T, (double)0.005,
		delta_x, delta_y, delta_z, 46.5, 46.5);

	return 0;
}