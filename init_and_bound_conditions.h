#pragma once
#include <omp.h>
#include "type_data.h"

template <class ftypePML>
ftypePML distanceH(int N, int delta, int i)
{
	ftypePML dist = 0.0;

	if (delta == 0)
		return 0.0;

	if (i < delta + 1)
		dist = (ftypePML)(delta + 1 - i);

	if (i > N + delta)
		dist = (ftypePML)(i - delta - N) - (ftypePML)0.5;

	return dist / (ftypePML)delta;
}

template <class ftypePML>
ftypePML distanceE(int N, int delta, int i)
{
	ftypePML dist = 0.0;

	if (delta == 0)
		return 0.0;

	if (i < delta + 1)
		dist = (ftypePML)(delta + 1 - i) - (ftypePML)0.5;

	if (i > N + delta)
		dist = (ftypePML)(i - delta - N);

	return dist / (ftypePML)delta;
}

template <class ftypePML>
void Init_Sigma_SoA(std::vector<data3d<ftypePML> >& arr_sigma, double n,
	int Nx, int Ny, int Nz, int delta_x, int delta_y, int delta_z,
	ftypePML sigma_x, ftypePML sigma_y)
{
	ftypePML var_Sigma_max_x = sigma_x;
	ftypePML var_Sigma_max_y = sigma_y;
	ftypePML var_Sigma_max_z = sigma_y;
	
	std::cout << var_Sigma_max_x << "  " << var_Sigma_max_y << "   " << var_Sigma_max_z << "    " << std::endl;

	for (int i = 0; i < 6;++i) {
		arr_sigma[i].Create(Nx, Ny, Nz, delta_x, delta_y, delta_z);
	}

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
				arr_sigma[0](i, j, k) = var_Sigma_max_x * pow(distanceE<ftypePML>(Nx, delta_x, i), n);
				arr_sigma[3](i, j, k) = var_Sigma_max_x * pow(distanceH<ftypePML>(Nx, delta_x, i), n);

				arr_sigma[1](i, j, k) = var_Sigma_max_y * pow(distanceE<ftypePML>(Ny, delta_y, j), n);
				arr_sigma[4](i, j, k) = var_Sigma_max_y * pow(distanceH<ftypePML>(Ny, delta_y, j), n);

				arr_sigma[2](i, j, k) = var_Sigma_max_z * pow(distanceE<ftypePML>(Nz, delta_z, k), n);
				arr_sigma[5](i, j, k) = var_Sigma_max_z * pow(distanceH<ftypePML>(Nz, delta_z, k), n);
			}
}

template <class ftypePML, class ftype = double>
void Init_Coeff_SoA(SplitFields<ftypePML>& split_field,
	std::vector<data3d<ftypePML> >& Sigma, ftype dt) {

	double Exy1, Exz1, Ezx1, Ezy1, Eyx1, Eyz1;
	double Bxy1, Bxz1, Bzx1, Bzy1, Byx1, Byz1;

	int Nx = split_field.Exy.Get_Nx();
	int Ny = split_field.Exy.Get_Ny();
	int Nz = split_field.Exy.Get_Nz();
	int delta_x = split_field.Exy.Get_deltaX();
	int delta_y = split_field.Exy.Get_deltaY();
	int delta_z = split_field.Exy.Get_deltaZ();

	for (int i = 1; i < Nx + 2 * delta_x + 1; i++)
		for (int j = 1; j < Ny + 2 * delta_y + 1; j++)
			for (int k = 1; k < Nz + 2 * delta_z + 1; k++) {
				if ((i >= delta_x + 1) && (i < Nx + delta_x + 1) &&
					(j >= delta_y + 1) && (j < Ny + delta_y + 1) &&
					(k >= delta_z + 1) && (k < Nz + delta_z + 1)){}
				else {
					Exy1 = exp(-dt * Sigma[1](i, j, k));
					Exz1 = exp(-dt * Sigma[2](i, j, k));
					Eyx1 = exp(-dt * Sigma[0](i, j, k));
					Eyz1 = exp(-dt * Sigma[2](i, j, k));
					Ezx1 = exp(-dt * Sigma[0](i, j, k));
					Ezy1 = exp(-dt * Sigma[1](i, j, k));

					Bxy1 = exp(-dt * Sigma[4](i, j, k));
					Bxz1 = exp(-dt * Sigma[5](i, j, k));
					Byx1 = exp(-dt * Sigma[3](i, j, k));
					Byz1 = exp(-dt * Sigma[5](i, j, k));
					Bzx1 = exp(-dt * Sigma[3](i, j, k));
					Bzy1 = exp(-dt * Sigma[4](i, j, k));

					split_field.coeff[0](i, j, k) = (ftypePML)Exy1;
					split_field.coeff[2](i, j, k) = (ftypePML)Exz1;
					split_field.coeff[4](i, j, k) = (ftypePML)Eyx1;
					split_field.coeff[6](i, j, k) = (ftypePML)Eyz1;
					split_field.coeff[8](i, j, k) = (ftypePML)Ezx1;
					split_field.coeff[10](i, j, k) = (ftypePML)Ezy1;

					split_field.coeff[12](i, j, k) = (ftypePML)Bxy1;
					split_field.coeff[14](i, j, k) = (ftypePML)Bxz1;
					split_field.coeff[16](i, j, k) = (ftypePML)Byx1;
					split_field.coeff[18](i, j, k) = (ftypePML)Byz1;
					split_field.coeff[20](i, j, k) = (ftypePML)Bzx1;
					split_field.coeff[22](i, j, k) = (ftypePML)Bzy1;

					if (Sigma[0](i, j, k) != (ftypePML)0.0) {
						split_field.coeff[5](i, j, k) = 1.0 / Sigma[0](i, j, k) - Eyx1 / Sigma[0](i, j, k);
						split_field.coeff[9](i, j, k) = 1.0 / Sigma[0](i, j, k) - Ezx1 / Sigma[0](i, j, k);
					} else {
						split_field.coeff[5](i, j, k) = dt;
						split_field.coeff[9](i, j, k) = dt;
					}
					if (Sigma[1](i, j, k) != (ftypePML)0.0) {
						split_field.coeff[1](i, j, k) = 1.0 / Sigma[1](i, j, k) - Exy1 / Sigma[1](i, j, k);
						split_field.coeff[11](i, j, k) = 1.0 / Sigma[1](i, j, k) - Ezy1 / Sigma[1](i, j, k);
					} else {
						split_field.coeff[1](i, j, k) = dt;
						split_field.coeff[11](i, j, k) = dt;
					}
					if (Sigma[2](i, j, k) != (ftypePML)0.0) {
						split_field.coeff[3](i, j, k) = 1.0 / Sigma[2](i, j, k) - Exz1 / Sigma[2](i, j, k);
						split_field.coeff[7](i, j, k) = 1.0 / Sigma[2](i, j, k) - Eyz1 / Sigma[2](i, j, k);
					} else {
						split_field.coeff[3](i, j, k) = dt;
						split_field.coeff[7](i, j, k) = dt;
					}
					if (Sigma[3](i, j, k) != (ftypePML)0.0) {
						split_field.coeff[17](i, j, k) = 1.0 / Sigma[3](i, j, k) - Byx1 / Sigma[3](i, j, k);
						split_field.coeff[21](i, j, k) = 1.0 / Sigma[3](i, j, k) - Bzx1 / Sigma[3](i, j, k);
					} else {
						split_field.coeff[17](i, j, k) = dt;
						split_field.coeff[21](i, j, k) = dt;
					}
					if (Sigma[4](i, j, k) != (ftypePML)0.0) {
						split_field.coeff[13](i, j, k) = 1.0 / Sigma[4](i, j, k) - Bxy1 / Sigma[4](i, j, k);
						split_field.coeff[23](i, j, k) = 1.0 / Sigma[4](i, j, k) - Bzy1 / Sigma[4](i, j, k);
					} else {
						split_field.coeff[13](i, j, k) = dt;
						split_field.coeff[23](i, j, k) = dt;
					}
					if (Sigma[5](i, j, k) != (ftypePML)0.0) {
						split_field.coeff[15](i, j, k) = 1.0 / Sigma[5](i, j, k) - Bxz1 / Sigma[5](i, j, k);
						split_field.coeff[19](i, j, k) = 1.0 / Sigma[5](i, j, k) - Byz1 / Sigma[5](i, j, k);
					} else {
						split_field.coeff[15](i, j, k) = dt;
						split_field.coeff[19](i, j, k) = dt;
					}
				}
			}
}

template <class ftype>
void Update_boundary_condit(Fields<ftype>& field, ftype dx, ftype dy, ftype dz, const ftype dt, int it,
	std::pair <ftype, ftype> ab, std::pair <ftype, ftype> cd, std::pair <ftype, ftype> fg)
{
	int Nx = field.Ex.Get_Nx();
	int Ny = field.Ex.Get_Ny();
	int Nz = field.Ex.Get_Nz();
	int delta_x = field.Ex.Get_deltaX();
	int delta_y = field.Ex.Get_deltaY();
	int delta_z = field.Ex.Get_deltaZ();
	ftype x1, x2, y1, y2, z1, z2, t1, t2;

	ftype ax_transverse = ab.second * (ftype)3. / (ftype)4.; // поперечное
	ftype ay_transverse = cd.second * (ftype)1. / (ftype)6.;
	ftype az_transverse = fg.second * (ftype)1. / (ftype)6.;

	//ftype ax_transverse = ab.second * (ftype)1. / (ftype)2.; // поперечное
	//ftype ay_transverse = cd.second * (ftype)1. / (ftype)2.;
	//ftype az_transverse = fg.second * (ftype)1. / (ftype)2.;

	ftype tp_x = ab.second / (ftype)6.;
	ftype tp_y = cd.second / (ftype)6.;
	ftype tp_z = fg.second / (ftype)6.;

	ftype lambda_ = (ftype)2. * (ftype)M_PI / (ftype)4.; // wvelength               //?
	ftype k_ = (ftype)2. * (ftype)M_PI / lambda_; // wavenumber = 2 pi/ wavelength
	ftype omega_ = (ftype)2. * (ftype)M_PI / lambda_; // 2 pi c /wavelength 
	ftype A = (ftype)1.; //amplitude

	t1 = dt * (ftype)it;
	t2 = t1 + (ftype)0.5 * dt;

	ftype x0 = (ab.second - ab.first) / (ftype)2.;
	ftype y0 = (cd.second - cd.first) / (ftype)2.;
	ftype z0 = (fg.second - fg.first) / (ftype)2.;
	ftype t0_x = (ftype)3. * tp_x;
	ftype t0_y = (ftype)3. * tp_y;
	ftype t0_z = (ftype)3. * tp_z;

	int index_start_x = delta_x + 1, index_start_y = delta_y + 1, index_start_z = delta_z + 1;

	int offset = 5;
	x1 = (ftype)(ab.first + dx * (ftype)offset + (ftype)0.5 * dx);
	x2 = (ftype)(ab.first + dx * (ftype)offset);

#pragma omp parallel for collapse(2)
	for (int j = 0; j < Ny; j++)
		for (int k = 0; k < Nz; k++) {
			y1 = (ftype)(cd.first + dy * (ftype)j + (ftype)0.5 * dy);
			y2 = (ftype)(cd.first + dy * (ftype)j);

			z1 = (ftype)(fg.first + dz * (ftype)k + (ftype)0.5 * dz);
			z2 = (ftype)(fg.first + dz * (ftype)k);

			// Ez and By has the same y coordinate
			field.Ey(offset + index_start_x, j + index_start_y, k + index_start_z) += A * exp(-(t1 - t0_x) * (t1 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t1 - k_ * x1); // sin(phase), phase = omega * t - k * x
			field.Bz(offset + index_start_x, j + index_start_y, k + index_start_z) += A * exp(-(t2 - t0_x) * (t2 - t0_x) / (tp_x * tp_x))
				* exp(-(y2 - y0) * (y2 - y0) / (ay_transverse * ay_transverse)) *
				sin(omega_ * t2 - k_ * x2); // sin(phase), phase = omega * t - k * x
		}
}

