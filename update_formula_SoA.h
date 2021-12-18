#pragma once
#include "type_data.h"

template <class ftype>
void Update_electric_field_SoA(Fields<ftype>& field,
	ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tEx, tEy, tEz;

	tEx = (field.Bz(i, j + 1, k) - field.Bz(i, j, k)) * (dt_y)
		-(field.By(i, j, k + 1) - field.By(i, j, k)) * (dt_z) + field.Ex(i, j, k);

	tEy = (field.Bx(i, j, k + 1) - field.Bx(i, j, k)) * (dt_z)
		-(field.Bz(i + 1, j, k) - field.Bz(i, j, k)) * (dt_x) + field.Ey(i, j, k);

	tEz = (field.By(i + 1, j, k) - field.By(i, j, k)) * (dt_x)
		-(field.Bx(i, j + 1, k) - field.Bx(i, j, k)) * (dt_y) + field.Ez(i, j, k);

	field.Ex(i, j, k) = tEx;
	field.Ey(i, j, k) = tEy;
	field.Ez(i, j, k) = tEz;
}

template <class ftype>
void Update_magnetic_field_SoA(Fields<ftype>& field,
	ftype dt_x, ftype dt_y, ftype dt_z, int i, int j, int k)
{
	ftype tBx, tBy, tBz;

	tBx = -(field.Ez(i, j, k) - field.Ez(i, j - 1, k)) * (dt_y)
		+(field.Ey(i, j, k) - field.Ey(i, j, k - 1)) * (dt_z) + field.Bx(i, j, k);

	tBy = -(field.Ex(i, j, k) - field.Ex(i, j, k - 1)) * (dt_z)
		+(field.Ez(i, j, k) - field.Ez(i - 1, j, k)) * (dt_x) + field.By(i, j, k);

	tBz = -(field.Ey(i, j, k) - field.Ey(i - 1, j, k)) * (dt_x)
		+(field.Ex(i, j, k) - field.Ex(i, j - 1, k)) * (dt_y) + field.Bz(i, j, k);

	field.Bx(i, j, k) = tBx;
	field.By(i, j, k) = tBy;
	field.Bz(i, j, k) = tBz;
}

template <class ftype, class ftypePML>
void Update_electric_field_PML_SoA(Fields<ftype>& field, SplitFields<ftypePML>& split_field,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tExy, tExz, tEyx, tEyz, tEzx, tEzy;
 
	tExy = split_field.Exy(i, j, k) * split_field.coeff[0](i, j, k) +
		((ftypePML)field.Bz(i, j + 1, k) - (ftypePML)field.Bz(i, j, k)) * split_field.coeff[1](i, j, k) * (_1dy);

	tExz = split_field.Exz(i, j, k) * split_field.coeff[2](i, j, k) -
		((ftypePML)field.By(i, j, k + 1) - (ftypePML)field.By(i, j, k)) * split_field.coeff[3](i, j, k) * (_1dz);

	tEyx = split_field.Eyx(i, j, k) * split_field.coeff[4](i, j, k) -
		((ftypePML)field.Bz(i + 1, j, k) - (ftypePML)field.Bz(i, j, k)) * split_field.coeff[5](i, j, k) * (_1dx);

	tEyz = split_field.Eyz(i, j, k) * split_field.coeff[6](i, j, k) +
		((ftypePML)field.Bx(i, j, k + 1) - (ftypePML)field.Bx(i, j, k)) * split_field.coeff[7](i, j, k) * (_1dz);

	tEzx = split_field.Ezx(i, j, k) * split_field.coeff[8](i, j, k) +
		((ftypePML)field.By(i + 1, j, k) - (ftypePML)field.By(i, j, k)) * split_field.coeff[9](i, j, k) * (_1dx);

	tEzy = split_field.Ezy(i, j, k) * split_field.coeff[10](i, j, k) -
		((ftypePML)field.Bx(i, j + 1, k) - (ftypePML)field.Bx(i, j, k)) * split_field.coeff[11](i, j, k) * (_1dy);

	split_field.Exy(i, j, k) = tExy;
	split_field.Exz(i, j, k) = tExz;
	split_field.Eyx(i, j, k) = tEyx;
	split_field.Eyz(i, j, k) = tEyz;
	split_field.Ezx(i, j, k) = tEzx;
	split_field.Ezy(i, j, k) = tEzy;

	field.Ex(i, j, k) = (ftype)(tExy + tExz);
	field.Ey(i, j, k) = (ftype)(tEyx + tEyz);
	field.Ez(i, j, k) = (ftype)(tEzx + tEzy);
}

template <class ftype, class ftypePML>
void Update_magnetic_field_PML_SoA(Fields<ftype>& field, SplitFields<ftypePML>& split_field,
	ftypePML _1dx, ftypePML _1dy, ftypePML _1dz, int i, int j, int k)
{
	ftypePML tBxy, tBxz, tByx, tByz, tBzx, tBzy;

	tBxy = split_field.Bxy(i, j, k) * split_field.coeff[12](i, j, k) -
		((ftypePML)field.Ez(i, j, k) - (ftypePML)field.Ez(i, j - 1, k)) * split_field.coeff[13](i, j, k) * (_1dy);

	tBxz = split_field.Bxz(i, j, k) * split_field.coeff[14](i, j, k) +
		((ftypePML)field.Ey(i, j, k) - (ftypePML)field.Ey(i, j, k - 1)) * split_field.coeff[15](i, j, k) * (_1dz);

	tByx = split_field.Byx(i, j, k) * split_field.coeff[16](i, j, k) +
		((ftypePML)field.Ez(i, j, k) - (ftypePML)field.Ez(i - 1, j, k)) * split_field.coeff[17](i, j, k) * (_1dx);

	tByz = split_field.Byz(i, j, k) * split_field.coeff[18](i, j, k) -
		((ftypePML)field.Ex(i, j, k) - (ftypePML)field.Ex(i, j, k - 1)) * split_field.coeff[19](i, j, k) * (_1dz);

	tBzx = split_field.Bzx(i, j, k) * split_field.coeff[20](i, j, k) -
		((ftypePML)field.Ey(i, j, k) - (ftypePML)field.Ey(i - 1, j, k)) * split_field.coeff[21](i, j, k) * (_1dx);

	tBzy = split_field.Bzy(i, j, k) * split_field.coeff[22](i, j, k) +
		((ftypePML)field.Ex(i, j, k) - (ftypePML)field.Ex(i, j - 1, k)) * split_field.coeff[23](i, j, k) * (_1dy);

	split_field.Bxy(i, j, k) = tBxy;
	split_field.Bxz(i, j, k) = tBxz;
	split_field.Byx(i, j, k) = tByx;
	split_field.Byz(i, j, k) = tByz;
	split_field.Bzx(i, j, k) = tBzx;
	split_field.Bzy(i, j, k) = tBzy;

	field.Bx(i, j, k) = (ftype)(tBxy + tBxz);
	field.By(i, j, k) = (ftype)(tByx + tByz);
	field.Bz(i, j, k) = (ftype)(tBzx + tBzy);
}
