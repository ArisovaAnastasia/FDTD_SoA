#pragma once
#include <vector>

template <typename type_data>
class data3d {
	std::vector<type_data> data;
	int n, m, k, delta_x, delta_y, delta_z;
public:
	void Create(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		data.resize((n + 2 * delta_x + 2) * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2));

		this->n = n; this->m = m; this->k = k;
		this->delta_x = delta_x; this->delta_y = delta_y; this->delta_z = delta_z;
	}
	inline type_data& operator()(int i, int j, int l) {
		return data[i * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2) + j * (k + 2 * delta_z + 2) + l];
	}
	int Get_Nx() {
		return n;
	}
	int Get_Ny() {
		return m;
	}
	int Get_Nz() {
		return k;
	}
	int Get_deltaX() {
		return delta_x;
	}
	int Get_deltaY() {
		return delta_y;
	}
	int Get_deltaZ() {
		return delta_z;
	}
};

template <typename type_data>
struct Fields {
	data3d<type_data> Ex, Ey, Ez, Bx, By, Bz;
	Fields(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		Ex.Create(n, m, k, delta_x, delta_y, delta_z);
		Ey.Create(n, m, k, delta_x, delta_y, delta_z);
		Ez.Create(n, m, k, delta_x, delta_y, delta_z);

		Bx.Create(n, m, k, delta_x, delta_y, delta_z);
		By.Create(n, m, k, delta_x, delta_y, delta_z);
		Bz.Create(n, m, k, delta_x, delta_y, delta_z);
	}
};

template <typename type_data>
struct SplitFields {
	data3d<type_data> Exy, Exz, Eyx, Eyz, Ezx, Ezy,
						Bxy, Bxz, Byx, Byz, Bzx, Bzy;
	std::vector<data3d<type_data> > coeff;
	SplitFields(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		Exy.Create(n, m, k, delta_x, delta_y, delta_z);
		Exz.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyx.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyz.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezx.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezy.Create(n, m, k, delta_x, delta_y, delta_z);

		Bxy.Create(n, m, k, delta_x, delta_y, delta_z);
		Bxz.Create(n, m, k, delta_x, delta_y, delta_z);
		Byx.Create(n, m, k, delta_x, delta_y, delta_z);
		Byz.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzx.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzy.Create(n, m, k, delta_x, delta_y, delta_z);

		coeff.resize(24);
		for (int i = 0; i < 24;++i) {
			coeff[i].Create(n, m, k, delta_x, delta_y, delta_z);
		}
	}
};
