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
		/*
		if (i<0 || j<0 || l<0 || i>=(n + 2 * delta_x + 2) || j>= (m + 2 * delta_y + 2) || l>= (k + 2 * delta_z + 2)){
			std::cout << "ERROR" << std::endl;
		}*/
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
	/*~data3d() {
		data.~vector();
	}*/
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
	/*~Fields() {
		Ex.~data3d();
		Ey.~data3d();
		Ez.~data3d();
		Bx.~data3d();
		By.~data3d();
		Bz.~data3d();
	}*/
};

template <typename type_data>
struct SplitFields {
	data3d<type_data> Exy, Exz, Eyx, Eyz, Ezx, Ezy,
						Bxy, Bxz, Byx, Byz, Bzx, Bzy;
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
	}
	/*~SplitFields() {
		Exy.~data3d();
		Exz.~data3d();
		Eyx.~data3d();
		Eyz.~data3d();
		Ezx.~data3d();
		Ezy.~data3d();

		Bxy.~data3d();
		Bxz.~data3d();
		Byx.~data3d();
		Byz.~data3d();
		Bzx.~data3d();
		Bzy.~data3d();
	}*/
};

template <typename type_data>
struct Sigma {
	data3d<type_data> sigmaEx, sigmaEy, sigmaEz, sigmaBx, sigmaBy, sigmaBz;
	Sigma(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		sigmaEx.Create(n, m, k, delta_x, delta_y, delta_z);
		sigmaEy.Create(n, m, k, delta_x, delta_y, delta_z);
		sigmaEz.Create(n, m, k, delta_x, delta_y, delta_z);

		sigmaBx.Create(n, m, k, delta_x, delta_y, delta_z);
		sigmaBy.Create(n, m, k, delta_x, delta_y, delta_z);
		sigmaBz.Create(n, m, k, delta_x, delta_y, delta_z);
	}
	/*~Sigma() {
		sigmaEx.~data3d();
		sigmaEy.~data3d();
		sigmaEz.~data3d();
		sigmaBx.~data3d();
		sigmaBy.~data3d();
		sigmaBz.~data3d();
	}*/
};
template <typename type_data>
struct Coefficient {
	data3d<type_data> Exy1, Exz1, Eyx1, Eyz1, Ezx1, Ezy1, Bxy1, Bxz1, Byx1, Byz1, Bzx1, Bzy1,
		Exy2, Exz2, Eyx2, Eyz2, Ezx2, Ezy2, Bxy2, Bxz2, Byx2, Byz2, Bzx2, Bzy2;
	Coefficient(int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		Exy1.Create(n, m, k, delta_x, delta_y, delta_z);
		Exz1.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyx1.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyz1.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezx1.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezy1.Create(n, m, k, delta_x, delta_y, delta_z);

		Bxy1.Create(n, m, k, delta_x, delta_y, delta_z);
		Bxz1.Create(n, m, k, delta_x, delta_y, delta_z);
		Byx1.Create(n, m, k, delta_x, delta_y, delta_z);
		Byz1.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzx1.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzy1.Create(n, m, k, delta_x, delta_y, delta_z);

		Exy2.Create(n, m, k, delta_x, delta_y, delta_z);
		Exz2.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyx2.Create(n, m, k, delta_x, delta_y, delta_z);
		Eyz2.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezx2.Create(n, m, k, delta_x, delta_y, delta_z);
		Ezy2.Create(n, m, k, delta_x, delta_y, delta_z);

		Bxy2.Create(n, m, k, delta_x, delta_y, delta_z);
		Bxz2.Create(n, m, k, delta_x, delta_y, delta_z);
		Byx2.Create(n, m, k, delta_x, delta_y, delta_z);
		Byz2.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzx2.Create(n, m, k, delta_x, delta_y, delta_z);
		Bzy2.Create(n, m, k, delta_x, delta_y, delta_z);
	}
	/*~Coefficient() {
		Exy1.~data3d();
		Exz1.~data3d();
		Eyx1.~data3d();
		Eyz1.~data3d();
		Ezx1.~data3d();
		Ezy1.~data3d();

		Bxy1.~data3d();
		Bxz1.~data3d();
		Byx1.~data3d();
		Byz1.~data3d();
		Bzx1.~data3d();
		Bzy1.~data3d();

		Exy2.~data3d();
		Exz2.~data3d();
		Eyx2.~data3d();
		Eyz2.~data3d();
		Ezx2.~data3d();
		Ezy2.~data3d();

		Bxy2.~data3d();
		Bxz2.~data3d();
		Byx2.~data3d();
		Byz2.~data3d();
		Bzx2.~data3d();
		Bzy2.~data3d();
	}*/

};