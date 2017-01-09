#include "Matrix3x4.h"
#include "Strides.h"

#include <cmath>
#include <algorithm>


namespace MathHelpers
{

	// class Mat3x4, хранение по строкам

	Mat3x4::Mat3x4()
	{
		for (int i = 0; i < 12; i++)
			_data[i] = 0;
		_mat = &_data[0];
	}

	Mat3x4::Mat3x4(Vec3 row1, Vec3 row2, Vec3 row3)
	{
		for (int i = 0; i < 3; i++)
		{
			_data[i * STRIDE4] = row1[i];
			_data[i * STRIDE4 + 1] = row2[i];
			_data[i * STRIDE4 + 2] = row3[i];
		}
		_mat = &_data[0];
	}

	Mat3x4::Mat3x4(double* data)
		:
			_mat(_data),
			_data()
	{
		for (int i = 0; i < 12; i++)
			_data[i] = data[i];
		_mat = &_data[0];
	}

	Mat3x4::Mat3x4(const Mat3x4& copy)
	{
		for (int i = 0; i < 12; i++)
		{
			_data[i] = copy._data[i];
			_mat = &_data[0];
		}
	}

	Mat3x4 Mat3x4::Identity()
	{
		Mat3x4 res;
		for (int i = 0; i < 3; i++)
			res.E(i, i) = 1.;
		return res;
	}

	void Mat3x4::ExportRow(unsigned int rowId, double* data) const
	{
		for (int i = 0; i < 3; i++)
		{
			data[i] = _mat[rowId * 4 + i];
		}
		data[3] = 0;
	}

	Mat3x4& Mat3x4::operator = (const Mat3x4& copy)
	{
		for (int i = 0; i < 12; i++)
		{
			_data[i] = copy[i];
			_mat = &_data[0];
		}

		return *this;
	}

	double Mat3x4::operator [] (unsigned int id) const
	{
		return _mat[id];
	}

	double& Mat3x4::operator [] (unsigned int id)
	{
		return _mat[id];
	}

	Vec3 Mat3x4::operator * (const Vec3& vec) const
	{
		Vec3 res;

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				res[i] += E(i,j) * vec[j];
			}
		}

		return res;
	}

	Vec3 Mat3x4::Tmul(const Vec3& vec) const
	{
		Vec3 res;

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				res[i] += E(j, i) * vec[j];
			}
		}
		
		return res;
	}

	//// скалярное произведение столбца матрицы на вектор
	//double Mat3x4::Cdot(int colId, const Vec3& vec) const // умножение колонки матрицы
	//{
	//	double res = 0.0;

	//	for (int i = 0; i < 3; i++)
	//	{
	//		res += _mat[colId * 4 + i] * vec[i];
	//	}

	//	return res;
	//}

	//// векторное произведение строки матрицы на вектор
	//Vec3 Mat3x4::Rcross(int rowId, const Vec3& vec) const // умножение колонки матрицы
	//{
	//	Vec3 row;

	//	row.Init(_mat[rowId], _mat[rowId + 4], _mat[rowId + 8]);

	//	return row.Cross(vec);
	//}

	Mat3x4 Mat3x4::Tmul(const Mat3x4& mat) const
	{
		Mat3x4 res;

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					res.E(i,j) += E(k,i) * mat.E(k,j);
				}
			}
		}
		
		return res;
	}

	Mat3x4 Mat3x4::operator * (const Mat3x4& mat) const
	{
		Mat3x4 res;

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					res.E(i,j) += E(i,k) * mat.E(k,j);
				}
			}
		}
		//
		//for (int j = 0; j < 3; j++)
		//{
		//	res.E(j, 3) = E(j, 3);
		//	for (int k = 0; k < 3; k++)
		//	{
		//		res.E(j,3) += E(j,k) * mat.E(k,3);
		//	}
		//}


		return res;
	}

	void Mat3x4::Export(double* data) const
	{
		std::copy(_mat, _mat + STRIDE3X4, data);
	}


	void Mat3x4::FImport(float* array)
	{
		for (int i = 0; i < 3; i++)
		{
			E(i, 3) = array[i];
			for (int j = 0; j < 3; j++)
				E(j, i) = array[(i + 1) * 3 + j];
		}
	}

	void Mat3x4::FExport(float* array) const
	{
		for (int i = 0; i < 3; i++)
		{
			array[i] = static_cast<float>(E(i, 3));
			for (int j = 0; j < 3; j++)
				array[(i + 1) * 3 + j] = static_cast<float>(E(j, i));
		}
	}

	void Mat3x4::SetTranslation(const Vec3& translation)
	{
		for (int i = 0; i < 3; i++)
			E(i, 3) = translation[i];
	}

	void Mat3x4::SetRotation(const Vec3& axis, double angle)
	{
		double cosA = cos(angle);
		double sinA = sin(angle);
		double ncoA = 1 - cos(angle);

		E(0, 0) = axis[0] * axis[0] * ncoA + cosA;
		E(0, 1) = axis[0] * axis[1] * ncoA - axis[2] * sinA;
		E(0, 2) = axis[0] * axis[2] * ncoA + axis[1] * sinA;
		E(1, 0) = axis[0] * axis[1] * ncoA + axis[2] * sinA;
		E(1, 1) = axis[1] * axis[1] * ncoA + cosA;
		E(1, 2) = axis[1] * axis[2] * ncoA - axis[0] * sinA;
		E(2, 0) = axis[0] * axis[2] * ncoA - axis[1] * sinA;
		E(2, 1) = axis[0] * axis[1] * ncoA + axis[0] * sinA;
		E(2, 2) = axis[2] * axis[2] * ncoA + cosA;
	}

		
	Mat3x4 Mat3x4::MakeXYZRotationMtx10(double* angles)
	{
		return MakeXYZRotationMtx01(angles).Tr();
	}

	Mat3x4 Mat3x4::MakeXYZRotationMtx01(double* angles)
	{
		double cosX = cos(angles[0]);
		double sinX = sin(angles[0]);
		double cosY = cos(angles[1]);
		double sinY = sin(angles[1]);
		double cosZ = cos(angles[2]);
		double sinZ = sin(angles[2]);

		Mat3x4 res;

		// A01

		res.E(0, 0) = cosY*cosZ;
		res.E(0, 1) = -cosY*sinZ;
		res.E(0, 2) = sinY;
		res.E(1, 0) = sinX*sinY*cosZ + cosX*sinZ;
		res.E(1, 1) = -sinX*sinY*sinZ + cosX*cosZ;
		res.E(1, 2) = -sinX*cosY;
		res.E(2, 0) = -cosX*sinY*cosZ + sinX*sinZ;
		res.E(2, 1) = cosX*sinY*sinZ + sinX*cosZ;
		res.E(2, 2) = cosX*cosY;

		return res;
	}

	Mat3x4 Mat3x4::MakeXYZRotationMtx01(double* sinarr, double* cosarr)
	{
		double cosX = cosarr[0];
		double sinX = sinarr[0];
		double cosY = cosarr[1];
		double sinY = sinarr[1];
		double cosZ = cosarr[2];
		double sinZ = sinarr[2];

		Mat3x4 res;

		// A01

		res.E(0, 0) = cosY*cosZ;
		res.E(0, 1) = -cosY*sinZ;
		res.E(0, 2) = sinY;
		res.E(1, 0) = sinX*sinY*cosZ + cosX*sinZ;
		res.E(1, 1) = -sinX*sinY*sinZ + cosX*cosZ;
		res.E(1, 2) = -sinX*cosY;
		res.E(2, 0) = -cosX*sinY*cosZ + sinX*sinZ;
		res.E(2, 1) = cosX*sinY*sinZ + sinX*cosZ;
		res.E(2, 2) = cosX*cosY;

		return res;
	}


	Mat3x4 Mat3x4::Tr() const
	{
		Mat3x4 mtx(*this);
		for (int i = 0; i < 3; i++)
			for (int j = i + 1; j < 3; j++)
				std::swap(mtx.E(i, j), mtx.E(j, i));
		return mtx;
	}

}
