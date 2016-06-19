#include "Matrix3x3.h"
#include <cmath>

#include <algorithm>


namespace MathHelpers
{

	Mat3::Mat3()
		:
			_mat(_data)
	{
		for (int i = 0; i < STRIDE3X3; i++)
		{
			_data[i] = 0;
		}
	}

	Mat3::Mat3
		(
			Vec3 row1,
			Vec3 row2,
			Vec3 row3
		)
		:
			_mat(_data)
	{
		for (int i = 0; i < STRIDE3; i++)
		{
			_data[i * STRIDE3] = row1[i];
			_data[i * STRIDE3 + 1] = row2[i];
			_data[i * STRIDE3 + 2] = row3[i];
		}
	}

	Mat3::Mat3
		(
			double* data
		)
		:
		_mat(_data)
	{
		std::copy(&data[0], &data[9], _data);
	}

	Mat3::Mat3
		(
			const Mat3& copy
		)
		:
			_mat(_data)
	{
		for (int i = 0; i < STRIDE3X3; i++)
		{
			_data[i] = copy[i];
		}
	}

	void Mat3::ExportRow(unsigned int rowId, double* data) const
	{
		for(int i = 0; i < 3; i++)
			data[i] = _mat[rowId*3+i];
	}

	void Mat3::Export(double* data) const
	{
		std::copy(_mat, _mat + STRIDE3X3, data);
	}

	void Mat3::Export(float* data) const
	{
		std::transform
			(
				_mat,
				_mat + STRIDE3X3,
				data,
				[](double val)
				{
					return static_cast<float>(val);
				}
			);
	}

	Mat3& Mat3::operator = (const Mat3& copy)
	{
		for(int i = 0; i < 9; i++)
		{
			_data[i] = copy[i];
			_mat = &_data[0];
		}

		return *this;
	}

	double Mat3::operator [] (unsigned int id) const
	{
		return _mat[id];
	}

	double& Mat3::operator [] (unsigned int id)
	{
		return _mat[id];
	}

	double& Mat3::E(int i, int j)
	{
		return _data[i*STRIDE3 + j];
	}

	double Mat3::E(int i, int j) const
	{
		return _data[i*STRIDE3 + j];
	}



	Vec3Ref Mat3::Row(unsigned int rowId)
	{
		return MakeVec3(_mat + rowId * STRIDE3);
	}

	Vec3CRef Mat3::Row(unsigned int rowId) const
	{
		return MakeVec3(_mat + rowId * STRIDE3);
	}

	//Vec3 Mat3::operator * (const Vec3& vec) const
	//{
	//	Vec3 res;
	//	for(int i = 0; i < 3; i++)
	//		for(int j = 0; j < 3; j++)
	//		{
	//			res[i] += _mat[j*3+i]*vec[j];
	//		}
	//		return res;
	//}

	//Mat3 Mat3::operator * (const Mat3& mat) const
	//{
	//	Mat3 res;
	//	for(int i = 0; i < 3; i++)
	//		for(int j = 0; j < 3; j++)
	//			for(int k = 0; k < 3; k++)
	//			{
	//				res[j*3+i] += _mat[k*3+i]*mat[j*3+k];
	//			}
	//			return res;
	//}

	//Vec3 Mat3::Tmul(const Vec3& vec) const
	//{
	//	Vec3 res;
	//	for(int i = 0; i < 3; i++)
	//		for(int j = 0; j < 3; j++)
	//		{
	//			res[i] += _mat[i*3+j]*vec[j];
	//		}
	//		return res;
	//}

	//Mat3 Mat3::Tmul(const Mat3& mat) const
	//{
	//	Mat3 res;
	//	for (int i = 0; i < 3; i++)
	//		for (int j = 0; j < 3; j++)
	//			for (int k = 0; k < 3; k++)
	//			{
	//				res[j * 3 + i] += _mat[i * 3 + k] * mat[j * 3 + k];
	//			}
	//	return res;
	//}


	Vec3 Mat3::operator * (const Vec3& vec) const
	{
		Vec3 res;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				res[i] += _mat[i * 3 + j] * vec[j];
			}
		return res;
	}

	Mat3 Mat3::operator * (const Mat3& mat) const
	{
		Mat3 res;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
				{
					res[i * 3 + j] += _mat[i * 3 + k] * mat[k * 3 + j];
				}
		return res;
	}

	Vec3 Mat3::Tmul(const Vec3& vec) const
	{
		Vec3 res;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				res[i] += _mat[j * 3 + i] * vec[j];
			}
		return res;
	}

	Mat3 Mat3::Tmul(const Mat3& mat) const
	{
		Mat3 res;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
				{
					res[i * 3 + j] += _mat[k * 3 + i] * mat[k * 3 + j];
				}
		return res;
	}


	// скалярное произведение столбца матрицы на вектор
	double Mat3::Cdot(int colId, const Vec3& vec) const // умножение колонки матрицы
	{
		double res = 0.;
		for(int i = 0; i < 3; i++)
			res += _mat[colId*3+i]*vec[i];
		return res;
	}

	// векторное произведение строки матрицы на вектор
	Vec3 Mat3::Rcross(int rowId, const Vec3& vec) const // умножение колонки матрицы
	{
		Vec3 row;
		row.Init(_mat[rowId], _mat[rowId+3], _mat[rowId+6]);
		return row.Cross(vec);
	}

	Mat3 Mat3::Identity()
	{
		Mat3 mtxI;
		std::fill(mtxI._data,mtxI._data+STRIDE3X3, 0.);
		mtxI[0] = 1.;
		mtxI[4] = 1.;
		mtxI[8] = 1.;
		return mtxI;
	}

	Mat3 Mat3::SetRotation(const Vec3& axis, double angle)
	{
		Mat3 mtx;
		double cosA = cos(angle);
		double sinA = sin(angle);
		double ncoA = 1-cos(angle);

		mtx.Row(0)[0] = axis[0]*axis[0]*ncoA + cosA;
		mtx.Row(0)[1] = axis[0]*axis[1]*ncoA - axis[2]*sinA;
		mtx.Row(0)[2] = axis[0]*axis[2]*ncoA + axis[1]*sinA;

		mtx.Row(1)[0] = axis[0]*axis[1]*ncoA + axis[2]*sinA;
		mtx.Row(1)[1] = axis[1]*axis[1]*ncoA + cosA;
		mtx.Row(1)[2] = axis[1]*axis[2]*ncoA - axis[0]*sinA;
		
		mtx.Row(2)[0] = axis[0]*axis[2]*ncoA - axis[1]*sinA;
		mtx.Row(2)[1] = axis[0]*axis[1]*ncoA + axis[0]*sinA;
		mtx.Row(2)[2] = axis[2]*axis[2]*ncoA + cosA;

		return mtx;
	}

	Mat3 Mat3::MakeXYZRotationMtx10(double* angles)
	{
		return MakeXYZRotationMtx01(angles).Tr();
	}

	Mat3 Mat3::MakeXYZRotationMtx01(double* angles)
	{
		double cosX = cos(angles[0]);
		double sinX = sin(angles[0]);
		double cosY = cos(angles[1]);
		double sinY = sin(angles[1]);
		double cosZ = cos(angles[2]);
		double sinZ = sin(angles[2]);

		Mat3 res;

		// A01

		res.E(0, 0) = cosY*cosZ;
		res.E(1, 0) = -cosY*sinZ;
		res.E(2, 0) = sinY;
		res.E(0, 1) = sinX*sinY*cosZ + cosX*sinZ;
		res.E(1, 1) = -sinX*sinY*sinZ + cosX*cosZ;
		res.E(2, 1) = -sinX*cosY;
		res.E(0, 2) = -cosX*sinY*cosZ + sinX*sinZ;
		res.E(1, 2) = cosX*sinY*sinZ + sinX*cosZ;
		res.E(2, 2) = cosX*cosY;

		return res;
	}



	Mat3 Mat3::Tr() const
	{
		Mat3 mtx(*this);
		std::swap(mtx.Row(0)[1], mtx.Row(1)[0]);
		std::swap(mtx.Row(0)[2], mtx.Row(2)[0]);
		std::swap(mtx.Row(1)[2], mtx.Row(2)[1]);
		return mtx;
	}

}
