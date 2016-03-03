#ifndef MATRIX_3_X_3_H

#define MATRIX_3_X_3_H


#include "Vector3.h"


namespace MathHelpers
{

	/**
	* Класс для матрицы 3x3 (хранение по столбцам)
	*/
	class Mat3
	{
	public:

		Mat3();
		
		Mat3
			(
				Vec3 row1,
				Vec3 row2,
				Vec3 row3
			);

		Mat3
			(
				double* data
			);

		Mat3
			(
				const Mat3& copy
			);


		void ExportRow(unsigned int rowId, double* data) const;

		void Export(double* data) const;

		void Export(float* data) const;


		Mat3& operator = (const Mat3& copy);

		double operator [] (unsigned int id) const;

		double& operator [] (unsigned int id);

		Vec3Ref Row (unsigned int rowId);

		Vec3CRef Row (unsigned int rowId) const;

		Vec3 operator * (const Vec3& vec) const;

		Mat3 operator * (const Mat3& mat) const;

		Mat3 Tr() const;

		Vec3 Tmul(const Vec3& vec) const;

		// скалярное произведение столбца матрицы на вектор
		double Cdot(int colId, const Vec3& vec) const; // умножение колонки матрицы

		// векторное произведение строки матрицы на вектор
		Vec3 Rcross(int rowId, const Vec3& vec) const; // умножение колонки матрицы

		Mat3 Tmul(const Mat3& mat) const;


		static
		Mat3 Identity();

		static
		Mat3 SetRotation(const Vec3& axis, double angle);
		double& E(int i, int j);
		double E(int i, int j) const;

		
		static Mat3 MakeXYZRotationMtx01(double* angles);
		static Mat3 MakeXYZRotationMtx10(double* angles);


	private:

		double* _mat;
		double _data[STRIDE3X3];

	};

}


#endif // MATRIX_3_X_3_H
