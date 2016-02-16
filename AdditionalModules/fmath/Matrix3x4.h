#ifndef MATRIX_3_X_4_H

#define MATRIX_3_X_4_H


#include "Vector3.h"


namespace MathHelpers
{

	/**
	* Класс для матрицы 3x4 (хранение по столбцам)
	*/
	class Mat3x4
	{
	public:

		Mat3x4();

		Mat3x4(Vec3 row1, Vec3 row2, Vec3 row3);

		Mat3x4(double* _data);

		Mat3x4(const Mat3x4& copy);


		void FImport(float* array);

		void FExport(float* array) const;

		void ExportRow(unsigned int rowId, double* data) const;


		Mat3x4& operator = (const Mat3x4& copy);

		double& E(int i, int j)
		{
			return _data[i*STRIDE4 + j];
		}

		double E(int i, int j) const
		{
			return _data[i*STRIDE4 + j];
		}

		double operator [] (unsigned int id) const;

		double& operator [] (unsigned int id);

		Vec3 operator * (const Vec3& vec) const;


		Vec3 Tmul(const Vec3& vec) const;

		//// скалярное произведение столбца матрицы на вектор
		//double Cdot(int colId, const Vec3& vec) const; // умножение колонки матрицы

		//// векторное произведение строки матрицы на вектор
		//Vec3 Rcross(int rowId, const Vec3& vec) const; // умножение колонки матрицы

		Mat3x4 Tmul(const Mat3x4& mat) const;

		Mat3x4 operator * (const Mat3x4& mat) const;
	
		static Mat3x4 Identity();

		void SetRotation(const Vec3& axis, double angle);
		void SetTranslation(const Vec3& translation);
	private:

		double *_mat;
		double _data[STRIDE3X4];

	};

}


#endif // MATRIX_3_X_4_H