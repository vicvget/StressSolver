#ifndef STIFFNESS_RHS_H

#define STIFFNESS_RHS_H


#include "AuxiliaryStressStuff.h"

#include <vector>
#include <string>


// вектор правых частей СЛАУ на основе матрицы жесткости (RHS - right hand side)
using StiffnessRHSVector = std::vector<double>;


/**
* Применить список наборов флагов закрепленных степеней свободы для элементов сеточного представления тела
* к вектору правых частей СЛАУ на основе матрицы жесткости (соответствующие элементы вектора будут обнулены)
* @param stiffnessRHSVector - вектор правых частей СЛАУ на основе матрицы жесткости
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
*/
void ApplySealedElementsToStiffnessRHSVector
	(
		StiffnessRHSVector& stiffnessRHSVector,
		const IsSealedFlagsList& isSealedFlagsList
	);

/**
* Инвертировать вектор правых частей СЛАУ на основе матрицы жесткости
* @param stiffnessRHSVector - вектор правых частей СЛАУ на основе матрицы жесткости
*/
void InvertStiffnessRHSVector
	(
		StiffnessRHSVector& stiffnessRHSVector
	);

/**
* Записать вектор правых частей СЛАУ на основе матрицы жесткости в текстовый файл
* @param stiffnessRHSTextFileName - наименование файла для записи вектора правых частей СЛАУ
* на основе матрицы жесткости
* @param stiffnessRHSVector - вектор правых частей СЛАУ на основе матрицы жесткости
*/
void WriteStiffnessRHSVectorToTextFile
	(
		const std::string& stiffnessRHSTextFileName,
		const StiffnessRHSVector& stiffnessRHSVector
	);

/**
* Записать вектор правых частей СЛАУ на основе матрицы жесткости в двоичный файл
* @param stiffnessRHSTextFileName - наименование файла для записи вектора правых частей СЛАУ
* на основе матрицы жесткости
* @param stiffnessRHSVector - вектор правых частей СЛАУ на основе матрицы жесткости
*/
void WriteStiffnessRHSVectorToBinaryFile
	(
		const std::string& stiffnessRHSBinaryFileName,
		const StiffnessRHSVector& stiffnessRHSVector
	);

/**
* Считать вектор правых частей СЛАУ на основе матрицы жесткости из двоичного файла
* @param stiffnessRHSTextFileName - наименование файла для чтения вектора правых частей СЛАУ
* на основе матрицы жесткости
* @return считанный вектор правых частей СЛАУ на основе матрицы жесткости
*/
StiffnessRHSVector ReadStiffnessRHSVectorFromBinaryFile
	(
		const std::string& stiffnessRHSBinaryFileName
	);


#endif // STIFFNESS_RHS_H