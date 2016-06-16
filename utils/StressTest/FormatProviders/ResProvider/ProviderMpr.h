#ifndef PROVIDER_MPR_H

#define PROVIDER_MPR_H


#include "ResultsWriter.h"

#define SOLVER_TYPES_NUMBER 2
enum SolverTypes 
{
	ST_Thermal = 0,
	ST_StressStrain = 1
};

#include "ResultsReader.h"
#include "ResultsRawReader.h"


/** Заголовок для файла результатов междисциплинарного расчета
*
* @author Getmanskiy Victor
*/
struct MultiphysicsResultsHeader
{

	// символ форматирования
	char FormatSymbol;

	// тип
	char Type;

	// количество скалярных значений на узел
	char ScalarsCount;

	// количество векторных значений на узел
	char VectorsCount;

	// количество узлов
	int PointsCount;

	// модельное время начала расчета
	float InitialTime;

	// модельное время окончания расчета
	float FinalTime;

	MultiphysicsResultsHeader(){}

	MultiphysicsResultsHeader
		(
			SolverTypes type,
			size_t pointsCount,
			float initialTime,
			float totalTime
		);

};


/** Заголовок для кадра междисциплинарного расчета
*
* @author Getmanskiy Victor
*/
struct MultiphysicsFrameHeader
{
	// текущее время
	float CurrentTime;

	// размер кадра
	int FrameSize;

	MultiphysicsFrameHeader(){};

	MultiphysicsFrameHeader
		(
			int frameSize,
			float currentTime
		)
		:
			CurrentTime(currentTime),
			FrameSize(frameSize)
	{
	}

};


/**
* Провайдер файла результатов расчета сетки
*
* @author Getmanskiy VIctor
*/
class ProviderMpr
{
public:
	/**
	* Записывается заголовок // TODO: WTF!
	* @param mprFileName - имя файла с результатами
	* @param mprHeader - заголовок
	*/
	void InitWriter(const string& mprFileName, const MultiphysicsResultsHeader* mprHeader);
	void InitReader(const string& mprFileName, MultiphysicsResultsHeader* mprHeader);
	
	/** Запись кадра
	* @param data - данные
	* @param dataSize - размер данных в байтах
	* @params currentTime - текущее время кадра
	*/
	void WriteFrame
		(
			const void* data,
			size_t dataSize,
			float currentTime
		);

	/** Чтение кадра
	* @param data - данные
	* @param dataSize - размер данных в байтах
	* @params frameNumber - номер кадра, возвращает номер последнего считанного кадра, если < 0, то читает до последнего кадра
	*/
	bool ReadFrame(void* data, size_t& dataSize, int& frameNumber);

	ProviderMpr() = default;

private:

	std::unique_ptr<ResultsWriter> _resultsWriter;
	std::unique_ptr<ResultsRawReader> _resultsReader;

	// запрещаем использование конструктора копирования
	ProviderMpr(const ProviderMpr&);

};


#endif // PROVIDER_MPR_H