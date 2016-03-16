#ifndef RESULTS_READER_H

#define RESULTS_READER_H


#include "AbstractProviderRezr.h"

#include <fstream>


using std::ifstream;


/**
* Класс для чтения файла результатов ФРУНД
* записанных в формате С++
* @author Sergeev Efim
*/
class ResultsReader
	:
		public AbstractProviderRezr
{
public:

	/** Конструктор
	* @param rezrName - имя файла с результатами
	* @param fadresName - имя файла со структурой результатов
	*/
	ResultsReader
		(
			const string& rezrName,
			const string& fadresName
		);

	virtual
	~ResultsReader() override;

	/** Открытие файла результатов
	* @return false в случае ошибки
	*/
	virtual
	bool OpenStream() override;

	/** Закрыть поток результатов
	*/
	virtual
	void CloseStream() override;

	/** Считывание заголовка
	*/
	virtual
	void ReadHeader() override;

	/** Чтение кадра
	* @return false, если кадр не прочитан целиком или конец потока
	*/
	virtual
	bool TryReadFrame() override;

	/**
	* Перейти к кадру с данным номером
	* @param frameNumber - номер кадра
	* @return признак успешного (true) или неуспешного (false) перехода
	*/
	virtual
	bool GoToFrame
		(
			int frameNumber
		)	override;

	
private:

	//файл с результатами
	ifstream _resultsFile;

};


#endif // RESULTS_READER_H