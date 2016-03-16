#ifndef ABSTRACT_PROVIDER_REZR_H

#define ABSTRACT_PROVIDER_REZR_H


#include "OldAddresses.h"

#include <vector>
#include <memory>


// Список номеров тел
typedef std::vector<int> NumbersOfBodies;


/**
* Заголовок файла результатов
*/
struct RezrHeader
{
	float simv, v, hz, tk;
	Integer nvar, kgr;
	float so;
};


/**
* Данные о положении тела
*/
struct RezrBodyPosition
{

	/// Массив координат x, y, z, rx, ry, rx;
	float Coordinates[6];

	/// Матрица поворота
	float Rotation[9];


	/**
	* Оператор сравнения двух положений тел
	* @param other - другое положение тела
	* @return признак, идентичны ли два положения тел
	*/
	bool operator ==
		(
			const RezrBodyPosition& other
		)	const;

};


/**
* Абстрактный провайдер доступа к данным результатов расчета MBS-решателя
*/
class AbstractProviderRezr
{
public:

	// Конструкторы/деструктор

	/**
	* Конструктор
	* @param rezrName - имя файла с результатами
	* @param fadresName - имя файла со структурой результатов
	*/
	AbstractProviderRezr
		(
			const string& rezrName,
			const string& fadresName
		);

	// запрещаем использование конструктора копирования
	AbstractProviderRezr(const AbstractProviderRezr&) = delete;
	
	/** Деструктор
	*/
	virtual
	~AbstractProviderRezr() = default;


	// Работа с файлом результатов

	/**
	* Открытие файла результатов
	* @return false в случае ошибки
	*/
	virtual
	bool OpenStream() = 0;

	/**
	* Закрыть поток результатов
	*/
	virtual
	void CloseStream() = 0;

	/**
	* Считывание заголовка
	*/
	virtual
	void ReadHeader() = 0;

	/**
	* Чтение кадра
	* @return false, если кадр не прочитан целиком или конец потока
	*/
	virtual
	bool TryReadFrame() = 0;

	/**
	* Перейти к кадру с данным номером
	* @param frameNumber - номер кадра
	* @return признак успешного (true) или неуспешного (false) перехода
	*/
	virtual
	bool GoToFrame
		(
			int frameNumber
		);


	// Селекторы

	/** Получить данные заголовка
	*/
	RezrHeader GetHeader() const;

	/** Получить по номеру тела данные о его положении
	* @param bodyNumber - номер тела
	* @return данные о положении тела в текущем кадре
	*/
	RezrBodyPosition GetRezrBodyPosition
		(
			size_t bodyNumber
		)	const;

	/** Изменить данные о положении тела по его номеру
	* @param bodyNumber - номер тела
	* @param bodyPosition - данные о положении тела в текущем кадре
	* @return признак того, что тело с данным номером было найдено
	*/
	bool GetRezrBodyPosition
		(
			size_t bodyNumber,
			RezrBodyPosition& rezrBodyPosition
		)	const;

	/** Получить список номеров тел
	* @return список номеров тел
	*/
	NumbersOfBodies GetNumbersOfBodies() const;

	/**
	* Получить текущее время (время, соответствующее текущему кадру)
	* @return текущее время
	*/
	float GetCurrentTime() const;

	/**
	* Получить номер текущего кадра (нумерация начианется с нуля)
	* @return номер текущего кадра
	*/
	int GetCurrentFrameNumber() const;

	/**
	* Получить списки смещений различных объектов
	* @return списки смещений различных объектов
	*/
	const OldAddresses& GetOldAddresses() const;


protected:

	// наименование файла с результатами решения
	string _rezrName;

	// списки смещений для различных объектов
	std::unique_ptr<OldAddresses> _oldAddresses;
	
	// буфер кадра (размерностью kgr)
	std::vector<float> _frameBuffer;

	// номер текущего кадра
	int _currentFrameNumber{};

	// количество кадров
	int _framesCount{};

	// текущее время
	float _currentTime{};

	// данные заголовка
	RezrHeader _header{};


	/** Вычислить количество кадров и выделить память
	*/
	void CalcFramesCount();

};


// умный указатель на провайдер доступа к данным результатов расчета MBS-решателя
using AbstractProviderRezrPointer = std::unique_ptr<AbstractProviderRezr>;


#endif // ABSTRACT_PROVIDER_REZR_H