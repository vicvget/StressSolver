#ifndef OLD_ADDRESSES_H

#define OLD_ADDRESSES_H


#include "Types.h"
#include "MechanicalObjects.h"

#include <string>


using std::string;


const string FadresName = "fadres.dat";


/**
* Класс списков смещений различных объектов
*/
class OldAddresses
{
public:

	// Конструкторы и деструктор

	/**
	* Конструктор по умолчанию
	*/
	OldAddresses();


	// Чтение и запись

	/**
	* Заполнение списков из файла fadres.dat
	* @param fadresName - имя файла fadres (возможно, с путем к нему)
	*/
	void ReadFadres
		(
			const string& fadresName = FadresName
		);


	// Селекторы

	/**
	* Получить значение числа степеней свободы для тел
	* @return ksts - число степеней свободы для тел
	*/
	int GetKsts() const;

	/**
	* Получить значение общего числа степеней свободы
	* @return kgr - общее число степеней свободы
	*/
	int GetKgr() const;

	/**
	* Получение списка смещений для степеней свободы тел
	* @return список смещений для степеней свободы тел
	*/
	const Bodies& GetBodies() const;

	/**
	* Получение списка смещений для соединительных элементов
	* @return список смещений для соединительных элементов
	*/
	const ConnectionElements& GetConnectionElements() const;

	/**
	* Получение списка смещений для матриц поворота
	* @return список смещений для матриц поворота
	*/
	const RotationMatrixes& GetRotationMatrixes() const;

	/**
	* Получение списка кэшированных степеней свободы для тел
	* @return кэш смещений для тел
	*/
	const BodiesDofOffsetsCache& GetBodiesDofOffsetsCache() const;

	/**
	* Проверить, есть ли данная степень свободы у тела с данным номером
	* @param bodyNumber - номер тела
	* @param freedom - тип степени свободы
	* @return признак, есть ли данная степень свободы у тела с данным номером (true) или ее нет (false)
	*/
	bool CheckBodyFreedom
		(
			Integer bodyNumber,
			FreedomType freedom
		)	const;

	/**
	* Проверить, есть ли данная степень свободы у тела с данным номером
	* @param bodyNumber - номер тела
	* @param motion - тип движения
	* @param direction - тип направления движения
	* @return признак, есть ли данная степень свободы у тела с данным номером (true) или ее нет (false)
	*/
	bool CheckBodyFreedom
		(
			Integer bodyNumber,
			MotionType motion,
			DirectionType direction
		)	const;

	/**
	* Проверить, есть ли данная степень свободы у соединительного элемента с данным номером
	* @param connectionElementNumber - номер соединительного элемента
	* @param freedom - тип степени свободы
	* @return признак, есть ли данная степень свободы у соединительного элемента с данным номером (true)
	* или ее нет (false)
	*/
	bool CheckConnectionElementFreedom
		(
			Integer connectionElementNumber,
			FreedomType freedom
		)	const;

	/**
	* Проверить, есть ли данная степень свободы у соединительного элемента с данным номером
	* @param connectionElementNumber - номер соединительного элемента
	* @param motion - тип движения
	* @param direction - тип направления движения
	* @return признак, есть ли данная степень свободы у соединительного элемента с данным номером (true)
	* или ее нет (false)
	*/
	bool CheckConnectionElementFreedom
		(
			Integer connectionElementNumber,
			MotionType motion,
			DirectionType direction
		)	const;

	/**
	* Получение точного смещения для тела
	* @param bodyNumber - номер тела
	* @param freedom - тип степени свободы
	* @param value - величина
	* @return смещение
	*/
	Integer GetBodyAddress
		(
			Integer bodyNumber,
			FreedomType freedom,
			Value value
		)	const;

	/**
	* Получение точного смещения для тела
	* @param bodyNumber - номер тела
	* @param motion - тип движения
	* @param direction - тип направления движения
	* @param value - величина
	* @return смещение
	*/
	Integer GetBodyAddress
		(
			Integer bodyNumber,
			MotionType motion,
			DirectionType direction,
			Value value
		)	const;

	/**
	* Получение точного смещения для соединительного элемента
	* @param connectionElementNumber - номер соединительного элемента
	* @param freedom - тип степени свободы
	* @param component - компонента характеристики силы
	* @return смещение
	*/
	Integer GetConnectionElementAddress
		(
			Integer connectionElementNumber,
			FreedomType freedom,
			CharacteristicComponent component
		)	const;

	/**
	* Получение точного смещения для соединительного элемента
	* @param connectionElementNumber - номер соединительного элемента
	* @param motion - тип движения
	* @param direction - тип направления движения
	* @param component - компонента характеристики силы
	* @return смещение
	*/
	Integer GetConnectionElementAddress
		(
			Integer connectionElementNumber,
			MotionType motion,
			DirectionType direction,
			CharacteristicComponent component
		)	const;

	/**
	* Получение точного значения смещения для матрицы поворота
	* @param rotationMatrixNumber - номер матрицы поворота
	* @return смещение
	*/
	Integer GetRotationMatrixAddress
		(
			Integer rotationMatrixNumber
		)	const;


protected:

	//поля

	//число степеней свободы для тел
	Integer _ksts;

	//общее число степеней свободы
	Integer _kgr;

	//список тел
	Bodies _bodies;

	//список соединительных элементов
	ConnectionElements _connectionElements;

	//список матриц поворота
	RotationMatrixes _rotationMatrixes;

	// список кэшированных степеней свободы для тел
	BodiesDofOffsetsCache _bodiesDofOffsetsCache;


	/**
	* Заполнение кэша степеней свободы для тел
	*/
	void BuildBodiesDofOffsetsCache();

};


#endif // OLD_ADDRESSES_H