#ifndef FPARAMS_SET_H

#define FPARAMS_SET_H


#include <cstddef>


/** Интерфейсный класс для общих методов в наборах параметров
* @author Victor Getmasnkiy
*/
class FParamsSet
{
protected:
	std::size_t currentId;
public:
	FParamsSet();/*

	size_t CurrentId() const { return currentId; }
	void CurrentId(size_t val) { currentId = val; }*/

	virtual ~FParamsSet();
	
	/** Получение количества наборов параметров
	* return количество наборов параметров
	*/
	virtual std::size_t Size() const = 0;

	/** Выбор используемого набора
	* @param id - номер набора параметров
	*/
	void Select(std::size_t id);

	std::size_t CurrentId() const { return currentId; }
};


#endif // FPARAMS_SET_H