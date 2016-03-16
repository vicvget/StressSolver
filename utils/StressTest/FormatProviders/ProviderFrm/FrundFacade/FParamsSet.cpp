#include "FParamsSet.h"

#include "../../../fcore/Exceptions/fcExceptions.h"


FParamsSet::FParamsSet()
{
	currentId = 0;
}

FParamsSet::~FParamsSet()
{
}

/** Выбор используемого набора
* @param id - номер набора параметров
*/
void FParamsSet::Select(std::size_t id)
{
	if(id < Size())
		currentId = id;
	else
		exceptions::ThrowCollectionOutOfBounds("FParamsSet");
}
