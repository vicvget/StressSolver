#ifndef FCUNIT_H

#define FCUNIT_H


#include "FElement.h"
#include "FGeometry.h"

//TODO: нормализовать комментарии

/**
* Класс для считывания управляющих элементов модели ФРУНД
* cunit - это тела со спец. номерами с 501 - ...
*
* @author Getmanskiy Victor
*/
class FCunit : public FElement
{
	// экранные координаты на диаграмме
	int _x, _y;

public:

	FCunit();

	FCunit
		(
			int orderIndex
		);

	FCunit
		(
			const char* path
		)
		:
			_x(),
			_y()
	{
		FElement::Load(path);
	}

	FCunit
		(
			ifstream& ifs
		)
		:
			_x(),
			_y()
	{
		Load(ifs);
	}

	int X() const {return _x;}
	void X(int val) { _x = val; }

	int Y() const {return _y;}
	void Y(int val) { _y = val; }

// См. комментарии в базовом классе
#pragma region overriden

	virtual
	bool Load
		(
			ifstream& stream
		)	override;

	void Save(ofstream& ofs) const;

	virtual
	void Output
		(
			ofstream& stream,
			int stage
		)	const override;

#pragma endregion

};


#endif // FCUNIT_H