#ifndef BOUNDARYNORMAL_H
#define BOUNDARYNORMAL_H


#include "Vertex.h"


/*
Класс содержащий нормали для грничных элементов сетки, описанной 
объектом OccRectilinearGrid
*/
class BoundaryNormal
{
private:
	
	// номер точки, для которой указана нормаль
	unsigned long _pointNumber{};
	
	// единичный вектор, описывающий нормаль
	VertexPoint _unitaryNormal;

public:

	/**
	* Конструктор по умолчанию
	*/
	BoundaryNormal() = default;

	/**
	* Конструктор
	*/
	BoundaryNormal
		(
			unsigned long pointNumber,
			const VertexPoint& unitaryNormal
		);
	
	/*
	* @return номер точки, для которой указана нормаль 
	*/
	unsigned long GetPointNumber() const {return _pointNumber;}
	
	/*
	* Устанавливает номер точки, для которой указана нормаль 
	*/
	void SetPointNumber(unsigned long pointNumber) {_pointNumber = pointNumber;}

	/*
	* * @return единичный вектор, описывающий нормаль
	*/
	const VertexPoint& GetUnitaryNormal() const {return _unitaryNormal;}
	
	/*
	* Устанавливает единичный вектор, описывающий нормаль
	*/
	void SetUnitaryNormal(const VertexPoint& unitaryNormal){ _unitaryNormal = unitaryNormal;}
};

#endif