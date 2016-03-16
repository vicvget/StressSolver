#ifndef FFORCE_H

#define FFORCE_H


#include "FElement.h"

//TODO: нормализовать комментарии

/**
* Класс для считывания нагрузок модели ФРУНД
*
* @author Getmanskiy Victor
*/
class FForce : public FElement
{
private:
	// номер тела 
	BodyNodeNumber _bodyNodeNumber;

	// координаты узла приложения
	FGeometryPoint _node;
	// номер характеристики (==number)
	int _charNumber;

public:

	FForce();

	FForce
		(
			const FForce& src
		);

	FForce
		(
			int orderIndex
		);

	FForce
		(
			const char* path
		);

	FForce
		(
			ifstream& ifs
		);

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


#pragma region Properties
	
	int BodyNumber() const { return _bodyNodeNumber.bodyNumber; }
	void BodyNumber(int val) { _bodyNodeNumber.bodyNumber = val; }

	int NodeNumber() const { return _bodyNodeNumber.nodeNumber; }
	void NodeNumber(int val) { _bodyNodeNumber.nodeNumber = val; }

	int CharNumber() const { return _charNumber; }
	void CharNumber(int val) { _charNumber = val; }

	FGeometryPoint Node() const { return _node; }
	void Node(FGeometryPoint val) { _node = val; }
#pragma  endregion

};


#endif // FFORCE_H