#ifndef FJOINT_H

#define FJOINT_H


#include "FElement.h"
#include "../../../fcore/Calculator/Calculator.h"

//TODO: нормализовать комментарии

/**
* Класс для считывания шарнира модели ФРУНД
*
* @author Getmanskiy Victor
*/
class FJoint : public FElement
{
private:
	// номер первого тела и узла крепления
	BodyNodeNumber _bodyNodeNumber1;
	// номер второго тела и узла крепления
	BodyNodeNumber _bodyNodeNumber2;
	// координаты узлов крепления // TODO удалить
	FGeometryPoint node1, node2;
	vector<int> _excludedDirections;

	// номер характеристики
	int charNumber;

	// 0 - обычный шарнир, 1 - пружина, линия действия силы по линии узлов, 
	// -1 - только на первое тело
	int _springType;
	// длина пружины (строковые выражения)
	string _sspringLength;
	// длина пружины
	mutable double _springLength;
	// & ... во второй строке
	string _tag;
	// id для сортировки при выводе
	mutable int _sortId;

public:

	FJoint();

	FJoint
		(
			const FJoint &src
		);

	FJoint
		(
			int orderIndex
		);
	
	//FJoint(const char* path);
	//FJoint(ifstream& stream);

	FJoint
		(
			const char* path
		)
		:
			charNumber(),
			_springType(),
			_springLength(),
			_sortId()
	{
		FElement::Load(path);
	}

	FJoint
		(
			ifstream& ifs
		)
		:
			charNumber(),
			_springType(),
			_springLength(),
			_sortId()
	{
		Load(ifs);
	}

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

#pragma  region GettersSetters
	int SortId() const { return _sortId; }
	void SortId(int val) const { _sortId = val; }
	FGeometryPoint Node2() const { return node2; }
	void Node2(FGeometryPoint val) { node2 = val; }
	FGeometryPoint Node1() const { return node1; }
	void Node1(FGeometryPoint val) { node1 = val; }
	int CharNumber() const { return charNumber; }
	void CharNumber(int val) { charNumber = val; }
	string GetTag() const {return _tag;};
	void SetTag(string value) {_tag = value;}
	int NodeNumber1() const { return _bodyNodeNumber1.nodeNumber;}
	void NodeNumber1(int value) {_bodyNodeNumber1.nodeNumber = value;}
	int NodeNumber2() const	{ return _bodyNodeNumber2.nodeNumber;}
	void NodeNumber2(int value) {_bodyNodeNumber2.nodeNumber = value;}
	int BodyNumber1() const { return _bodyNodeNumber1.bodyNumber;}
	void BodyNumber1(int value) {_bodyNodeNumber1.bodyNumber = value;}
	int BodyNumber2() const	{ return _bodyNodeNumber2.bodyNumber;}
	void BodyNumber2(int value) {_bodyNodeNumber2.bodyNumber = value;}
	int SpringType() const { return _springType; }
	void SpringType(int val) { _springType = val; }
	string SspringLength() const { return _sspringLength; }
	void SspringLength(string val) { _sspringLength = val; }
	double SpringLength() const { return _springLength; }
	void SpringLength(double val) const { _springLength = val; }
	vector<int> ExcludedDirections() const { return _excludedDirections; }
	void ExcludedDirections(vector<int> val) { _excludedDirections = val; }
	const BodyNodeNumber& BodyNodeNumber1() const { return _bodyNodeNumber1; }
	void BodyNodeNumber1(const BodyNodeNumber& val) { _bodyNodeNumber1 = val; }
	const BodyNodeNumber& BodyNodeNumber2() const { return _bodyNodeNumber2; }
	void BodyNodeNumber2(const BodyNodeNumber& val) { _bodyNodeNumber2 = val; }

#pragma endregion

	/** Заменить численные параметры по переменным в калькуляторе
	* @param calc - калькулятор
	*/
	void Eval(Calc::Calculator &calc) const;
};


#endif // FJOINT_H