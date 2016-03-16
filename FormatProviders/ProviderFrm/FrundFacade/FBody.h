#ifndef FBODY_H

#define FBODY_H


#include "FElement.h"
#include "FGeometry.h"
#include "FMultiphysics.h"
#include "FMph.h"
#include "../../../fcore/Calculator/Calculator.h"


struct FBodyDof
{

	// степени свободы поступательные 0-1 (для упругого тела номер узла)
	int x, y, z;
	// степени свободы вращательные 0,1,3 (3 - для больших поворотов)
	int rx, ry, rz;


	FBodyDof();

	void Init(int x, int y, int z, int rx, int ry, int rz);

	bool IsSmallMovementsMode() const;

	bool IsFixed() const;

};


struct FBodyInertia
{

	// инерционные параметры (значения)
	mutable double m, jx, jy, jz;

	// инерционные параметры (строковые выражения)
	string strm, strjx, strjy, strjz;

	// ограничения по минимальной массе и моменту инерции
	mutable double minMass, minInertia;


	FBodyInertia();

	void Eval(Calc::Calculator &calc) const;

	void ToStrings(string& strm_out, string& strjx_out, string& strjy_out, string& strjz_out) const;


private:

	static
	string ToStringByMinValue(const string& str, double value, double minValue);

};


struct FFlexibleBodyParams
{
	// тело - упругое
	int isFlexible;

	// номер формы упругого тела
	int flexFormNumber;

	// путь к файлу с упругими формами
	string flexFormPath;

	// узлы упругого тела, задающие плоскость, по которой определяется н.у. и матрица поворота
	int node1, node2, node3;

	// признак заморозки ???
	int freeze;


	FFlexibleBodyParams();

};


/**
* Класс для считывания тел модели ФРУНД
* номера тел (number) идут не подряд и не все, но ссылки именно по ним
* номер 49 - это дорога, нумерация с 1, номеров удаленных из модели тел нет
* после 500 тела нумерация идет с 5000
*
* @author Getmanskiy Victor
*/
class FBody : public FElement
{
	// Степени свободы
	FBodyDof _dof;

	// Инерционные параметры
	FBodyInertia _inertia;

	// Параметры упругого тела
	FFlexibleBodyParams _flexible;

	// каркасная геометрия (GEO)
	int _geometryId;
	//FGeometry* _geometry;
	// параметры вспомогательного расчета (MPH) old
	// FMultiphysics* _multiphysics;
	// параметры вспомогательного расчета (MPH)
	FMph* _mph;
public:

#pragma region Properties
	bool HasGeometry() const {return _geometryId >= 0;}

	int Id() const {return _id;}
	void Id(int val ){_id = val;}
	int GeometryId() const { return _geometryId; }
	void GeometryId(int val) { _geometryId = val; }
	double M() const {return _inertia.m;}
	void M(double val) {_inertia.m= val;}
	double Jx() const {return _inertia.jx;}
	void Jx(double val){ _inertia.jx = val;}
	double Jy() const {return _inertia.jy;}
	void Jy(double val) {_inertia.jy = val;}
	double Jz() const {return _inertia.jz;}
	void Jz(double val) { _inertia.jz = val; }
	int Z() const { return _dof.z; }
	void Z(int val) { _dof.z = val; }
	int Y() const { return _dof.y; }
	void Y(int val) { _dof.y = val; }
	int X() const { return _dof.x; }
	void X(int val) { _dof.x = val; }
	int Rx() const { return _dof.rx; }
	void Rx(int val) { _dof.rx = val; }
	int Ry() const { return _dof.ry; }
	void Ry(int val) { _dof.ry = val; }
	int Rz() const { return _dof.rz; }
	void Rz(int val) { _dof.rz = val; }
	string Strjz() const { return _inertia.strjz; }
	void Strjz(string val) { _inertia.strjz = val; }
	string Strjy() const { return _inertia.strjy; }
	void Strjy(string val) { _inertia.strjy = val; }
	string Strjx() const { return _inertia.strjx; }
	void Strjx(string val) { _inertia.strjx = val; }
	string Strm() const { return _inertia.strm; }
	void Strm(string val) { _inertia.strm = val; }
	int Freeze() const { return _flexible.freeze; }
	void Freeze(int val) { _flexible.freeze = val; }
	double MinInertia() const { return _inertia.minInertia; }
	void MinInertia(double val) const { _inertia.minInertia = val; }
	double MinMass() const { return _inertia.minMass; }
	void MinMass(double val) const { _inertia.minMass = val; }
	int FlexFormNumber() const { return _flexible.flexFormNumber; }
	void FlexFormNumber(int val) { _flexible.flexFormNumber = val; }
	string FlexFormPath() const { return _flexible.flexFormPath; }
	void FlexFormPath(string val) { _flexible.flexFormPath = val; }
	int Node3() const { return _flexible.node3; }
	void Node3(int val) { _flexible.node3 = val; }
	int Node2() const { return _flexible.node2; }
	void Node2(int val) { _flexible.node2 = val; }
	int Node1() const { return _flexible.node1; }
	void Node1(int val) { _flexible.node1 = val; }

	bool IsFlexible() const {return _flexible.isFlexible > 0;}
	void IsFlexible(int value) {_flexible.isFlexible = value;}

	#pragma endregion
	
	//#pragma
	FBody();
	FBody(const FBody &src);
	FBody(const char* path){Load(path);};
	FBody(ifstream& ifs){Load(ifs);};

	//FGeometry* GetGeometry() const {return _geometry;}
	//FGeometry* Geometry() const { return _geometry; }
	//void Geometry(FGeometry* val) { _geometry = val; }
	
	FMph* GetMph() {return this->_mph;}
	FMph* Mph() const { return _mph; }
	void Mph(FMph* val) { _mph = val; }

	/**
	* Признак малых движений
	*/
	bool IsSmallMovementsMode() const;
	
	/** Добавить рештель для тела
	* @param solverType - тип решателя
	* @return решатель
	*/
	SolverParamsBase* AddNewSolver(SolverTypes solverType);

// См. комментарии в базовом классе
#pragma region overriden
	void Load(const char* path);

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

	void SaveByIndex() const;
#pragma endregion

	/** Возвращает параметры тела
	* @params - параметры тела
	*/
	void GetBodyDof(FBodyDof& dof) const {dof = _dof;}

	/** Вывод LBA, требуется для старой логики работы AVMODEL
	*/
	void OutputLba() const;
	
	int GetIdToOutput() const;

	/** Вывод LBA в поток, требуется для старой логики работы AVMODEL
	* @param stream - поток
	*/
	void OutputLba(ofstream& stream) const;

	// Сохранение геометрии в .mdl формат (0 = failed)
	//int SaveLbm(AvmodelInterface& avmodel);
	bool operator==(const FBody &other) const;
	
	/** Заменить численные инерционные параметры по переменным в калькуляторе
	* @param calc - калькулятор
	*/
	void Eval(Calc::Calculator &calc) const;
};


#endif // FBODY_H