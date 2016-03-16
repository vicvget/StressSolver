#ifndef BASIC_TYPES_H

#define BASIC_TYPES_H


//#include "../../../fcore/fcore.h"
//#include "../../../fcore/Exceptions/fcExceptions.h"
//#include "../../../fcore/wrappers/FileRoutines.h"
#include "FGeometryPoint.h"

#ifdef GNUCPP
#include "../../../fcore/wrappers/fcUnistd.h"
#else

#include <cstring>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>

#endif

// максимальное ограничение по максимальной жесткости податливой кинематической пары
const double MaxStiffness = 1e100;

typedef char StaticCharString[1024];

enum class Dof
{
	x = 1,
	y,
	z,
	rx,
	ry,
	rz
};

struct StaticCharStringStruct
{
	StaticCharString str;
};

struct FDebugParams
{
	double minMass{};
	double minInertia{};
	double maxStiffness{MaxStiffness};


	bool Load();

	bool Save() const;


private:

	bool isLoaded{};

};

/** Структура для хранения параметра 
*/
struct Parameter
{
	// имя параметра
	string name;
	//string name;
	// вычислимое выражение выражение (значение параметра)
	string expression;
	//string expression;
	// описание параметра
	string description;
	//string description;
};

/** Связь в каркасной геометрии
*/
struct FGeometryLink
{

	// номер первой точки
	int node1;

	// номер второй точки
	int node2;

	FGeometryLink();

	FGeometryLink
		(
			int n1,
			int n2
		);

};

struct FGeometryCircle
{
	int code; //3 или 5???
	int centerNode;
	int axisNode;
	int nPoints;
	float Radius;
	float angle1;
	float angle2;
};

/** Компонент характеристики обобщенной силы
*/
struct CharComponent
{
	// TODO: жесткость
	mutable double stiffness;
	// тип характеристики
	int type;
	// направление (3 поступательных и 3 вращательных)
	int direction; // 1-6
	// набор параметров для заданного типа характеристики
	vector<string> params;

	CharComponent();
};

/**
* Для хранения номера тела и номера узла
*/
struct BodyNodeNumber
{
	int bodyNumber;
	int nodeNumber;

	BodyNodeNumber();
	BodyNodeNumber(int bodyNumber, int nodeNumber);
};

bool operator == (const BodyNodeNumber& bodyNodeNumber1, const BodyNodeNumber& bodyNodeNumber2);


#endif // BASIC_TYPES_H