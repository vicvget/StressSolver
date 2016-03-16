 #ifndef FMODEL_BUILDER_H

#define FMODEL_BUILDER_H


#include "FForceCharBuilder.h"
#include "FGeometryBuilder.h"

#include "../../ProviderFrm/FrundFacade/FModel.h"
#include "../../../AdditionalModules/fmath/Vector3.h"
#include "../../../AdditionalModules/fmath/Matrix3x3.h"


using MathHelpers::Vec3;
using MathHelpers::Mat3;


class FModelBuilder
{
protected:
	FModel _model;
	int _charId;
	int _sphericalStiffId;
	int _cylindricalStiffId[3];
	int _inPlaneStiffId[3];
	int _groundBodyNumber;

	double GetDistance
		(
			const BodyNodeNumber& bodyNode1,
			const BodyNodeNumber& bodyNode2
		)	const;

	// TODO: создать узел в GCS
	int AddGroundNode(const BodyNodeNumber& bodyNode1);

public:
	FModelBuilder();
	const FModel& Get() const;
	
	static
	void CheckShiftDirection(int direction);


	int GetSphericalJointCharId();
	int GetCylindricalJointCharId(int direction);
	int GetInPlaneJointCharId(int direction);

	PGeometryBuilder GetGeometryBuilder(int geometryId) const;
	void SetGeometry(int geometryId, const FGeometry& geometry);

	// Изменение геометрии тела за счет добавления точки в глобальной системе координат (Global Coordinate System)
	// @param bodyNumber - номер тела
	// @param gcsPoint - точка в GCS
	// @returns номер добавленной точки
	int AddGCSPointToBodyGeometry(int bodyNumber, const Vec3& gcsPoint);

	// Изменение геометрии тела за счет добавления точки в локальной системе координат (Local Coordinate System)
	// @param bodyNumber - номер тела
	// @param lcsPoint - точка в LCS
	// @returns номер добавленной точки
	int AddLCSPointToBodyGeometry(int bodyNumber, const Vec3& lcsPoint);

	void Setup(const string& name);

	// Добавление гравитации
	int AddGravityForce(int direction, int bodyNumber);

	// Добавление зафиксированного тела (0 степеней свободы)
	int AddFixedBodyByBoundTransform(const string& name, Vec3 minPoint, Vec3 maxPoint, Vec3 cmPoint, Mat3 mtxTransform, double m, Vec3 inertia);

	// Добавление земли "GROUND" (0 степеней свободы)
	int AddGroundBody();

	// Добавление свободного тела (6 степеней свободы)
	int AddFreeBodyByBoundTransform(const string& name, Vec3 minPoint, Vec3 maxPoint, Vec3 cmPoint, Mat3 mtxTransform, double m, Vec3 inertia);

	//Добавление свободного тела (6 степеней свободы) по точкам из маркеров
	int AddFreeBodyByMarkersPoints(const string& name, vector<Vec3> markersPoints, Vec3 cmPoint, Mat3 mtxTransform, double mass, Vec3 inertia);

	// Добавление цилиндрического шарнира в СК первого тела (direction - направление вектора оси вращения)
	int AddRevolvedJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2);
	int AddRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2);

	// Добавление сферического шарнира в СК первого тела
	int AddSphericalJoint(const string& name, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2);
	int AddSphericalJoint(const string& name, BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2);

	// Добавление шарнира, ограничивающего 2 поступательные степени в СК первого тела (direction - направление вектора оси вращения)
	int AddCylindricalJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2);
	int AddCylindricalJoint(const string& name, int direction, BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2);

	// Добавление шарнира, ограничивающего 1 поступательную степень в СК первого тела (direction - направление ограничения)
	int AddInPlaneJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2);
	int AddInPlaneJoint(const string& name, int direction, BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2);

	// Добавление 2 шарниров для моделирования зафиксированного (2 поступательные степени) цилиндрического шарнира с заданным расстоянием вдоль оси
	int AddCylindricalRevolvedJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance);

	// Добавление 2 шарниров для моделирования зафиксированного (2 поступательные степени) цилиндрического шарнира с заданным расстоянием вдоль оси
	int AddCylindricalRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance);

	// Добавление 2 шарниров для моделирования зафиксированного цилиндрического шарнира с заданным расстоянием вдоль оси
	int AddRevolvedJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance);

	// Добавление 2 шарниров для моделирования зафиксированного цилиндрического шарнира с заданным расстоянием вдоль оси
	int AddRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance);


	// Добавление 2 шарниров для моделирования зафиксированного цилиндрического шарнира с заданным расстоянием вдоль оси
	int AddTranslationalJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance);

	// Добавление 2 шарниров для моделирования зафиксированного цилиндрического шарнира с заданным расстоянием вдоль оси
	int AddTranslationalJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance);

	// добавление 3 шарниров: сферического,цилиндрического и inplane
	int AddFixedJoint(const string& name, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance);

	// добавление 3 шарниров: сферического,цилиндрического и inplane
	int AddFixedJoint(const string& name, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance);

	// Добавление упруго-демпфирующей связи
	int AddSpringDamper(const string& name, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double stiffness, double damping);

#pragma region Single Reference Joints
	int AddRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1);

	// Добавление сферического шарнира в СК первого тела
	int AddSphericalJoint(const string& name, BodyNodeNumber& bodyNode1);

	// Добавление шарнира, ограничивающего 2 поступательные степени в СК первого тела (direction - направление вектора оси вращения)
	int AddCylindricalJoint(const string& name, int direction, BodyNodeNumber& bodyNode1);

	// Добавление шарнира, ограничивающего 1 поступательную степень в СК первого тела (direction - направление ограничения)
	int AddInPlaneJoint(const string& name, int direction, BodyNodeNumber& bodyNode1);

	// Добавление 2 шарниров для моделирования зафиксированного (-3 translational DOF) цилиндрического шарнира с заданным расстоянием вдоль оси
	int AddRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, double distance);

	// Добавление 2 шарниров для моделирования зафиксированного (-2 translational DOF) цилиндрического шарнира с заданным расстоянием вдоль оси
	int AddCylindricalRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, double distance);


	// Добавление 2 шарниров для моделирования зафиксированного цилиндрического шарнира с заданным расстоянием вдоль оси
	int AddTranslationalJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, double distance);
	
	//добавление 3 шарниров: сферического,цилиндрического и inplane
	int AddFixedJoint(const string& name, const BodyNodeNumber& bodyNode1, double distance);

	// Добавление упруго-демпфирующей связи
	int AddSpringDamper(const string& name, const BodyNodeNumber& bodyNode1, double stiffness, double damping);

	// Добавление упруго-демпфирующей связи
	int AddSpringType2Damper(const string& name, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double stiffness1, double stiffness2, double lInterval, double damping);

#pragma endregion

	// Добавление постоянной силы
	int AddConstantForce(const string& name, int direction, const BodyNodeNumber& bodyNode, double param);

	// Добавление постоянной силы
	int AddConstantForce(const string& name, const BodyNodeNumber& bodyNode, const Vec3& force);

	// Добавление постоянной силы
	int AddConstantForce(const string& name, const BodyNodeNumber& bodyNodeNumber1, const BodyNodeNumber& bodyNodeNumber2, double force);

	// Добавление набора именованного параметров расчета и установка текущим
	int AddSolverParamsSet(const string& name);

	// Добавление начальных условий из текущего положения тел в текущий набор параметров расчета
	int AddCurrentIco(const string& name);

	int GetGroundBody();

	int AddConstantForceSpring(const string& name, const BodyNodeNumber& bodyNodeNumber1, const BodyNodeNumber& bodyNodeNumber2, double force);
/* TODO:
	// Добавление постоянной силы
	void AddConstantForce(const string& name, int direction, int bodyNumber, int nodeNumber);

	// Добавление земли (тело без степеней свободы)
	void AddGroundBody(const string& name);

	// Добавление свободного тела (6 степеней свободы)
	void AddFreeBody(const string& name);

	// Добавление цилиндрического шарнира в СК первого тела (direction - направление вектора оси вращения)
	void AddRevolvedJoint(const string& name, int direction, int bodyNumber1, int bodyNumber2);

	// Добавление шарнира, ограничивающего перемещение в плоскости
	void AddPlanarJoint(const string& name, int direction, int bodyNumber1, int bodyNumber2);
	*/
};


#endif // FMODEL_BUILDER_H