#include "FModelBuilder.h"

#include "FBodyBuilder.h"
#include "FGeometryBuilder.h"
#include "FForceBuilder.h"
#include "FForceCharBuilder.h"
#include "FJointBuilder.h"
#include "FJointCharBuilder.h"

#include "../../../fcore/wrappers/StringRoutines.h"

#include <cmath>


FModelBuilder::FModelBuilder(): 
	_charId(0),
	_sphericalStiffId(0),
	_groundBodyNumber(0)
{
	std::fill(_cylindricalStiffId, _cylindricalStiffId+3, 0);
	std::fill(_inPlaneStiffId, _inPlaneStiffId+3, 0);
}

const FModel& FModelBuilder::Get() const {return _model;}

void FModelBuilder::CheckShiftDirection(int direction)
{
	if(direction > 3 || direction < 1)
		exceptions::ThrowMessage("Shift direction out of range (1-3)");
}

PGeometryBuilder FModelBuilder::GetGeometryBuilder(int geometryId) const
{
	return PGeometryBuilder(new FGeometryBuilder(_model.GetGeometry(geometryId)));
}

void FModelBuilder::SetGeometry(int geometryId, const FGeometry& geometry)
{
	_model.GetGeometry(geometryId) = geometry;
}

double FModelBuilder::GetDistance
	(
		const BodyNodeNumber& bodyNode1,
		const BodyNodeNumber& bodyNode2
	)	const
{
	const FGeometryPoint& point1 = _model.GetTransformedNode(bodyNode1);
	const FGeometryPoint& point2 = _model.GetTransformedNode(bodyNode2);

	return (point1 - point2).Magnitude() / 1000.0; //to meters
}

int FModelBuilder::GetSphericalJointCharId()
{
	const string sphericalCharName("SphericalChar");
	if(!_sphericalStiffId)
	{
		_sphericalStiffId = _model.GetJointCharsCount(); // 0-based
		FJointSphericalCharBuilder sphericalCharBuilder;
		sphericalCharBuilder.Setup(sphericalCharName, _sphericalStiffId);
		_sphericalStiffId++; // store 1- based
		_model.AddJointChar(sphericalCharBuilder.Get());
	}
	return _sphericalStiffId;
}

int FModelBuilder::GetCylindricalJointCharId(int direction)
{
	CheckShiftDirection(direction);
	int directionID = direction-1;

	const string XYZ = "XYZ";
	const string charName = string("CylindricalChar") + XYZ[directionID];
	if(!_cylindricalStiffId[directionID])
	{
		_cylindricalStiffId[directionID] = _model.GetJointCharsCount(); // 0-based
		FJointCylindricalCharBuilder charBuilder;
		charBuilder.Setup(direction, charName, _cylindricalStiffId[directionID]);
		_cylindricalStiffId[directionID]++;
		_model.AddJointChar(charBuilder.Get());
	}
	return _cylindricalStiffId[directionID];
}

int FModelBuilder::GetInPlaneJointCharId(int direction)
{
	CheckShiftDirection(direction);
	int directionID = direction-1;

	const string XYZ = "XYZ";
	const string charName = string("InPLaneChar") + XYZ[directionID];
	if(!_inPlaneStiffId[directionID])
	{
		_inPlaneStiffId[directionID] = _model.GetJointCharsCount(); // 0-based
		FJointInPlaneCharBuilder charBuilder;
		charBuilder.Setup(direction, charName, _inPlaneStiffId[directionID]);
		_inPlaneStiffId[directionID]++;
		_model.AddJointChar(charBuilder.Get());
	}
	return _inPlaneStiffId[directionID];
}

void FModelBuilder::Setup(const string& name) 
{
	_model.Name(name);
}

int FModelBuilder::AddGravityForce(int direction, int bodyNumber)
{
	CheckShiftDirection(direction);
	const string gravityCharName("GravityChar");
	const string gravityName("Gravity");

	//if(!_gravityCharId)
	{
		_charId = _model.GetForceCharsCount(); // 0-based
		FForceGravityCharBuilder gravityCharBuilder;
		gravityCharBuilder.Setup(direction, gravityCharName, _charId);
		_model.AddForceChar(gravityCharBuilder.Get());
	}
	int gravityId = _model.GetForcesCount(); // 0-based
	FForceBuilder forceGravityBuilder;
	auto body = _model.GetBodyByNumber(bodyNumber);
	int nodeNumber = _model.GetBodyGeometry(body.Id()).CmNodeNumber() + 1; // 1-based
	forceGravityBuilder.Setup(gravityName, gravityId, bodyNumber, nodeNumber, _charId+1);
	return _model.AddForce(forceGravityBuilder.Get());
}

// Добавление свободного тела (6 степеней свободы)
int FModelBuilder::AddFixedBodyByBoundTransform(const string& name, Vec3 minPoint, Vec3 maxPoint, Vec3 cmPoint, Mat3 mtxTransform, double mass, Vec3 inertia)
{
	int bodyId = _model.GetBodiesCount(); // 0-based
	FFixedBodyBuilder bodyBuilder;
	FGeometryBoundBuilder geometryBuilder;
	int geometryId = _model.GetGeometriesCount();
	bodyBuilder.Setup(name, bodyId, mass, inertia, geometryId);
	geometryBuilder.SetupTransformedBound(name, bodyId, bodyBuilder.Get().Number(), minPoint, maxPoint, mtxTransform, cmPoint);
	_model.AddGeometry(geometryBuilder.Get());
	return _model.AddBody(bodyBuilder.Get());
}

// Добавление свободного тела (6 степеней свободы)
int FModelBuilder::AddGroundBody()
{
	int bodyId = _model.GetBodiesCount(); // 0-based
	FFixedBodyBuilder bodyBuilder;
	FGeometryBoundBuilder geometryBuilder;
	Vec3 minPoint, maxPoint, cmPoint, inertia;
	minPoint.Init(-0.005, -0.005, -0.005);
	maxPoint.Init(0.005, 0.005, 0.005);
	cmPoint.Init(0., 0., 0.);
	inertia.Init(1., 1., 1.);
	double data[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
	double mass = 1.;
	const string name = "GROUND";
	Mat3 mtxTransform(data);
	int geometryId = _model.GetGeometriesCount();
	
	bodyBuilder.Setup(name, bodyId, mass, inertia, geometryId);
	geometryBuilder.SetupTransformedBound(name, bodyId, bodyBuilder.Get().Number(), minPoint, maxPoint, mtxTransform, cmPoint);
	_model.AddGeometry(geometryBuilder.Get());
	_groundBodyNumber = bodyBuilder.Get().Number();
	return _model.AddBody(bodyBuilder.Get());
}

// Добавление свободного тела (6 степеней свободы)
int FModelBuilder::AddFreeBodyByBoundTransform(const string& name, Vec3 minPoint, Vec3 maxPoint, Vec3 cmPoint, Mat3 mtxTransform, double mass, Vec3 inertia)
{
	int bodyId = _model.GetBodiesCount(); // 0-based
	FFreeBodyBuilder bodyBuilder;
	FGeometryBoundBuilder geometryBuilder;
	int geometryId = _model.GetGeometriesCount();
	bodyBuilder.Setup(name, bodyId, mass, inertia, geometryId);
	
	geometryBuilder.SetupTransformedBound(name, bodyId, bodyBuilder.Get().Number(), minPoint, maxPoint, mtxTransform, cmPoint);

	_model.AddGeometry(geometryBuilder.Get());
	return _model.AddBody(bodyBuilder.Get());
}

int FModelBuilder::AddFreeBodyByMarkersPoints(const string& name, vector<Vec3> markersPoints, Vec3 cmPoint, Mat3 mtxTransform, double mass, Vec3 inertia)
{
	int bodyId = _model.GetBodiesCount(); // 0-based
	FFreeBodyBuilder bodyBuilder;
	FGeometryBoundBuilder geometryBuilder;
	int geometryId = _model.GetGeometriesCount();
	bodyBuilder.Setup(name, bodyId, mass, inertia, geometryId);
	/*if (bodyId == 48)
		int debugvar = 0;*/
	geometryBuilder.SetupMarkersBound(name, bodyId, bodyBuilder.Get().Number(), markersPoints, mtxTransform, cmPoint);
	_model.AddGeometry(geometryBuilder.Get());
	return _model.AddBody(bodyBuilder.Get());
}

// Добавление сферического шарнира в СК первого тела
int FModelBuilder::AddSphericalJoint(const string& name, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2)
{
	int jointId = _model.GetJointsCount(); // 0-based
	FJointBuilder jointBuilder;
	jointBuilder.Setup(name, jointId,  bodyNumber1, nodeNumber1, bodyNumber2, nodeNumber2, GetSphericalJointCharId());
	return _model.AddJoint(jointBuilder.Get());
}

int FModelBuilder::AddSphericalJoint(const string& name, BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2)
{
	int jointId = _model.GetJointsCount(); // 0-based
	FJointBuilder jointBuilder;
	jointBuilder.Setup(name, jointId,  bodyNode1, bodyNode2, GetSphericalJointCharId());
	return _model.AddJoint(jointBuilder.Get());
}

int FModelBuilder::AddSphericalJoint(const string& name, BodyNodeNumber& bodyNode1)
{
	BodyNodeNumber groundNode(GetGroundBody(), AddGroundNode(bodyNode1));
	return AddSphericalJoint(name, bodyNode1, groundNode);
}

// Добавление шарнира, ограничивающего 2 поступательные степени в СК первого тела (direction - направление вектора оси вращения)
int FModelBuilder::AddCylindricalJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2)
{
	int jointId = _model.GetJointsCount(); // 0-based
	FJointBuilder jointBuilder;
	jointBuilder.Setup(name, jointId,  bodyNumber1, nodeNumber1, bodyNumber2, nodeNumber2, GetCylindricalJointCharId(direction));
	return _model.AddJoint(jointBuilder.Get());
}

int FModelBuilder::AddCylindricalJoint(const string& name, int direction, BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2)
{
	int jointId = _model.GetJointsCount(); // 0-based
	FJointBuilder jointBuilder;
	jointBuilder.Setup(name, jointId,  bodyNode1, bodyNode2, GetCylindricalJointCharId(direction));
	return _model.AddJoint(jointBuilder.Get());
}

int FModelBuilder::AddCylindricalJoint(const string& name, int direction, BodyNodeNumber& bodyNode1)
{
	BodyNodeNumber groundNode(GetGroundBody(), AddGroundNode(bodyNode1));
	return AddCylindricalJoint(name, direction, bodyNode1, groundNode);
}

// Добавление шарнира, ограничивающего 1 поступательную степень в СК первого тела (direction - направление ограничения)
int FModelBuilder::AddInPlaneJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2)
{
	int jointId = _model.GetJointsCount(); // 0-based
	FJointBuilder jointBuilder;
	jointBuilder.Setup(name, jointId,  bodyNumber1, nodeNumber1, bodyNumber2, nodeNumber2, GetInPlaneJointCharId(direction));
	return _model.AddJoint(jointBuilder.Get());
}

int FModelBuilder::AddInPlaneJoint(const string& name, int direction, BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2)
{
	int jointId = _model.GetJointsCount(); // 0-based
	FJointBuilder jointBuilder;
	jointBuilder.Setup(name, jointId,  bodyNode1, bodyNode2, GetInPlaneJointCharId(direction));
	return _model.AddJoint(jointBuilder.Get());
}

int FModelBuilder::AddInPlaneJoint(const string& name, int direction, BodyNodeNumber& bodyNode1)
{
	BodyNodeNumber groundNode(GetGroundBody(), AddGroundNode(bodyNode1));
	return AddInPlaneJoint(name, direction, bodyNode1, groundNode);
}
// Добавление 2 шарниров для моделирования зафиксированного цилиндрического шарнира с заданным расстоянием вдоль оси
int FModelBuilder::AddTranslationalJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance)
{
	return AddTranslationalJoint(name, direction, bodyNode1.bodyNumber, bodyNode1.nodeNumber, bodyNode2.bodyNumber, bodyNode2.nodeNumber, distance);
}


// Добавление 2 шарниров для моделирования зафиксированного цилиндрического шарнира с заданным расстоянием вдоль оси
int FModelBuilder::AddTranslationalJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance)
{
	FBody& body1 = _model.GetBodyByNumber(bodyNumber1);
	FBody& body2 = _model.GetBodyByNumber(bodyNumber2);
	// TODO: логика создания 2 шарниров, требуется модификация геометрии тел

	const FGeometry& geom1 = _model.GetBodyGeometry(body1.Id());
	const FGeometry& geom2 = _model.GetBodyGeometry(body2.Id());

	
	FTransform transform1 = geom1.GetTransform();	
	Vec3 axis;// = transform1..GetDirectionAxis(direction-1);
	axis[direction-1]=1.;
	Vec3 point11 = geom1.GetGeometryPoint(nodeNumber1-1);
	
	Vec3 point12 = point11+axis*distance;

	FTransform transform2 = geom2.GetTransform();
	// TODO: переменная point21 не используется
	// Vec3 point21 = geom2.GetGeometryPoint(nodeNumber2-1);
	Vec3 point22 = transform1.ToLCS(point12, transform2); // перевод в lcs второго тела

	// вставка новых точек
	int nodeNumber12 = AddLCSPointToBodyGeometry(bodyNumber1, point12);
	int nodeNumber22 = AddLCSPointToBodyGeometry(bodyNumber2, point22);

	// вставка шарниров
	AddCylindricalJoint(name+"Cylindrical1",direction, bodyNumber1, nodeNumber1, bodyNumber2, nodeNumber2);
	return AddCylindricalJoint(name+"Cylindrical2",direction,bodyNumber1, nodeNumber12, bodyNumber2, nodeNumber22);
}

int FModelBuilder::AddTranslationalJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, double distance)
{
	BodyNodeNumber groundNode(GetGroundBody(), AddGroundNode(bodyNode1));
	return AddTranslationalJoint(name, direction, bodyNode1, groundNode, distance);
}

// Добавление 2 шарниров для моделирования зафиксированного цилиндрического шарнира с заданным расстоянием вдоль оси
int FModelBuilder::AddRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance)
{
	return AddRevolvedJoint(name, direction, bodyNode1.bodyNumber, bodyNode1.nodeNumber, bodyNode2.bodyNumber, bodyNode2.nodeNumber, distance);
}


// Добавление 2 шарниров для моделирования зафиксированного цилиндрического шарнира с заданным расстоянием вдоль оси
int FModelBuilder::AddRevolvedJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance)
{
	FBody& body1 = _model.GetBodyByNumber(bodyNumber1);
	FBody& body2 = _model.GetBodyByNumber(bodyNumber2);
	// TODO: логика создания 2 шарниров, требуется модификация геометрии тел

	const FGeometry& geom1 = _model.GetBodyGeometry(body1.Id());
	const FGeometry& geom2 = _model.GetBodyGeometry(body2.Id());

	
	FTransform transform1 = geom1.GetTransform();	
	Vec3 axis;// = transform1.GetDirectionAxis(direction-1);
	axis[direction-1]=1.;
	Vec3 point11 = geom1.GetGeometryPoint(nodeNumber1-1);
	Vec3 point12 = point11+axis*distance;

	FTransform transform2 = geom2.GetTransform();
	// TODO: переменная point21 не используется
	// Vec3 point21 = geom2.GetGeometryPoint(nodeNumber2-1);
	Vec3 point22 = transform1.ToLCS(point12, transform2); // перевод в lcs второго тела

	// вставка новых точек
	int nodeNumber12 = AddLCSPointToBodyGeometry(bodyNumber1, point12);
	int nodeNumber22 = AddLCSPointToBodyGeometry(bodyNumber2, point22);

	// вставка шарниров
	AddSphericalJoint(name+"Spherical",bodyNumber1, nodeNumber1, bodyNumber2, nodeNumber2);
	return AddCylindricalJoint(name+"Cylindrical",direction,bodyNumber1, nodeNumber12, bodyNumber2, nodeNumber22);
}

// Добавление 2 шарниров для моделирования зафиксированного цилиндрического шарнира с заданным расстоянием вдоль оси
int FModelBuilder::AddCylindricalRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance)
{
	return AddCylindricalRevolvedJoint(name, direction, bodyNode1.bodyNumber, bodyNode1.nodeNumber, bodyNode2.bodyNumber, bodyNode2.nodeNumber, distance);
}

// Добавление 2 шарниров для моделирования зафиксированного (ограничены 2 поступательные степени) цилиндрического шарнира с заданным расстоянием вдоль оси
int FModelBuilder::AddCylindricalRevolvedJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance)
{
	FBody& body1 = _model.GetBodyByNumber(bodyNumber1);
	FBody& body2 = _model.GetBodyByNumber(bodyNumber2);

	// TODO: убрать или переписать закомментированный код
	// if (bodyNumber2 == 50)
	//	int debugvar = 1;

	// TODO: логика создания 2 шарниров, требуется модификация геометрии тел

	const FGeometry& geom1 = _model.GetBodyGeometry(body1.Id());
	const FGeometry& geom2 = _model.GetBodyGeometry(body2.Id());


	FTransform transform1 = geom1.GetTransform();	
	Vec3 axis;// = transform1.GetDirectionAxis(direction-1);
	axis[direction-1]=1.;
	Vec3 point11 = geom1.GetGeometryPoint(nodeNumber1-1);
	Vec3 point12 = point11+axis*distance;

	FTransform transform2 = geom2.GetTransform();
	// TODO: неиспользуемая переменная point21
	// Vec3 point21 = geom2.GetGeometryPoint(nodeNumber2-1);
	Vec3 point22 = transform1.ToLCS(point12, transform2); // перевод в lcs второго тела

	// вставка новых точек
	int nodeNumber12 = AddLCSPointToBodyGeometry(bodyNumber1, point12);
	int nodeNumber22 = AddLCSPointToBodyGeometry(bodyNumber2, point22);

	// вставка шарниров
	AddCylindricalJoint(name+"CylindricalMain",direction,bodyNumber1, nodeNumber1, bodyNumber2, nodeNumber2);
	return AddCylindricalJoint(name+"CylindricalSub",direction,bodyNumber1, nodeNumber12, bodyNumber2, nodeNumber22);
}

int FModelBuilder::AddCylindricalRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, double distance)
{
	BodyNodeNumber groundNode(GetGroundBody(), AddGroundNode(bodyNode1));
	return AddCylindricalRevolvedJoint(name, direction, bodyNode1, groundNode, distance);
}

int FModelBuilder::AddRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, double distance)
{
	BodyNodeNumber groundNode(GetGroundBody(), AddGroundNode(bodyNode1));
	return AddRevolvedJoint(name, direction, bodyNode1, groundNode, distance);
}

// Добавление 3 шарниров 
int FModelBuilder::AddFixedJoint(const string& name, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance)
{
	return AddFixedJoint(name, bodyNode1.bodyNumber, bodyNode1.nodeNumber, bodyNode2.bodyNumber, bodyNode2.nodeNumber, distance);
}

// Добавление 3 шарниров 
int FModelBuilder::AddFixedJoint(const string& name, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance)
{
	int directionX = 1; 
	int directionY = 2;
	int directionZ = 3;

	FBody& body1 = _model.GetBodyByNumber(bodyNumber1);
	FBody& body2 = _model.GetBodyByNumber(bodyNumber2);
	// TODO: логика создания 2 шарниров, требуется модификация геометрии тел

	const FGeometry& geom1 = _model.GetBodyGeometry(body1.Id());
	const FGeometry& geom2 = _model.GetBodyGeometry(body2.Id());

	
	FTransform transform1 = geom1.GetTransform();	
	Vec3 axis;// = transform1.GetDirectionAxis(direction-1);
	axis[directionY-1]=1.;
	Vec3 point11 = geom1.GetGeometryPoint(nodeNumber1-1);
	Vec3 point12 = point11+axis*distance;

	FTransform transform2 = geom2.GetTransform();
	// TODO: неиспользуемая переменная
	// Vec3 point21 = geom2.GetGeometryPoint(nodeNumber2-1);
	Vec3 point22 = transform1.ToLCS(point12, transform2); // перевод в lcs второго тела

	// вставка новых точек
	int nodeNumber12 = AddLCSPointToBodyGeometry(bodyNumber1, point12);
	int nodeNumber22 = AddLCSPointToBodyGeometry(bodyNumber2, point22);

	FTransform transform13 = geom1.GetTransform();	
	Vec3 axis3;// = transform1.GetDirectionAxis(direction-1);
	axis3[directionX-1]=1.;
	Vec3 point_11 = geom1.GetGeometryPoint(nodeNumber1-1);
	Vec3 point13 = point_11+axis3*distance;

	FTransform transform23 = geom2.GetTransform();
	// TODO: неиспользуемая переменная
	// Vec3 point31 = geom2.GetGeometryPoint(nodeNumber2-1);
	Vec3 point32 = transform13.ToLCS(point13, transform23); // перевод в lcs второго тела

	// вставка новых точек
	int nodeNumber13 = AddLCSPointToBodyGeometry(bodyNumber1, point13);
	int nodeNumber32 = AddLCSPointToBodyGeometry(bodyNumber2, point32);
	// вставка шарниров
	AddSphericalJoint(name+"Spherical", bodyNumber1, nodeNumber1, bodyNumber2, nodeNumber2);
	AddCylindricalJoint(name+"Cylindrical", directionY, bodyNumber1, nodeNumber12, bodyNumber2, nodeNumber22);
	return AddInPlaneJoint(name+"Inplane", directionZ, bodyNumber1, nodeNumber13, bodyNumber2, nodeNumber32);
}

int FModelBuilder::AddFixedJoint(const string& name, const BodyNodeNumber& bodyNode1, double distance)
{
	BodyNodeNumber groundNode(GetGroundBody(), AddGroundNode(bodyNode1));
	return AddFixedJoint(name, bodyNode1, groundNode, distance);
}

int FModelBuilder::AddSpringDamper(const string& name, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double stiffness, double damping)
{
	int jointId = _model.GetJointsCount(); // 0-based
	FJointBuilder jointBuilder;

	int jointCharId = _model.GetJointCharsCount(); // 0-based
	FJointSpringDamperCharBuilder charBuilder;
	charBuilder.Setup(name, jointCharId, stiffness, damping);
	_model.AddJointChar(charBuilder.Get());
	jointBuilder.Setup(name, jointId,  bodyNode1, bodyNode2, jointCharId+1); // 1-based
	jointBuilder.SetSpringType(1);
	jointBuilder.SetLength(this->GetDistance(bodyNode1, bodyNode2));
	return _model.AddJoint(jointBuilder.Get());

}

int FModelBuilder::AddSpringDamper(const string& name, const BodyNodeNumber& bodyNode1, double stiffness, double damping)
{
	BodyNodeNumber groundNode(GetGroundBody(), AddGroundNode(bodyNode1));
	return AddSpringDamper(name, bodyNode1, groundNode, stiffness, damping);
}

int FModelBuilder::AddConstantForce(const string& name, int direction, const BodyNodeNumber& bodyNode, double param)
{
	_charId = _model.GetForceCharsCount(); // 0-based
	FForceConstantCharBuilder charBuilder;
	charBuilder.Setup(direction, name, _charId, param);
	_model.AddForceChar(charBuilder.Get());

	int forceId = _model.GetForcesCount(); // 0-based
	FForceBuilder forceBuilder;
	forceBuilder.Setup(name, forceId, bodyNode, _charId+1); // 1-based
	return _model.AddForce(forceBuilder.Get());	
}

int FModelBuilder::AddConstantForce(const string& name, const BodyNodeNumber& bodyNode, const Vec3& force)
{
	_charId = _model.GetForceCharsCount(); // 0-based
	FForceConstant3DCharBuilder charBuilder;
	charBuilder.Setup(name, _charId, force);
	_model.AddForceChar(charBuilder.Get());

	int forceId = _model.GetForcesCount(); // 0-based
	FForceBuilder forceBuilder;
	forceBuilder.Setup(name, forceId, bodyNode, _charId+1); // 1-based
	return _model.AddForce(forceBuilder.Get());	
}

int FModelBuilder::AddConstantForce(const string& name, const BodyNodeNumber& bodyNodeNumber1, const BodyNodeNumber& bodyNodeNumber2, double force)
{
	int geomId1 = _model.GetBodyByNumber(bodyNodeNumber1.bodyNumber).GeometryId();
	int geomId2 = _model.GetBodyByNumber(bodyNodeNumber2.bodyNumber).GeometryId();

	FGeometryPoint point1 = _model.GetGeometry(geomId1).GetNode(bodyNodeNumber1.nodeNumber);
	FGeometryPoint point2 = _model.GetGeometry(geomId2).GetNode(bodyNodeNumber2.nodeNumber);
	FGeometryPoint  delta = point2 - point1;
	Vec3 direction = MathHelpers::MakeVec3(delta.x,delta.y,delta.z);
	if(std::abs(direction.Magnitude()) > 1e-6)
	{
		direction = direction * (force / direction.Magnitude());
		return AddConstantForce(name, bodyNodeNumber2, direction);
	}
	return -1;
}

int FModelBuilder::AddConstantForceSpring(const string& name, const BodyNodeNumber& bodyNodeNumber1, const BodyNodeNumber& bodyNodeNumber2, double force)
{
	const double delta = 1e-6;
	double c1 = force / delta;
	int springId = AddSpringType2Damper(name, bodyNodeNumber1, bodyNodeNumber2, c1, 1e-6, delta, 50); 
	double springLength = _model.GetJoint(springId).SpringLength() + 1.;
	_model.GetJoint(springId).SpringLength(springLength);
	_model.GetJoint(springId).SspringLength(NumberToString(springLength));
	
	return springId;
}


int FModelBuilder::AddGCSPointToBodyGeometry(int bodyNumber, const Vec3& gcsPoint)
{
	int geomId = _model.GetBodyByNumber(bodyNumber).GeometryId();
	PGeometryBuilder geomBuilder = GetGeometryBuilder(geomId);
	const FTransform& transform = geomBuilder->Get().GetTransform();
	int pointNumber = geomBuilder->AddPoint(transform.FromGCS(gcsPoint));
	SetGeometry(geomId, geomBuilder->Get());

	return pointNumber;
}

int FModelBuilder::AddLCSPointToBodyGeometry(int bodyNumber, const Vec3& lcsPoint)
{
	int geomId = _model.GetBodyByNumber(bodyNumber).GeometryId();
	PGeometryBuilder geomBuilder = GetGeometryBuilder(geomId);
	int pointNumber = geomBuilder->AddPoint(lcsPoint);
	SetGeometry(geomId, geomBuilder->Get());

	return pointNumber;
}

int FModelBuilder::AddGroundNode(const BodyNodeNumber& bodyNode1)
{
	const FBody& body1 = _model.GetBodyByNumber(bodyNode1.bodyNumber);
	// TODO: логика создания 2 шарниров, требуется модификация геометрии тел

	const FGeometry& geom1 = _model.GetBodyGeometry(body1.Id());

	FTransform transform1 = geom1.GetTransform();
	Vec3 point11 = geom1.GetGeometryPoint(bodyNode1.nodeNumber-1);
	Vec3 point12 = transform1.ToGCS(point11);

	// вставка новых точек
	return AddLCSPointToBodyGeometry(GetGroundBody(), point12);
}

int FModelBuilder::GetGroundBody()
{
	if(!_groundBodyNumber)
		return AddGroundBody();
	return _groundBodyNumber;
}

int FModelBuilder::AddSpringType2Damper(const string& name, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double stiffness1, double stiffness2, double lInterval, double damping)
{
	int jointId = _model.GetJointsCount(); // 0-based
	FJointBuilder jointBuilder;

	int jointCharId = _model.GetJointCharsCount(); // 0-based
	FJointSpringType2DamperCharBuilder charBuilder;
	charBuilder.Setup(name, jointCharId, stiffness1, stiffness2, lInterval, damping);
	_model.AddJointChar(charBuilder.Get());
	jointBuilder.Setup(name, jointId,  bodyNode1, bodyNode2, jointCharId+1); // 1-based
	jointBuilder.SetSpringType(1);
	jointBuilder.SetLength(this->GetDistance(bodyNode1, bodyNode2));
	return _model.AddJoint(jointBuilder.Get());
}
