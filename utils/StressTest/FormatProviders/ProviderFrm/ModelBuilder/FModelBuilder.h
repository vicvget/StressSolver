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

	// TODO: ������� ���� � GCS
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

	// ��������� ��������� ���� �� ���� ���������� ����� � ���������� ������� ��������� (Global Coordinate System)
	// @param bodyNumber - ����� ����
	// @param gcsPoint - ����� � GCS
	// @returns ����� ����������� �����
	int AddGCSPointToBodyGeometry(int bodyNumber, const Vec3& gcsPoint);

	// ��������� ��������� ���� �� ���� ���������� ����� � ��������� ������� ��������� (Local Coordinate System)
	// @param bodyNumber - ����� ����
	// @param lcsPoint - ����� � LCS
	// @returns ����� ����������� �����
	int AddLCSPointToBodyGeometry(int bodyNumber, const Vec3& lcsPoint);

	void Setup(const string& name);

	// ���������� ����������
	int AddGravityForce(int direction, int bodyNumber);

	// ���������� ���������������� ���� (0 �������� �������)
	int AddFixedBodyByBoundTransform(const string& name, Vec3 minPoint, Vec3 maxPoint, Vec3 cmPoint, Mat3 mtxTransform, double m, Vec3 inertia);

	// ���������� ����� "GROUND" (0 �������� �������)
	int AddGroundBody();

	// ���������� ���������� ���� (6 �������� �������)
	int AddFreeBodyByBoundTransform(const string& name, Vec3 minPoint, Vec3 maxPoint, Vec3 cmPoint, Mat3 mtxTransform, double m, Vec3 inertia);

	//���������� ���������� ���� (6 �������� �������) �� ������ �� ��������
	int AddFreeBodyByMarkersPoints(const string& name, vector<Vec3> markersPoints, Vec3 cmPoint, Mat3 mtxTransform, double mass, Vec3 inertia);

	// ���������� ��������������� ������� � �� ������� ���� (direction - ����������� ������� ��� ��������)
	int AddRevolvedJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2);
	int AddRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2);

	// ���������� ������������ ������� � �� ������� ����
	int AddSphericalJoint(const string& name, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2);
	int AddSphericalJoint(const string& name, BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2);

	// ���������� �������, ��������������� 2 �������������� ������� � �� ������� ���� (direction - ����������� ������� ��� ��������)
	int AddCylindricalJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2);
	int AddCylindricalJoint(const string& name, int direction, BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2);

	// ���������� �������, ��������������� 1 �������������� ������� � �� ������� ���� (direction - ����������� �����������)
	int AddInPlaneJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2);
	int AddInPlaneJoint(const string& name, int direction, BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2);

	// ���������� 2 �������� ��� ������������� ���������������� (2 �������������� �������) ��������������� ������� � �������� ����������� ����� ���
	int AddCylindricalRevolvedJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance);

	// ���������� 2 �������� ��� ������������� ���������������� (2 �������������� �������) ��������������� ������� � �������� ����������� ����� ���
	int AddCylindricalRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance);

	// ���������� 2 �������� ��� ������������� ���������������� ��������������� ������� � �������� ����������� ����� ���
	int AddRevolvedJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance);

	// ���������� 2 �������� ��� ������������� ���������������� ��������������� ������� � �������� ����������� ����� ���
	int AddRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance);


	// ���������� 2 �������� ��� ������������� ���������������� ��������������� ������� � �������� ����������� ����� ���
	int AddTranslationalJoint(const string& name, int direction, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance);

	// ���������� 2 �������� ��� ������������� ���������������� ��������������� ������� � �������� ����������� ����� ���
	int AddTranslationalJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance);

	// ���������� 3 ��������: ������������,��������������� � inplane
	int AddFixedJoint(const string& name, int bodyNumber1, int nodeNumber1, int bodyNumber2, int nodeNumber2, double distance);

	// ���������� 3 ��������: ������������,��������������� � inplane
	int AddFixedJoint(const string& name, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double distance);

	// ���������� ������-������������ �����
	int AddSpringDamper(const string& name, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double stiffness, double damping);

#pragma region Single Reference Joints
	int AddRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1);

	// ���������� ������������ ������� � �� ������� ����
	int AddSphericalJoint(const string& name, BodyNodeNumber& bodyNode1);

	// ���������� �������, ��������������� 2 �������������� ������� � �� ������� ���� (direction - ����������� ������� ��� ��������)
	int AddCylindricalJoint(const string& name, int direction, BodyNodeNumber& bodyNode1);

	// ���������� �������, ��������������� 1 �������������� ������� � �� ������� ���� (direction - ����������� �����������)
	int AddInPlaneJoint(const string& name, int direction, BodyNodeNumber& bodyNode1);

	// ���������� 2 �������� ��� ������������� ���������������� (-3 translational DOF) ��������������� ������� � �������� ����������� ����� ���
	int AddRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, double distance);

	// ���������� 2 �������� ��� ������������� ���������������� (-2 translational DOF) ��������������� ������� � �������� ����������� ����� ���
	int AddCylindricalRevolvedJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, double distance);


	// ���������� 2 �������� ��� ������������� ���������������� ��������������� ������� � �������� ����������� ����� ���
	int AddTranslationalJoint(const string& name, int direction, const BodyNodeNumber& bodyNode1, double distance);
	
	//���������� 3 ��������: ������������,��������������� � inplane
	int AddFixedJoint(const string& name, const BodyNodeNumber& bodyNode1, double distance);

	// ���������� ������-������������ �����
	int AddSpringDamper(const string& name, const BodyNodeNumber& bodyNode1, double stiffness, double damping);

	// ���������� ������-������������ �����
	int AddSpringType2Damper(const string& name, const BodyNodeNumber& bodyNode1, const BodyNodeNumber& bodyNode2, double stiffness1, double stiffness2, double lInterval, double damping);

#pragma endregion

	// ���������� ���������� ����
	int AddConstantForce(const string& name, int direction, const BodyNodeNumber& bodyNode, double param);

	// ���������� ���������� ����
	int AddConstantForce(const string& name, const BodyNodeNumber& bodyNode, const Vec3& force);

	// ���������� ���������� ����
	int AddConstantForce(const string& name, const BodyNodeNumber& bodyNodeNumber1, const BodyNodeNumber& bodyNodeNumber2, double force);

	// ���������� ������ ������������ ���������� ������� � ��������� �������
	int AddSolverParamsSet(const string& name);

	// ���������� ��������� ������� �� �������� ��������� ��� � ������� ����� ���������� �������
	int AddCurrentIco(const string& name);

	int GetGroundBody();

	int AddConstantForceSpring(const string& name, const BodyNodeNumber& bodyNodeNumber1, const BodyNodeNumber& bodyNodeNumber2, double force);
/* TODO:
	// ���������� ���������� ����
	void AddConstantForce(const string& name, int direction, int bodyNumber, int nodeNumber);

	// ���������� ����� (���� ��� �������� �������)
	void AddGroundBody(const string& name);

	// ���������� ���������� ���� (6 �������� �������)
	void AddFreeBody(const string& name);

	// ���������� ��������������� ������� � �� ������� ���� (direction - ����������� ������� ��� ��������)
	void AddRevolvedJoint(const string& name, int direction, int bodyNumber1, int bodyNumber2);

	// ���������� �������, ��������������� ����������� � ���������
	void AddPlanarJoint(const string& name, int direction, int bodyNumber1, int bodyNumber2);
	*/
};


#endif // FMODEL_BUILDER_H