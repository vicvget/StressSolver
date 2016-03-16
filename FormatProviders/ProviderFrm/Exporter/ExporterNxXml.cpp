#include "ExporterNxXml.h"

#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <limits>
#include <algorithm>

#include "../../../Fcore/fcore.h"
#include "../../../Fcore/wrappers/StringRoutines.h"
#include "../../../AdditionalModules/RealsComparing/RealComparing.h"
#include "../../../AdditionalModules/Workarounds/GNU_gcc/string/to_string.h"
#include "../../../AdditionalModules/Workarounds/GNU_gcc/iomanip/put_time.h"


#define TOLERANCE 1e-8


enum class NxJointType
{
	Unsupported = 0,
	InPlane = 1,
	InLine = 2,
	Spherical = 3
};

enum class NxJointForceType
{
	Unsupported = 0,
	DirectedForce = 1,
	SpringDamper = 2
};

enum class NxForceType
{
	Unsupported = 0,
	Gravity = 1,
	Constant = 2
};


const string GetDofStringValue(int direction)
{
	switch (direction)
	{
		case 0: return "x";
		case 1: return "y";
		case 2: return "z";
	}
	return "";
}

/** Определение значения направленной силы, моделируемой пружиной
* @param jointChar - характеристика шарнира
* @return - значение силы
*/
double GetJointForceValue(const FJointChar& jointChar)
{
	const auto& params = 
		jointChar
		.GetSpringParams()
		.front()
		.params;
	
	double forceValue = 
		StringToNumber<double>(params.at(0)) *
		StringToNumber<double>(params.at(1));

	return forceValue;
}

/** Определение направленной силы
* @param forceChar - характеристика
* @return - вектор силы
*/
bool GetForceValue(const FForceChar& forceChar, Vec3Data& vectorForce)
{
	bool result = false;
	for(auto const& param: forceChar.Params())
	{
		int directionId = param.direction-1;
		if(directionId < 3) 
		{
			vectorForce[directionId] = StringToNumber<double>(param.params.at(0));
			result = true;
		}
	}
	return result;
}

/** Определение направленной силы
* @param forceChar - характеристика
* @return - вектор силы
*/
bool GetTorqueValue(const FForceChar& forceChar, Vec3Data& vectorTorque)
{
	bool result = false;
	for(auto const& param: forceChar.Params())
	{
		int directionId = param.direction-1;
		if(directionId > 2 && directionId < 6)
		{
			directionId -= 3;
			vectorTorque[directionId] = StringToNumber<double>(param.params.at(0));
		}
	}
	return result;
}

/** Определение значения направленной силы, моделируемой пружиной
* @param jointChar - характеристика шарнира
* @return - значение силы
*/
void GetSpringDamperParams(const FJointChar& jointChar, double& stiffnessFactor, double& dampingFactor)
{
	stiffnessFactor = StringToNumber<double>(
		jointChar
		.GetSpringParams()
		.front()
		.params
		.at(0));
		
	dampingFactor = StringToNumber<double>(
		jointChar
		.GetDampingParams()
		.front()
		.params
		.at(0));
}


/** Определение типа силы по характеристике
* @param jointChar - характеристика силы
*/
NxForceType GetForceType(const FForceChar& forceChar)
{
	const auto& params = forceChar.Params();
	std::size_t nDofs = params.size();
	if(nDofs == 1 && forceChar.Params().front().type == 3)
		return NxForceType::Gravity;

	for(const auto& param:params)
	{
		if(param.type != 2)
			return NxForceType::Unsupported;
	}
	return NxForceType::Constant;
}


/** Определение типа направленной силы по характеристике (сила или пружина)
* @param jointChar - характеристика шарнира
*/
NxJointForceType GetJointForceType(const FJointChar& jointChar)
{
	const auto& springParams = jointChar.GetSpringParams();
	std::size_t nRestrictedDofs = springParams.size();
	if(nRestrictedDofs != 1)
		return NxJointForceType::Unsupported;

	CharComponent dampingParamsComponent = jointChar
		.GetDampingParams()
		.front();

	CharComponent springParamsComponent = springParams.front();
	if(springParamsComponent.direction != 1)
		return NxJointForceType::Unsupported;
	if(springParamsComponent.type == 2)
	{
		return NxJointForceType::DirectedForce;
	}
	else if(springParamsComponent.type == 1 && dampingParamsComponent.type == 50)
	{
		return NxJointForceType::SpringDamper;
	}
	return NxJointForceType::Unsupported;
}

/** Определение типа кинематического шарнира по характеристике (InLine, InPlane, Spherical)
* @param jointChar - характеристика шарнира
* @param direction - направление нормали к плоскости в случае InPlane или  разрешенное направление в случае InLine
*/
NxJointType GetJointType(const FJointChar& jointChar, int& direction)
{
	auto const& springParams = jointChar.GetSpringParams();
	std::size_t nRestrictedDofs = springParams.size();
	if(nRestrictedDofs > 3)
		return NxJointType::Unsupported;
	
	bool dofs[3];
	std::fill(dofs, dofs+3, false);
	for(auto const& springParam:springParams)
	{
		if(springParam.type != 15 || springParam.direction > 3)
			return NxJointType::Unsupported;
		dofs[springParam.direction-1] = true;
	}
	direction = 1;
	for(std::size_t i = 0; i < 3; i++)
	{
		if((nRestrictedDofs == 1) == dofs[i]) //(true+true & false+false)
		{
			direction = i+1;
			break;
		}
	}
	return (NxJointType)nRestrictedDofs;
}

string GetNextId(int& prevId)
{
	stringstream ss;
	ss << "id" << prevId++;
	return ss.str();
}

string GetRefId(int id)
{
	stringstream ss;
	ss << "#id" << id;
	return ss.str();
}

double ToTolerance(double src)
{
	const double tolerance = TOLERANCE;
	if (std::abs(src) < tolerance)
		return 0.;
	else if (std::abs(1-src) < tolerance)
		return 1.;
	else if (std::abs(1+src) < tolerance)
		return -1.;
	return src;
}

string VecToString(const Vec3& vec)
{
	stringstream ss;
	ss.precision(15);
	ss << ToTolerance(vec[0]) << ' ' << ToTolerance(vec[1]) << ' ' << ToTolerance(vec[2]);
	return ss.str();
}

Mat3 MakeZBasis(const Vec3& zAxis)
{
	Vec3 destZAxis = MathHelpers::MakeVec3(0., 0., 1.);
	Vec3 normalVector = zAxis.Cross(destZAxis);
	if(normalVector.Magnitude() < std::numeric_limits<double>::epsilon())
	{
		return Mat3::Identity();
	}

	Vec3 xAxis, yAxis;

	if(DoubleComparing::IsZero(zAxis.X()))
	{
		xAxis = MathHelpers::MakeVec3(1.,0,0);
		yAxis = zAxis.Cross(xAxis);
	} 
	else if(DoubleComparing::IsZero(zAxis.Y()))
	{
		yAxis = MathHelpers::MakeVec3(0,1.,0);
		xAxis = yAxis.Cross(zAxis);
	}
	else if(DoubleComparing::IsZero(zAxis.Z()))
	{
		xAxis = MathHelpers::MakeVec3(0,0,1.);
		yAxis = zAxis.Cross(xAxis);
	}
	else
	{
		double minComponent = std::abs(zAxis[0]);
		int id = 0;
		for(int i = 1; i < 3; i++)
		{
			if(minComponent > std::abs(zAxis[i]))
			{
				minComponent = std::abs(zAxis[i]);
				id = i;
			}
		}
		if(id == 0)
		{
			xAxis = MathHelpers::MakeVec3(0., zAxis[2], -zAxis[1]).Normalize();
		}
		else if(id == 1)
		{
			xAxis = MathHelpers::MakeVec3(zAxis[2], 0., -zAxis[0]).Normalize();
		}
		else
		{
			xAxis = MathHelpers::MakeVec3(zAxis[1], -zAxis[0], 0.).Normalize();		
		}
		yAxis = zAxis.Cross(xAxis);
	}

	Mat3 matrix(xAxis, yAxis, zAxis);
	return matrix.Tr();
}

void GetBound(const std::vector<MathHelpers::Vec3>& points, MathHelpers::Vec3& minBound, MathHelpers::Vec3& maxBound)
{
	if(points.size() == 0)
		return;

	double minX, minY, minZ, maxX, maxY, maxZ;
	minX = maxX = points[0].X();
	minY = maxY = points[0].Y();
	minZ = maxZ = points[0].Z();

	for(int i = 1; i < points.size(); i++)
	{
		if(points[i].X() < minX) minX = points[i].X();
		else if(points[i].X() > maxX) maxX = points[i].X();

		if(points[i].Y() < minY) minY = points[i].Y();
		else if(points[i].Y() > maxY) maxY = points[i].Y();	

		if(points[i].Z() < minZ) minZ = points[i].Z();
		else if(points[i].Z() > maxZ) maxZ = points[i].Z();	
	}

	minBound = MathHelpers::MakeVec3(minX, minY, minZ);
	maxBound = MathHelpers::MakeVec3(maxX, maxY, maxZ);
}

string MatVecToString(const Mat3& mat, const Vec3& vec)
{
	stringstream ss;
	ss.precision(15);
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			ss << ToTolerance(mat[i*3+j]) << ' ';
		}
		ss << "0 ";
	}
	ss << VecToString(vec) << " 1";
	
	return ss.str();
}


ExporterNxXml::ExporterNxXml
	(
		const FModel* model
	)
	:
		ExporterFModel(model)
{
	_model->Eval();
	if(MechanismOccurrence::Load("MechanismOccurrence.xml", _mechanismOccurrenceData))
		std::cout << "MechanismOccurrence\n";
}

TiXmlElement* ExporterNxXml::ExportBody(int& id, int bodyId, TiXmlElement* parent) const
{
	const FBody& body = _model->GetBody(bodyId);
	const FGeometry& geo = _model->GetBodyGeometry(bodyId);
	FTransform tr = geo.GetTransform();

	BodyInfo& bodyInfo = _bodyInfos[bodyId];
	FBodyDof dof;
	body.GetBodyDof(dof);
	bodyInfo._isGround = dof.IsFixed();

	if(bodyInfo._isGround)
	{
		return nullptr;
	}

	TiXmlElement* xmlBody = AddMechanismRevisionView(id, body.Name(), "link", parent);
	
	Vec3 boundMinPoint, boundMaxPoint;
	std::vector<Vec3> points;
	geo.GetPoints(points, true);
	//GetBound(points, boundMinPoint, boundMaxPoint);

	geo.GetBound(boundMinPoint, boundMaxPoint);
	boundMinPoint = tr.ToGCS(boundMinPoint);
	boundMaxPoint = tr.ToGCS(boundMaxPoint);
	int cmNodeFModelId = geo.CmNodeNumber();

	Vec3 cmPoint = points[cmNodeFModelId];
	
	int cmNodeXmlId = -1;
	AddInertia(id, body.M(), 
		MakeVec3(body.Jx(), body.Jy(), body.Jz()),
		cmPoint,
		tr.GetTranslation(),
		tr.GetRotation(),
		xmlBody,
		cmNodeXmlId);
	bodyInfo._cmNodeXmlId = cmNodeXmlId;

	TiXmlElement* xmlBound = AddElement(id, "Bound", xmlBody);  
	string strBound = VecToString(boundMinPoint) + " " + VecToString(boundMaxPoint);
	xmlBound->SetAttribute("values", strBound.c_str());

	TiXmlElement* xmlUserData = AddElement(id, "UserData", xmlBound);  
	AddUserValue("BOUNDING_TYPE", 0, xmlUserData);
	
	// Mechanism occurrence if required
	auto mechanismOccurrenceId = _mechanismOccurrenceData.find(body.Name());
	if(mechanismOccurrenceId!=_mechanismOccurrenceData.end())
	{
		TiXmlElement* xmlMechanismOccurrence = AddElement(id, "MechanismOccurrence", xmlBody);  
		TiXmlElement* xmlReference = AddElement(id, "Reference", xmlMechanismOccurrence);  
		TiXmlElement* xmlLocatorRef = AddElement(id, "LocatorRef", xmlReference);  
		string locationRefValue = "#PLMXML(UGPrt-doc/occ[@id=" + to_string(mechanismOccurrenceId->second) +"])";
		xmlLocatorRef->SetAttribute("locationRef", locationRefValue.c_str());
	}

	bodyInfo.nodeIds.resize(points.size());

	// Insert markers
	for(int i = 0; i < points.size(); i++)
	{
		if(i != cmNodeFModelId) // cm is already inserted
		{
			bodyInfo.nodeIds[i] = id;
			AddMarker(id, "", points[i], Mat3::Identity(), xmlBody);
		}
		else
		{
			bodyInfo.nodeIds[cmNodeFModelId] = cmNodeXmlId;
		}

	}
	std::cout << body.Name() << std::endl;
	return xmlBody;
}

void ExporterNxXml::ExportTo(const string& file) const
{
	TiXmlDocument doc;  
	doc.LinkEndChild(new TiXmlDeclaration("1.0", "utf-8", ""));  

	TiXmlElement * root = new TiXmlElement("PLMXML");  
	doc.LinkEndChild(root);  
	root->SetAttribute("xmlns", "http://www.plmxml.org/Schemas/PLMXMLSchema");
	root->SetAttribute("schemaVersion", 6);

	time_t time = std::time(nullptr);
	auto localTime = *std::localtime(&time);
	stringstream ssDate, ssTime;
	ssDate << put_time(&localTime, "%Y-%m-%d");
	root->SetAttribute("date", ssDate.str().c_str());	

	ssTime << put_time(&localTime, "%H:%M:%S");
	root->SetAttribute("time", ssTime.str().c_str());	
	
	root->SetAttribute("author", "NX 8.5.3.3");

	int id = 1;

	// Add <MechanismRevisionView>
	TiXmlElement *xmlModelHeader = AddMechanismRevisionView(id, _model->Header().Name().c_str(), "mechanism", root);

	// Add <PropertyGroup>
	TiXmlElement *xmlPropertyGroup = AddElement(id, "PropertyGroup", xmlModelHeader);
	
	// Add <Gravity>
	TiXmlElement *xmlGravity = AddElement(id, "Gravity", xmlPropertyGroup);

	// Setup <Gravity>
	Vec3 gravityVector = MathHelpers::MakeVec3(0, 0, -9.81);
	xmlGravity->SetAttribute("value", VecToString(gravityVector).c_str());

	// Add <MechanismParameters>
	TiXmlElement *xmlMechanismParameters = AddElement(id, "MechanismParameters", xmlModelHeader);

	// Setup <MechanismParameters>
	xmlMechanismParameters->SetAttribute("maxIntegratorStep", "0.01");
	xmlMechanismParameters->SetAttribute("integratorError", "0.001");

	// Add <SolverParameter>
	TiXmlElement *xmlSolverParameter = AddElement(id, "SolverParameter", xmlMechanismParameters);

	// Setup <SolverParameter>
	xmlSolverParameter->SetAttribute("value", "10");
	xmlSolverParameter->SetAttribute("title", "maxIntegratorIterations");

	// Add <SolverParameter>
	TiXmlElement *xmlSolverParameter2 = AddElement(id, "SolverParameter", xmlMechanismParameters);

	// Setup <SolverParameter>
	xmlSolverParameter2->SetAttribute("value", "50");
	xmlSolverParameter2->SetAttribute("title", "maxKinematicsIterations");

	// Add <SolverParameter>
	TiXmlElement *xmlSolverParameter3 = AddElement(id, "SolverParameter", xmlMechanismParameters);

	// Setup <SolverParameter>
	xmlSolverParameter3->SetAttribute("value", "100");
	xmlSolverParameter3->SetAttribute("title", "maxEquilibriumIterations");

	std::cout << "//// Bodies:\n";
	_bodyInfos.resize(_model->CountBodies());
	for(int i = 0; i<_model->CountBodies(); i++)
	{
		ExportBody(id, i, root);
	}

	//<Constraint id="id141" type="joint" />
	int refId = id;

	std::cout << "//// Joints:\n";
	AddConstraintSeparator(id, "joint", root);
	for(int i = 0; i<_model->GetJointsCount(); i++)
	{
		ExportJoint(id, refId, i, root);
	}
	//<Constraint id="id241" type="force" /> 
	refId = id;

	std::cout << "//// Joint Forces:\n";
	AddConstraintSeparator(id, "force", root);
	for(int i = 0; i<_model->GetJointsCount(); i++)
	{
		ExportJointForces(id, refId, i, root); // joint forces + spring-dampers
	}
	// export forces
	for(int i = 0; i<_model->GetForcesCount(); i++)
	{
		ExportForces(id, refId, i, root, xmlModelHeader); // vector forces
	}
	// export momentums
	for(int i = 0; i<_model->GetForcesCount(); i++)
	{
		ExportTorques(id, refId, i, root, xmlModelHeader); // vector momentums
	}

	// Export SpringDampers
	if(_springDamperInfos.size() > 0)
	{
		std::cout << "//// Dampers:\n";
		refId = id;
		AddConstraintSeparator(id, "damper", root);
		ExportDampers(id, refId, root);
		
		std::cout << "//// Springs:\n";
		refId = id;
		AddConstraintSeparator(id, "spring", root);
		ExportSprings(id, refId, root);
	}

	// Коррекция матриц преобразования маркеров
	for (const auto& marker : _markerInfos)
	{
		marker.second.AddTransformToXML();
	}

	doc.SaveFile(file.c_str());
}

string ExporterNxXml::GetExt() const
{
	return EXT_XML;
}

TiXmlElement* ExporterNxXml::AddTransform(int& id, const Vec3& translate, TiXmlElement* parent) const
{
	TiXmlElement *xmlTransform = AddElement(id, "Transform", parent);
	xmlTransform->SetAttribute("type", "translate");
	
	xmlTransform->LinkEndChild(new TiXmlText(VecToString(translate).c_str()));
	return xmlTransform;
}

TiXmlElement* ExporterNxXml::AddTransform(int& id, const Vec3& translate, const Mat3& rotation, TiXmlElement* parent) const
{
	TiXmlElement *xmlTransform = AddElement(id, "Transform", parent);
	xmlTransform->LinkEndChild(new TiXmlText(MatVecToString(rotation, translate).c_str()));
	return xmlTransform;
}

TiXmlElement* ExporterNxXml::AddMarker(int& id, const string& name, TiXmlElement* parent) const
{
	TiXmlElement *xmlMarker = AddElement(id, "Marker", parent);
	if(!name.empty())
		xmlMarker->SetAttribute("name", name.c_str());
	return xmlMarker;
}

TiXmlElement* ExporterNxXml::AddMarker(int& id, const string& name, const Vec3& translate, const Mat3& rotation, TiXmlElement* parent) const
{
	int prevId = id;
	TiXmlElement *xmlMarker = AddElement(id, "Marker", parent);
	if(!name.empty())
		xmlMarker->SetAttribute("name", name.c_str());
	TiXmlElement *xmlTransform = AddTransform(id, translate, rotation, xmlMarker);
	_markerInfos[prevId] = MarkerInfo(translate, rotation, xmlTransform);
	return xmlMarker;
}

TiXmlElement* ExporterNxXml::AddInertia(int& id, double mass, const Vec3& inertia, const Vec3& cmPoint, const Vec3& translate, const Mat3& rotation, TiXmlElement* parent, int& cmNodeId) const
{
	//<CentreOfMass>
	//<MassProperty>
	//   <UserData>
	//      <UserValue>
	//<MechanismInertia>
	//<InitialVelocity>
	//<Marker>
	TiXmlElement *xmlPropertyGroup = AddElement(id, "PropertyGroup", parent);
	TiXmlElement *xmlCentreOfMass = AddElement(id, "CentreOfMass", xmlPropertyGroup);
	xmlCentreOfMass->SetAttribute("value", VecToString(cmPoint).c_str());

	TiXmlElement *xmlMassProperty = AddElement(id, "MassProperty", xmlPropertyGroup);
	TiXmlElement *xmlValueWithUnit = AddElement(id, "ValueWithUnit", xmlMassProperty);

	TiXmlElement *xmlUserData = AddElement(id, "UserData", xmlMassProperty);
	AddUserValue("MassType", 0, xmlUserData);
	AddUserValue("UserDefined", 0, xmlUserData);

	const string& strMass = to_string(mass);
	xmlValueWithUnit->SetAttribute("value", strMass.c_str());

	TiXmlElement *xmlMechanismInertia = AddElement(id, "MechanismInertia", xmlPropertyGroup);
	string inertiaString = VecToString(inertia) + " 0 0 0";
	xmlMechanismInertia->SetAttribute("value", inertiaString.c_str());
	xmlMechanismInertia->SetAttribute("markerRef", GetRefId(id).c_str());

	cmNodeId = id;
	AddMarker(id, "", translate, rotation, parent);
	AddElement(id, "InitialVelocity", xmlPropertyGroup);

	return xmlPropertyGroup;
}

TiXmlElement* ExporterNxXml::AddElement(int& id, const char* element, TiXmlElement* parent) const
{
	TiXmlElement *xmlElement = new TiXmlElement(element);
	parent->LinkEndChild(xmlElement);
	xmlElement->SetAttribute("id", GetNextId(id).c_str());
	return xmlElement;
}

TiXmlElement* ExporterNxXml::AddElement(int& id, const char* element, const char* type, TiXmlElement* parent) const
{
	TiXmlElement* res = AddElement(id, element, parent);
	res->SetAttribute("type", type);
	return res;
}

TiXmlElement* ExporterNxXml::AddMechanismRevisionView(int& id, const string& name, const string subtype, TiXmlElement* parent) const
{
	TiXmlElement *xmlElement = AddElement(id, "MechanismRevisionView", parent);
	xmlElement->SetAttribute("name", name.c_str());
	xmlElement->SetAttribute("subType", subtype.c_str());
	return xmlElement;
}

TiXmlElement* ExporterNxXml::AddUserValue(const char* title, int value, TiXmlElement *parent) const
{
	TiXmlElement *xmlUserValue = new TiXmlElement("UserValue");
	parent->LinkEndChild(xmlUserValue);
	xmlUserValue->SetAttribute("value", value);
	xmlUserValue->SetAttribute("title", title);
	return xmlUserValue;
}

TiXmlElement* ExporterNxXml::AddUserValue(const char* title, const char* value, TiXmlElement *parent) const
{
	TiXmlElement *xmlUserValue = new TiXmlElement("UserValue");
	parent->LinkEndChild(xmlUserValue);
	xmlUserValue->SetAttribute("value", value);
	xmlUserValue->SetAttribute("title", title);
	return xmlUserValue;
}


TiXmlElement* ExporterNxXml::ExportJoint(int& id, int instancedRef, int jointId, TiXmlElement* parent) const
{
	const FJoint& joint = _model->GetJoint(jointId);
	const FJointChar& jointChar = _model->GetJointCharByNumber(joint.CharNumber());
	
	int direction;
	NxJointType type = GetJointType(jointChar, direction);
	if(type == NxJointType::Unsupported)
		return nullptr;

	string refIdConstraint = GetRefId(instancedRef);
	//TiXmlElement* xmlConstraint = AddElement(id, "Constraint", "joint", parent);
	TiXmlElement* xmlConstraintInstance = AddElement(id, "ConstraintInstance", parent);
	xmlConstraintInstance->SetAttribute("name", joint.Name().c_str());
	xmlConstraintInstance->SetAttribute("instancedRef", refIdConstraint.c_str());

	int bodyIds[2] = {joint.BodyNumber1()-1, joint.BodyNumber2()-1};
	int nodeIds[2] = {joint.NodeNumber1()-1, joint.NodeNumber2()-1};
	//int markerIds[2];

// получаем вектор в системе глобальной координат
	MathHelpers::Vec3 directionVec = MathHelpers::MakeVec3(0,0,0);
	directionVec[direction-1] = 1.;
	const FGeometry& geo = _model->GetBodyGeometry(bodyIds[0]);
	FTransform tr = geo.GetTransform();
	directionVec = tr.ToGCSRotation(directionVec);
// строим по вектору Z-базис
	MathHelpers::Mat3 zBasis;
	zBasis = MakeZBasis(directionVec);

	Vec3 markerTranslation; // if Ground it is zero (default)
	bool err = false;
	for(int currentRefId = 0; currentRefId < 2; currentRefId++)
	{
		const BodyInfo& bodyInfo = _bodyInfos.at(bodyIds[currentRefId]);
		if(!(bodyInfo._isGround))
		{
			int nodeId = nodeIds[currentRefId];
			if(nodeId >= 0 && nodeId < bodyInfo.nodeIds.size())
			{			
				TiXmlElement* xmlConstraintTargetRef = AddElement(id, "ConstraintTargetRef", xmlConstraintInstance);
				int markerId = bodyInfo.nodeIds.at(nodeId);
				if(type == NxJointType::InLine)
				{
					auto markerIt = _markerInfos.find(markerId);
					if(markerIt != _markerInfos.end())
					{
						markerIt->second._rotation = zBasis;
						if(currentRefId > 0) // not Ground if second RefId!
							markerTranslation = markerIt->second._translation;
					}
					else
					{
						std::cerr << "Invalid marker: " << markerId << " bodyId: " << bodyIds[currentRefId] << std::endl;
						err = true;
					}
				}
				xmlConstraintTargetRef->SetAttribute("targetRef", GetRefId(markerId).c_str());
			}
			else
			{
				std::cerr << "Invalid nodeId: " << nodeId << " bodyId: " << bodyIds[currentRefId] << std::endl;
			}
		}
	}

	if(!err)
	{
		string types[] = {"", "inplane", "inline", "spherical"};
	
		TiXmlElement* xmlJointData = AddElement(id, "JointData", xmlConstraintInstance);
		xmlJointData->SetAttribute("type", types[(int)type].c_str());
		TiXmlElement* xmlUserData = AddElement(id, "UserData", xmlJointData);

		AddUserValue("I_DIRECTION", (VecToString(markerTranslation) + string(" ") + VecToString(directionVec)).c_str(), xmlUserData);
		AddUserValue("J_DIRECTION", (VecToString(markerTranslation) + string(" ") + VecToString(directionVec)).c_str(), xmlUserData);
		std::cout << joint.Name() << std::endl;
		return xmlConstraintInstance;
	}

	return nullptr;
}


TiXmlElement* ExporterNxXml::ExportJointForces(int& id, int instancedRef, int jointId, TiXmlElement* parent) const
{
	const FJoint& joint = _model->GetJoint(jointId);
	if(joint.SpringType() == 0)
	{
		return nullptr;
	}

	const FJointChar& jointChar = _model->GetJointCharByNumber(joint.CharNumber());
	NxJointForceType type  = GetJointForceType(jointChar);
	TiXmlElement* xmlConstraintInstance = nullptr;

	switch(type)
	{
		case NxJointForceType::DirectedForce:
			{
				string refIdConstraint = GetRefId(instancedRef);
				xmlConstraintInstance = AddElement(id, "ConstraintInstance", parent);
				xmlConstraintInstance->SetAttribute("name", joint.Name().c_str());
				xmlConstraintInstance->SetAttribute("instancedRef", refIdConstraint.c_str());

				int bodyIds[2] = {joint.BodyNumber1()-1, joint.BodyNumber2()-1};
				int nodeIds[2] = {joint.NodeNumber1()-1, joint.NodeNumber2()-1};
				bool err = false;
				for(int i = 0; i < 2; i++)
				{
					AddConstraintTargetRef(id, xmlConstraintInstance, bodyIds[i], nodeIds[i]);
					//const BodyInfo& bodyInfo = _bodyInfos.at(bodyIds[i]);
					//if(!(bodyInfo._isGround))
					//{
					//	int nodeId = nodeIds[i];
					//	if(nodeId >= 0 && nodeId < bodyInfo.nodeIds.size())
					//	{			
					//		TiXmlElement* xmlConstraintTargetRef = AddElement(id, "ConstraintTargetRef", xmlConstraintInstance);
					//		int markerId = bodyInfo.nodeIds.at(nodeId);
					//		xmlConstraintTargetRef->SetAttribute("targetRef", GetRefId(markerId).c_str());
					//	}
					//	else
					//	{
					//		std::cerr << "Invalid nodeId: " << nodeId << " bodyId: " << bodyIds[i] << std::endl;
					//		err = true;
					//	}
					//}
				}
				if(err) return nullptr;

				TiXmlElement* xmlForceData = AddElement(id, "ForceData", xmlConstraintInstance);
				TiXmlElement* xmlForceComponent = AddElement(id, "ForceComponent", xmlForceData);
				xmlForceComponent->SetAttribute("force", NumberToString(GetJointForceValue(jointChar)).c_str());
				std::cout << joint.Name() << std::endl;
			}
			break;

		case NxJointForceType::SpringDamper:
			{
				string refIdConstraint = GetRefId(instancedRef);

				//TiXmlElement* xmlConstraintInstance = AddElement(id, "ConstraintInstance", parent);
				//xmlConstraintInstance->SetAttribute("name", joint.Name().c_str());
				//xmlConstraintInstance->SetAttribute("instancedRef", refIdConstraint.c_str());

				int bodyIds[2] = {joint.BodyNumber1()-1, joint.BodyNumber2()-1};
				int nodeIds[2] = {joint.NodeNumber1()-1, joint.NodeNumber2()-1};
				SpringDamperInfo springDamper;
				springDamper._name = joint.Name();
				
				GetSpringDamperParams(jointChar, springDamper._stiffnessFactor, springDamper._dampingFactor);
				springDamper._length = joint.SpringLength();
				springDamper._preload = 0;
				bool err = false;
				for(int i = 0; i < 2; i++)
				{
					const BodyInfo& bodyInfo = _bodyInfos.at(bodyIds[i]);
					if(!(bodyInfo._isGround))
					{
						int nodeId = nodeIds[i];
						if(nodeId >= 0 && nodeId < bodyInfo.nodeIds.size())
						{			
							//TiXmlElement* xmlConstraintTargetRef = AddElement(id, "ConstraintTargetRef", xmlConstraintInstance);
							int markerId = bodyInfo.nodeIds.at(nodeId);
							springDamper._markerIds.push_back(markerId);
							//xmlConstraintTargetRef->SetAttribute("targetRef", GetRefId(markerId).c_str());
						}
						else
						{
							std::cerr << "Invalid nodeId: " << nodeId << " bodyId: " << bodyIds[i] << std::endl;
							err = true;
						}
					}
				}
				if(!err)
				{
					_springDamperInfos.push_back(springDamper);
				}
					
				//TiXmlElement* xmlForceData = AddElement(id, "ForceData", xmlConstraintInstance);
				//TiXmlElement* xmlForceComponent = AddElement(id, "ForceComponent", xmlForceData);
				//xmlForceComponent->SetAttribute("force", NumberToString(GetJointForceValue(jointChar)).c_str());
			}
			break;
		case NxJointForceType::Unsupported:
		default:
			return nullptr;
	}

	return xmlConstraintInstance;
}

TiXmlElement* ExporterNxXml::ExportForces(int& id, int instancedRef, int forceId, TiXmlElement* parent, TiXmlElement* header) const
{
	const FForce& force = _model->GetForce(forceId);
	const FForceChar& forceChar = _model->GetForceCharByNumber(force.CharNumber());
	NxForceType type  = GetForceType(forceChar);
	TiXmlElement* xmlConstraintInstance = nullptr;

	switch(type)
	{
	case NxForceType::Constant:
		{
			Vec3Data vectorForce;
			if(GetForceValue(forceChar, vectorForce))
			{
				string refIdConstraint = GetRefId(instancedRef);
				xmlConstraintInstance = AddElement(id, "ConstraintInstance", parent);
				xmlConstraintInstance->SetAttribute("name", force.Name().c_str());
				xmlConstraintInstance->SetAttribute("instancedRef", refIdConstraint.c_str());

				int bodyId = force.BodyNumber()-1;
				int nodeId = force.NodeNumber()-1;
				AddConstraintTargetRef(id, xmlConstraintInstance, bodyId, nodeId);
				int markerId = id;
				AddMarker(id, "default_name", header);
				AddConstraintTargetRefToMarker(id, xmlConstraintInstance, markerId);

				TiXmlElement* xmlForceData = AddElement(id, "ForceData", xmlConstraintInstance);
				xmlForceData->SetAttribute("scalar", "false");

				for(int i = 0; i < 3; i++)
				{
					TiXmlElement* xmlForceComponent = AddElement(id, "ForceComponent", xmlForceData);
					xmlForceComponent->SetAttribute("degreeOfFreedom", GetDofStringValue(i).c_str());
					xmlForceComponent->SetAttribute("force", to_string(vectorForce[i]).c_str());
				}
				std::cout << force.Name() << std::endl;
			}
		}
		break;

	case NxForceType::Unsupported:
	default:
		return nullptr;
	}

	return xmlConstraintInstance;
}


TiXmlElement* ExporterNxXml::ExportTorques(int& id, int instancedRef, int forceId, TiXmlElement* parent, TiXmlElement* header) const
{
	const FForce& force = _model->GetForce(forceId);
	const FForceChar& forceChar = _model->GetForceCharByNumber(force.CharNumber());
	NxForceType type  = GetForceType(forceChar);
	TiXmlElement* xmlConstraintInstance = nullptr;

	switch(type)
	{
	case NxForceType::Constant:
		{
			Vec3Data vectorTorque;
			if(GetForceValue(forceChar, vectorTorque))
			{
				string refIdConstraint = GetRefId(instancedRef);
				xmlConstraintInstance = AddElement(id, "ConstraintInstance", parent);
				xmlConstraintInstance->SetAttribute("name", force.Name().c_str());
				xmlConstraintInstance->SetAttribute("instancedRef", refIdConstraint.c_str());

				int bodyId = force.BodyNumber()-1;
				int nodeId = force.NodeNumber()-1;
				AddConstraintTargetRef(id, xmlConstraintInstance, bodyId, nodeId);
				int markerId = id;
				AddMarker(id, "default_name", header);
				AddConstraintTargetRefToMarker(id, xmlConstraintInstance, markerId);

				TiXmlElement* xmlMomentumData = AddElement(id, "ForceData", xmlConstraintInstance);
				xmlMomentumData->SetAttribute("translation", "false");
				xmlMomentumData->SetAttribute("scalar", "false");
				for(int i = 0; i < 3; i++)
				{
					TiXmlElement* xmlMomentumComponent = AddElement(id, "ForceComponent", xmlMomentumData);
					xmlMomentumComponent->SetAttribute("degreeOfFreedom", GetDofStringValue(i).c_str());
					xmlMomentumComponent->SetAttribute("force", to_string(vectorTorque[i]).c_str());
				}
				std::cout << force.Name() << std::endl;
			}
		}
		break;

	case NxForceType::Unsupported:
	default:
		return nullptr;
	}

	return xmlConstraintInstance;
}



TiXmlElement* ExporterNxXml::AddConstraintSeparator(int& id, const char* type, TiXmlElement* parent) const
{
	TiXmlElement *xmlConstraint = AddElement(id, "Constraint", parent);
	xmlConstraint->SetAttribute("type", type);
	return xmlConstraint;
}

void ExporterNxXml::ExportDampers(int& id, int instancedRef, TiXmlElement* parent) const
{
	for(const SpringDamperInfo& springDamper: _springDamperInfos)
	{
		string refIdConstraint = GetRefId(instancedRef);
		
		TiXmlElement* xmlConstraintInstance = AddElement(id, "ConstraintInstance", parent);
		string damperName = "D" + springDamper._name;
		
		xmlConstraintInstance->SetAttribute("name", damperName.c_str());
		xmlConstraintInstance->SetAttribute("instancedRef", refIdConstraint.c_str());
		for(int i=0; i < springDamper._markerIds.size(); i++)
		{
			TiXmlElement* xmlConstraintTargetRef = AddElement(id, "ConstraintTargetRef", xmlConstraintInstance);
			xmlConstraintTargetRef->SetAttribute("targetRef", GetRefId(springDamper._markerIds[i]).c_str());
		}
		TiXmlElement* xmlDampingData = AddElement(id, "DampingData", xmlConstraintInstance);
		xmlDampingData->SetAttribute("coefficient", NumberToString(springDamper._dampingFactor).c_str());
		
		std::cout << "Damper " << springDamper._name << std::endl;
	}
}

void ExporterNxXml::ExportSprings(int& id, int instancedRef, TiXmlElement* parent) const
{
	for(const SpringDamperInfo& springDamper: _springDamperInfos)
	{
		string refIdConstraint = GetRefId(instancedRef);

		TiXmlElement* xmlConstraintInstance = AddElement(id, "ConstraintInstance", parent);
		xmlConstraintInstance->SetAttribute("name", springDamper._name.c_str());
		xmlConstraintInstance->SetAttribute("instancedRef", refIdConstraint.c_str());
		for(int i=0; i < springDamper._markerIds.size(); i++)
		{
			TiXmlElement* xmlConstraintTargetRef = AddElement(id, "ConstraintTargetRef", xmlConstraintInstance);
			xmlConstraintTargetRef->SetAttribute("targetRef", GetRefId(springDamper._markerIds[i]).c_str());
		}
		TiXmlElement* xmlSpringData = AddElement(id, "SpringData", xmlConstraintInstance);
		xmlSpringData->SetAttribute("stiffness", NumberToString(springDamper._stiffnessFactor).c_str());
		xmlSpringData->SetAttribute("preloadForce", NumberToString(springDamper._preload).c_str());
		xmlSpringData->SetAttribute("free", NumberToString(springDamper._length).c_str());

		std::cout << "Spring " << springDamper._name << std::endl;
	}
}

TiXmlElement* ExporterNxXml::AddConstraintTargetRef(int& id, TiXmlElement* xmlConstraintInstance, int bodyId, int nodeId) const
{
	const BodyInfo& bodyInfo = _bodyInfos.at(bodyId);
	if(!(bodyInfo._isGround))
	{
		if(!(nodeId >= 0 && nodeId < bodyInfo.nodeIds.size()))
		{			
			std::cerr << "Invalid nodeId: " << nodeId << " bodyId: " << bodyId << std::endl;
			return nullptr;
		}
		int markerId = bodyInfo.nodeIds.at(nodeId);
		return AddConstraintTargetRefToMarker(id, xmlConstraintInstance, markerId);
	}
	return nullptr;
}

TiXmlElement* ExporterNxXml::AddConstraintTargetRefToMarker(int& id, TiXmlElement* xmlConstraintInstance, int markerId) const
{
	TiXmlElement* xmlConstraintTargetRef = AddElement(id, "ConstraintTargetRef", xmlConstraintInstance);
	xmlConstraintTargetRef->SetAttribute("targetRef", GetRefId(markerId).c_str());
	return xmlConstraintTargetRef;
}


void ExporterNxXml::MarkerInfo::AddTransformToXML() const
{
	_xmlTransform->Clear();
	_xmlTransform->LinkEndChild(new TiXmlText(MatVecToString(_rotation, _translation).c_str()));
}
