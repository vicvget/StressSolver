#include "Assignment.h"

#include "../../../fcore/wrappers/DirectoryLister.h"
#include "../../ResProvider/OldAddresses.h"

#include <cmath>


using std::endl;


Assignment::Assignment()
	:
		newGeometry(nullptr)
{
	
}

Assignment::~Assignment()
{
	delete newGeometry;
}

int Assignment::Assign
	(
		const string& geoFileName,
		int bodyNumber
	)
{
	ifstream istream(geoFileName);

	if (!istream.is_open())
	{
		exceptions::ThrowFileNotOpened(geoFileName);
		
		// невозможно открыть файл с новой геометрией тела
		return 1;
	}
	newGeometry = new FGeometry(istream);
	istream.close();
	
	string modelFileName;

	if (FindModelFile(modelFileName) != 0)
	{
		//файл модели "*.frm" не был найден
		return 2;
	}

	JointChars jointChars;
	Joints joints;

	if (Init(modelFileName, bodyNumber, jointChars, joints) != 0)
	{
		//проблемы при инициализации данных о теле и соединительных элементах
		return 3;
	}

	/*
	// преобразование в координат геометрии в соответствии с начальными условиями
	newGeometry->GlobalTransform
		(
			body.geometry->GetTransformMtx()
		);
	
	// преобразование в координат геометрии в соответствии с начальными условиями
	body.geometry->GlobalTransform
		(
			body.geometry->GetTransformMtx()
		);
	*/
	
	oldAddresses.ReadFadres();

	if (FormSForces(jointChars, joints) != 0)
	{
		//проблемы при формировании файла "sforces.dat"
		return 4;
	}

	return 0;
}

int Assignment::FindModelFile
	(
		string& modelFileName
	)
{
	fs::DirectoryLister dl(".", "*.frm");
	string curPath;
	string path;

	dl.Init();
	while (!dl.IsLast())
	{
		dl.Next(curPath);
		if (curPath.find("backup") == string::npos)
		{
			path = curPath;
			break;
		}
	}
	if (path.empty())
	{
		return 1;
	}
	modelFileName = path;

	return 0;
}

int Assignment::Init
	(
		const string& modelFileName,
		int bodyNumber,
		JointChars& jointChars,
		Joints& joints
	)
{
	if (!fs::FileExist(modelFileName.c_str()))
	{
		exceptions::ThrowFileNotOpened(modelFileName);

		// невозможно открыть файл с моделью
		return 1;
	}
	
	FModel model(modelFileName.c_str());

	_body = model.GetBodyByNumber(bodyNumber);
	_geometry = model.GetBodyGeometry(bodyNumber);
	joints = model.Joints();
	for (const FJoint& joint: joints)
	{
		jointChars.emplace(joint.CharNumber(), model.GetJointCharByNumber(joint.CharNumber()));
	}
/*
	int type = 0;
	string index;
	bool bodyFound = false;

	istream >> type >> index;
	while(!istream.eof())
	{
		index = index.append(".elf");
		switch(type)
		{
		case 1: // body
			{
				if (!bodyFound)
				{
					FBody tempBody(index.c_str());

					if (tempBody.Number() == bodyNumber)
					{
						body = tempBody;
						bodyFound = true;
					}
				}
			}
			break;

		case 2: // joint char
			{
				FJointChar jointChar(index.c_str());
				jointChars.insert(JointCharPair(jointChar.Number(), jointChar));
			}
			break;

		case 5: // joint
			{
				FJoint joint(index.c_str());
				joints.push_back(joint);
			}
			break;
			
		default:
			break;
		}
		istream >> type >> index;		
	}
	istream.close();
	*/
	return 0;
}

float FindDistance
	(
		const FGeometryPoint& point1,
		const FGeometryPoint& point2
	)
{
	return
		pow
		(
			pow(point1.x - point2.x, 2) +
			pow(point1.y - point2.y, 2) +
			pow(point1.z - point2.z, 2),
			0.5f
		);
}

void Assignment::FindNearestPoint
	(
		const FGeometryPoint& point,
		int& nearestPointNumber,
		FGeometryPoint& nearestPoint
	)	const
{
	float minDistance;

	minDistance = pow(10.0f, 7.0f);
	for
		(
			vector<FGeometryPoint>::const_iterator it = newGeometry->Nodes().begin();
			it != newGeometry->Nodes().end();
			++it
		)
	{
		float distance = FindDistance(point, *it);

		if (distance < minDistance)
		{
			nearestPointNumber = (int)(it - newGeometry->Nodes().begin()) + 1;
			nearestPoint = *it;
			minDistance = distance;
		}
	}
}

void Assignment::TransformCoordinates
	(
		const FGeometryPoint& cmNode,
		const FGeometryPoint& point,
		const float (&transformMatrix)[12],
		FGeometryPoint& transformedPoint
	)
{
	const float* shiftVector;
	const float* rotationMatrix;

	shiftVector = transformMatrix;
	rotationMatrix = transformMatrix + 3;
	for(int i = 0; i < 3; i++)
	{
		transformedPoint.coords[i] = 0;
		for(int j = 0; j < 3; j++)
		{
			transformedPoint.coords[i] += (point.coords[j] - cmNode.coords[j]) * rotationMatrix[j * 3 + i];
		}
		transformedPoint.coords[i] += 1000 * shiftVector[i];
	}
}

int Assignment::FormSForces
	(
		const JointChars& jointChars,
		const Joints& joints
	)
{
	char fileName[100];

	sprintf(fileName, "sforces_%d.dat", _body.Number());

	ofstream ostream(fileName);

	if(!ostream.is_open())
	{
		exceptions::ThrowFileNotOpened(fileName);
		
		// невозможно открыть файл с моделью
		return 1;
	}

	FGeometryPoint node;
	int nearestPointNumber;
	FGeometryPoint nearestPoint;
	FGeometryPoint transformedPoint;

	JointChars::const_iterator jt;
	Joints connectedJoints;

	for (Joints::const_iterator it = joints.begin(); it != joints.end(); ++it)
	{
		if
			(
				(it->BodyNumber1() == _body.Number()) ||
				(it->BodyNumber2() == _body.Number())
			)
		{
			connectedJoints.push_back(*it);
		}
	}
	ostream << connectedJoints.size() << endl;
	for
		(
			Joints::const_iterator it = connectedJoints.begin();
			it != connectedJoints.end();
			++it
		)
	{
		if (it->BodyNumber1() == _body.Number())
		{
			node = (*it).Node1();
		}
		else
		{
			node = (*it).Node2();
		}
		FindNearestPoint
			(
				node,
				nearestPointNumber,
				nearestPoint
			);
		TransformCoordinates
			(
				_geometry.CmNode(),
				nearestPoint,
				_geometry.GetTransformArray(),
				transformedPoint
			);
		
		int charNumber = (*it).CharNumber();

		jt = jointChars.find(charNumber);
		if (jt == jointChars.end())
		{
			// ошибка при поиске характеристики с.э.
			return 2;
		}

		ostream << (*it).Number() << endl;
		
		const vector<CharComponent>& springParams = jt->second.GetSpringParams();
		
		ostream
			<< transformedPoint.x << ' '
			<< transformedPoint.y << ' '
			<< transformedPoint.z << endl;
		ostream << springParams.size() << endl;
		for
			(
				vector<CharComponent>::const_iterator kt = springParams.begin();
				kt != springParams.end();
				++kt
			)
		{
			ostream
				<< (*kt).direction << ' '
				<< oldAddresses.GetConnectionElementAddress
					(
						(*it).Number(),
						(FreedomType)(*kt).direction,
						CharacteristicComponent::Strain
					)
					+ 1
				<< endl;
		}
		ostream << endl;
	}
	ostream.close();

	return 0;
}
