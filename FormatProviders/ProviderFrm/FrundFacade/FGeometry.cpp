#include "FGeometry.h"

#include "TransformMatrix.h"

#define _USE_MATH_DEFINES

#include <cmath>


using std::setw;


FGeometry::FGeometry(ifstream& stream)
{
	this->Load(stream);
}

FGeometry::FGeometry(ifstream& stream, float scale)
{
	this->Load(stream, scale);
}


// Генерация окружностей (в массивы points и links)
void FGeometry::GenerateAllNodes()
{
	_transformedNodes.assign(_nodes.begin(), _nodes.end());
	_transformedLinks.assign(_links.begin(), _links.end());
	for (int k = 0; k < _circles.size(); k++)
	{
		FGeometryCircle& circle = _circles[k];
		circle.nPoints = circle.nPoints > 100 ? 100 : circle.nPoints;
		circle.nPoints = circle.nPoints < 3 ? 3 : circle.nPoints;

		FGeometryPoint center;
		FGeometryPoint axis;
		FGeometryPoint sqAxis;

		if (circle.code == 3)
		{
			center = _nodes[circle.axisNode - 1];
			axis = center;
			switch (circle.centerNode)
			{
			case 1:
				axis.z += 1000;
				break;

			case 2:
				axis.x += 1000;
				break;

			case 3:
				axis.y += 1000;
				break;
			}
		}
		else
		{
			center = _nodes[circle.centerNode - 1];
			axis = _nodes[circle.axisNode - 1];
		}

		for (int i = 0; i < 3; i++)
		{
			axis.coords[i] = -center.coords[i] + axis.coords[i];
		}

		float axisMod = sqrt(axis.x*axis.x + axis.y*axis.y + axis.z*axis.z);
		if (axisMod < 1E-10)
		{
			// TODO:
			//throw("нулевая длина нормали к плоскости окружности");
			axis.coords[0] = axis.coords[1] = 0;
			axis.coords[2] = 1;
			axisMod = 1;
		}
		for (int i = 0; i < 3; i++)
		{
			axis.coords[i] /= axisMod;
			sqAxis.coords[i] = axis.coords[i] * axis.coords[i];
		}

		float rotationMtx[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
		if (fabs(fabs(axis.x) - 1) > 1E-10)
		{
			rotationMtx[0] = axis.x;
			rotationMtx[1] = axis.y;
			rotationMtx[2] = axis.z;
			rotationMtx[3] = -rotationMtx[1];
			rotationMtx[6] = -rotationMtx[2];

			rotationMtx[8] = (sqAxis.y + sqAxis.z * axis.x) / (1 - sqAxis.x);
			rotationMtx[4] = rotationMtx[8] * axis.x + sqAxis.z;
			rotationMtx[5] = -axis.y*axis.z / (axis.x + 1);
			rotationMtx[7] = rotationMtx[5];
		}

		float deltaAngle = (float)(360. / circle.nPoints);
		for (int k = 0; k < circle.nPoints; k++)
		{
			FGeometryPoint oldPoint;
			float angle = (float)(k*deltaAngle*asin(1.) / 90.);
			oldPoint.x = 0;
			oldPoint.y = circle.Radius*cos(angle);
			oldPoint.z = circle.Radius*sin(angle);

			FGeometryPoint p;
			for (int i = 0; i < 3; i++)
			{
				p.coords[i] = 0;
				for (int j = 0; j < 3; j++)
				{
					p.coords[i] += oldPoint.coords[j] * rotationMtx[j * 3 + i];
				}
				p.coords[i] += center.coords[i];
			}
			_transformedNodes.push_back(p);
			if (k > 0)
			{
				_transformedLinks.push_back(*(new FGeometryLink(
					_transformedNodes.size() - 1, _transformedNodes.size())));
			}
		}
		_transformedLinks.push_back(*(new FGeometryLink(
			_transformedNodes.size(), _transformedNodes.size() - circle.nPoints + 1)));
	}
}

void FGeometry::GlobalTransform
	(
		const float (&transform)[12]
	)
{
	// global transformation
	const float* rotationMtx = transform + 3;
	for(int k = 0; k < _transformedNodes.size(); k++)
	{
		FGeometryPoint oldPoint = _transformedNodes[k];
		FGeometryPoint p;
		for(int i = 0; i < 3; i++)
		{
			p.coords[i] = 0;
			for(int j = 0; j < 3; j++)
			{
				p.coords[i] += oldPoint.coords[j]*rotationMtx[j*3+i];
			}
			_transformedNodes[k].coords[i] = p.coords[i] + transform[i];
		}		
	}
}

// Сохранение геометрии в .lba формат
void FGeometry::OutputLba(ofstream& stream) const
{
	stream << _transformedNodes.size() + 1 << ' ' <<  
		_transformedNodes.size() - 1 << " 0 " <<
		_transformedLinks.size() << " 0 0\n";
	int l = 12;
	stream << setw(4) << _cmNodeNumber+1 << ' ' <<
		setw(l) << _transformedNodes[_cmNodeNumber].x << ' ' <<
		setw(l) << _transformedNodes[_cmNodeNumber].y << ' ' <<
		setw(l) << _transformedNodes[_cmNodeNumber].z << '\n';

	for(int i = 0; i < _transformedNodes.size(); i++)
	{
		if(i != _cmNodeNumber)
		{
			stream << setw(4) << i+1 << ' ' <<
				setw(l) << _transformedNodes[i].x << ' ' <<
				setw(l) << _transformedNodes[i].y << ' ' <<
				setw(l) << _transformedNodes[i].z << '\n';		
		}
	}
	stream  << setw(4) << _transformedNodes.size()+1 << ' ' << 
		setw(l) << 0.0 << ' ' <<
		setw(l) << 0.0 << ' ' <<
		setw(l) << 0.0 << '\n';

	for(int i = 0; i < _transformedLinks.size(); i++)
	{
		stream << setw(4) << _transformedLinks[i].node1 << ' ' <<
			setw(4) << _transformedLinks[i].node2 << '\n';
	}
}

void FGeometry::Load(ifstream& stream, float scale)
{
	SetExt(EXT_BODY_GEOMETRY);
	for (int i = 0; i < 12; i++)
	{
		_transformMtx[i] = 0.0;
	}
	_transformMtx[3] = _transformMtx[7] = _transformMtx[11] = 1.0;
	_isTransformed = false;

	int nNodes;

	stream >> nNodes;
	_code60.clear();
	_nodes.clear();
	for (int i = 0; i < nNodes; i++)
	{
		FGeometryPoint node;

		stream >> node.x >> node.y >> node.z;
		node.x *= scale;
		node.y *= scale;
		node.z *= scale;
		_nodes.push_back(node);
	}

	int code;

	stream >> code;
	_links.clear();
	while (!(code == 300 || code == 999) && !stream.eof())
	{
		switch(code)
		{
		case 2:
			{
				FGeometryLink link;
				int linkColorId;

				stream >> link.node1 >> link.node2 >> linkColorId;
				_links.push_back(link);
			}
			break;

		case 3:
		case 4:
		case 5:
			{
				FGeometryCircle circle;

				stream
					>> circle.centerNode >> circle.axisNode
					>> circle.Radius >> circle.nPoints
					>> circle.angle1 >> circle.angle2;
				circle.code = code;
				_circles.push_back(circle);
			}
			break;

		case 100:
			{
				stream >> _cmNodeNumber;
				_cmNodeNumber--;
			}
			break;

		case 30:
			{
				float buf[12];

				for (int i = 0; i < 12; i++)
				{
					stream >> buf[i];
				}
				SetTransformArray(buf);
			}
			break;

		case 60:
			fs::ReadLineString(stream, _code60);
			break;

		default:
			{
				string tmp;

				fs::ReadLineString(stream, tmp);
			}
		}
		stream >> code;
	}
	GenerateAllNodes();
	//GlobalTransform(transformMtx); в LBA этого нет!
}

// Загрузка геометрии из .elg формата
bool FGeometry::Load(ifstream& stream)
{
	Load(stream, 1.0);

	return true;
}

void FGeometry::Load(const char* path)
{
	DefaultLoad(path);
}

void FGeometry::Save(ofstream& ofs) const
{
	ofs << Nodes().size() << std::endl;
	ofs.setf(ios_base::fixed);
	ofs.precision(6);
	for(int i = 0; i < Nodes().size(); i++)
		ofs << GetNode(i).x << ' ' 
			<< GetNode(i).y << ' ' 
			<< GetNode(i).z << std::endl;
	ofs.unsetf(ios_base::fixed);
	for(int i = 0; i < Links().size(); i++)
		ofs << "2 " << GetLink(i).node1
		<<	' ' << GetLink(i).node2 << " 4" << std::endl;
	for(int i = 0; i < Circles().size(); i++)
		ofs << GetCircle(i).code << ' '
		<< GetCircle(i).centerNode << ' ' 
		<< GetCircle(i).axisNode << ' ' 
		<< setiosflags(ios_base::fixed) << GetCircle(i).Radius << ' ' 
		<< resetiosflags(ios_base::fixed) << GetCircle(i).nPoints << ' ' 
		<< setiosflags(ios_base::fixed) << GetCircle(i).angle1 << ' ' 
		<< setiosflags(ios_base::fixed) << GetCircle(i).angle2 << std::endl;
	ofs << "100 " << _cmNodeNumber+1 << std::endl;
	if(_code60.size() > 0)	
	{
		ofs << "60" << _code60 << std::endl;
	}
	if(_isTransformed)
	{
		ofs << "30";
		for(int i = 0; i < 3; i++)
		{
			ofs << ' ' << _transformMtx[i];
		}
		ofs << std::endl << _transformMtx[3];
		for(int i = 4; i < 12; i++)
		{
			ofs << ' ' << _transformMtx[i];
		}
		ofs << std::endl;
	}	
	ofs << "999" << std::endl;
}

bool FGeometry::Combine(int* ids, FGeometryPoint* points, int count)	
{
	GenerateAllNodes(); // build accessable nodes
	if(count == 0)
	{
		return false;
	}

	TransformMatrix mTransform;
	TransformMatrix mRotate1;
	TransformMatrix mRotate2;

	if(count > 1)
	{
		FGeometryPoint vt01(points[0], points[1]);
		FGeometryPoint v01(_transformedNodes[ids[0]], _transformedNodes[ids[1]]);

		FGeometryPoint axis = v01.Cross(vt01);
		double angle = v01.Angle(vt01);
		if(IS_ZERO(axis.Magnitude()))
			axis = v01.Normal();
		mRotate2.SetRotation(axis, angle);
	}
	if(count > 2)
	{
		FGeometryPoint vt01(points[0], points[1]);
		FGeometryPoint vt02(points[0], points[2]);
		FGeometryPoint pt2 = _transformedNodes[ids[2]];
		FGeometryPoint pt1 = _transformedNodes[ids[1]];
		FGeometryPoint pt0 = _transformedNodes[ids[0]];
		pt0.Transform(mRotate2,CmNode());
		pt1.Transform(mRotate2,CmNode());
		pt2.Transform(mRotate2,CmNode());
		FGeometryPoint v02(pt0, pt2);

		FGeometryPoint normal = vt01.Cross(v02);
		FGeometryPoint normalTarget = vt01.Cross(vt02);

		FGeometryPoint axis = normal.Cross(normalTarget);
		double angle = normal.Angle(normalTarget);
		if(IS_ZERO(axis.Magnitude()))
			axis = vt01;
		
		mRotate1.SetRotation(axis, angle);
	}
	mTransform = mRotate1*mRotate2;
	FGeometryPoint point0Transformed = _transformedNodes[ids[0]];
	point0Transformed.Transform(mTransform, CmNode());
	FGeometryPoint vtr(point0Transformed,points[0]);
	mTransform.SetTranslation(vtr);

	mTransform.FillRepresentation(_transformMtx);
	_isTransformed = true;
	return true;
}

FGeometryPoint FGeometry::GetLocalTransformedNode(size_t id) const
{
	FGeometryPoint point = _transformedNodes[id];
	float mtx[12];

	GetTransformArray(mtx);
	point.Transform(mtx, CmNode());

	return point;
}
FGeometry::FGeometry()
{
	_isTransformed = false;
	TransformMatrix mtx;
	mtx.FillRepresentation(_transformMtx);
	_cmNodeNumber = -1;
	SetExt(EXT_BODY_GEOMETRY);
}


//static 
void FGeometry::AnglesToMtx(double a, double b, double c, double* mtx)
{
	double sinA = sin(a);
	double sinB = sin(b);
	double sinC = sin(c);
	double cosA = cos(a);
	double cosB = cos(b);
	double cosC = cos(c);
	
	mtx[0] =  cosB * cosC;
	mtx[1] = -sinA * sinB * cosC + cosA * sinC;
	mtx[2] = -cosA*sinB*cosC + sinA*sinC;
	mtx[3] =  -cosB*sinC;
	mtx[4] = -sinA * sinB*sinC + cosA * cosC;
	mtx[5] = cosA * sinB * sinC + sinA*cosC;
	mtx[6] =  sinB;
	mtx[7] =  -sinA*cosB;
	mtx[8] =  cosA*cosB;
}

void FGeometry::AddNode(const FGeometryPoint& node)
{
	_nodes.push_back(node);
}

void FGeometry::AddLink(const FGeometryLink& link)
{
	_links.push_back(link);
}

FTransform FGeometry::GetTransform() const
{
	double transform[12];
	for(int i = 0; i < 12; i++)
		transform[i] = _transformMtx[i];
	return FTransform(transform+3, transform);
}

Vec3 FGeometry::GetGeometryPoint(int pointId) const
{
	return MathHelpers::MakeVec3(_nodes[pointId].x*1e-3,_nodes[pointId].y*1e-3,_nodes[pointId].z*1e-3);
}

void FGeometry::GetBound(MathHelpers::Vec3& minPoint, MathHelpers::Vec3& maxPoint) const
{
	if(_nodes.size() == 0)
		return;
	
	double minX, minY, minZ, maxX, maxY, maxZ;
	minX = maxX = _nodes[0].x;
	minY = maxY = _nodes[0].y;
	minZ = maxZ = _nodes[0].z;

	for(int i = 1; i < _nodes.size(); i++)
	{
		if(_nodes[i].x < minX) minX = _nodes[i].x;
		else if(_nodes[i].x > maxX) maxX = _nodes[i].x;

		if(_nodes[i].y < minY) minY = _nodes[i].y;
		else if(_nodes[i].y > maxY) maxY = _nodes[i].y;	
	
		if(_nodes[i].z < minZ) minZ = _nodes[i].z;
		else if(_nodes[i].z > maxZ) maxZ = _nodes[i].z;	
	}
	
	minPoint = MathHelpers::MakeVec3(minX, minY, minZ)*1e-3;
	maxPoint = MathHelpers::MakeVec3(maxX, maxY, maxZ)*1e-3;
}

void FGeometry::GetPoints(std::vector<MathHelpers::Vec3>& points, bool toGCS/*=false*/) const
{
	const FTransform& tr = GetTransform();
	points.clear();
	for(int i = 0; i < CountNodes(); i++)
	{
		Vec3 point = MathHelpers::MakeVec3(_nodes[i].x, _nodes[i].y, _nodes[i].z);
		point = point * 0.001;
		if(toGCS)
			point = tr.ToGCS(point);
		points.push_back(point);
	}
}
