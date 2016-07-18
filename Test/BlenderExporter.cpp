#include "BlenderExporter.h"

bool BlenderExporter::Init(const string& filename, size_t nelements, size_t stride, double* coordinates, double size)
{
	_nelements = nelements;
	_stride = stride;
	_stride2 = _stride * 2;
	_size = size;
	_matStride = stride * 3;
	_frame = 0;
	_coordinates = coordinates;

	_ofs.open(filename);
	return _ofs.is_open();
}

BlenderExporter::~BlenderExporter()
{
	Close();
}

void BlenderExporter::Close()
{
	if (_ofs.is_open())
		_ofs.close();
}

void BlenderExporter::Init()
{
	Init("blender.py", _solver->_nElements, _solver->vecStride, _solver->GetElementShift(0), _solver->_gridStep);
	WriteHeader();
	WriteBody();
}

void BlenderExporter::WriteFrame()
{
	AddFrame(_solver->GetElementShift(0), _solver->GetDataRotaionMtx());
}

void BlenderExporter::Finalize()
{
	WriteFooter();
	Close();
}

void BlenderExporter::AddFrame(double* coordinates, double* rotationMatrices)
{
	_frame++;
	_ofs
		<< std::endl << "# frame " << _frame << std::endl;
	for (size_t elementId = 0; elementId < _nelements; elementId++)
	{
		_ofs
			<< "obj = objects[" << elementId << "]" << std::endl
			<< GenerateFrameForElement(elementId, coordinates + _stride2 * elementId, rotationMatrices + _matStride * elementId)
			<< "obj.matrix_world = Matrix(m)" << std::endl
			<< "obj.keyframe_insert(data_path = \"rotation_euler\", frame = " << _frame << ")" << std::endl
			<< "obj.keyframe_insert(data_path = \"location\", frame = " << _frame << ")" << std::endl;
	}
	//_frame++;
}

string BlenderExporter::PrintVertices() const
{
	stringstream ss;
	double h = _size*0.3;
	for (size_t vertexId = 0; vertexId < 8; vertexId++)
	{
		ss << '(';
		switch (vertexId)
		{
		case 0:
			ss << -h << ',' << -h << ',' << -h;
			break;
		case 1:
			ss << -h << ',' << h << ',' << -h;
			break;
		case 2:
			ss << h << ',' << h << ',' << -h;
			break;
		case 3:
			ss << h << ',' << -h << ',' << -h;
			break;
		case 4:
			ss << -h << ',' << -h << ',' << h;
			break;
		case 5:
			ss << -h << ',' << h << ',' << h;
			break;
		case 6:
			ss << h << ',' << h << ',' << h;
			break;
		case 7:
			ss << h << ',' << -h << ',' << h;
			break;
		}
		ss << ')';
		if (vertexId < 7) ss << ',';
	}
	return ss.str();
}

void BlenderExporter::WriteHeader()
{
	_ofs
		<< "import bpy" << std::endl
		<< "import math" << std::endl
		<< "from mathutils import Matrix" << std::endl
		<< std::endl
		<< "def remove_all_mesh_objects() :" << std::endl
		<< "   # gather list of items of interest." << std::endl
		<< "   candidate_list = [item.name for item in bpy.data.objects if item.type ==\"MESH\"]" << std::endl
		<< "   # select them only." << std::endl
		<< "   for object_name in candidate_list:" << std::endl
		<< "     bpy.data.objects[object_name].select = True" << std::endl
		<< "   # remove all selected." << std::endl
		<< "   bpy.ops.object.delete()" << std::endl
		<< "   # remove the meshes, they have no users anymore." << std::endl
		<< "   for item in bpy.data.meshes:" << std::endl
		<< "     bpy.data.meshes.remove(item)" << std::endl
		<< std::endl
		<< "def make_cube(x, y, z, name) :" << std::endl
		<< "   verts = [" << PrintVertices() << "]" << std::endl
		<< "   faces = [(0,1,2,3), (0,1,5,4), (1,2,6,5), (2,3,7,6), (3,0,4,7), (4,5,6,7)]" << std::endl
		<< std::endl
		<< "   me = bpy.data.meshes.new(name)" << std::endl
		<< "   obj = bpy.data.objects.new(name, me)" << std::endl
		<< "   bpy.context.scene.objects.link(obj)" << std::endl
		<< "   me.from_pydata(verts, [], faces)" << std::endl
		<< "   me.update(calc_edges=True)" << std::endl
		<< "   obj.location = (x,y,z)" << std::endl
		<< "   return obj" << std::endl
		<< std::endl
		<< "remove_all_mesh_objects()" << std::endl;
}

void BlenderExporter::WriteBody()
{
	_ofs << "objects = []" << std::endl;
	for (size_t elementId = 0; elementId < _nelements; elementId++)
	{
		double* shift = _coordinates + elementId * _stride2;
		_ofs << "objects.append(make_cube("
			<< shift[0] << ", "
			<< shift[1] << ", "
			<< shift[2] << ", \"e"
			<< elementId << "\"))" << std::endl;
	}
}

void BlenderExporter::WriteFooter()
{
	_ofs << "scene = bpy.context.scene" << std::endl
		<< "scene.frame_end = 1" << std::endl
		<< "scene.frame_start = 1" << std::endl
		<< "scene.frame_end = " << _frame+1 << std::endl;
}

string BlenderExporter::GenerateFrameForElement(size_t elementId, double* shift, double* rotationMatrix) const
{
	stringstream ss;
	ss << "m=[";
	for (size_t row = 0; row < 3; row++)
	{
		double* matrixRow = rotationMatrix + row*_stride;
		ss << "(";
		for (size_t element = 0; element < 3; element++)
		{
			ss << ((float)matrixRow[element]) << ", ";
		}
		ss << ((float)shift[row]) << ")," << std::endl;
	}
	ss << "(0, 0, 0, 1)]" << std::endl;
	return ss.str();
}
