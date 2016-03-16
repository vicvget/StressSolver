#pragma once
#include <fstream>
#include <iostream>

using namespace std;

class BlenderExporter
{
private:
	ofstream ofs;

public:
	bool Init(const string& filename)
	{
		ofs.open(filename);
		return ofs.is_open();
	}

	void Close()
	{
		if (ofs.is_open())
			ofs.close();
	}

	void WriteHeader()
	{
		ofs
			<< "import bpy" << std::endl
			<< "import math" << std::endl
			<< "from mathutils import Matrix" << std::endl
			<< std::endl
			<< "def remove_all_mesh_objects() :" << std::endl
			<< "   # gather list of items of interest." << std::endl
			<< "   candidate_list = [item.name for item in bpy.data.objects if item.type ==" << std::endl"MESH\"]" << std::endl
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
			<< "   verts = [(-1, -1, -1), (-1, 1, -1), (1, 1, -1), (1, -1, -1), (-1, -1, 1), (-1, 1, 1), (1, 1, 1), (1, -1, 1)]" << std::endl
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

	void WriteBody()
	{
	}

	void WriteFooter()
	{
	}

obj = make_cube(1,1,1, "cube1")
obj1 = make_cube(3,1,1, "cube2")
#obj.select = True
#obj.rotation_mode = 'XYZ'
scene = bpy.context.scene
scene.frame_start = 1
scene.frame_end = 100

#obj.animation_data_clear()
scene.animation_data_clear()
#matrix = Matrix([1,1,1],[0,1,0],[1,0,0])

for x in range(scene.frame_start,scene.frame_end+1,1):
    angle = 2*math.pi/(scene.frame_end+1 - scene.frame_start)*x
    obj.matrix_world = Matrix.Translation(obj.location) * Matrix.Rotation(angle, 4, 'Z')    
    obj.keyframe_insert(data_path = "rotation_euler", frame = x)
    obj.keyframe_insert(data_path = "location", frame = x)
	


};