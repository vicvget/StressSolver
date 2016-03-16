#include "MeshDataGraphFormatter.h"
#include <algorithm>
MeshDataGraphFormatter::MeshDataGraphFormatter(OccRectilinearGrid* newgrid):MeshDataProvider(newgrid)
{
}

MeshDataGraphFormatter::MeshDataGraphFormatter(OccRectilinearGrid *newgrid,vector<BoundaryCondition *> Surfaces):MeshDataProvider(newgrid,Surfaces)
{
}

MeshDataGraphFormatter::~MeshDataGraphFormatter(void)
{
}

/**
* ������� ��������������� ����� � ������ ���������� �������� ���� Vertex
* @return ��������� 1 ���� �������������� �������
*/
int MeshDataGraphFormatter::Vectorize()
{
double pntx,pnty,pntz;

	int nx,ny,nz;
	double xs,ys,zs;

	EnumerationMesh();

	_grid->GetGridSize(&nx,&ny,&nz);
	_grid->GetStartPoint(&xs,&ys,&zs);

	graph.reserve(nx*ny*nz);
	//����������� ������ ����� � ������ ������ 
	for(int z=0;z<nz;z++)
	{
		for(int y=0;y<ny;y++)
		{
			for(int x=0;x<nx;x++) 
			{
				int tmp=(*_grid)[x][y][z];
				if(tmp>0)
				{
					//вычислить координаты
					pntx = xs+x*_dx;
					pnty = ys+y*_dy;
					pntz = zs+z*_dz;
					//������� �������
					Vertex *tmpVert = new Vertex();
					tmpVert->pt.x = pntx;
					tmpVert->pt.y = pnty;
					tmpVert->pt.z = pntz;
					tmpVert->_globalNumber = tmp;
					try
					{
						if(z<nz-1)
						if((*_grid)[x][y][z+1]>0)
						{
							tmpVert->adjVert.push_back((*_grid)[x][y][z+1]);
						}
						if(y<ny-1)
						if((*_grid)[x][y+1][z]>0)
						{
							tmpVert->adjVert.push_back((*_grid)[x][y+1][z]);
						}
						if(x<nx-1)
						if((*_grid)[x+1][y][z]>0)
						{
							tmpVert->adjVert.push_back((*_grid)[x+1][y][z] );
						}
						if(z>0)
						if((*_grid)[x][y][z-1]>0)
						{
							tmpVert->adjVert.push_back((*_grid)[x][y][z-1]);
						}
						if(y>0)
						if((*_grid)[x][y-1][z]>0)
						{
							tmpVert->adjVert.push_back((*_grid)[x][y-1][z]);
						}
						if(x>0)
						if((*_grid)[x-1][y][z]>0)
						{
							tmpVert->adjVert.push_back((*_grid)[x-1][y][z]);
						}
					}
					catch (...)
					{
						return -1;
					}
					std::sort(tmpVert->adjVert.begin(),tmpVert->adjVert.end());
					graph.push_back(tmpVert);
				}
			}
		}
	}
	//��������� �������
	//bc = _bc;
	return 1;
}

/**
* ������� ���������� ���������� �������� ���� Vertex
* @return
*/
vector<Vertex *> MeshDataGraphFormatter::GetArrayOfVertex()
{
	return graph;
}
/**
* ������� ���������� ��������� ������� ��� �����
* @return ������ ���������� �� ������� ���� BoundaryCondition
*/
vector<BoundaryCondition *> MeshDataGraphFormatter::GetArrayOfBC()
{
	return bc;
}