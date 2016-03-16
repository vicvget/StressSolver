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
* Функция преобразовывает сетку в массив указателей объектов типа Vertex
* @return возвращет 1 если преобразование успешно
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
	//превращение сжатой сетки в массив вершин 
	for(int z=0;z<nz;z++)
	{
		for(int y=0;y<ny;y++)
		{
			for(int x=0;x<nx;x++) 
			{
				int tmp=(*_grid)[x][y][z];
				if(tmp>0)
				{
					//РІС‹С‡РёСЃР»РёС‚СЊ РєРѕРѕСЂРґРёРЅР°С‚С‹
					pntx = xs+x*_dx;
					pnty = ys+y*_dy;
					pntz = zs+z*_dz;
					//создать вершину
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
	//граничные условия
	//bc = _bc;
	return 1;
}

/**
* Функция возвращает указателей объектов типа Vertex
* @return
*/
vector<Vertex *> MeshDataGraphFormatter::GetArrayOfVertex()
{
	return graph;
}
/**
* Функция возвращает граничные условия для сетки
* @return вектор указателей на объекты типа BoundaryCondition
*/
vector<BoundaryCondition *> MeshDataGraphFormatter::GetArrayOfBC()
{
	return bc;
}