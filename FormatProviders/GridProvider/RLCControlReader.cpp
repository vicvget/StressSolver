#include "RLCControlReader.h"

#include <stdexcept>


using std::cout;


RLCControlReader::RLCControlReader(void)
{
}

RLCControlReader::~RLCControlReader(void)
{
}

OccRectilinearGrid* RLCControlReader::fillGrid(ifstream *fs)
{
	try
	{
		//чтение точек
		int nx, ny, nz;
		//чтение количества узлов
		double xs, ys, zs;
		fs->read(reinterpret_cast<char*>(&_nodesCount), sizeof(int));
		fs->read(reinterpret_cast<char*>(&xs), sizeof(double));
		fs->read(reinterpret_cast<char*>(&ys), sizeof(double));
		fs->read(reinterpret_cast<char*>(&zs), sizeof(double));
		//чтение заголовка
		//чтение количества ячеек по осям
		fs->read(reinterpret_cast<char*>(&nx), sizeof(int));
		fs->read(reinterpret_cast<char*>(&ny), sizeof(int));
		fs->read(reinterpret_cast<char*>(&nz), sizeof(int));
		//чтение размера ячеек
		fs->read(reinterpret_cast<char*>(&_dx), sizeof(double));
		fs->read(reinterpret_cast<char*>(&_dy), sizeof(double));
		fs->read(reinterpret_cast<char*>(&_dz), sizeof(double));

		// если nz и ny бостаточно большие, а nx равно нулю, то получаем оооочень долгий цикл
		// если кол-во элементов вдоль какой-либо из осей равно 0, то нет смысла дальше продолжать работу
		if (nx != 0 && ny != 0 && nz != 0)
		{
			_grid = new OccRectilinearGrid(nx, ny, nz);
			_grid->SetStartPoint(xs, ys, zs);
			_grid->SetElementSize(_dx, _dy, _dz);
			//порядок обхода x y z
			bool tmp;
			for (int z = 0; z < nz; z++)
			{
				for (int y = 0; y < ny; y++)
				{
					for (int x = 0; x < nx; x++)
					{
						fs->read(reinterpret_cast<char*>(&tmp), sizeof(bool));
						(*_grid)[x][y][z] = tmp ? 1 : 0;
					}
				}
			}

			//чтение граничных условий
			int sizes = 0;
			fs->read(reinterpret_cast<char*>(&sizes), sizeof(int));

			for (int i = 0; i < sizes; i++)
			{
				BoundaryCondition *face = new BoundaryCondition();
				int tmp;
				char buf[256];
				string id;
				fs->read((char*)&tmp, sizeof(int));
				fs->read((char*)&buf, tmp);
				buf[tmp] = 0;
				face->SetID(buf);

				string var;
				fs->read((char*)&tmp, sizeof(int));
				fs->read((char*)&buf, tmp);
				buf[tmp] = 0;
				face->SetVar(buf);

				string name;
				fs->read((char*)&tmp, sizeof(int));
				fs->read((char*)&buf, tmp);
				buf[tmp] = 0;
				face->SetName(buf);

				int numBP = 0;
				fs->read(reinterpret_cast<char*>(&numBP), sizeof(int));
				for (int j = 0; j < numBP; j++)
				{
					int pnt;
					fs->read(reinterpret_cast<char*>(&pnt), sizeof(int));
					face->AddPoint(pnt);
				}
				_surfaces.push_back(face);
			}
		}
	}
	catch (const std::exception&)
	{
		cout << "Cant read file rlc" << std::endl;
		delete _grid;

		return NULL;
	}
	catch (...)
	{
		cout << "Cant read file rlc" << std::endl;
		delete _grid;

		return NULL;
	}
	//cout << "Read file rlc completed" << std::endl;
	return _grid;
}