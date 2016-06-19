#include "CommonSolversTest.h"

#include "CommonSolvers.h"
#include "../FormatProviders/GridProvider/RLCControlWriter.h"


using std::string;


namespace SpecialSolversTest
{

	using namespace SpecialSolvers;


	// Вспомогательные функции создания сетки

	void CreateTestGrid
		(
			double*& nodes,
			int*& links,
			int& nLinks,
			int nx,
			int ny,
			int nz,
			double step,
			const string& fileName
		)
	{
		nodes = new double [nx * ny * nz * 3];
		nLinks = (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1);
		links = new int [2 * nLinks];

		int nyz = nz * ny;

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int k = 0; k < nz; k++)
				{
					int offset = (nyz * i + nz * j + k) * 3;

					nodes[offset] = step * i;
					nodes[offset + 1] = step * j;
					nodes[offset + 2] = step * k;
				}
			}
		}

		int linksId = 0;

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int k = 0; k < nz; k++)
				{
					int nNode = nyz * i + nz * j + k;

					if (k < nz - 1)
					{
						links[linksId++] = nNode;
						links[linksId++] = nNode + 1;
					}
					if (j < ny - 1)
					{
						links[linksId++] = nNode;
						links[linksId++] = nNode + nz;
					}
					if (i < nx - 1)
					{
						links[linksId++] = nNode;
						links[linksId++] = nNode + nyz;
					}
				}
			}
		}

		// make RLC
		OccRectilinearGrid occRectilinearGrid(nx, ny, nz); //nx,ny,nz

		occRectilinearGrid.SetStartPoint(0., 0., 0.);
		step *= 1000.0;
		occRectilinearGrid.SetElementSize(step, step, step);
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int k = 0; k < nz; k++)
				{
					occRectilinearGrid[i][j][k] = 1;
				}
			}
		}

		RLCControlWriter rlcControlWriter(&occRectilinearGrid);

		rlcControlWriter.DumpMeshToFile(fileName.c_str()); 
	}

	void CreateTestGrid
		(
			const GridParams& gridParams,
			const string& rlcGridFileName
		)
	{
		double* nodes = nullptr;
		int* links = nullptr;
		int nLinks;

		CreateTestGrid
			(
				nodes,
				links,
				nLinks,
				gridParams._nx,
				gridParams._ny,
				gridParams._nz,
				gridParams._gridStep,
				rlcGridFileName
			);

		delete[] nodes;
		delete[] links;
	}

}