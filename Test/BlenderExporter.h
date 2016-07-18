#pragma once

#include <fstream>
#include <sstream>
#include "BaseExporter.h"

using namespace std;

class BlenderExporter: public BaseExporter
{
private:
	ofstream _ofs;
	double _size;
	size_t _nelements;
	size_t _stride;
	size_t _stride2;
	size_t _matStride;
	size_t _frame;
	double* _coordinates;

public:
	BlenderExporter(Stress::StressStrainCppIterativeSolver* solver) : BaseExporter(solver) {};
	~BlenderExporter();

	bool Init(
		const string& filename,
		size_t nelements,
		size_t stride,
		double* coordinates,
		double size);

	void WriteHeader();
	void WriteBody();
	void AddFrame(double* coordinates, double* rotationMatrices);
	void WriteFooter();
	void Close();

#pragma region overriden
	virtual void Init();
	virtual void WriteFrame();
	virtual void Finalize();
#pragma endregion overriden

private:
	string GenerateFrameForElement(size_t elementId, double* shift, double* rotationMatrix) const;
	string PrintVertices() const;
};

