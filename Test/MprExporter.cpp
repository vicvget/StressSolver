#include "MprExporter.h"

void MprExporter::Init()
{
	MultiphysicsResultsHeader mprHeader
		(
			SolverTypes::ST_StressStrain,
			_solver->_nElements,
			_integrationParams._initialTime,
			_integrationParams.TotalTime()
		);

	string fileResults = _solver->_uid + ".mpr";
	_writer.InitWriter(fileResults, &mprHeader);
	SetInitialized();
}

void MprExporter::WriteFrame(float time)
{
	float* data = _solver->GetMemoryPointer();
	int dataSize = _solver->GetMemorySize();
	_writer.WriteFrame(data, dataSize, time);

}

void MprExporter::Finalize()
{
}
