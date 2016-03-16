#include "ProviderMpr.h"


void ProviderMpr::InitWriter(const string& mprFileName,
						const MultiphysicsResultsHeader* mprHeader)
{
	_resultsWriter = std::make_unique<ResultsWriter>(mprFileName);
	_resultsWriter->WriteBuffer(mprHeader, sizeof(MultiphysicsResultsHeader));
}

void ProviderMpr::InitReader(const string& mprFileName,
	MultiphysicsResultsHeader* mprHeader)
{
	_resultsReader = std::make_unique<ResultsRawReader>(mprFileName);
	_resultsReader->ReadBuffer(mprHeader, sizeof(MultiphysicsResultsHeader));
}

void ProviderMpr::WriteFrame
	(
		const void* data,
		size_t frameSize,
		float currentTime
	)
{
	MultiphysicsFrameHeader mprFrameHeader(frameSize, currentTime);

	_resultsWriter->WriteBuffer(&mprFrameHeader, sizeof(MultiphysicsFrameHeader));
	_resultsWriter->WriteBuffer(data, frameSize);
}

bool ProviderMpr::ReadFrame(void* data, size_t& dataSize, int& frameNumber)
{
	int frameId = 0;

	MultiphysicsFrameHeader frameHeader;
	if (!_resultsReader->ReadBuffer(&frameHeader, sizeof(MultiphysicsFrameHeader)))
		return false;

		
	while (_resultsReader->ReadBuffer(data, frameHeader.FrameSize))
	{
		dataSize = frameHeader.FrameSize;
		if (!_resultsReader->ReadBuffer(&frameHeader, sizeof(MultiphysicsFrameHeader)))
			break;
		frameId++;
		if (!(frameNumber < 0))
			if (frameId == frameNumber)
				break;
	}
	frameNumber = frameId;
	return true;
}

MultiphysicsResultsHeader::MultiphysicsResultsHeader
	(
		SolverTypes type,
		size_t pointsCount,
		float initialTime,
		float totalTime
	)
{
	FormatSymbol = 'M';
	PointsCount = pointsCount;
	InitialTime = initialTime;
	FinalTime = totalTime;
	switch (type)
	{
		case ST_Thermal:
			Type = (char)type;
			ScalarsCount = 1;
			VectorsCount = 0;
			break;

		case ST_StressStrain:
			Type = (char)type;
			ScalarsCount = 1;
			VectorsCount = 1;
			break;
	}
}