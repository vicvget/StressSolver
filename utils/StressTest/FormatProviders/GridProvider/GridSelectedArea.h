#ifndef GRID_SELECTED_AREA_H

#define GRID_SELECTED_AREA_H


#include "MeshDataProvider.h"

#include <vector>
#include <array>


/**
* ���������� ������� ����� � ���� ���������������
*/
class GridSelectedArea final
{
public:

	// ���� ������

	// ������ �������� ����� �����
	using PointList = std::vector<int>;

	// ������ �� ���� ������� �������� �����
	using PointLists = std::array<PointList, 3>;


	// ������������ ���������� ������� �����

	/**
	* �������� ������� ����� � ���� �������������� � �������� ������� � ���������
	* @param meshDataProvider - ��������� ������� � ������ �����
	* @param centerPosition - ��������� ������ ���������� �������
	* @param areaSizes - ������� ���������� ������� (������ ���� ����� ���������������) � ����� �����
	*/
	void Form
		(
			const MeshDataProvider* meshDataProvider,
			const int centerPosition[3],
			const int areaSizes[3]
		);


	// ���������

	/**
	* �������� ������ �����, �������� � ���������� �������
	* @return ������ �����, �������� � ���������� �������
	*/
	const PointList& GetAreaPoints() const;


	// ������ � ���� / ������ �� �����

	/**
	* �������� ���� � ���������� �������� ����� (� ���� ���������������)
	* @param selectedAreaFileName - ������������ ����� � ���������� �������� �����
	* @return ������� �������� (true) ��� ���������� (false) ������ �����
	*/
	bool ToFile
		(
			const std::string& selectedAreaFileName
		)	const;

	/**
	* ������� ���� � ���������� �������� ����� (� ���� ���������������)
	* @param selectedAreaFileName - ������������ ����� � ���������� �������� �����
	* @return ������� ��������� (true) ��� ����������� (false) ������ �����
	*/
	bool FromFile
		(
			const std::string& selectedAreaFileName
		);

	/**
	* ������� ���� � ���������� �������� ����� (� ���� ���������������)
	* @param selectedAreaFileName - ������������ ����� � ���������� �������� �����
	* @param meshDataProvider - ��������� ������� � ������ �����
	* @param centerPosition - ��������� ������ ���������� ������� � �����
	* @param areaSizes - ������� ���������� ������� (������ ���� ����� ���������������) � ����� �����
	* @return ������� ��������� (true) ��� ����������� (false) �������� �����
	*/
	static
	bool CreateFile
		(
			const std::string& selectedAreaFileName,
			const MeshDataProvider* meshDataProvider,
			const int centerPosition[3],
			const int areaSizes[3]
		);


private:

	// ������ �����, �������� � ���������� �������
	PointList _areaPoints;

	// ������ �����, �������� � ���������� ���������� �������,
	// ���������� ������������ ������������ ���������� � ������ ��������
	PointLists _subAreasPoints;


	// ��������������� ������� ������������ ���������� �������

	/**
	* �������� ������� ����� � ���� �������������� � �������� ������� � ���������
	* @param meshDataProvider - ��������� ������� � ������ �����
	* @param centerPosition - ��������� ������ ���������� �������
	* @param areaSizes - ������� ���������� ������� (������ ���� ����� ���������������) � ����� �����
	* @param selectedArea - ������ ������� ����� �����, ������������ ���������� �������
	* @param selectedSubAreas - ������ ������� ������� �����, ���������� ��������� ���������� �������,
	* ����������������� ������������ ���� � ����������� ����� ����� ���������� �������
	*/
	static
	void SelectArea
		(
			const MeshDataProvider* meshDataProvider,
			const int centerPosition[3],
			const int areaSizes[3],
			PointList& selectedArea,
			PointLists& selectedSubAreas
		);


	// ��������������� ������� ������ � ���� / ������ �� �����

	/**
	* �������� � ���� ������ ������� ����� �����
	* @return ������� �������� (true) ��� ���������� (false) ������ �����
	*/
	static
	bool WritePointList
		(
			ofstream& outputFile,
			const PointList& pointList
		);

	/**
	* ������� � ���� ������ ������� ����� �����
	* @return ������� ��������� (true) ��� ����������� (false) ������ �����
	*/
	static
	bool ReadPointList
		(
			ifstream& inputFile,
			PointList& pointList
		);

	/**
	* �������� ���� � ���������� �������� ����� (� ���� ���������������)
	* @param selectedAreaFileName - ������������ ����� � ���������� �������� �����
	* @param selectedArea - ������ ������� ����� �����, ������������ ���������� �������
	* @param selectedSubAreas - ������ ������� ������� �����, ���������� ��������� ���������� �������,
	* ����������������� ������������ ���� � ����������� ����� ����� ���������� �������
	* @return ������� �������� (true) ��� ���������� (false) ������ �����
	*/
	static
	bool WriteFile
		(
			const std::string& selectedAreaFileName,
			const PointList& selectedArea,
			const PointLists& selectedSubAreas
		);

	/**
	* ������� ���� � ���������� �������� ����� (� ���� ���������������)
	* @param selectedAreaFileName - ������������ ����� � ���������� �������� �����
	* @param selectedArea - ������ ������� ����� �����, ������������ ���������� �������
	* @param selectedSubAreas - ������ ������� ������� �����, ���������� ��������� ���������� �������,
	* ����������������� ������������ ���� � ����������� ����� ����� ���������� �������
	* @return ������� ��������� (true) ��� ����������� (false) ������ �����
	*/
	static
	bool ReadFile
		(
			const std::string& selectedAreaFileName,
			PointList& selectedArea,
			PointLists& selectedSubAreas
		);

};


#endif // GRID_SELECTED_AREA_H