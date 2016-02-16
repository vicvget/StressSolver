#ifndef AUXILIARY_STRESS_STUFF_H

#define AUXILIARY_STRESS_STUFF_H


#include <bitset>
#include <vector>


// ���������� �������� ������� � ������ ��������, �� ������� ����������� ����
// ��� ������� ����������-���������������� ���������
const size_t FreedomsCount = 6;


// ����� ������, ���������� �� ��������������� ������� �������
// � ��������� �������� ��������� ������������� ����
using IsSealedFlags = std::bitset<FreedomsCount>;


// ������ ������� ������ ������������ �������� ������� ��� ��������� ��������� ������������� ����
using IsSealedFlagsList = std::vector<IsSealedFlags>;


#endif // AUXILIARY_STRESS_STUFF_H