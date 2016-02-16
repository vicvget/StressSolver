#include "CsrSymmetricMatrix.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

void CsrSymmetricMatrix::Sort()
{
	for(int k = 0; k < _size; k++)
	{
		for(int i = ia[k]; i < ia[k+1];i++)
		{
			for(int ii = i+1; ii < ia[k+1]; ii++)
			{
				if(ja[ii-1]< ja[i-1])
				{
					std::swap(ja[ii-1],ja[i-1]);
					std::swap(a[ii-1],a[i-1]);
				}
			}				
		}	
	}
}

void CsrSymmetricMatrix::RemoveNearestToZero()
{
	for(int k = 0; k < _size; k++)
	{
		for(int i = ia[k]; i < ia[k+1];i++)
		{
			if(fabs(a[i-1]) < 1e-5)
			{
				a[i-1] = 0;
			}
		}
	}		
}

void CsrSymmetricMatrix::Print() const
{
	for (int k = 0; k < _size; k++)
	{
		for (int i = ia[k]; i < ia[k + 1]; i++)
		{
			std::cout << "a(" << ja[i - 1] << ")=" << a[i - 1] << " ";
		}
		std::cout << std::endl;
	}

}

void CsrSymmetricMatrix::PrintPortrait() const
{
	const std::string zero = "  ";
	const std::string nonzero = "()";
	for(int k = 0; k < _size; k++)
	{
		std::cout << std::setw(2) << k+1 << ' ';
		//			if(k+1 == 14)
		//				int dbg=0;
		int jprev = 1;
		for(int i = ia[k]; i < ia[k+1];i++)
		{
			for(int j = jprev; j < ja[i-1]; j++)
			{
				std::cout << zero;
			}
			//if(fabs(a[i-1])<1e-9)
			//	std::cout << zero;
			//else
			std::cout << nonzero;
			jprev = ja[i-1]+1;
		}
		for(int j = jprev-1; j < _size; j++)
		{
			std::cout << zero;
		}
		std::cout << std::endl;
	}
}

void CsrSymmetricMatrix::Allocate()
{
	ia	= new int[_size+1];
	rhs = new double[_size];
	x	= new double[_size];
	ja	= new int[_memsize];
	a	= new double[_memsize];
	memset(ia,0,sizeof(int)*_size);
	memset(x,0,sizeof(double)*_size);
	memset(rhs,0,sizeof(double)*_size);
}

void CsrSymmetricMatrix::Allocate1d(int nNodes)
{
	_memsize = 6*nNodes-1;
	_size = nNodes;
	Allocate();
}

void CsrSymmetricMatrix::Allocate3d(int nNodes)
{
	_memsize = 18*6*nNodes;
	_size = 6*nNodes;
	Allocate();
}


void CsrSymmetricMatrix::ApplyRhs(int varId, double value)
{
	rhs[varId]=value;
}

void CsrSymmetricMatrix::Fix3dNode(int nodeNumber)
{
	for(int i = 0; i < 6; i++)
	{
		int id = nodeNumber*6+i;
		FixVar(id);
	}
}

void CsrSymmetricMatrix::Fix3dNode(int nodeNumber, int direction)
{
	int id = nodeNumber * 6 + direction;

	FixVar(id);
}

void CsrSymmetricMatrix::Fix3dNodeRotation(int nodeNumber)
{
	for(int i = 5; i < 6; i++)
	{
		int id = nodeNumber*6+i;
		FixVar(id);
	}
}

void CsrSymmetricMatrix::ApplyRhsTo3dNode(int nodeNumber, int direction, double value)
{
	int id = nodeNumber*6+direction;
	ApplyRhs(id, value);
}

void CsrSymmetricMatrix::FixVar(int varId)
{		
	if(varId < _size)
	{
		if(varId == _size-1)
		{
			a[ia[varId]-1] = 1;
		}
		else
		{
			int delta = ia[varId+1]-ia[varId]-1;			
			a[ia[varId]-1]=1;
			for(int i = ia[varId]; i < ia[varId+1]-1; i++)
			{
				a[i]=0;
			}

			for(int i = ia[varId]; i < ia[_size-1]; i++)
			{
				a[i] = a[i+delta];
				ja[i] = ja[i+delta];
			}
			for(int i = varId+1; i < _size+1; i++)
			{
				ia[i]-=delta;
			}
		}
		rhs[varId]=0;
	}
}

CsrSymmetricMatrix::CsrSymmetricMatrix()
	:
		ia(nullptr),
		ja(nullptr),
		a(nullptr),
		rhs(nullptr),
		x(nullptr),
		_size(),
		_memsize()
{

}

CsrSymmetricMatrix::CsrSymmetricMatrix
	(
		int size
	)
	:
		ia(nullptr),
		ja(nullptr),
		a(nullptr),
		_size(size),
		_memsize()
{
	x = new double[size];
	rhs = new double[size];
}

CsrSymmetricMatrix::~CsrSymmetricMatrix()
{
	if(ia != NULL)
		delete [] ia;
	if(ja != NULL)
		delete [] ja;
	if(a != NULL)
		delete [] a;
	if(rhs != NULL)
		delete [] rhs;
	if(x != NULL)
		delete [] x;
}

