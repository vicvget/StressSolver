#pragma OPENCL EXTENSION cl_khr_fp64:enable

//#define DEBUG__

typedef struct LocalSolveParams
{
	double elasticModulus;
	double dampingFactorAngular;
	double dampingFactorLinear;
	double elasticFactorLinear;
	double elasticFactorAngular;
	int nElements;
	double3* DataInternal;/////
	double3* NeighbourDataInternal;//////
	double* RotationMatrix;/////////////
	double* NeighbourRotationMatrixes;//////////
	int* LinkedElements;//////
	double3* RadiusVectors;//////
	double3* boundaryForce;
	double3* elementStressFactorCache;
	double3* stressScalingFactor;
	double3* strainsLinear;
	double3* strainsAngular;
	double3* coordinates;
	uint* isFree;
	int* existLinks;
	uint uid;
	double _timeStep;
	double3 _hDR1;
	double3 _hDR2;
	double3 _hDR3;
	double3 _initR;
	double3 _hDDR1;
	double3 _hDDR2;
	double3 _hDDR3;
	double3 _initDR;
	double FrameMatrix[12];
} SolveParams;

void clCalculateStrains
(
	double3 * shiftStrains,		// выход деформаций
	double3 * velocityStrains,	// выход изм. скоростей
	double3 * DataInternal,
	double3 * NeighbourDataInternal,
	double3 * RadiusVector,
	double  * RotationMatrix,
	double  * NeighbourRotationMatrix
);

void CalculateForces   (SolveParams* params);
void RotationSolve1    (SolveParams* params);
void RotationSolve2    (SolveParams* params);
void RotationSolve3    (SolveParams* params);
void RotationSolve4    (SolveParams* params);
uint IsSingularityAngle(SolveParams* params);
void UpdateMatrix      (SolveParams* params);
void UpdateRHS         (SolveParams* params);
void mulPrivateMM(double* a, double* b, double* res);
void mulPrivateMV(double* a,  double* b,   double* res);

double Get(double3* vector, size_t index)
{
	return *((double*)vector + index);
}

void Set(double3* vector, size_t index,double value)
{
	*(((double*)vector) + index ) = value;
}

//служебные константы, честно спёртые из решателя на С++
constant uint vecStride = 4;				// смещение векторов
constant uint vecStride2 = 8;
constant uint matStride = 12; // смещения матриц поворота (3x3, 3x4)

//инициализирующая итерация производится на хосте, потом скидываю буферы на девайс

//всё стараться хранить в локальной памяти.
//хранить копии обрабатываемых элементов в локальной памяти
//хранить индексы линковки в локальной памяти. индексы отчекрыжить и хранить участок, соответствующий нашей рабочей группе
//хранить копии линкованных элементов в локальной памяти. в том же порядке, в котором они идут в массиве линкованных по индексам
//хранить матрицы и радиус-векторы

__kernel void clCalculate(
	global double * _dataInternal, //[Vx][Va][Vdx][Vda][Vddx][Vdda]
	global double * _dataRotationMatrix, //matrix[12]
	global uint * _linkedElements,//int[nElements]
	global double* radiusVectors,//6*vector3
	global double* boundaryForces,//vector3[nElements]
	global uint* boundaryFree,//6*int[nElements]
	global double* _elementStressFactorCache, //Vector3[nElements]
	global double* _stressScalingFactor,//vector3
	global double* _stress,//vector3 linear vector3 angle [nElements]
	global double* _coordinates,//
	double _elasticModulus,			// модуль упругости (используется только для вычисления напряжений)
	double _dampingFactorAngular,	// приведенный коэффициент углового демпфирования
	double _dampingFactorLinear,	// приведенный коэффициент линейного демпфирования
	double _elasticFactorLinear,	// приведенный коэффициент линейной жесткости
	double _elasticFactorAngular,	// приведенный коэффициент угловой жесткости*/
	int _nElements,
	double _timeStep
)
{
	private SolveParams params;
	params._timeStep=_timeStep;
	params.elasticModulus = _elasticModulus;		// модуль упругости (используется только для вычисления напряжений)
	params.dampingFactorAngular=_dampingFactorAngular;	// приведенный коэффициент углового демпфирования
	params.dampingFactorLinear=_dampingFactorLinear;	// приведенный коэффициент линейного демпфирования
	params.elasticFactorLinear=_elasticFactorLinear;	// приведенный коэффициент линейной жесткости
	params.elasticFactorAngular=_elasticFactorAngular;	// приведенный коэффициент угловой жесткости
	params.nElements = _nElements;
	//[Vx][Va][Vdx][Vda][Vddx][Vdda]
	private double3 DataInternal[6];/////
	params.DataInternal=DataInternal;
	private double3 NeighbourDataInternal[36];//////
	params.NeighbourDataInternal = NeighbourDataInternal;

	private double RotationMatrix[12];/////////////
	params.RotationMatrix=RotationMatrix;
	private double NeighbourRotationMatrixes[72];//////////
	params.NeighbourRotationMatrixes=NeighbourRotationMatrixes;

	private int LinkedElements[6];//////
	params.LinkedElements=LinkedElements;

	private double3 RadiusVectors[6];//////
	params.RadiusVectors=RadiusVectors;
	private double3 boundaryForce[2];
	params.boundaryForce=boundaryForce;
	private double3 elementStressFactorCache;
	params.elementStressFactorCache=&elementStressFactorCache;
	private double3 stressScalingFactor;
	params.stressScalingFactor=&stressScalingFactor;
	private double3 strainsLinear;
	params.strainsLinear=&strainsLinear;
	private double3 strainsAngular;
	params.strainsAngular=&strainsAngular;
	private double3 coordinates;
	params.coordinates=&coordinates;

	private uint isFree[6];
	params.isFree=isFree;

	uint uid = get_global_id(0);
	params.uid=uid;
	uint nElements = params.nElements;
	if (uid >= params.nElements)
	{
		return;
	}
		int existLinks[6] = {1,1,1,1,1,1};
		params.existLinks = existLinks;
		//Копирую индексы прилинкованных элементов из глобальной памяти в приватную
		for (int i = 0; i < 6; i++)
		{
			params.LinkedElements[i] = _linkedElements[i + uid*6];
			if (params.LinkedElements[i]==0)
			{
				params.existLinks[i] = 0;
			}
			else
			{
				params.LinkedElements[i]--;
			}
			//if(existLinks[i] != 0)
				//printf("%d -> %d\n",uid,LinkedElements[i]);
		}
		//printf("========\n");

		//граничные условия
		params.boundaryForce[0] = vload3(0, boundaryForces+vecStride2*uid);
		params.boundaryForce[1] = vload3(0, boundaryForces+vecStride2*uid+vecStride);


		for (int i = 0; i < 6; i++)
		{
			params.isFree[i] = boundaryFree[uid * 6 + i];
		}

		*params.elementStressFactorCache = vload3(0, _elementStressFactorCache + vecStride*uid);

		*params.coordinates = vload3(0, _coordinates + 3 * uid);

		*params.stressScalingFactor = vload3(0, _stressScalingFactor);

		*params.strainsLinear = vload3(0, _stress + vecStride2*uid);
		*params.strainsAngular = vload3(0, _stress + vecStride2*uid + vecStride);

		//загружаем в приватную память параметры для интегрирования из глобальной
		params.DataInternal[0] = vload3(0, (_dataInternal + (uid * vecStride2)));
		params.DataInternal[1] = vload3(0, (_dataInternal + (uid * vecStride2 + vecStride)));
		params.DataInternal[2] = vload3(0, (_dataInternal + (nElements * vecStride2 + uid * vecStride2)));
		params.DataInternal[3] = vload3(0, (_dataInternal + (nElements * vecStride2 + uid * vecStride2 + vecStride)));
		params.DataInternal[4] = vload3(0, (_dataInternal + (nElements * vecStride2 * 2 + uid * vecStride2)));
		params.DataInternal[5] = vload3(0, (_dataInternal + (nElements * vecStride2 * 2 + uid * vecStride2 + vecStride)));

		/*printf("%d ->\n", uid);
		printf("%v3f\n", DataInternal[0]);
		printf("%v3f\n", DataInternal[1]);
		printf("%v3f\n", DataInternal[2]);
		printf("%v3f\n", DataInternal[3]);
		printf("%v3f\n", DataInternal[4]);
		printf("%v3f\n\n", DataInternal[5]);*/

		//таже и для соседних элементов, чтобы иметь быстрый доступ
		for (int i = 0; i < 6; i++)
		{
			params.NeighbourDataInternal[i * 6 + 0] = vload3(0, (_dataInternal + (LinkedElements[i] * vecStride2)));
			params.NeighbourDataInternal[i * 6 + 1] = vload3(0, (_dataInternal + (LinkedElements[i] * vecStride2 + vecStride)));
			params.NeighbourDataInternal[i * 6 + 2] = vload3(0, (_dataInternal + (nElements * vecStride2 + LinkedElements[i] * vecStride2)));
			params.NeighbourDataInternal[i * 6 + 3] = vload3(0, (_dataInternal + (nElements * vecStride2 + LinkedElements[i] * vecStride2 + vecStride)));
			params.NeighbourDataInternal[i * 6 + 4] = vload3(0, (_dataInternal + (nElements * vecStride2 * 2 + LinkedElements[i] * vecStride2)));
			params.NeighbourDataInternal[i * 6 + 5] = vload3(0, (_dataInternal + (nElements * vecStride2 * 2 + LinkedElements[i] * vecStride2 + vecStride)));


		}

		for (int i = 0; i < 6; i++)
		{
			params.RadiusVectors[i] = vload3(0, radiusVectors + i*vecStride);// double3(radiusVectors[i * vecStride], radiusVectors[i * vecStride + 1], radiusVectors[i * vecStride + 2]);
		}

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				params.RotationMatrix[i * 4 + j] = _dataRotationMatrix[matStride * uid + i * 4 + j];
			}
		}

		for (int nbh=0;nbh<6;nbh++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					double t = _dataRotationMatrix[matStride * LinkedElements[nbh] + i * 4 + j];
					params.NeighbourRotationMatrixes[matStride * nbh + i * 4 + j] = t;
				}
			}
		}
		//предполагается, что инициализирующая итерация прошла на хосте!
		//служебные массивы для интегрирования
		//координаты перед началом итерации
		double3 _initX;
		//скорости перед началом итерации
		double3 _initDX;
		//получивниеся приращения второй производной на втором этапе метода
		double3 _hDDX1=(double3)(0,0,0);
		//третьего
		double3 _hDDX2=(double3)(0,0,0);
		//четвертого
		double3 _hDDX3=(double3)(0,0,0);



		_initX = params.DataInternal[0];
		_initDX = (double3)(0,0,0);
		params._initR = params.DataInternal[1];
		params._initDR = params.DataInternal[3];
		params._hDR1 = params.DataInternal[3];
		double* elementMtx = params.RotationMatrix;
		double* rframeMtx = params.FrameMatrix;
		for(int i=0;i<matStride;i++)
		{
			rframeMtx[i]=elementMtx[i];
		}
		params._hDR2=(double3)(0,0,0);
		params._hDR3=(double3)(0,0,0);
		//phase 2
		//интегрирование линейное

		_hDDX1 = params.DataInternal[4] * params._timeStep;
		params.DataInternal[0] += params.DataInternal[2] * params._timeStep * 0.5;
		params.DataInternal[2] += _hDDX1 * 0.5;


		//интегрирование вращений
		RotationSolve2(&params);
		CalculateForces(&params);

		//phase 3
		_hDDX2 = params.DataInternal[4] * params._timeStep;
		params.DataInternal[0] += _hDDX1 * params._timeStep * 0.25;
		params.DataInternal[2] = _initDX + _hDDX2 * 0.5;

		//	_stageRK = 3;
		RotationSolve3(&params);
		CalculateForces(&params);

		//phase 4

		_hDDX3 = params.DataInternal[4] * params._timeStep;
		params.DataInternal[0] = _initX + (_initDX + _hDDX2 * 0.5) * params._timeStep;
		params.DataInternal[2] = _initDX + _hDDX3;


		//	_stageRK = 4;
		RotationSolve4(&params);
 		CalculateForces(&params);

		//phase 5

		double3 sDDX = _hDDX2 + _hDDX3;
		params.DataInternal[0] = _initX + (_initDX + (_hDDX1 + sDDX) / 6.0) * params._timeStep;
		params.DataInternal[2] = _initDX + (_hDDX1 + sDDX + sDDX + params.DataInternal[4] * params._timeStep) / 6.0;

		//printf("pos %v3f\n\n",params.DataInternal[0]);

		RotationSolve1(&params);
		CalculateForces(&params);
/*
#ifdef DEBUG__
		printf("%v3f\n",params.DataInternal[0]);
		printf("%v3f\n",params.DataInternal[1]);
		printf("%v3f\n",params.DataInternal[2]);
		printf("%v3f\n",params.DataInternal[3]);
		printf("%v3f\n",params.DataInternal[4]);
		printf("%v3f\n",params.DataInternal[5]);
		printf("===============================\n");
		printf("%v3f\n",params.strainsLinear[0]);
		printf("%v3f\n",params.strainsAngular[0]);
		printf("===============================\n\n");
#endif*/
		////
		////Загрузка обратно в глобальную память//////////////////////////////////////////////////////////////////////////////////
		////

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				/*#ifdef DEBUG__
				printf("%f ",params.RotationMatrix[i * 4 + j]);
				#endif*/
				_dataRotationMatrix[matStride * uid + i * 4 + j] = params.RotationMatrix[i * 4 + j];
			}
			/*#ifdef DEBUG__
			printf("\n");
			#endif*/
		}
		/*#ifdef DEBUG__
		printf("===============================\n\n");
		#endif*/
		vstore3(params.DataInternal[0], 0, (_dataInternal +						       uid * vecStride2             ));
		vstore3(params.DataInternal[1], 0, (_dataInternal +						       uid * vecStride2 + vecStride ));
		vstore3(params.DataInternal[2], 0, (_dataInternal + (nElements * vecStride2     + uid * vecStride2			   )));
		vstore3(params.DataInternal[3], 0, (_dataInternal + (nElements * vecStride2     + uid * vecStride2 + vecStride)));
		vstore3(params.DataInternal[4], 0, (_dataInternal + (nElements * vecStride2 * 2 + uid * vecStride2            )));
		vstore3(params.DataInternal[5], 0, (_dataInternal + (nElements * vecStride2 * 2 + uid * vecStride2 + vecStride)));
		_coordinates[3*uid]=Get(&DataInternal[0],0);
		_coordinates[3*uid+1]=Get(&DataInternal[0],1);
		_coordinates[3*uid+2]=Get(&DataInternal[0],2);

		vstore3(params.strainsLinear[0], 0, _stress + vecStride2*uid);
		vstore3(params.strainsAngular[0], 0, _stress + vecStride2*uid + vecStride);

}

void RotationSolve1(SolveParams* params)
{/*
	double3 sDDX = params->_hDDR2 + params->_hDDR3;
	params->DataInternal[1] = params-> _initR + (params->_initDR + (params->_hDDR1 + sDDX) / 6.0) * params->_timeStep;
	params->DataInternal[3] = params->_initDR + (params->_hDDR1 + sDDX + sDDX + params->DataInternal[5] * params->_timeStep) / 6.0;*/

	params->DataInternal[1] = params->_initR + params->_hDR1 * 0.5;
	UpdateRHS(params);
	params->_hDR2=params->DataInternal[3];
	UpdateMatrix(params);
}

void RotationSolve2(SolveParams* params)
{/*
	params->_hDDR1 = params->DataInternal[5] * params->_timeStep;
	params->DataInternal[1] += params->DataInternal[3] * params->_timeStep * 0.5;
	params->DataInternal[3] += params->_hDDR1 * 0.5;*/

	params->DataInternal[1] = params->_initR + params->_hDR2 * 0.5;
	UpdateRHS(params);
	params->_hDR3=params->DataInternal[3];
	UpdateMatrix(params);
}

void RotationSolve3(SolveParams* params)
{/*
	params->_hDDR2 = params->DataInternal[5] * params->_timeStep;
	params->DataInternal[1] += params->_hDDR1 * params->_timeStep * 0.25;
	params->DataInternal[3] = params->_initDR + params->_hDDR2 * 0.5;*/

	params->DataInternal[1] = params->_initR + params->_hDR3;
	UpdateRHS(params);
	UpdateMatrix(params);
}

void RotationSolve4(SolveParams* params)
{/*
	params->_hDDR3 = params->DataInternal[5] * params->_timeStep;
	params->DataInternal[1] = params->_initR + (params->_initDR + params->_hDDR2 * 0.5) * params->_timeStep;
	params->DataInternal[3] = params->_initDR + params->_hDDR3;*/

	params->DataInternal[1] = params->_initR + (params->_hDR1 + 2 * (params->_hDR2 + params->_hDR3) + params->DataInternal[3]) / 6.0;

	UpdateRHS(params);
	// проверка сингулярности

	double* elementMtx = params->RotationMatrix;
	double* rframeMtx = params->FrameMatrix;
	if (IsSingularityAngle(params))
	{
		params->DataInternal[3]=(double3)(0,0,0);
		params->_hDR1=(double3)(0,0,0);
		params->_hDR2=(double3)(0,0,0);
		params->_hDR3=(double3)(0,0,0);
		params->DataInternal[1]=(double3)(0,0,0);
		params->_initR=(double3)(0,0,0);
		for(int i=0;i<matStride;i++)
		{
			rframeMtx[i]=elementMtx[i];
		}
	}
	UpdateMatrix(params);
}

uint IsSingularityAngle(SolveParams* params)
{
	double angle = params->DataInternal[1].y;
	return ((M_PI_2 - (angle - ((int)(angle / M_PI)) * M_PI)) < 0.1 && (M_PI_2 - (angle - ((int)(angle / M_PI)) * M_PI))>0) || ((M_PI_2 - (angle - ((int)(angle / M_PI)) * M_PI))>-0.1&&(M_PI_2 - (angle - ((int)(angle / M_PI)) * M_PI))<0);
}

void UpdateMatrix(SolveParams* params)
{
	double* rotationMtx = params->RotationMatrix;
	double* rframeMtx = params->FrameMatrix;
	double* angles = &params->DataInternal[1];
	double cosX = cos(angles[0]);
	double sinX = sin(angles[0]);
	double cosY = cos(angles[1]);
	double sinY = sin(angles[1]);
	double cosZ = cos(angles[2]);
	double sinZ = sin(angles[2]);
	double newMtx[12];
	newMtx[0*4 + 0] = cosY*cosZ;
	newMtx[0*4 + 1] = -cosY*sinZ;
	newMtx[0*4 + 2] = sinY;
	newMtx[0*4 + 3] = 0;
	newMtx[1*4 + 0] = sinX*sinY*cosZ + cosX*sinZ;
	newMtx[1*4 + 1] = -sinX*sinY*sinZ + cosX*cosZ;
	newMtx[1*4 + 2] = -sinX*cosY;
	newMtx[1*4 + 3] = 0;
	newMtx[2*4 + 0] = -cosX*sinY*cosZ + sinX*sinZ;
	newMtx[2*4 + 1] = cosX*sinY*sinZ + sinX*cosZ;
	newMtx[2*4 + 2] = cosX*cosY;
	newMtx[2*4 + 3] = 0;
	mulPrivateMM(rframeMtx,newMtx,rotationMtx);//(rframe*newMtx).Export(rotationMtx);

}

void UpdateRHS(SolveParams* params)
{
	double3 angles = params->DataInternal[1];
	double3 elementW = *params->strainsAngular;

	double cosY = cos(Get(&angles,1));
	double tanY = tan(Get(&angles,1));
	double sinZ = sin(Get(&angles,2));
	double cosZ = cos(Get(&angles,2));

	double wx = Get(&elementW,0);
	double wy = Get(&elementW,1);
	double wz = Get(&elementW,2);

/*
	// Проверки сходимости
	double controlAngleX = (-sinZ*wx - cosZ*wy) / cosZ;
	double controlAngleY = (wy*cosZ + wx*sinZ)*tanY;
	double maxValue = 400.0 / (_timeStep*_timeStep);
	bool result = true;
	if (std::abs(controlAngleX) > maxValue)
	{
		std::cout << "small step for 1 euler eq\n";
		result = false;
	}
	if (std::abs(controlAngleY) > maxValue)
	{
		std::cout << "small step for 3 euler eq\n";
		result = false;
	}
*/

	double* derivatives = &params->DataInternal[3];

	derivatives[0] = (-wy*sinZ + wx*cosZ) / cosY * params->_timeStep;
	derivatives[1] = ( wx*sinZ + wy*cosZ) *params->_timeStep;
	derivatives[2] = ((wy*sinZ - wx*cosZ) * tanY + wz) * params->_timeStep;
}

//Рассчет напряжений и сил
void CalculateForces(SolveParams* params)
{

	int exclusive_dofs[3][2] = { { 1, 2 },{ 0, 2 },{ 1, 3 } };
	double3 shiftStrains[2];
	shiftStrains[0] = (double3)(0, 0, 0);
	shiftStrains[1] = (double3)(0, 0, 0);
	double3 velocityStrains[2];
	velocityStrains[0] = (double3)(0, 0, 0);
	velocityStrains[1] = (double3)(0, 0, 0);
	params->DataInternal[4] = (double3)(0, 0, 0);
	params->DataInternal[5] = (double3)(0, 0, 0);
	*params->strainsLinear = (double3)(0,0,0);
	*params->strainsAngular = (double3)(0, 0, 0);
	for (int i = 0; i < 6; i++)
	{
			clCalculateStrains(shiftStrains, velocityStrains, params->DataInternal, &params->NeighbourDataInternal[i * 6], &params->RadiusVectors[i], params->RotationMatrix, &params->NeighbourRotationMatrixes[matStride*i]);
			double3 linear_strains = shiftStrains[0];
			double3 angular_strains = shiftStrains[1];
			double3 linear_vstrains = velocityStrains[0];
			double3 angular_vstrains = velocityStrains[1];
			//текущая компонента напряженности
			double temp0 = Get(&linear_strains, i % 3) * Get(params->elementStressFactorCache, i % 3) * Get(params->stressScalingFactor, i % 3);
			//в зависимости от стороны, которая сейчас рассчитывается устанавливается знак
			#ifdef DEBUG__
			if (temp0>=0.000000000001 && params->existLinks[i])
			{
				printf("shift strains %10.20v3f\n", shiftStrains[0]);
				printf("velocity strains %10.20v3f\n", velocityStrains[0]);
				printf("temp 0 %10.20f element %d side %d\n\n", temp0,params->uid,i);
			}
			#endif
			if(i<3)
			{
				Set(params->strainsLinear, i % 3, temp0*params->existLinks[i] + Get(params->strainsLinear,i%3));
			}
			else
			{
				Set(params->strainsLinear, i % 3, -temp0*params->existLinks[i] + Get(params->strainsLinear,i%3));
			}
			int dof0 = exclusive_dofs[i % 3][0];
			int dof1 = exclusive_dofs[i % 3][1];
			temp0 = Get(&linear_strains, dof1) * Get(params->elementStressFactorCache, i % 3) * Get(params->stressScalingFactor, i % 3);
			double temp1 = Get(&linear_strains, dof0) * Get(params->elementStressFactorCache, i % 3) *  Get(params->stressScalingFactor, i % 3);
			if(i<3)
			{
				Set(params->strainsAngular, dof0, temp0*params->existLinks[i] + Get(params->strainsAngular, dof0));
				Set(params->strainsAngular, dof1, temp1*params->existLinks[i] + Get(params->strainsAngular, dof1));
			}
			else
			{
				Set(params->strainsAngular, dof0, -temp0*params->existLinks[i] + Get(params->strainsAngular, dof0));
				Set(params->strainsAngular, dof1, -temp1*params->existLinks[i] + Get(params->strainsAngular, dof1));
			}
			double3 vForce1 = -linear_vstrains * params->dampingFactorLinear - linear_strains * params->elasticFactorLinear;
			double3 vTorque = -angular_vstrains * params->dampingFactorAngular - angular_strains * params->elasticFactorAngular;
			double3 vForce0;
			mulPrivateMV(params->RotationMatrix, &vForce1, &vForce0);
			double3 vForce1Torque = cross(params->RadiusVectors[i], vForce1);
			double3 vTorque1 = vForce1Torque + vTorque;
			params->DataInternal[4] += (vForce0*params->existLinks[i])*params->isFree[i%3];
			if (i == 0)
			{
				params->DataInternal[4] += params->boundaryForce[0] * params->isFree[i%3];
			}

			params->DataInternal[5] += (vTorque1*params->existLinks[i])*params->isFree[i%3+3];
			if(i == 0)
			{
				params->DataInternal[5] += params->boundaryForce[1] * params->isFree[i%3+3];
			}
	}
	#ifdef DEBUG__
	if (length(params->DataInternal[4])>=0.000000000001 && params->existLinks[0])
	{
		printf("element %d\ntotal force = %10.20v3f\n\n", params->uid,params->DataInternal[4]);
	}
	#endif
}

void clCalculateStrains
(
	double3 * shiftStrains,		// выход деформаций
	double3 * velocityStrains,	// выход изм. скоростей
	double3 * DataInternal,
	double3 * NeighbourDataInternal,
	double3 * RadiusVector,
	double  * RotationMatrix,
	double  * NeighbourRotationMatrix
)
{
	double* matA01 = RotationMatrix;
		//matA01 = matA01.Tr();
		//matA02 = matA02.Tr();
		//Mat3 matA12 = matA01.Tmul(matA02);
		double matA21[12]; //= matA02.Tmul(matA01);
		mulPrivateMM(NeighbourRotationMatrix, matA01, &matA21);

		double3 vecC1 = *RadiusVector;
		double3 vecC2 = -vecC1;

		double3 vecP1 = DataInternal[0];
		double3 vecP2 = NeighbourDataInternal[0];
		double3 vecR1 = DataInternal[1];
		double3 vecR2 = NeighbourDataInternal[1];
		double3 vecV1 = DataInternal[2];
		double3 vecV2 = NeighbourDataInternal[2];
		double3 vecW1 = DataInternal[3];
		double3 vecW2 = NeighbourDataInternal[3];

		/*

		printf("%v3f\n", NeighbourDataInternal[ 0]);
		printf("%v3f\n", NeighbourDataInternal[ 1]);
		printf("%v3f\n", NeighbourDataInternal[ 2]);
		printf("%v3f\n", NeighbourDataInternal[ 3]);
		printf("%v3f\n", NeighbourDataInternal[ 4]);
		printf("%v3f\n\n", NeighbourDataInternal[ 5]);



		printf("%v3f\n", NeighbourDataInternal[0]);
		printf("%v3f\n", NeighbourDataInternal[1]);
		printf("%v3f\n", NeighbourDataInternal[2]);
		printf("%v3f\n", NeighbourDataInternal[3]);
		printf("%v3f\n", NeighbourDataInternal[4]);
		printf("%v3f\n\n", NeighbourDataInternal[5]);

		*/

		/*printf("%v3f\n",   vecP1);
		printf("%v3f\n\n", vecP2);
		*/


		double3 temp0, temp1;
		// переводим вектор линии точек связи С2-С1 в СК1

		mulPrivateMV(matA21,&vecC2,&temp0);
		double3 temp2 = (vecP2 - vecP1);
		mulPrivateMV(matA01, &temp2, &temp1);
		double3 lShiftStrains = vecC1 - temp0 - temp1;

		shiftStrains[0] = lShiftStrains;
		/*
		if (get_global_id(0) == 99)
		{
			printf("Radius vector %v3f\n", vecC1);
			printf("temp0 %v3f\n", temp0);
			printf("temp1 %v3f\n", temp1);
			printf("temp2 %v3f\n", temp2);
			printf("linear shift %v3f\n", shiftStrains[0]);
			printf("=========\n");
		}*/
		// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1

		temp1 = (vecV1 - vecV2);
		mulPrivateMV(matA01, &temp1, &temp0);
		temp0 = cross(vecW2, vecC2);
		double3 temp3;
		mulPrivateMV(matA21, &temp0, &temp3);

		double3 lVelocityStrains = temp1 + cross(vecW1,vecC1) - temp3;

		velocityStrains[0] = lVelocityStrains;


		// для угловых степеней свободы - просто берутся разница углов и угловых скоростей
		lShiftStrains = (vecR1 - vecR2);
		lVelocityStrains = (vecW1 - vecW2);
		shiftStrains[1] = lShiftStrains;
		velocityStrains[1] = lVelocityStrains;
}



void mulPrivateMM(double* a, double* b, double * res)
{
	for (int i = 0; i < 3; i++)
	{

		for (int j = 0; j < 3; j++)
		{
			res[i * 4 + j] = 0;
			for (int k = 0; k < 3; k++)
			{
				res[i*4+j] += a[k*4+i] * b[k*4+j];

			}
			//printf("%f ",b[i*4+j]);
		}
		//printf("\n");
	}
}

void mulPrivateMV(double * a, double* b, double* res)
{
	for (int i = 0; i < 3; i++) {
		res[i] = 0;
		for (int j = 0; j < 3; j++) {
			res[i] += a[i * 4 + j] * b[j];
		}
	}
}
