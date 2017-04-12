#include "StressStrainCppIterativeSolverKNC2.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"

using MathHelpers::Vec3;
using MathHelpers::Vec3Ref;
using MathHelpers::MakeVec3;
using MathHelpers::Mat3;
using MathHelpers::Mat3x4;


//__m512i izmm00 = _mm512_load_epi64(pmtx);    //A2A1
//__m512i izmm01 = _mm512_load_epi64(pmtx+8);	 //B1A3
//__m512i izmm02 = _mm512_load_epi64(pmtx+16); //B3B2
//__m512i izmm03 = _mm512_alignr_epi32(izmm00, izmm01, 8); //A1B1
//__m512i izmm04 = _mm512_alignr_epi32(izmm01, izmm02, 8); //A3B3
//
//__m512d zmm00 = (__m512d)_mm512_alignr_epi32(izmm03, izmm03, 8); //B1A1
//__m512d zmm01 = (__m512d)_mm512_alignr_epi32(izmm02, izmm00, 8); //B2A2
//__m512d zmm02 = (__m512d)_mm512_alignr_epi32(izmm04, izmm04, 8); //B3A3
//
////Вместо D1C1,D2C2,D3C3
////нам надо получить С1B1,C2B2,C3B3
//
//izmm03 = _mm512_load_epi64(pmtr);    //C2C1
//izmm04 = _mm512_load_epi64(pmtr+8);  //D1C3 // D1 не используется
//izmm05 = _mm512_alignr_epi32(izmm02, izmm03, 8); //B2C2
//
//zmm03 = (__m512d)_mm512_alignr_epi32(izmm03, izmm01, 8); //C1B1
//zmm04 = (__m512d)_mm512_alignr_epi32(izmm02, izmm02, 8); //C2B2
//zmm05 = (__m512d)_mm512_alignr_epi32(izmm04, izmm02, 8); //C3B3
//
//// остальная часть кода такая же как KNC

namespace Stress
{
	
		StressStrainCppIterativeSolverKNC2::StressStrainCppIterativeSolverKNC2
		(
			double* params,
			int* links,
			int nLinks,
			double *gridElements,
			int nElements,
			double gridStep,
			double timeStep,
			int numThreads,
			int stride
		)
			:
			StressStrainCppIterativeSolverKNC
			(
				params,
				links,
				nLinks,
				gridElements,
				nElements,
				gridStep,
				timeStep,
				numThreads,
				stride
			)
		{
			std::cout << "KNC2 SOLVER" << std::endl << std::flush;
		}
		// virtual
		StressStrainCppIterativeSolverKNC2::~StressStrainCppIterativeSolverKNC2()
		{
		}
		void StressStrainCppIterativeSolverKNC2::CalculateForces()
		{
			_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
			static int it = 0;

			__declspec(align(64)) double strains[16], velocityStrains[16];

			for (int elementId1 = 0; elementId1 < _nElements; elementId1++)
			{
				memset(GetElementAcceleration(elementId1), 0u, sizeof(double)*vecStride2);
				memset(GetElementStress(elementId1), 0u, sizeof(double)*vecStride2);
			}

			const int exclusive_dofs[][2] = { { 1, 2 },{ 0, 2 },{ 1, 3 } };

			//#pragma omp parallel for private (strains, velocityStrains) num_threads(_numThreads)
			for (int elementId1 = 0; elementId1 + 1 < _nElements; elementId1 += 2)
			{
				// обход x-,y-,z-
				for (int dof = 0; dof < 3; dof++)
				{
					size_t elementId2 = GetLinkedElement(elementId1, dof);
					//cout<<"privet!  sosed"<<elementId1<<std::endl;
					if (elementId2)
					{
						elementId2--;
						//cout<<"privet!  "<<dof<<elementId2<<elementId1<<_nElements <<std::endl;

						CalculateStrains(dof, strains, velocityStrains, elementId1, elementId2);

						//cout<< "ghbdtn"<<strains[0]<<" "<<strains[1]<<" "<< strains[2]<<" "<< strains[4]<<" "<< strains[5]<<" "<< strains[6]<<std::endl<<std::flush;
						for (int rep = 0; rep<2; rep++)
						{
							Vec3Ref linear_strains = MakeVec3(&strains[0] + rep*vecStride);
							Vec3Ref angular_strains = MakeVec3(&strains[0] + rep*vecStride + 8);
							Vec3Ref linear_vstrains = MakeVec3(&velocityStrains[0] + rep*vecStride);
							Vec3Ref angular_vstrains = MakeVec3(&velocityStrains[0] + rep*vecStride + 8);

							// нормальные напряжения
							GetElementStress(elementId1 + rep)[dof] += linear_strains[dof] * GetElementStressFactors(elementId1 + rep)[dof] * _stressScalingFactors[dof];
							GetElementStress(elementId2 + rep)[dof] += linear_strains[dof] * GetElementStressFactors(elementId2 + rep)[dof] * _stressScalingFactors[dof];

							// степени свободы смещений, участвующих в создании касательных напряжений
							int dof0 = exclusive_dofs[dof][0];
							int dof1 = exclusive_dofs[dof][1];

							GetElementStressAngular(elementId1 + rep)[dof0] += linear_strains[dof1] * GetElementStressFactors(elementId1)[dof] * _stressScalingFactors[dof];
							GetElementStressAngular(elementId1 + rep)[dof1] += linear_strains[dof0] * GetElementStressFactors(elementId1)[dof] * _stressScalingFactors[dof];

							GetElementStressAngular(elementId2 + rep)[dof0] += linear_strains[dof1] * GetElementStressFactors(elementId2)[dof] * _stressScalingFactors[dof];
							GetElementStressAngular(elementId2 + rep)[dof1] += linear_strains[dof0] * GetElementStressFactors(elementId2)[dof] * _stressScalingFactors[dof];

							// сила и момент из полученных деформаций
							Vec3 vForce1 = -linear_vstrains * _dampingFactorLinear - linear_strains * _elasticFactorLinear;
							Vec3 vTorque = -angular_vstrains * _dampingFactorAngular - angular_strains * _elasticFactorAngular;

							Vec3 vForce0; // сила в СК0
							Vec3 vForce2; // сила в СК второго тела

							if (vecStride == 4)
							{
								Mat3x4 matA01(GetRotationMatrix(elementId1 + rep));
								Mat3x4 matA02(GetRotationMatrix(elementId2 + rep));
								vForce0 = matA01*vForce1;
								vForce2 = matA02.Tmul(vForce0);
							}
							else
							{
								Mat3 matA01(GetRotationMatrix(elementId1 + rep));
								Mat3 matA02(GetRotationMatrix(elementId2 + rep));
								vForce0 = matA01*vForce1;
								vForce2 = matA02.Tmul(vForce0);
							}

							Vec3Ref vR = MakeVec3(GetRadiusVector(dof));
							Vec3 vForce1Torque = vR.Cross(vForce1);
							Vec3 vForce2Torque = vR.Cross(vForce2); //(-R and -vForce2 gives +vForce2Torque)

																	// Full torque
							Vec3 vTorque1 = vForce1Torque + vTorque;
							Vec3 vTorque2 = vForce2Torque - vTorque;

							MakeVec3(GetElementAcceleration(elementId1 + rep)) += vForce0;
							MakeVec3(GetElementAccelerationAngular(elementId1 + rep)) += vTorque1;
							MakeVec3(GetElementAcceleration(elementId2 + rep)) -= vForce0;
							MakeVec3(GetElementAccelerationAngular(elementId2 + rep)) += vTorque2;
						}
					}
				}
			}
			ApplyBoundary(); // модифицирует силы и моменты
			ApplyMass();	 // вычисляет ускорения делением сил на массы и моментов на моменты инерции
		}
		void StressStrainCppIterativeSolverKNC2::CrossProductTwice
		(
			double* v1,
			double* v2,
			double* v3,
			double* res
		)	const
		{
			res[0] = v1[1] * v2[2] - v1[2] * v2[1];
			res[1] = -v1[0] * v2[2] + v1[2] * v2[0];
			res[2] = v1[0] * v2[1] - v1[1] * v2[0];

			res[4] = v3[1] * v2[2] - v3[2] * v2[1];
			res[5] = -v3[0] * v2[2] + v3[2] * v2[0];
			res[6] = v3[0] * v2[1] - v3[1] * v2[0];
		}
		void StressStrainCppIterativeSolverKNC2::CalculateStrains
		(
			size_t side,			// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
			double *shiftStrains,		// выход деформаций
			double *velocityStrains,	// выход изм. скоростей
			size_t nodeId1,				// номер узла 1
			size_t nodeId2					// номер узла 2
		) const
		{
#if defined(USE_KNC) || defined(USE_KNL)

			// Start AVX code
			double* pmtx = GetRotationMatrix(nodeId1);
			double* pmtr = GetRotationMatrix(nodeId2);
			__m512i izmm00, izmm01, izmm02, izmm03, izmm04, izmm05;
			__m512d zmm00, zmm01, zmm02, zmm03, zmm04, zmm05;
			if (side == 1)
			{

				izmm00 = _mm512_load_epi64(pmtx);    //A2A1
				izmm01 = _mm512_load_epi64(pmtx + 8);	 //B1A3
				izmm02 = _mm512_load_epi64(pmtx + 16); //B3B2
				izmm03 = _mm512_alignr_epi32(izmm00, izmm01, 8); //A1B1
				izmm04 = _mm512_alignr_epi32(izmm01, izmm02, 8); //A3B3

				zmm00 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm03, izmm03, 8)); //B1A1
				zmm01 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm02, izmm00, 8)); //B2A2
				zmm02 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm04, izmm04, 8)); //B3A3

																		 //Вместо D1C1,D2C2,D3C3
																		 //нам надо получить С1B1,C2B2,C3B3

				izmm03 = _mm512_castpd_si512(_mm512_loadu_pd(pmtr));    //C2C1
				izmm04 = _mm512_castpd_si512(_mm512_loadu_pd(pmtr + 8));  //D1C3 // D1 не используется
				izmm05 = _mm512_alignr_epi32(izmm04, izmm00, 8); 
				izmm02 = _mm512_alignr_epi32(izmm03, izmm03, 8);

				zmm03 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm00, izmm02, 8)); //C1B1
				zmm04 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm05, izmm03, 8)); //C2B2
				zmm05 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm01, izmm05, 8)); //C3B3
			}
			else
			{
				izmm00 = _mm512_load_epi64(pmtx);    //A2A1
				izmm01 = _mm512_load_epi64(pmtx + 8);	 //B1A3
				izmm02 = _mm512_load_epi64(pmtx + 16); //B3B2
				izmm03 = _mm512_alignr_epi32(izmm00, izmm01, 8); //A1B1
				izmm04 = _mm512_alignr_epi32(izmm01, izmm02, 8); //A3B3

				zmm00 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm03, izmm03, 8)); //B1A1
				zmm01 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm02, izmm00, 8)); //B2A2
				zmm02 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm04, izmm04, 8)); //B3A3

				izmm00 = _mm512_load_epi64(pmtr);    //C2C1
				izmm01 = _mm512_load_epi64(pmtr + 8);  //D1C3
				izmm02 = _mm512_load_epi64(pmtr + 16); //D3D2
				izmm03 = _mm512_alignr_epi32(izmm00, izmm01, 8); //C1D1
				izmm04 = _mm512_alignr_epi32(izmm01, izmm02, 8); //C3D3

				zmm03 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm03, izmm03, 8)); //D1C1
				zmm04 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm02, izmm00, 8)); //D2C2
				zmm05 = _mm512_castsi512_pd(_mm512_alignr_epi32(izmm04, izmm04, 8)); //D3C3
			}

			__m512d zmm06 = _mm512_swizzle_pd(zmm00, _MM_SWIZ_REG_AAAA); //A11-B11
			__m512d zmm07 = _mm512_swizzle_pd(zmm00, _MM_SWIZ_REG_BBBB); //A21-B21
			__m512d zmm08 = _mm512_swizzle_pd(zmm00, _MM_SWIZ_REG_CCCC); //A31-B31

			__m512d zmm09 = _mm512_mul_pd(zmm06, zmm03);
			zmm09 = _mm512_fmadd_pd(zmm07, zmm04, zmm09);
			zmm09 = _mm512_fmadd_pd(zmm08, zmm05, zmm09);


			zmm06 = _mm512_swizzle_pd(zmm01, _MM_SWIZ_REG_AAAA); //A12-B12
			zmm07 = _mm512_swizzle_pd(zmm01, _MM_SWIZ_REG_BBBB); //A22-B22
			zmm08 = _mm512_swizzle_pd(zmm01, _MM_SWIZ_REG_CCCC); //A32-B32

			__m512d zmm10 = _mm512_mul_pd(zmm06, zmm03);
			zmm10 = _mm512_fmadd_pd(zmm07, zmm04, zmm10);
			zmm10 = _mm512_fmadd_pd(zmm08, zmm05, zmm10);


			zmm06 = _mm512_swizzle_pd(zmm02, _MM_SWIZ_REG_AAAA); //A13-B13
			zmm07 = _mm512_swizzle_pd(zmm02, _MM_SWIZ_REG_BBBB); //A23-B23
			zmm08 = _mm512_swizzle_pd(zmm02, _MM_SWIZ_REG_CCCC); //A33-B33

			__m512d zmm11 = _mm512_mul_pd(zmm06, zmm03);
			zmm11 = _mm512_fmadd_pd(zmm07, zmm04, zmm11);
			zmm11 = _mm512_fmadd_pd(zmm08, zmm05, zmm11);
			// Матрица A_{21} сформирована

			//std::cout << "Multiply add completed" << std::endl << std::flush;

			double* vecC1 = GetRadiusVector(side); // ? ?
			__m512d ivecC1 = _mm512_extload_pd(vecC1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE);
			__m512d vecDP1 = _mm512_sub_pd(
				_mm512_load_pd(GetElementShift(nodeId1)),
				_mm512_load_pd(GetElementShift(nodeId2))); // P1-P2
			zmm06 = _mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_castpd_si512(vecDP1), _mm512_castpd_si512(vecDP1), 8));
			__m512d vecDP2 = _mm512_sub_pd(
				_mm512_load_pd(GetElementShift(nodeId1 + 1)),
				_mm512_load_pd(GetElementShift(nodeId2 + 1))); // P1-P2
			__m512d vecDP = _mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_castpd_si512(vecDP2), _mm512_castpd_si512(zmm06), 8)); //DP1 DP2 ?
			 //std::cout << "Radius vector loaded" << std::endl << std::flush;


			 //__declspec(align(64)) double tmp[8]
			 //double tmp[8]  __attribute__((aligned(64)));
			//double* tmp = _buffer;
			//_mm512_store_pd(tmp, ivecC1);
			//std::cout << "Stored to tmp" << std::endl << std::flush;

			__m512d matA02el1 = _mm512_set1_pd(vecC1[0]);
			__m512d matA02el2 = _mm512_set1_pd(vecC1[1]);
			__m512d matA02el3 = _mm512_set1_pd(vecC1[2]);
			
			__declspec(align(64)) double tmp[8] = { 0 };

			matA02el1 = _mm512_mul_pd(zmm09, matA02el1);
			matA02el2 = _mm512_mul_pd(zmm10, matA02el2);
			matA02el3 = _mm512_mul_pd(zmm11, matA02el3);

			__m512d mul1 = _mm512_add_pd(_mm512_add_pd(matA02el1, matA02el2), matA02el3);

			_mm512_store_pd(tmp, vecDP);
			matA02el1 = _mm512_swizzle_pd(vecDP, _MM_SWIZ_REG_AAAA); //_mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_castpd_si512(_mm512_set1_pd(tmp[4])), _mm512_castpd_si512(_mm512_set1_pd(tmp[0])), 8));
			matA02el2 = _mm512_swizzle_pd(vecDP, _MM_SWIZ_REG_BBBB); //_mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_castpd_si512(_mm512_set1_pd(tmp[5])), _mm512_castpd_si512(_mm512_set1_pd(tmp[1])), 8));
			matA02el3 = _mm512_swizzle_pd(vecDP, _MM_SWIZ_REG_CCCC); //_mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_castpd_si512(_mm512_set1_pd(tmp[6])), _mm512_castpd_si512(_mm512_set1_pd(tmp[2])), 8));


			//????
			matA02el1 = _mm512_mul_pd(zmm09, matA02el1);
			matA02el2 = _mm512_mul_pd(zmm10, matA02el2);
			matA02el3 = _mm512_mul_pd(zmm11, matA02el3);

			__m512d mul2 = _mm512_add_pd(_mm512_add_pd(matA02el1, matA02el2), matA02el3);
			__m512d res = _mm512_add_pd(_mm512_add_pd(ivecC1, mul1), mul2);

			_mm512_store_pd(shiftStrains, res); // получено SL, линейные компоненты

												// Расчет VL
												// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
												// Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) + vecW1.Cross(vecC1) - matA12*(vecW2.Cross(vecC2));
												// vecC2 = -vecC1
												// vecC2.Cross(vecW1) = -vecC1.Cross(vecW1)

			__declspec(align(64)) double cp1[8] = { 0 };
			__declspec(align(64)) double cp2[8] = { 0 };

												//double cp1[8] __attribute__((aligned(64)));
												//double cp2[8] __attribute__((aligned(64)));


			CrossProductTwice(GetElementVelocityAngular(nodeId1), GetRadiusVector(side), GetElementVelocityAngular(nodeId1 + 1), cp1);	// [w1 x c1]
			CrossProductTwice(GetElementVelocityAngular(nodeId2), GetRadiusVector(side), GetElementVelocityAngular(nodeId2 + 1), cp2);	// -[w2 x c2] = [w2 x c1]


																																		//
			__m512d cp1r = _mm512_load_pd(&cp1[0]);
			__m512d vecDV = _mm512_sub_pd(
				_mm512_load_pd(GetElementVelocity(nodeId1)),
				_mm512_load_pd(GetElementVelocity(nodeId2))); // V1-V2

			__m512d vecDV2 = _mm512_sub_pd(
				_mm512_load_pd(GetElementVelocity(nodeId1 + 1)),
				_mm512_load_pd(GetElementVelocity(nodeId2 + 1))); // V1-V2
			__m512d vecCP2 = _mm512_load_pd(cp2);												  //
			matA02el1 = _mm512_swizzle_pd(vecCP2, _MM_SWIZ_REG_AAAA); //_mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_castpd_si512(_mm512_set1_pd(cp2[4])), _mm512_castpd_si512(_mm512_set1_pd(cp2[0])), 8));
			matA02el2 = _mm512_swizzle_pd(vecCP2, _MM_SWIZ_REG_BBBB); //_mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_castpd_si512(_mm512_set1_pd(cp2[5])), _mm512_castpd_si512(_mm512_set1_pd(cp2[1])), 8));
			matA02el3 = _mm512_swizzle_pd(vecCP2, _MM_SWIZ_REG_CCCC); //_mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_castpd_si512(_mm512_set1_pd(cp2[6])), _mm512_castpd_si512(_mm512_set1_pd(cp2[2])), 8));

			matA02el1 = _mm512_mul_pd(zmm09, matA02el1);
			matA02el2 = _mm512_mul_pd(zmm10, matA02el2);
			matA02el3 = _mm512_mul_pd(zmm11, matA02el3);

			mul1 = _mm512_add_pd(_mm512_add_pd(matA02el1, matA02el2), matA02el3);

			//_mm512_store_pd(tmp, vecDV);
			//_mm512_store_pd(tmp + 8, vecDV2);
			vecDV = _mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_alignr_epi32(_mm512_castpd_si512(vecDV2), _mm512_castpd_si512(vecDV2), 8), _mm512_castpd_si512(vecDV), 8));

			matA02el1 = _mm512_swizzle_pd(vecDV, _MM_SWIZ_REG_AAAA);// _mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_castpd_si512(_mm512_set1_pd(tmp[0])), _mm512_castpd_si512(_mm512_set1_pd(tmp[8])), 8));
			matA02el2 = _mm512_swizzle_pd(vecDV, _MM_SWIZ_REG_BBBB);// _mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_castpd_si512(_mm512_set1_pd(tmp[1])), _mm512_castpd_si512(_mm512_set1_pd(tmp[9])), 8));
			matA02el3 = _mm512_swizzle_pd(vecDV, _MM_SWIZ_REG_CCCC);// _mm512_castsi512_pd(_mm512_alignr_epi32(_mm512_castpd_si512(_mm512_set1_pd(tmp[2])), _mm512_castpd_si512(_mm512_set1_pd(tmp[10])), 8));


			matA02el1 = _mm512_mul_pd(zmm09, matA02el1);
			matA02el2 = _mm512_mul_pd(zmm10, matA02el2);
			matA02el3 = _mm512_mul_pd(zmm11, matA02el3);

			mul2 = _mm512_add_pd(_mm512_add_pd(matA02el1, matA02el2), matA02el3);

			res = _mm512_add_pd(_mm512_add_pd(cp1r, mul1), mul2);

			_mm512_store_pd(velocityStrains, res); // получено VL, линейные компоненты

			double* sp1 = GetElementShiftAngular(nodeId1);
			double* sp2 = GetElementShiftAngular(nodeId2);
			double* vp1 = GetElementVelocityAngular(nodeId1);
			double* vp2 = GetElementVelocityAngular(nodeId2);
			for (size_t i = 0; i < 3; i++)
			{
				shiftStrains[i + 8] = sp1[i] - sp2[i];
				velocityStrains[i + 8] = vp1[i] - vp2[i];
			}
			sp1 = GetElementShiftAngular(nodeId1+1);
			sp2 = GetElementShiftAngular(nodeId2+1);
			vp1 = GetElementVelocityAngular(nodeId1+1);
			vp2 = GetElementVelocityAngular(nodeId2+1);
			for (size_t i = 0; i < 3; i++)
			{
				shiftStrains[i + vecStride+8] = sp1[i] - sp2[i];
				velocityStrains[i + vecStride+8] = vp1[i] - vp2[i];
			}
#endif
		}


}