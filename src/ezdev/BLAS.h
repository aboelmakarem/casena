// BLAS.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 08/01/2018

#ifndef BLAS_H_
#define BLAS_H_

namespace EZ
{
	namespace Math
	{
		class BLAS
		{
		public:
			static void Copy(const unsigned int& size,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y);
			static void Swap(const unsigned int& size,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y);
			static void Scale(const unsigned int& size,const unsigned int& x_inc,double* x,const double& factor);
			static void ScaleAndAdd(const unsigned int& size,const unsigned int& x_inc,double* x,const double& factor,const unsigned int& y_inc,double* y);
			static void ComputeGivensRotation(const double& a,const double& b,double& angle_cosine,double& angle_sine);
			static void ApplyGivensRotation(const unsigned int& size,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y,const double& angle_cosine,const double& angle_sine);
			static double AbsoluteSum(const unsigned int& size,const unsigned int& x_inc,double* x);
			static double DotProduct(const unsigned int& size,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y);
			static double SquaredNorm(const unsigned int& size,const unsigned int& x_inc,double* x);
			static double Norm(const unsigned int& size,const unsigned int& x_inc,double* x);
			static unsigned int MaxAbsoluteIndex(const unsigned int& size,const unsigned int& x_inc,double* x);
			static void MatrixVectorProduct(const unsigned int& row_count,const unsigned int& column_count,const double& alpha,double* matrix,const unsigned int& x_inc,double* x,const double& beta,const unsigned int& y_inc,double* y,const int& operation = 1);
			/*static void TriangularMatrixVectorProduct(const unsigned int& size,double* matrix,const unsigned int& x_inc,double* x,const int& matrix_type = 1);
			static void TriangularMatrixVectorProduct(const unsigned int& size,double* matrix,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y,const int& matrix_type = 1);
			static void SolveTriangularSystem(const unsigned int& size,double* matrix,const unsigned int& x_inc,double* x,const int& matrix_type = 1);*/
			static void DiagonalMatrixVectorProduct(const unsigned int& size,double* matrix,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y,const int& matrix_type = 1);
			static bool Test(const unsigned int& size,const bool& comprehensive,const unsigned int& paranoia = 1);

		private:
			static bool Copy();
			static bool Swap();
			static bool Scale();
			static bool ScaleAndAdd();
			static bool ComputeGivensRotation();
			static bool ApplyGivensRotation();
			static bool AbsoluteSum();
			static bool DotProduct();
			static bool SquaredNorm();
			static bool Norm();
			static bool MaxAbsoluteIndex();
			static bool MatrixVectorProduct();
			/*static bool TriangularMatrixVectorProduct();
			static bool SolveTriangularSystem();*/
			static bool DiagonalMatrixVectorProduct();
			static void VectorRandomPopulate(const unsigned int& size,double* x);
			static void MatrixRandomPopulate(const unsigned int& row_count,const unsigned int& column_count,double* matrix,const bool& non_singular = false);
			static void MakeUpperTriangular(const unsigned int& size,double* matrix);
			static void MakeLowerTriangular(const unsigned int& size,double* matrix);
			static void MakeUnitDiagonal(const unsigned int& size,double* matrix);
			static void SetValue(const unsigned int& size,double* x,const double& value = 0.0);
			static bool CompareVectors(const unsigned int& size,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y,const double& y_factor,const unsigned int& y_shift_inc,double* y_shift = 0);
			static unsigned int test_size;
			static unsigned int test_paranoia;
			static double tolerance;
		};
	}
}

#endif

