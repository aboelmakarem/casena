// BLAS.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 08/01/2018

#include "BLAS.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"

namespace EZ
{
	namespace Math
	{
		unsigned int BLAS::test_size = 10;
		unsigned int BLAS::test_paranoia = 1;
		double BLAS::tolerance = 1.0e-9;
		void BLAS::Copy(const unsigned int& size,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y)
		{
			// This function copies N elements of the array X to the array Y. The increments x_inc and y_inc 
			// are used so that the copied elements can be spaced at different intervals. This enables copying 
			// matrix rows into columns and the like regardless whether the matrix stores its data 
			// row-wise or column wise. 
			if(size == 0)				return;
			if((x_inc == 1) && (y_inc == 1))
			{
				// if we are doing a continuous, non strided, copy, do a fast block memory copy
				memcpy(y,x,size*sizeof(double));
			}
			else
			{
				unsigned int x_index = 0;
				unsigned int y_index = 0;
				for(unsigned int i = 0 ; i < size ; i++)
				{
					y[y_index] = x[x_index];
					x_index += x_inc;
					y_index += y_inc;
				}
			}
		}
		void BLAS::Swap(const unsigned int& size,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y)
		{
			// This function swaps N elements of the array X and Y arrays. The increments x_inc and y_inc 
			// are used so that the swapped elements can be spaced at different intervals. This enables swapping 
			// matrix rows into columns and the like regardless whether the matrix stores its data 
			// row-wise or column wise. 
			if(size == 0)				return;
			double temp = 0.0;
			if((x_inc == 1) && (y_inc == 1))
			{
				// unroll the loops 3 at a time
				unsigned int step = 3;
				unsigned int i = 0;
				unsigned int m = size%step;
				for(i = 0 ; i < m ; i++)
				{
					temp = x[i];
					x[i] = y[i];
					y[i] = temp;
				}
				i = m;
				while(i < size)
				{
					temp = x[i];
					x[i] = y[i];
					y[i] = temp;
					i++;
					temp = x[i];
					x[i] = y[i];
					y[i] = temp;
					i++;
					temp = x[i];
					x[i] = y[i];
					y[i] = temp;
					i++;
				}
			}
			else
			{
				unsigned int x_index = 0;
				unsigned int y_index = 0;
				for(unsigned int i = 0 ; i < size ; i++)
				{
					temp = x[x_index];
					x[x_index] = y[y_index];
					y[y_index] = temp;
					x_index += x_inc;
					y_index += y_inc;
				}
			}
		}
		void BLAS::Scale(const unsigned int& size,const unsigned int& x_inc,double* x,const double& factor)
		{
			// This function scales N elements of the array X spaced x_inc apart by the factor Factor.
			if(size == 0)					return;
			if(x_inc == 0)					return;
			if(fabs(factor - 1.0) < tolerance)		return;
			if(x_inc == 1)
			{
				// unroll the loops 5 at a time
				unsigned int step = 5;
				unsigned int i = 0;
				unsigned int m = size%step;
				for(i = 0 ; i < m ; i++)
				{
					x[i] = x[i]*factor;
				}
				i = m;
				while(i < size)
				{
					x[i] = x[i]*factor;
					x[i + 1] = x[i + 1]*factor;
					x[i + 2] = x[i + 2]*factor;
					x[i + 3] = x[i + 3]*factor;
					x[i + 4] = x[i + 4]*factor;
					i += step;
				}
			}
			else
			{
				unsigned int x_index = 0;
				for(unsigned int i = 0 ; i < size ; i++)
				{
					x[x_index] = x[x_index]*factor;
					x_index += x_inc;
				}
			}
		}
		void BLAS::ScaleAndAdd(const unsigned int& size,const unsigned int& x_inc,double* x,const double& factor,const unsigned int& y_inc,double* y)
		{
			// This function scales N elements of the array X spaced x_inc apart by the factor Factor then adds 
			// the scaled array to another array Y with elements spaced y_inc apart and stores the result in Y.
			if(size == 0)					return;
			if((x_inc == 1) && (y_inc == 1))
			{
				// unroll the loops 4 at a time
				unsigned int step = 4;
				unsigned int i = 0;
				unsigned int m = size%step;
				for(i = 0 ; i < m ; i++)
				{
					y[i] += x[i]*factor;
				}
				i = m;
				while(i < size)
				{
					y[i] += x[i]*factor;
					y[i + 1] += x[i + 1]*factor;
					y[i + 2] += x[i + 2]*factor;
					y[i + 3] += x[i + 3]*factor;
					i += step;
				}
			}
			else
			{
				unsigned int x_index = 0;
				unsigned int y_index = 0;
				for(unsigned int i = 0 ; i < size ; i++)
				{
					y[y_index] += x[x_index]*factor;
					x_index += x_inc;
					y_index += y_inc;
				}
			}
		}
		void BLAS::ComputeGivensRotation(const double& a,const double& b,double& angle_cosine,double& angle_sine)
		{
			// This function computes the Givens rotation that would transform the 
			// vector [A B]^T to the vector [R 0]^T. This can be done by premultiplying 
			// the vector by G where G is an 2x2 orthonormal matrix of the form [c s ; -s c]
			// and its effect is to rotate the vector by some angle theta in its plane.
			// The cosine and sine of the rotation angle are stored in angle_cosine and angle_sine respectively.
			// The solution for C and S are C = A/R and S = B/R where R = sqrt(A^2 + B^2). 
			if(fabs(b) < tolerance)
			{
				// if B is already zero, no rotation is needed
				angle_cosine = 1.0;
				angle_sine = 0.0;
				return;
			}
			// for stability, divide by the larger of the two vector entries
			double r = 0.0;
			double t = 0.0;
			if(fabs(b) > fabs(a))
			{
				t = a/b;
				r = sqrt(1.0 + t*t);
				angle_sine = 1.0/r;
				angle_cosine = angle_sine*t;
			}
			else
			{
				t = b/a;
				r = sqrt(1.0 + t*t);
				angle_cosine = 1.0/r;
				angle_sine = angle_cosine*t;
			}
		}
		void BLAS::ApplyGivensRotation(const unsigned int& size,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y,const double& angle_cosine,const double& angle_sine)
		{
			// This function applies a Givens rotation described by the cosine and the sine of an angle 
			// (C and S) on the two arrays X and Y both of size N with elements spaced x_inc and y_inc apart. 
			// Since applying Givens rotation to a matrix amounts to modifying two of its rows, a typical 
			// use case for this function is to pass the arrays of the two rows that will be modified
			// and they will be modified in place by this function. This is why the indices of the two 
			// rows to be modified is immaterial because the actual arrays, not their indices, will be passed.
			// The specification of x_inc and y_inc allows passing any rows directly even if the matrix is 
			// stored column-wise. 
			if(size == 0)				return;
			double temp = 0.0;
			if((x_inc == 1) && (y_inc == 1))
			{
				for(unsigned int i = 0 ; i < size ; i++)
				{
					temp = angle_cosine*x[i] + angle_sine*y[i];
					y[i] = angle_cosine*y[i] - angle_sine*x[i];
					x[i] = temp;
				}
			}
			else
			{
				unsigned int x_index = 0;
				unsigned int y_index = 0;
				for(unsigned int i = 0 ; i < size ; i++)
				{
					temp = angle_cosine*x[x_index] + angle_sine*y[y_index];
					y[y_index] = angle_cosine*y[y_index] - angle_sine*x[x_index];
					x[x_index] = temp;
					x_index += x_inc;
					y_index += y_inc;
				}
			}
		}
		double BLAS::AbsoluteSum(const unsigned int& size,const unsigned int& x_inc,double* x)
		{
			// This function computes the sum of N absolute values of the array X entries spaced x_inc 
			// apart.
			double sum = 0.0;
			if(size == 0)				return sum;
			if(x_inc == 0)				return sum;
			if(x_inc == 1)
			{
				// unroll the loops 6 at a time
				unsigned int step = 6;
				unsigned int i = 0;
				unsigned int m = size%step;
				for(i = 0 ; i < m ; i++)
				{
					sum += fabs(x[i]);
				}
				i = m;
				while(i < size)
				{
					sum += fabs(x[i]) + fabs(x[i + 1]) + fabs(x[i + 2]) + fabs(x[i + 3]) + fabs(x[i + 4]) + fabs(x[i + 5]);
					i += step;
				}
			}
			else
			{
				unsigned int x_index = 0;
				for(unsigned int i = 0 ; i < size ; i++)
				{
					sum += fabs(x[x_index]);
					x_index += x_inc;
				}
			}
			return sum;
		}
		double BLAS::DotProduct(const unsigned int& size,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y)
		{
			// This function computes the dot product of the arrays X and Y such that the dot product is the 
			// sum of the products of the corresponding array elements spaced x_inc and y_inc elements apart and 
			// the sum runs up to N product pairs. In case of x_inc = y_inc = 1, this is the regular dot product 
			// between two arrays. The specification of x_inc and y_inc is useful so that this function can be 
			// used to compute the dot product between a row and a column vector for example in two matrices. 
			// In this case, x_inc and y_inc can be used to specify the column elements if its matrix is stored 
			// row wise. 
			double product = 0.0;
			if(size == 0)					return product;
			if((x_inc == 1) && (y_inc == 1))
			{
				// unroll the loops 5 at a time
				unsigned int step = 5;
				unsigned int i = 0;
				unsigned int m = size%step;
				for(i = 0 ; i < m ; i++)
				{
					product += x[i]*y[i];
				}
				i = m;
				while(i < size)
				{
					product += x[i]*y[i] + x[i + 1]*y[i + 1] + x[i + 2]*y[i + 2] + x[i + 3]*y[i + 3] + x[i + 4]*y[i + 4];
					i += step;
				}
			}
			else
			{
				unsigned int x_index = 0;
				unsigned int y_index = 0;
				for(unsigned int i = 0 ; i < size ; i++)
				{
					product += x[x_index]*y[y_index];
					x_index += x_inc;
					y_index += y_inc;
				}
			}
			return product;
		}
		double BLAS::SquaredNorm(const unsigned int& size,const unsigned int& x_inc,double* x)
		{
			// This function computes the square of the Euclidean norm of N values of the array X 
			// entries spaced x_inc apart.
			double norm = 0.0;
			if(size == 0)					return norm;
			if(x_inc == 0)					return norm;
			if(size == 1)					return (x[0]*x[0]);
			unsigned int x_index = 0;
			for(unsigned int i = 0 ; i < size ; i++)
			{
				norm += x[x_index]*x[x_index];
				x_index += x_inc;
			}
			return norm;
		}
		double BLAS::Norm(const unsigned int& size,const unsigned int& x_inc,double* x)
		{
			// This function computes the Euclidean norm of N values of the array X 
			// entries spaced x_inc apart.
			return sqrt(SquaredNorm(size,x_inc,x));
		}
		unsigned int BLAS::MaxAbsoluteIndex(const unsigned int& size,const unsigned int& x_inc,double* x)
		{
			// This function finds the index of the first element in the array X that has the maximum 
			// absolute value out of N elements for the elements that are spaced x_inc apart.
			unsigned int index = 0;
			if(size < 2)				return index;
			if(x_inc == 0)				return index;
			double max = 0.0;
			if(x_inc == 1)
			{
				max = fabs(x[0]);
				for(unsigned int i = 1 ; i < size ; i++)
				{
					if(fabs(x[i]) > max)
					{
						max = fabs(x[i]);
						index = i;
					}
				}
			}
			else
			{
				unsigned int x_index = 0;
				max = fabs(x[0]);
				for(unsigned int i = 0 ; i < size ; i++)
				{
					if(fabs(x[x_index]) > max)
					{
						max = fabs(x[x_index]);
						index = i;
					}
					x_index += x_inc;
				}
			}
			return index;
		}
		void BLAS::MatrixVectorProduct(const unsigned int& row_count,const unsigned int& column_count,const double& alpha,double* matrix,const unsigned int& x_inc,double* x,const double& beta,const unsigned int& y_inc,double* y,const int& operation)
		{
			// This function performs the operation 
			// y = alpha*A*x + beta*y if Operation = 1
			// or 
			// y = alpha*A^T*x + beta*y if Operation = 2
			// in both cases, the results are stored in Y. 
			// A is an MxN matrix, x_inc and y_inc are used to set the entry spacing in the X and Y arrays 
			// respectively, alpha and beta are scalars
			if(x_inc == 0)			return;
			if(y_inc == 0)			return;
			if(row_count == 0)		return;
			if(column_count == 0)	return;
			if(operation < 1)		return;
			if(operation > 2)		return;
			// start by forming the product beta*y if needed and store in Y
			if(operation == 2)		Scale(column_count,y_inc,y,beta);
			else					Scale(row_count,y_inc,y,beta);
			// then form the product alpha*A*x or alpha*A^T*x as instructed if needed and add it to Y
			if(fabs(alpha) > tolerance)
			{
				if(operation == 2)
				{
					// form alpha*A^T*x and add it to Y
					if(y_inc == 1)
					{
						for(unsigned int i = 0 ; i < column_count ; i++)
						{
							y[i] += alpha*DotProduct(row_count,column_count,matrix + i,x_inc,x);
						}
					}
					else
					{
						unsigned int y_index = 0;
						for(unsigned int i = 0 ; i < column_count ; i++)
						{
							y[y_index] += alpha*DotProduct(row_count,column_count,matrix + i,x_inc,x);
							y_index += y_inc;
						}
					}
				}
				else
				{
					// form alpha*A*x and add it to Y
					if(y_inc == 1)
					{
						for(unsigned int i = 0 ; i < row_count ; i++)
						{
							y[i] += alpha*DotProduct(column_count,1,matrix + i*column_count,x_inc,x);
						}
					}
					else
					{
						unsigned int y_index = 0;
						for(unsigned int i = 0 ; i < row_count ; i++)
						{
							y[y_index] += alpha*DotProduct(column_count,1,matrix + i*column_count,x_inc,x);
							y_index += y_inc;
						}
					}
				}
			}
		}
		/*void BLAS::TriangularMatrixVectorProduct(const unsigned int& size,double* matrix,const unsigned int& x_inc,double* x,const int& matrix_type)
		{
			// this function is a wrapper for the one that takes in separate input and output 
			// arrays. 
			TriangularMatrixVectorProduct(size,matrix,x_inc,x,x_inc,x,matrix_type);
		}
		void BLAS::TriangularMatrixVectorProduct(const unsigned int& size,double* matrix,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y,const int& matrix_type)
		{
			// This function performs the operation 
			// y <-- A*x 
			// or 
			// y <-- A^T*x
			// depending on the matrix type where A is an nxn triangular matrix and X is a vector
			// x_inc and y_inc are used to set the entry spacing in the X and Y arrays
			// matrix types:
			// 1: upper triangular, don't transpose, with a non-unit diagonal
			// 2: upper triangular, don't transpose, with a unit diagonal
			// 3: lower triangular, don't transpose, with a non-unit diagonal
			// 4: lower triangular, don't transpose, with a unit diagonal
			// 5: upper triangular, transpose, with a non-unit diagonal
			// 6: upper triangular, transpose, with a unit diagonal
			// 7: lower triangular, transpose, with a non-unit diagonal
			// 8: lower triangular, transpose, with a unit diagonal
			if(matrix_type < 1)		return;
			if(matrix_type > 8)		return;
			if(x_inc == 0)			return;
			if(y_inc == 0)			return;
			if(size == 0)			return;
			double temp = 0.0;
			// The multiplication order is as follows, instead of performing the regular row-wise 
			// matrix vector product, we will go over all the columns of the A matrix, multiply the jth 
			// column by the jth scalar in the X vector and accumulate the product in the result. 
			// While this does not save us any operations, we are able to skip computations if the jth 
			// vector element is zero, in which case the entire scalar-column product is skipped. This 
			// works only for triangular matrices though because otherwise, the accumulation would pollute 
			// the output array (which is the same as the input array) before future iterations get to 
			// use the original values needed for their execution. 
			// use the matrix as is
			if(matrix_type < 5)
			{
				// use the matrix as is
				if(matrix_type < 3)
				{
					// this is an upper triangular matrix
					if(x_inc == 1)
					{
						for(unsigned int j = 0 ; j < size ; j++)
						{
							if(fabs(x[j]) < tolerance)		continue;
							// don't we need to clear x[j] before we start accumulation ? 
							// no actually, its old value, after scaling if needed, will be 
							// used as the initial accumulation value since it is in this 
							// iteration that the accumulation will begin by setting the 
							// initial value (or just leaving it as is) and then adding other 
							// values to it in future iterations.
							temp = x[j];
							for(unsigned int i = 0 ; i < j ; i++)
							{
								y[i] += temp*matrix[i*size + j];
							}
							// why aren't we incrementing here ? because this is the first 
							// time this vector element is modified, if the main diagonal 
							// is one, then leave it as is because it will be the initial 
							// value for future accumulations which will begin next iteration, 
							// if it the main diagonal element is not one, scale it and use 
							// it as the start value for future accumulation.
							if(matrix_type < 2)		y[j] = temp*matrix[j*size + j];
							else 					y[j] = temp;
						}
					}
					else
					{
						unsigned int y_row_index = 0;
						unsigned int y_column_index = 0;
						unsigned int x_column_index = 0;
						for(unsigned int j = 0 ; j < size ; j++)
						{
							if(fabs(x[x_column_index]) < tolerance)		continue;
							temp = x[x_column_index];
							y_row_index = 0;
							for(unsigned int i = 0 ; i < j ; i++)
							{
								y[y_row_index] += temp*matrix[i*size + j];
								y_row_index += y_inc;
							}
							if(matrix_type < 2)		y[y_column_index] = temp*matrix[j*size + j];
							else 					y[y_column_index] = temp;
							x_column_index += x_inc;
							y_column_index += y_inc;
						}
					}
				}
				else
				{
					// this is a lower triangular matrix
					if(x_inc == 1)
					{
						for(unsigned int j = (size - 1) ; j >= 0 ; j--)
						{
							if(fabs(x[j]) < tolerance)
							{
								if(j == 0)		break;
								continue;
							}
							temp = x[j];
							for(unsigned int i = (size - 1) ; i > j ; i--)
							{
								y[i] += temp*matrix[i*size + j];
							}
							if(matrix_type < 4)		y[j] = temp*matrix[j*size + j];
							else 					y[j] = temp;
							if(j == 0)		break;
						}
					}
					else
					{
						unsigned int y_row_index = 0;
						unsigned int y_column_index = (size - 1)*y_inc;
						unsigned int x_column_index = (size - 1)*x_inc;
						for(unsigned int j = (size - 1) ; j >= 0 ; j--)
						{
							if(fabs(x[x_column_index]) < tolerance)
							{
								if(j == 0)		break;
								continue;
							}
							temp = x[x_column_index];
							y_row_index = (size - 1)*y_inc;
							for(unsigned int i = (size - 1) ; i > j ; i--)
							{
								y[y_row_index] += temp*matrix[i*size + j];
								y_row_index -= y_inc;
							}
							if(matrix_type < 4)		y[y_column_index] = temp*matrix[j*size + j];
							else 					y[y_column_index] = temp;
							if(j == 0)		break;
							x_column_index -= x_inc;
							y_column_index -= y_inc;
						}
					}
				}
			}
			else
			{
				// use the matrix transpose
				if(matrix_type < 7)
				{
					// this is an upper triangular matrix
					if(x_inc == 1)
					{
						for(unsigned int j = (size - 1) ; j >= 0 ; j--)
						{
							if(fabs(x[j]) < tolerance)
							{
								if(j == 0)		break;
								continue;
							}
							temp = x[j];
							if(matrix_type < 6)		y[j] = temp*matrix[j*size + j];
							else 					y[j] = temp;
							for(unsigned int i = (size - 1) ; i > j ; i--)
							{
								y[i] += temp*matrix[j*size + i];
							}
							if(j == 0)				break;
						}
					}
					else
					{
						unsigned int y_row_index = 0;
						unsigned int y_column_index = (size - 1)*y_inc;
						unsigned int x_column_index = (size - 1)*x_inc;
						for(unsigned int j = (size - 1) ; j >= 0 ; j--)
						{
							if(fabs(x[x_column_index]) < tolerance)
							{
								if(j == 0)		break;
								continue;
							}
							temp = x[x_column_index];
							if(matrix_type < 6)		y[y_column_index] = temp*matrix[j*size + j];
							else 					y[y_column_index] = temp;
							y_row_index = (size - 1)*y_inc;
							for(unsigned int i = (size - 1) ; i > j ; i--)
							{
								y[y_row_index] += temp*matrix[j*size + i];
								y_row_index -= y_inc;
							}
							if(j == 0)				break;
							x_column_index -= x_inc;
							y_column_index -= y_inc;
						}
					}
				}
				else
				{
					// this is a lower triangular matrix
					if(x_inc == 1)
					{
						for(unsigned int j = 0 ; j < size ; j++)
						{
							if(fabs(x[j]) < tolerance)		continue;
							temp = x[j];
							for(unsigned int i = 0 ; i < j ; i++)
							{
								y[i] += temp*matrix[j*size + i];
							}
							if(matrix_type < 8)		y[j] = temp*matrix[j*size + j];
							else 					y[j] = temp;
						}
					}
					else
					{
						unsigned int y_row_index = 0;
						unsigned int y_column_index = 0;
						unsigned int x_column_index = 0;
						for(unsigned int j = 0 ; j < size ; j++)
						{
							if(fabs(x[x_column_index]) < tolerance)		continue;
							temp = x[x_column_index];
							y_row_index = 0;
							for(unsigned int i = 0 ; i < j ; i++)
							{
								y[y_row_index] += temp*matrix[j*size + i];
								y_row_index += y_inc;
							}
							if(matrix_type < 8)		y[y_column_index] = temp*matrix[j*size + j];
							else 					y[y_column_index] = temp;
							x_column_index += x_inc;
							y_column_index += y_inc;
						}
					}
				}
			}
		}
		void BLAS::SolveTriangularSystem(const unsigned int& size,double* matrix,const unsigned int& x_inc,double* x,const int& matrix_type)
		{
			// This function solves the system of equations of the form
			// A*x = b
			// or 
			// A^T*x = b
			// depending on the matrix type where A is an nxn triangular matrix and b is a vector. 
			// The X input array holds the values for the b vector on entry and is overwritten with 
			// the solution on exit. x_inc is used to set the entry spacing in the X array.
			// matrix types:
			// 1: upper triangular, don't transpose, with a non-unit diagonal
			// 2: upper triangular, don't transpose, with a unit diagonal
			// 3: lower triangular, don't transpose, with a non-unit diagonal
			// 4: lower triangular, don't transpose, with a unit diagonal
			// 5: upper triangular, transpose, with a non-unit diagonal
			// 6: upper triangular, transpose, with a unit diagonal
			// 7: lower triangular, transpose, with a non-unit diagonal
			// 8: lower triangular, transpose, with a unit diagonal
			if(matrix_type < 1)		return;
			if(matrix_type > 8)		return;
			if(x_inc == 0)			return;
			if(size == 0)			return;
			if(matrix_type < 5)
			{
				// use the matrix as is
				if(matrix_type < 3)
				{
					// this is an upper triangular matrix
					if(x_inc == 1)
					{
						for(int j = (size - 1) ; j >= 0 ; j--)
						{
							if(fabs(x[j]) < tolerance)		continue;
							if(matrix_type < 2)			x[j] = x[j]/matrix[j*size + j];
							for(int i = (j - 1) ; i >= 0 ; i--)
							{
								x[i] -= matrix[i*size + j]*x[j];
							}
						}
					}
					else
					{
						unsigned int x_row_index = 0;
						unsigned int x_column_index = (size - 1)*x_inc;
						for(int j = (size - 1) ; j >= 0 ; j--)
						{
							if(fabs(x[x_column_index]) < tolerance)		continue;
							if(matrix_type < 2)			x[x_column_index] = x[x_column_index]/matrix[j*size + j];
							x_row_index = (j - 1)*x_inc;
							for(int i = (j - 1) ; i >= 0 ; i--)
							{
								x[x_row_index] -= matrix[i*size + j]*x[x_column_index];
								x_row_index -= x_inc;
							}
							x_column_index -= x_inc;
						}
					}
				}
				else
				{
					// this is a lower triangular matrix
					if(x_inc == 1)
					{
						for(unsigned int j = 0 ; j < size ; j++)
						{
							if(fabs(x[j]) < tolerance)		continue;
							if(matrix_type < 4)			x[j] = x[j]/matrix[j*size + j];
							for(unsigned int i = (j + 1) ; i < size ; i++)
							{
								x[i] -= matrix[i*size + j]*x[j];
							}
						}
					}
					else
					{
						unsigned int x_row_index = 0;
						unsigned int x_column_index = 0;
						for(unsigned int j = 0 ; j < size ; j++)
						{
							if(fabs(x[x_column_index]) < tolerance)	continue;
							if(matrix_type < 4)			x[x_column_index] = x[x_column_index]/matrix[j*size + j];
							x_row_index = (j + 1)*x_inc;
							for(unsigned int i = (j + 1) ; i < size ; i++)
							{
								x[x_row_index] -= matrix[i*size + j]*x[x_column_index];
								x_row_index += x_inc;
							}
							x_column_index += x_inc;
						}
					}
				}
			}
			else
			{
				// use the matrix transpose
				if(matrix_type < 7)
				{
					// this is an upper triangular matrix
					if(x_inc == 1)
					{
						for(unsigned int j = 0 ; j < size ; j++)
						{
							if(fabs(x[j]) < tolerance)			continue;
							if(matrix_type < 6)			x[j] = x[j]/matrix[j*size + j];
							for(unsigned int i = (j + 1) ; i < size ; i++)
							{
								x[i] -= matrix[j*size + i]*x[j];
							}
						}
					}
					else
					{
						unsigned int x_row_index = 0;
						unsigned int x_column_index = 0;
						for(unsigned int j = 0 ; j < size ; j++)
						{
							if(fabs(x[x_column_index]) < tolerance)		continue;
							if(matrix_type < 6)			x[x_column_index] = x[x_column_index]/matrix[j*size + j];
							x_row_index = (j + 1)*x_inc;
							for(unsigned int i = (j + 1) ; i < size ; i++)
							{
								x[x_row_index] -= matrix[j*size + i]*x[x_column_index];
								x_row_index += x_inc;
							}
							x_column_index += x_inc;
						}
					}
				}
				else
				{
					// this is a lower triangular matrix
					if(x_inc == 1)
					{
						for(int j = (size - 1) ; j >= 0 ; j--)
						{
							if(fabs(x[j]) < tolerance)		continue;
							if(matrix_type < 8)			x[j] = x[j]/matrix[j*size + j];
							for(int i = (j - 1) ; i >= 0 ; i--)
							{
								x[i] -= matrix[j*size + i]*x[j];
							}
						}
					}
					else
					{
						unsigned int x_row_index = 0;
						unsigned int x_column_index = (size - 1)*x_inc;
						for(int j = (size - 1) ; j >= 0 ; j--)
						{
							if(fabs(x[x_column_index]) < tolerance)		continue;
							if(matrix_type < 8)			x[x_column_index] = x[x_column_index]/matrix[j*size + j];
							x_row_index = (j - 1)*x_inc;
							for(int i = (j - 1) ; i >= 0 ; i--)
							{
								x[x_row_index] -= matrix[j*size + i]*x[x_column_index];
								x_row_index -= x_inc;
							}
							x_column_index -= x_inc;
						}
					}
				}
			}
		}*/
		void BLAS::DiagonalMatrixVectorProduct(const unsigned int& size,double* matrix,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y,const int& matrix_type)
		{
			// This function multiplies the diagonal matrix A with the column  vector X and stores the 
			// result in the column vector Y where A is an nxn triangular matrix and X is a vector
			// x_inc and y_inc are used to set the entry spacing in the X and Y arrays. If the matrix 
			// type is 1, then the full matrix is passed in A, if it is 2, then A is just an array of 
			// the diagonal entries of the matrix A
			if(x_inc == 0)			return;
			if(y_inc == 0)			return;
			if(size == 0)			return;
			if(matrix_type < 1)		return;
			if(matrix_type > 2)		return;
			if(matrix_type == 1)
			{
				// A is the full matrix
				if(x_inc == 1)
				{
					for(unsigned int j = 0 ; j < size ; j++)
					{
						y[j] = x[j]*matrix[j*size + j];
					}
				}
				else
				{
					unsigned int y_column_index = 0;
					unsigned int x_column_index = 0;
					for(unsigned int j = 0 ; j < size ; j++)
					{
						y[y_column_index] = x[x_column_index]*matrix[j*size + j];
						x_column_index += x_inc;
						y_column_index += y_inc;
					}
				}
			}
			else if(matrix_type == 2)
			{
				// A is the array of diagonal entries
				if(x_inc == 1)
				{
					for(unsigned int j = 0 ; j < size ; j++)
					{
						y[j] = x[j]*matrix[j];
					}
				}
				else
				{
					unsigned int y_column_index = 0;
					unsigned int x_column_index = 0;
					for(unsigned int j = 0 ; j < size ; j++)
					{
						y[y_column_index] = x[x_column_index]*matrix[j];
						x_column_index += x_inc;
						y_column_index += y_inc;
					}
				}
			}
		}
		bool BLAS::Test(const unsigned int& size,const bool& comprehensive,const unsigned int& paranoia)
		{
			test_size = size;
			test_paranoia = paranoia;
			bool pass = true;
			if(comprehensive)
			{
				if(!Copy())
				{
					printf("BLAS Copy failed test\n");
					pass = false;
				}
				if(!Swap())
				{
					printf("BLAS Swap failed test\n");
					pass = false;
				}
				if(!Scale())
				{
					printf("BLAS Scale failed test\n");
					pass = false;
				}
				if(!ScaleAndAdd())
				{
					printf("BLAS ScaleAndAdd failed test\n");
					pass = false;
				}
				if(!ComputeGivensRotation())
				{
					printf("BLAS ComputeGivensRotation failed test\n");
					pass = false;
				}
				if(!ApplyGivensRotation())
				{
					printf("BLAS ApplyGivensRotation failed test\n");
					pass = false;
				}
				if(!AbsoluteSum())
				{
					printf("BLAS AbsoluteSum failed test\n");
					pass = false;
				}
				if(!DotProduct())
				{
					printf("BLAS DotProduct failed test\n");
					pass = false;
				}
				if(!SquaredNorm())
				{
					printf("BLAS SquaredNorm failed test\n");
					pass = false;
				}
				if(!Norm())
				{
					printf("BLAS Norm failed test\n");
					pass = false;
				}
				if(!MaxAbsoluteIndex())
				{
					printf("BLAS MaxAbsoluteIndex failed test\n");
					pass = false;
				}
				if(!MatrixVectorProduct())
				{
					printf("BLAS MatrixVectorProduct failed test\n");
					pass = false;
				}
				/*if(!TriangularMatrixVectorProduct())
				{
					printf("BLAS TriangularMatrixVectorProduct failed test\n");
					pass = false;
				}
				if(!SolveTriangularSystem())
				{
					printf("BLAS SolveTriangularSystem failed test\n");
					pass = false;
				}*/
				if(!DiagonalMatrixVectorProduct())
				{
					printf("BLAS DiagonalMatrixVectorProduct failed test\n");
					pass = false;
				}
			}
			if(!DiagonalMatrixVectorProduct())
			{
				printf("BLAS DiagonalMatrixVectorProduct failed test\n");
				pass = false;
			}
			return pass;
		}
		bool BLAS::Copy()
		{
			double* source = 0;
			double* destination = 0;
			bool pass = true;
			// test unit increment copy in source and destination
			unsigned int source_stride = 1;
			unsigned int destination_stride = 1;
			source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,destination_stride,destination);
				pass = CompareVectors(test_size,source_stride,source,destination_stride,destination,1.0,0,0);
				if(!pass)			break;
			}
			delete [] source;
			delete [] destination;
			if(!pass)				return false;
			// test unit increment copy in source but not in destination, generate a random number 
			// for strides but keep it under 10.
			source_stride = 1;
			while(destination_stride <= 1)
			{
				destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,destination_stride,destination);
				pass = CompareVectors(test_size,source_stride,source,destination_stride,destination,1.0,0,0);
				if(!pass)			break;
			}
			delete [] source;
			delete [] destination;
			if(!pass)				return false;
			// test strided increment copy in source but unit increment in destination, generate a random number 
			// for strides but keep it under 10.
			while(source_stride <= 1)
			{
				source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			destination_stride = 1;
			source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,destination_stride,destination);
				pass = CompareVectors(test_size,source_stride,source,destination_stride,destination,1.0,0,0);
				if(!pass)			break;
			}
			delete [] source;
			delete [] destination;
			if(!pass)				return false;
			// test strided increment copy in both source and destination, generate a random number 
			// for strides but keep it under 10.
			source_stride = 1;
			destination_stride = 1;
			while(source_stride <= 1)
			{
				source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			while(destination_stride <= 1)
			{
				destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,destination_stride,destination);
				pass = CompareVectors(test_size,source_stride,source,destination_stride,destination,1.0,0,0);
				if(!pass)			break;
			}
			delete [] source;
			delete [] destination;
			return pass;
		}
		bool BLAS::Swap()
		{
			// This function uses Copy(). It is assumed that it passed the test. 
			double* source = 0;
			double* destination = 0;
			double* original_source = 0;
			double* original_destination = 0;
			bool pass = true;
			// test unit increment swap in source and destination
			unsigned int source_stride = 1;
			unsigned int destination_stride = 1;
			source = new double[test_size*source_stride];
			original_source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			original_destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,source_stride,original_source);
				VectorRandomPopulate(test_size*destination_stride,destination);
				Copy(test_size,destination_stride,destination,destination_stride,original_destination);
				Swap(test_size,source_stride,source,destination_stride,destination);
				pass = CompareVectors(test_size,source_stride,source,destination_stride,original_destination,1.0,0,0);
				pass = CompareVectors(test_size,destination_stride,destination,source_stride,original_source,1.0,0,0);
				if(!pass)			break;
			}
			delete [] source;
			delete [] original_source;
			delete [] destination;
			delete [] original_destination;
			if(!pass)				return false;
			// test unit increment swap in source but not in destination, generate a random number 
			// for strides but keep it under 10.
			source_stride = 1;
			while(destination_stride <= 1)
			{
				destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			source = new double[test_size*source_stride];
			original_source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			original_destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,source_stride,original_source);
				VectorRandomPopulate(test_size*destination_stride,destination);
				Copy(test_size,destination_stride,destination,destination_stride,original_destination);
				Swap(test_size,source_stride,source,destination_stride,destination);
				pass = CompareVectors(test_size,source_stride,source,destination_stride,original_destination,1.0,0,0);
				pass = CompareVectors(test_size,destination_stride,destination,source_stride,original_source,1.0,0,0);
				if(!pass)			break;
			}
			delete [] source;
			delete [] original_source;
			delete [] destination;
			delete [] original_destination;
			if(!pass)				return false;
			// test strided increment swap in source but unit increment in destination, generate a random number 
			// for strides but keep it under 10.
			while(source_stride <= 1)
			{
				source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			destination_stride = 1;
			source = new double[test_size*source_stride];
			original_source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			original_destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,source_stride,original_source);
				VectorRandomPopulate(test_size*destination_stride,destination);
				Copy(test_size,destination_stride,destination,destination_stride,original_destination);
				Swap(test_size,source_stride,source,destination_stride,destination);
				pass = CompareVectors(test_size,source_stride,source,destination_stride,original_destination,1.0,0,0);
				pass = CompareVectors(test_size,destination_stride,destination,source_stride,original_source,1.0,0,0);
				if(!pass)			break;
			}
			delete [] source;
			delete [] original_source;
			delete [] destination;
			delete [] original_destination;
			if(!pass)				return false;
			// test strided increment swap in both source and destination, generate a random number 
			// for strides but keep it under 10.
			source_stride = 1;
			destination_stride = 1;
			while(source_stride <= 1)
			{
				source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			while(destination_stride <= 1)
			{
				destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			source = new double[test_size*source_stride];
			original_source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			original_destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,source_stride,original_source);
				VectorRandomPopulate(test_size*destination_stride,destination);
				Copy(test_size,destination_stride,destination,destination_stride,original_destination);
				Swap(test_size,source_stride,source,destination_stride,destination);
				pass = CompareVectors(test_size,source_stride,source,destination_stride,original_destination,1.0,0,0);
				pass = CompareVectors(test_size,destination_stride,destination,source_stride,original_source,1.0,0,0);
				if(!pass)			break;
			}
			delete [] source;
			delete [] original_source;
			delete [] destination;
			delete [] original_destination;
			return pass;
		}
		bool BLAS::Scale()
		{
			// This function uses Copy(). It is assumed that it passed the test. 
			double* working_array = 0;
			double* original_array = 0;
			double factor = 1.0;
			bool pass = true;
			// factors will range from -100.0 to 100.0
			// test unit increment scaling
			unsigned int stride = 1;
			working_array = new double[test_size*stride];
			original_array = new double[test_size];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*stride,working_array);
				factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				Copy(test_size,stride,working_array,1,original_array);
				Scale(test_size,stride,working_array,factor);
				pass = CompareVectors(test_size,stride,working_array,1,original_array,factor,0,0);
				if(!pass)
				{
					printf("BLAS Scale test failed for factor %e\n",factor);
					break;
				}
			}
			delete [] working_array;
			if(!pass)
			{
				delete [] original_array;
				return false;
			}
			// test strided increment scaling, generate a random number for strides but keep it
			// under 10.
			while(stride <= 1)
			{
				stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			working_array = new double[test_size*stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*stride,working_array);
				factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				Copy(test_size,stride,working_array,1,original_array);
				Scale(test_size,stride,working_array,factor);
				pass = CompareVectors(test_size,stride,working_array,1,original_array,factor,0,0);
				if(!pass)
				{
					printf("BLAS Scale test failed for factor %e\n",factor);
					break;
				}
			}
			delete [] working_array;
			delete [] original_array;
			return pass;
		}
		bool BLAS::ScaleAndAdd()
		{
			// This function uses Copy(). It is assumed that it passed the test. 
			double* source = 0;
			double* destination = 0;
			double* original_source = 0;
			double* original_destination = 0;
			double factor = 0.0;
			bool pass = true;
			// factors will range from -100.0 to 100.0
			// test unit increment scale and add in source and destination
			unsigned int source_stride = 1;
			unsigned int destination_stride = 1;
			source = new double[test_size*source_stride];
			original_source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			original_destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,source_stride,original_source);
				VectorRandomPopulate(test_size*destination_stride,destination);
				Copy(test_size,destination_stride,destination,destination_stride,original_destination);
				factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				ScaleAndAdd(test_size,source_stride,source,factor,destination_stride,destination);
				pass = CompareVectors(test_size,destination_stride,destination,source_stride,original_source,factor,destination_stride,original_destination);
				if(!pass)			break;
			}
			delete [] source;
			delete [] original_source;
			delete [] destination;
			delete [] original_destination;
			if(!pass)				return false;
			// test unit increment scale and add in source but not in destination, generate a random number 
			// for strides but keep it under 10.
			source_stride = 1;
			while(destination_stride <= 1)
			{
				destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			source = new double[test_size*source_stride];
			original_source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			original_destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,source_stride,original_source);
				VectorRandomPopulate(test_size*destination_stride,destination);
				Copy(test_size,destination_stride,destination,destination_stride,original_destination);
				factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				ScaleAndAdd(test_size,source_stride,source,factor,destination_stride,destination);
				pass = CompareVectors(test_size,destination_stride,destination,source_stride,original_source,factor,destination_stride,original_destination);
				if(!pass)			break;
			}
			delete [] source;
			delete [] original_source;
			delete [] destination;
			delete [] original_destination;
			if(!pass)				return false;
			// test strided increment scale and add in source but unit increment in destination, generate a random number 
			// for strides but keep it under 10.
			while(source_stride <= 1)
			{
				source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			destination_stride = 1;
			source = new double[test_size*source_stride];
			original_source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			original_destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,source_stride,original_source);
				VectorRandomPopulate(test_size*destination_stride,destination);
				Copy(test_size,destination_stride,destination,destination_stride,original_destination);
				factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				ScaleAndAdd(test_size,source_stride,source,factor,destination_stride,destination);
				pass = CompareVectors(test_size,destination_stride,destination,source_stride,original_source,factor,destination_stride,original_destination);
				if(!pass)			break;
			}
			delete [] source;
			delete [] original_source;
			delete [] destination;
			delete [] original_destination;
			if(!pass)				return false;
			// test strided increment scale and add in both source and destination, generate a random number 
			// for strides but keep it under 10.
			source_stride = 1;
			destination_stride = 1;
			while(source_stride <= 1)
			{
				source_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			while(destination_stride <= 1)
			{
				destination_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			source = new double[test_size*source_stride];
			original_source = new double[test_size*source_stride];
			destination = new double[test_size*destination_stride];
			original_destination = new double[test_size*destination_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*source_stride,source);
				Copy(test_size,source_stride,source,source_stride,original_source);
				VectorRandomPopulate(test_size*destination_stride,destination);
				Copy(test_size,destination_stride,destination,destination_stride,original_destination);
				factor = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				ScaleAndAdd(test_size,source_stride,source,factor,destination_stride,destination);
				pass = CompareVectors(test_size,destination_stride,destination,source_stride,original_source,factor,destination_stride,original_destination);
				if(!pass)			break;
			}
			delete [] source;
			delete [] original_source;
			delete [] destination;
			delete [] original_destination;
			return pass;
		}
		bool BLAS::ComputeGivensRotation()
		{
			bool pass = true;
			double a = 0.0;
			double b = 0.0;
			double angle_cosine = 0.0;
			double angle_sine = 0.0;
			double x = 0.0;
			double y = 0.0;
			double r = 0.0;
			// generate A and B in the range -100.0 to 100.0
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				a = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				b = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				r = sqrt(a*a + b*b);
				ComputeGivensRotation(a,b,angle_cosine,angle_sine);
				x = angle_cosine*a + angle_sine*b;
				y = -angle_sine*a + angle_cosine*b;
				if(fabs(y) > tolerance)					pass = false;
				if(fabs(fabs(x) - r) > tolerance)				pass = false;
				r = angle_cosine*angle_cosine + angle_sine*angle_sine;
				if(fabs(r - 1.0) > tolerance)			pass = false;
				if(!pass)									break;
			}
			return pass;
		}
		bool BLAS::ApplyGivensRotation()
		{
			// ApplyGivensRotation() passes if ComputeGivensRotation() passes
			return ComputeGivensRotation();
		}
		bool BLAS::AbsoluteSum()
		{
			// The test for AbsoluteSum will be another implementation of AbsoluteSum
			return true;
		}
		bool BLAS::DotProduct()
		{
			double* x = 0;
			double* y = 0;
			double product = 0.0;
			bool pass = true;
			// test unit increment dot product in X and Y
			unsigned int x_stride = 1;
			unsigned int y_stride = 1;
			x = new double[test_size*x_stride];
			y = new double[test_size*y_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*x_stride,x);
				VectorRandomPopulate(test_size*y_stride,y);
				product = DotProduct(test_size,x_stride,x,y_stride,y);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					product -= x[j*x_stride]*y[j*y_stride];
				}
				pass = (fabs(product) < tolerance);
				if(!pass)
				{
					printf("BLAS DotProduct test failed with error : %e\n",product);
					break;
				}
			}
			delete [] x;
			delete [] y;
			if(!pass)				return false;
			// test unit increment dot product in X but not in Y, generate a random number 
			// for strides but keep it under 10.
			x_stride = 1;
			while(y_stride <= 1)
			{
				y_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			x = new double[test_size*x_stride];
			y = new double[test_size*y_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*x_stride,x);
				VectorRandomPopulate(test_size*y_stride,y);
				product = DotProduct(test_size,x_stride,x,y_stride,y);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					product -= x[j*x_stride]*y[j*y_stride];
				}
				pass = (fabs(product) < tolerance);
				if(!pass)
				{
					printf("BLAS DotProduct test failed with error : %e\n",product);
					break;
				}
			}
			delete [] x;
			delete [] y;
			if(!pass)				return false;
			// test strided increment dot product in X but unit increment in Y, generate a random number 
			// for strides but keep it under 10.
			while(x_stride <= 1)
			{
				x_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			y_stride = 1;
			x = new double[test_size*x_stride];
			y = new double[test_size*y_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*x_stride,x);
				VectorRandomPopulate(test_size*y_stride,y);
				product = DotProduct(test_size,x_stride,x,y_stride,y);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					product -= x[j*x_stride]*y[j*y_stride];
				}
				pass = (fabs(product) < tolerance);
				if(!pass)
				{
					printf("BLAS DotProduct test failed with error : %e\n",product);
					break;
				}
			}
			delete [] x;
			delete [] y;
			if(!pass)				return false;
			// test strided increment dot product in both X and Y, generate a random number 
			// for strides but keep it under 10.
			x_stride = 1;
			y_stride = 1;
			while(x_stride <= 1)
			{
				x_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			while(y_stride <= 1)
			{
				y_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			x = new double[test_size*x_stride];
			y = new double[test_size*y_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*x_stride,x);
				VectorRandomPopulate(test_size*y_stride,y);
				product = DotProduct(test_size,x_stride,x,y_stride,y);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					product -= x[j*x_stride]*y[j*y_stride];
				}
				pass = (fabs(product) < tolerance);
				if(!pass)
				{
					printf("BLAS DotProduct test failed with error : %e\n",product);
					break;
				}
			}
			delete [] x;
			delete [] y;
			return pass;
		}
		bool BLAS::SquaredNorm()
		{
			// This function uses DotProduct(). It is assumed that it passed the test. 
			double* x = 0;
			double squared_norm = 0.0;
			double product = 0.0;
			bool pass = true;
			// test unit increment array squared norm
			unsigned int stride = 1;
			x = new double[test_size*stride]; 
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*stride,x);
				squared_norm = SquaredNorm(test_size,stride,x);
				product = DotProduct(test_size,stride,x,stride,x);
				pass = fabs(squared_norm - product) < tolerance;
				if(!pass)				break;
			}
			delete [] x;
			if(!pass)				return false;
			// test strided array squared norm
			stride = 1;
			while(stride <= 1)
			{
				stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			x = new double[test_size*stride]; 
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*stride,x);
				squared_norm = SquaredNorm(test_size,stride,x);
				product = DotProduct(test_size,stride,x,stride,x);
				pass = fabs(squared_norm - product) < tolerance;
				if(!pass)				break;
			}
			delete [] x;
			return pass;
		}
		bool BLAS::Norm()
		{
			// This function uses SquaredNorm(). It is assumed that it passed the test. 
			double* x = 0;
			double squared_norm = 0.0;
			double norm = 0.0;
			bool pass = true;
			// test unit increment array squared norm
			unsigned int stride = 1;
			x = new double[test_size*stride]; 
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*stride,x);
				squared_norm = SquaredNorm(test_size,stride,x);
				norm = Norm(test_size,stride,x);
				pass = fabs(squared_norm - norm*norm) < tolerance;
				if(!pass)				break;
			}
			delete [] x;
			if(!pass)				return false;
			// test strided array squared norm
			stride = 1;
			while(stride <= 1)
			{
				stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			x = new double[test_size*stride]; 
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				VectorRandomPopulate(test_size*stride,x);
				squared_norm = SquaredNorm(test_size,stride,x);
				norm = Norm(test_size,stride,x);
				pass = fabs(squared_norm - norm*norm) < tolerance;
				if(!pass)				break;
			}
			delete [] x;
			return pass;
		}
		bool BLAS::MaxAbsoluteIndex()
		{
			// The test for MaxAbsoluteIndex will be another implementation of MaxAbsoluteIndex
			return true;
		}
		bool BLAS::MatrixVectorProduct()
		{
			// This function uses Copy() and DotProduct(). It is assumed that
			// they passed the test. 
			double* A = 0;
			double* x = 0;
			double* y = 0;
			double* original_y = 0;
			double alpha = 0.0;
			double beta = 0.0;
			double target = 0.0;
			bool pass = true;
			// alpha and beta factors will be chosen in the range -100.0 to 100.0
			// test unit increment matrix vector product in X and Y
			unsigned int x_stride = 1;
			unsigned int y_stride = 1;
			A = new double[test_size*test_size];
			x = new double[test_size*x_stride];
			y = new double[test_size*y_stride];
			original_y = new double[test_size*y_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				MatrixRandomPopulate(test_size,test_size,A);
				VectorRandomPopulate(test_size*x_stride,x);
				VectorRandomPopulate(test_size*y_stride,y);
				Copy(test_size,y_stride,y,y_stride,original_y);
				alpha = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				beta = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				// run both no-transpose and transpose tests
				MatrixVectorProduct(test_size,test_size,alpha,A,x_stride,x,beta,y_stride,y,1);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					target = alpha*DotProduct(test_size,1,A + j*test_size,x_stride,x) + beta*original_y[j*y_stride];
					if(fabs(target - y[j*y_stride]) > tolerance)
					{
						pass = false;
						break;
					}
				}
				if(!pass)			break;
				// generate a new Y vector because the old one was overwritten
				VectorRandomPopulate(test_size*y_stride,y);
				Copy(test_size,y_stride,y,y_stride,original_y);
				MatrixVectorProduct(test_size,test_size,alpha,A,x_stride,x,beta,y_stride,y,2);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					target = alpha*DotProduct(test_size,test_size,A + j,x_stride,x) + beta*original_y[j*y_stride];
					if(fabs(target - y[j*y_stride]) > tolerance)
					{
						pass = false;
						break;
					}
				}
				if(!pass)			break;
			}
			delete [] x;
			delete [] y;
			delete [] original_y;
			if(!pass)
			{
				delete [] A;
				return false;
			}
			// test unit increment matrix vector product in X but not in Y, generate a random number 
			// for strides but keep it under 10.
			x_stride = 1;
			while(y_stride <= 1)
			{
				y_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			x = new double[test_size*x_stride];
			y = new double[test_size*y_stride];
			original_y = new double[test_size*y_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				MatrixRandomPopulate(test_size,test_size,A);
				VectorRandomPopulate(test_size*x_stride,x);
				VectorRandomPopulate(test_size*y_stride,y);
				Copy(test_size,y_stride,y,y_stride,original_y);
				alpha = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				beta = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				// run both no-transpose and transpose tests
				MatrixVectorProduct(test_size,test_size,alpha,A,x_stride,x,beta,y_stride,y,1);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					target = alpha*DotProduct(test_size,1,A + j*test_size,x_stride,x) + beta*original_y[j*y_stride];
					if(fabs(target - y[j*y_stride]) > tolerance)
					{
						pass = false;
						break;
					}
				}
				if(!pass)			break;
				// generate a new Y vector because the old one was overwritten
				VectorRandomPopulate(test_size*y_stride,y);
				Copy(test_size,y_stride,y,y_stride,original_y);
				MatrixVectorProduct(test_size,test_size,alpha,A,x_stride,x,beta,y_stride,y,2);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					target = alpha*DotProduct(test_size,test_size,A + j,x_stride,x) + beta*original_y[j*y_stride];
					if(fabs(target - y[j*y_stride]) > tolerance)
					{
						pass = false;
						break;
					}
				}
				if(!pass)			break;
			}
			delete [] x;
			delete [] y;
			delete [] original_y;
			if(!pass)
			{
				delete [] A;
				return false;
			}
			// test strided increment matrix vector product in X but unit increment in Y, generate a random number 
			// for strides but keep it under 10.
			while(x_stride <= 1)
			{
				x_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			y_stride = 1;
			x = new double[test_size*x_stride];
			y = new double[test_size*y_stride];
			original_y = new double[test_size*y_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				MatrixRandomPopulate(test_size,test_size,A);
				VectorRandomPopulate(test_size*x_stride,x);
				VectorRandomPopulate(test_size*y_stride,y);
				Copy(test_size,y_stride,y,y_stride,original_y);
				alpha = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				beta = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				// run both no-transpose and transpose tests
				MatrixVectorProduct(test_size,test_size,alpha,A,x_stride,x,beta,y_stride,y,1);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					target = alpha*DotProduct(test_size,1,A + j*test_size,x_stride,x) + beta*original_y[j*y_stride];
					if(fabs(target - y[j*y_stride]) > tolerance)
					{
						pass = false;
						break;
					}
				}
				if(!pass)			break;
				// generate a new Y vector because the old one was overwritten
				VectorRandomPopulate(test_size*y_stride,y);
				Copy(test_size,y_stride,y,y_stride,original_y);
				MatrixVectorProduct(test_size,test_size,alpha,A,x_stride,x,beta,y_stride,y,2);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					target = alpha*DotProduct(test_size,test_size,A + j,x_stride,x) + beta*original_y[j*y_stride];
					if(fabs(target - y[j*y_stride]) > tolerance)
					{
						pass = false;
						break;
					}
				}
				if(!pass)			break;
			}
			delete [] x;
			delete [] y;
			delete [] original_y;
			if(!pass)
			{
				delete [] A;
				return false;
			}
			// test strided increment copy in both source and destination, generate a random number 
			// for strides but keep it under 10.
			x_stride = 1;
			y_stride = 1;
			while(x_stride <= 1)
			{
				x_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			while(y_stride <= 1)
			{
				y_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);
			}
			x = new double[test_size*x_stride];
			y = new double[test_size*y_stride];
			original_y = new double[test_size*y_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				MatrixRandomPopulate(test_size,test_size,A);
				VectorRandomPopulate(test_size*x_stride,x);
				VectorRandomPopulate(test_size*y_stride,y);
				Copy(test_size,y_stride,y,y_stride,original_y);
				alpha = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				beta = 200.0*(((double)rand())/((double)RAND_MAX) - 0.5);
				// run both no-transpose and transpose tests
				MatrixVectorProduct(test_size,test_size,alpha,A,x_stride,x,beta,y_stride,y,1);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					target = alpha*DotProduct(test_size,1,A + j*test_size,x_stride,x) + beta*original_y[j*y_stride];
					if(fabs(target - y[j*y_stride]) > tolerance)
					{
						pass = false;
						break;
					}
				}
				if(!pass)			break;
				// generate a new Y vector because the old one was overwritten
				VectorRandomPopulate(test_size*y_stride,y);
				Copy(test_size,y_stride,y,y_stride,original_y);
				MatrixVectorProduct(test_size,test_size,alpha,A,x_stride,x,beta,y_stride,y,2);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					target = alpha*DotProduct(test_size,test_size,A + j,x_stride,x) + beta*original_y[j*y_stride];
					if(fabs(target - y[j*y_stride]) > tolerance)
					{
						pass = false;
						break;
					}
				}
				if(!pass)			break;
			}
			delete [] x;
			delete [] y;
			delete [] original_y;
			delete [] A;
			return pass;
		}
		/*bool BLAS::TriangularMatrixVectorProduct()
		{
			// This function uses MatrixVectorProduct(). It is assumed that it passed 
			// the test.
			double* pdA = 0;
			double* pdX = 0;
			double* pdY = 0;
			double* pdFullX = 0;
			bool bAllPass = true;
			bool bPass = true;
			// run tests on all triangular matrix types
			// test unit increment matrix vector product in X
			unsigned int iXStride = 1;
			unsigned int iYStride = 1;
			pdA = new double[test_size*test_size];
			pdX = new double[test_size*iXStride];
			pdY = new double[test_size*iYStride];
			pdFullX = new double[test_size*iXStride];
			int iMatrixTypes[8] = {1,5,2,6,3,7,4,8};
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				// go over all matrix types
				for(unsigned int j = 0 ; j < 8 ; j++)
				{
					// XX calculation
					if(iMatrixTypes[j] == 1)
					{
						MatrixRandomPopulate(test_size,test_size,pdA);
						MakeUpperTriangular(test_size,pdA);
					}
					else if(iMatrixTypes[j] == 3)
					{
						MatrixRandomPopulate(test_size,test_size,pdA);
						MakeLowerTriangular(test_size,pdA);
					}
					if((iMatrixTypes[j]%2 == 0))				MakeUnitDiagonal(test_size,pdA);
					VectorRandomPopulate(test_size*iXStride,pdX);
					// better set this array to zero to save on scaling during full matrix vector product
					SetValue(test_size*iXStride,pdFullX,0.0);
					MatrixVectorProduct(test_size,test_size,1.0,pdA,iXStride,pdX,1.0,iXStride,pdFullX,(iMatrixTypes[j] - 1)/4 + 1);
					TriangularMatrixVectorProduct(test_size,pdA,iXStride,pdX,iMatrixTypes[j]);
					bPass = CompareVectors(test_size,iXStride,pdX,iXStride,pdFullX,1.0,0,0);
					if(!bPass)			printf("XX unit stride type %d TriangularMatrixVectorProduct() test failed\n",iMatrixTypes[j]);
					bAllPass &= bPass;

					// XY calculation
					VectorRandomPopulate(test_size*iXStride,pdX);
					VectorRandomPopulate(test_size*iYStride,pdY);
					SetValue(test_size*iXStride,pdFullX,0.0);
					MatrixVectorProduct(test_size,test_size,1.0,pdA,iXStride,pdX,1.0,iXStride,pdFullX,(iMatrixTypes[j] - 1)/4 + 1);
					TriangularMatrixVectorProduct(test_size,pdA,iXStride,pdX,iYStride,pdY,iMatrixTypes[j]);
					bPass = CompareVectors(test_size,iXStride,pdFullX,iYStride,pdY,1.0,0,0);
					if(!bPass)			printf("XY unit stride type %d TriangularMatrixVectorProduct() test failed\n",iMatrixTypes[j]);
					bAllPass &= bPass;
				}
			}
			delete [] pdX;
			delete [] pdY;
			delete [] pdFullX;
			// test strided increment matrix vector product in X and Y, generate a random number 
			// for strides but keep it under 10.
			while(iXStride <= 1)
			{
				iXStride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			while(iYStride <= 1)
			{
				iYStride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			pdX = new double[test_size*iXStride];
			pdY = new double[test_size*iYStride];
			pdFullX = new double[test_size*iXStride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				// go over all matrix types
				for(unsigned int j = 0 ; j < 8 ; j++)
				{
					// XX calculation
					if(iMatrixTypes[j] == 1)
					{
						MatrixRandomPopulate(test_size,test_size,pdA);
						MakeUpperTriangular(test_size,pdA);
					}
					else if(iMatrixTypes[j] == 3)
					{
						MatrixRandomPopulate(test_size,test_size,pdA);
						MakeLowerTriangular(test_size,pdA);
					}
					if((iMatrixTypes[j]%2 == 0))				MakeUnitDiagonal(test_size,pdA);
					VectorRandomPopulate(test_size*iXStride,pdX);
					// better set this array to zero to save on scaling during full matrix vector product
					SetValue(test_size*iXStride,pdFullX,0.0);
					MatrixVectorProduct(test_size,test_size,1.0,pdA,iXStride,pdX,1.0,iXStride,pdFullX,(iMatrixTypes[j] - 1)/4 + 1);
					TriangularMatrixVectorProduct(test_size,pdA,iXStride,pdX,iMatrixTypes[j]);
					bPass = CompareVectors(test_size,iXStride,pdX,iXStride,pdFullX,1.0,0,0);
					if(!bPass)			printf("XX strided type %d TriangularMatrixVectorProduct() test failed\n",iMatrixTypes[j]);
					bAllPass &= bPass;

					// XY calculation
					VectorRandomPopulate(test_size*iXStride,pdX);
					VectorRandomPopulate(test_size*iYStride,pdY);
					SetValue(test_size*iXStride,pdFullX,0.0);
					MatrixVectorProduct(test_size,test_size,1.0,pdA,iXStride,pdX,1.0,iXStride,pdFullX,(iMatrixTypes[j] - 1)/4 + 1);
					TriangularMatrixVectorProduct(test_size,pdA,iXStride,pdX,iYStride,pdY,iMatrixTypes[j]);
					bPass = CompareVectors(test_size,iXStride,pdFullX,iYStride,pdY,1.0,0,0);
					if(!bPass)			printf("XY strided type %d TriangularMatrixVectorProduct() test failed\n",iMatrixTypes[j]);
					bAllPass &= bPass;
				}
			}
			delete [] pdX;
			delete [] pdY;
			delete [] pdFullX;
			delete [] pdA;
			return bAllPass;
		}
		bool BLAS::SolveTriangularSystem()
		{
			// This function uses Copy() and TriangularMatrixVectorProduct(). It is
			// assumed that they passed. 
			// the test.
			double* pdA = 0;
			double* pdX = 0;
			double* pdB = 0;
			bool bAllPass = true;
			bool bPass = true;
			// run tests on all triangular matrix types
			// test unit increment matrix vector product in X
			unsigned int iStride = 1;
			pdA = new double[test_size*test_size];
			pdX = new double[test_size*iStride];
			pdB = new double[test_size*iStride];
			int iMatrixTypes[8] = {1,5,2,6,3,7,4,8};
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				for(unsigned int j = 0 ; j < 8 ; j++)
				{
					if(iMatrixTypes[j] == 1)
					{
						MatrixRandomPopulate(test_size,test_size,pdA,true);
						MakeUpperTriangular(test_size,pdA);
					}
					else if(iMatrixTypes[j] == 3)
					{
						MatrixRandomPopulate(test_size,test_size,pdA,true);
						MakeLowerTriangular(test_size,pdA);
					}
					if((iMatrixTypes[j]%2) == 0)			MakeUnitDiagonal(test_size,pdA);
					VectorRandomPopulate(test_size*iStride,pdX);
					Copy(test_size,iStride,pdX,iStride,pdB);
					SolveTriangularSystem(test_size,pdA,iStride,pdX,iMatrixTypes[j]);
					TriangularMatrixVectorProduct(test_size,pdA,iStride,pdX,iMatrixTypes[j]);
					bPass = CompareVectors(test_size,iStride,pdX,iStride,pdB,1.0,0,0);
					if(!bPass)			printf("unit stride type %d SolveTriangularSystem() test failed\n",iMatrixTypes[j]);
					bAllPass &= bPass;
				}
			}
			delete [] pdX;
			delete [] pdB;
			// test strided increment matrix vector product in X, generate a random number 
			// for strides but keep it under 10.
			while(iStride <= 1)
			{
				iStride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			pdX = new double[test_size*iStride];
			pdB = new double[test_size*iStride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				for(unsigned int j = 0 ; j < 8 ; j++)
				{
					if(iMatrixTypes[j] == 1)
					{
						MatrixRandomPopulate(test_size,test_size,pdA,true);
						MakeUpperTriangular(test_size,pdA);
					}
					else if(iMatrixTypes[j] == 3)
					{
						MatrixRandomPopulate(test_size,test_size,pdA,true);
						MakeLowerTriangular(test_size,pdA);
					}
					if((iMatrixTypes[j]%2) == 0)			MakeUnitDiagonal(test_size,pdA);
					VectorRandomPopulate(test_size*iStride,pdX);
					Copy(test_size,iStride,pdX,iStride,pdB);
					SolveTriangularSystem(test_size,pdA,iStride,pdX,iMatrixTypes[j]);
					TriangularMatrixVectorProduct(test_size,pdA,iStride,pdX,iMatrixTypes[j]);
					bPass = CompareVectors(test_size,iStride,pdX,iStride,pdB,1.0,0,0);
					if(!bPass)			printf("strided type %d SolveTriangularSystem() test failed\n",iMatrixTypes[j]);
					bAllPass &= bPass;
				}
			}
			delete [] pdX;
			delete [] pdB;
			delete [] pdA;
			return bAllPass;
		}*/
		bool BLAS::DiagonalMatrixVectorProduct()
		{
			// This function uses MatrixVectorProduct(). It is assumed that it passed 
			// the test.
			double* full_A = 0;
			double* diagonal_A = 0;
			double* x = 0;
			double* y = 0;
			double* full_x = 0;
			bool all_pass = true;
			bool pass = true;
			// run tests on all triangular matrix types
			// test unit increment matrix vector product in X
			unsigned int x_stride = 1;
			unsigned int y_stride = 1;
			full_A = new double[test_size*test_size];
			diagonal_A = new double[test_size];
			x = new double[test_size*x_stride];
			y = new double[test_size*y_stride];
			full_x = new double[test_size*x_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				MatrixRandomPopulate(test_size,test_size,full_A);
				MakeUpperTriangular(test_size,full_A);
				MakeLowerTriangular(test_size,full_A);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					diagonal_A[j] = full_A[j*test_size + j];
				}
				VectorRandomPopulate(test_size*x_stride,x);
				VectorRandomPopulate(test_size*y_stride,y);
				SetValue(test_size*x_stride,full_x,0.0);
				MatrixVectorProduct(test_size,test_size,1.0,full_A,x_stride,x,1.0,x_stride,full_x,1);
				// test matrix type 1
				DiagonalMatrixVectorProduct(test_size,full_A,x_stride,x,y_stride,y,1);
				pass = CompareVectors(test_size,x_stride,full_x,y_stride,y,1.0,0,0);
				if(!pass)			printf("unit stride type 1 DiagonalMatrixVectorProduct() test failed\n");
				all_pass &= pass;
				// test matrix type 2
				SetValue(test_size*y_stride,y,0.0);
				DiagonalMatrixVectorProduct(test_size,diagonal_A,x_stride,x,y_stride,y,2);
				pass = CompareVectors(test_size,x_stride,full_x,y_stride,y,1.0,0,0);
				if(!pass)			printf("unit stride type 2 DiagonalMatrixVectorProduct() test failed\n");
				all_pass &= pass;
			}
			delete [] x;
			delete [] y;
			delete [] full_x;
			// test strided increment matrix vector product in X and Y, generate a random number 
			// for strides but keep it under 10.
			while(x_stride <= 1)
			{
				x_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			while(y_stride <= 1)
			{
				y_stride = (unsigned int)floor(10.0*((double)rand())/((double)RAND_MAX) + 0.5);	
			}
			x = new double[test_size*x_stride];
			y = new double[test_size*y_stride];
			full_x = new double[test_size*x_stride];
			for(unsigned int i = 0 ; i < test_paranoia ; i++)
			{
				MatrixRandomPopulate(test_size,test_size,full_A);
				MakeUpperTriangular(test_size,full_A);
				MakeLowerTriangular(test_size,full_A);
				for(unsigned int j = 0 ; j < test_size ; j++)
				{
					diagonal_A[j] = full_A[j*test_size + j];
				}
				VectorRandomPopulate(test_size*x_stride,x);
				VectorRandomPopulate(test_size*y_stride,y);
				SetValue(test_size*x_stride,full_x,0.0);
				MatrixVectorProduct(test_size,test_size,1.0,full_A,x_stride,x,1.0,x_stride,full_x,1);
				// test matrix type 1
				DiagonalMatrixVectorProduct(test_size,full_A,x_stride,x,y_stride,y,1);
				pass = CompareVectors(test_size,x_stride,full_x,y_stride,y,1.0,0,0);
				if(!pass)			printf("strided type 1 DiagonalMatrixVectorProduct() test failed\n");
				all_pass &= pass;
				// test matrix type 2
				SetValue(test_size*y_stride,y,0.0);
				DiagonalMatrixVectorProduct(test_size,diagonal_A,x_stride,x,y_stride,y,2);
				pass = CompareVectors(test_size,x_stride,full_x,y_stride,y,1.0,0,0);
				if(!pass)			printf("strided type 2 DiagonalMatrixVectorProduct() test failed\n");
				all_pass &= pass;
			}
			delete [] x;
			delete [] y;
			delete [] full_x;
			delete [] full_A;
			delete [] diagonal_A;
			return all_pass;
		}
		void BLAS::VectorRandomPopulate(const unsigned int& size,double* x)
		{
			double factor = 1.0/(double)RAND_MAX;
			double range = 5.0;
			for(unsigned int i = 0 ; i < size ; i++)
			{
				x[i] = range*(((double)rand())*factor - 0.5);
			}
		}
		void BLAS::MatrixRandomPopulate(const unsigned int& row_count,const unsigned int& column_count,double* matrix,const bool& non_singular)
		{
			double factor = 1.0/(double)RAND_MAX;
			unsigned int size = row_count*column_count;
			double range = 5.0;
			for(unsigned int i = 0 ; i < size ; i++)
			{
				matrix[i] = range*(((double)rand())*factor - 0.5);
			}
			// if the matrix needs to be non-singular, modify its main diagonal elements so 
			// that each of them is larger than the first norm of its row
			if(non_singular)
			{
				double sum = 0.0;
				for(unsigned int i = 0 ; i < row_count ; i++)
				{
					// sum the absolute values of the row entries but skip the main diagonal 
					// element
					sum = AbsoluteSum(column_count,1,&matrix[i*row_count]);
					sum -= fabs(matrix[i*row_count + i]);
					if(sum > tolerance)
					{
						// Make the main diagonal element larger than the first norm of its row.
						// Here, we set the growth factor to a randomly generated real number.
						// The sign of the main diagonal element is maintained.
						range = 3.0;
						double base = 1.0;
						double growth_factor = range*(((double)rand())*factor) + base;
						if(matrix[i*row_count + i] < 0.0)			growth_factor = -growth_factor;
						matrix[i*row_count + i] = growth_factor*sum;
					}
					else
					{
						// this is a zero row, make sure that the main diagonal element is not 
						// zero
						while(fabs(matrix[i*row_count + i]) < tolerance)
						{
							matrix[i*row_count + i] = range*(((double)rand())*factor - 0.5);
						}
					}
				}
			}
		}
		void BLAS::MakeUpperTriangular(const unsigned int& size,double* matrix)
		{
			for(unsigned int i = 0 ; i < size ; i++)
			{
				for(unsigned int j = 0 ; j < i ; j++)
				{
					matrix[i*size + j] = 0.0;
				}
			}
		}
		void BLAS::MakeLowerTriangular(const unsigned int& size,double* matrix)
		{
			for(unsigned int i = 0 ; i < size ; i++)
			{
				for(unsigned int j = (i + 1) ; j < size ; j++)
				{
					matrix[i*size + j] = 0.0;
				}
			}
		}
		void BLAS::MakeUnitDiagonal(const unsigned int& size,double* matrix)
		{
			double diagonal = 0.0;
			for(unsigned int i = 0 ; i < size ; i++)
			{
				diagonal = matrix[i*size + i];
				for(unsigned int j = 0 ; j < size ; j++)
				{
					matrix[i*size + j] /= diagonal;
				}
			}
		}
		void BLAS::SetValue(const unsigned int& size,double* x,const double& value)
		{
			for(unsigned int i = 0 ; i < size ; i++)
			{
				x[i] = value;
			}
		}
		bool BLAS::CompareVectors(const unsigned int& size,const unsigned int& x_inc,double* x,const unsigned int& y_inc,double* y,const double& y_factor,const unsigned int& y_shift_inc,double* y_shift)
		{
			bool pass = true;
			if(y_shift != 0)
			{
				for(unsigned int i = 0 ; i < size ; i++)
				{
					if(fabs(x[i*x_inc] - (y[i*y_inc]*y_factor + y_shift[i*y_shift_inc])) > tolerance)
					{
						printf("vector values mismatch : %u : %u : %u : %e : %e : %e\n",size,x_inc,y_inc,x[i*x_inc],y[i*y_inc]*y_factor + y_shift[i*y_shift_inc],x[i*x_inc] - (y[i*y_inc]*y_factor + y_shift[i*y_shift_inc]));
						pass = false;
					}
				}
			}
			else
			{
				for(unsigned int i = 0 ; i < size ; i++)
				{
					if(fabs(x[i*x_inc] - (y[i*y_inc]*y_factor)) > tolerance)
					{
						printf("vector values mismatch : %u : %u : %u : %e : %e : %e\n",size,x_inc,y_inc,x[i*x_inc],y[i*y_inc]*y_factor,x[i*x_inc]- y[i*y_inc]*y_factor);
						pass = false;
					}
				}
			}
			return pass;
		}
	}
}

