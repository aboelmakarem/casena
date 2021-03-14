// Matrix.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 06/16/2019

#include "Matrix.h"
#include "string.h"
#include "BLAS.h"
#include "math.h"
#include "Random.h"
#include "stdio.h"

namespace EZ
{
	namespace Math
	{
		BaseMatrix::BaseMatrix(){Initialize();}
		BaseMatrix::BaseMatrix(const BaseMatrix& matrix)
		{
			Initialize();
			*this = matrix;
		}
		BaseMatrix::~BaseMatrix(){Reset();}
		BaseMatrix& BaseMatrix::operator=(const BaseMatrix& matrix)
		{
			valid = matrix.valid;
			return *this;
		}
		void BaseMatrix::Reset(){Initialize();}
		unsigned int BaseMatrix::RowCount() const{return row_count;}
		unsigned int BaseMatrix::ColumnCount() const{return column_count;}
		Matrix BaseMatrix::Solve(const Matrix& rhs) const{return SolveJacobi(rhs);}
		Matrix BaseMatrix::SolveJacobi(const Matrix& rhs) const{return SolveRelaxedJacobi(rhs,1.0);}
		Matrix BaseMatrix::SolveRelaxedJacobi(const Matrix& rhs,const double& relaxation_factor) const
		{
			// Solve the system of equations the matrix of which is this matrix and the right 
			// hand side of which is rhs using the relaxed Jacobi iteration method. The error 
			// is the root mean square error and the relaxation factor is given. The matrix 
			// needs to be diagonally dominant for the iteration to converge but this check 
			// is not made here. 
			double tolerance = 1.0e-6;
			double error = 100.0*tolerance;
			unsigned int rhs_count = rhs.column_count;
			Matrix x(row_count,rhs_count);
			Matrix x_new(row_count,rhs_count);
			Matrix y(row_count,rhs_count);
			// Copy the right hand side matrix to y.
			y = rhs;
			// Create the diagonal matrix and its inverse, store the main diagonal only.
			Matrix diagonal(row_count,1);
			Matrix inverse_diagonal(row_count,1);
			double temp = 0.0;
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				temp = (*this)(i,i);
				diagonal(i,0,temp);
				inverse_diagonal(i,0,1.0/temp);
			}
			while(error > tolerance)
			{
				// multiply the diagonal matrix by the current solution
				y = x.DiagonalPreMultiply(diagonal);
				// multiply the current solution x by the matrix A and subtract the product 
				// from the right hand side while taking out the product of the diagonal matrix 
				// and the current solution. 
				y = rhs - (*this)*x + y;
				// multiply the inverse of the diagonal matrix by y
				x_new = y.DiagonalPreMultiply(inverse_diagonal);
				// relax the solution if needed
				if(fabs(relaxation_factor - 1.0) > 1.0e-3)		x_new.Blend(x,relaxation_factor);
				// compute error
				error = sqrt((x_new - x).SquaredNorm());
				// update solution
				x = x_new;
			}
			return x;
		}
		Matrix BaseMatrix::Invert() const
		{
			Matrix rhs(row_count,column_count);
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				rhs(i,i,1.0);
			}
			Matrix inverse = Solve(rhs);
			return inverse;
		}
		bool BaseMatrix::Valid() const{return valid;}
		void BaseMatrix::Initialize()
		{
			row_count = 0;
			column_count = 0;
			valid = true;
		}

		Matrix::Matrix(){Initialize();}
		Matrix::Matrix(const Matrix& matrix) : BaseMatrix(matrix)
		{
			Initialize();
			*this = matrix;
		}
		Matrix::Matrix(const unsigned int& size)
		{
			Initialize();
			Allocate(size,size);
		}
		Matrix::Matrix(const unsigned int& target_row_count,const unsigned int& target_column_count)
		{
			Initialize();
			Allocate(target_row_count,target_column_count);
		}
		Matrix::~Matrix(){Reset();}
		Matrix& Matrix::operator=(const Matrix& matrix)
		{
			BaseMatrix::operator=(matrix);
			Allocate(matrix.row_count,matrix.column_count);
			BLAS::Copy(row_count*column_count,1,matrix.entries,1,entries);
			return *this;
		}
		void Matrix::Reset()
		{
			if(entries != 0)			delete [] entries;
			BaseMatrix::Reset();
			Initialize();
		}
		void Matrix::Allocate(const unsigned int& target_row_count,const unsigned int& target_column_count)
		{
			if(row_count == target_row_count)
			{
				if(column_count == target_column_count)		return;
			}
			Reset();
			row_count = target_row_count;
			column_count = target_column_count;
			unsigned int size = row_count*column_count;
			entries = new double[size];
			memset(entries,0,size*sizeof(double));
		}
		void Matrix::Randomize()
		{
			// This function populates the matrix with random values
			unsigned int size = row_count*column_count;
			for(unsigned int i = 0 ; i < size ; i++)
			{
				entries[i] = Random::Uniform(-10.0,10.0);
			}
		}
		void Matrix::RandomizeDiagonallyDominant()
		{
			// This function populates the matrix with random values while guaranteeing 
			// that the matrix is diagonally dominant.
			Randomize();
			double sum = 0.0;
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				// sum all absolute values of row entries
				sum = BLAS::AbsoluteSum(column_count,1,entries + i*column_count);
				// subtract the main diagonal entry
				sum -= entries[i*column_count + i];
				// update the main diagonal term by scaling the sum with a factor that 
				// is greater than one
				entries[i*column_count + i] = Random::Uniform(1.0,10.0)*sum;
			}
		}
		void Matrix::SetRow(const unsigned int& row_index,double* row_entries)
		{
			// Set the row of index row_index with the entries in row_entries
			// The array row_entries is assumed to have been properly allocated 
			// and populated. No further checks are made here. 
			BLAS::Copy(column_count,1,row_entries,1,&entries[row_index*column_count]);
		}
		double Matrix::operator()(const unsigned int& row_index) const{return entries[row_index*column_count];}
		double Matrix::operator()(const unsigned int& row_index,const unsigned int& column_index) const{return entries[row_index*column_count + column_index];}
		void Matrix::operator()(const unsigned int& row_index,const unsigned int& column_index,const double& value){entries[row_index*column_count + column_index] = value;}
		void Matrix::Increment(const unsigned int& row_index,const unsigned int& column_index,const double& value){entries[row_index*column_count + column_index] += value;}
		void Matrix::Decrement(const unsigned int& row_index,const unsigned int& column_index,const double& value){entries[row_index*column_count + column_index] -= value;}
		Matrix Matrix::operator+(const Matrix& matrix) const
		{
			// The function computes the sum of this matrix and the input matrix
			// The matrix dimensions must be consistent, no checks are made to ensure that. 
			Matrix result(row_count,column_count);
			BLAS::Copy(row_count*column_count,1,entries,1,result.entries);
			BLAS::ScaleAndAdd(row_count*column_count,1,matrix.entries,1.0,1,result.entries);
			return result;
		}
		void Matrix::operator+=(const Matrix& matrix) const
		{
			// The function computes the sum of this matrix and the input matrix
			// The matrix dimensions must be consistent, no checks are made to ensure that. 
			BLAS::ScaleAndAdd(row_count*column_count,1,matrix.entries,1.0,1,entries);
		}
		Matrix Matrix::operator-(const Matrix& matrix) const
		{
			// The function computes the difference between this matrix and the input matrix
			// The matrix dimensions must be consistent, no checks are made to ensure that. 
			Matrix result(row_count,column_count);
			BLAS::Copy(row_count*column_count,1,entries,1,result.entries);
			BLAS::ScaleAndAdd(row_count*column_count,1,matrix.entries,-1.0,1,result.entries);
			return result;
		}
		void Matrix::operator-=(const Matrix& matrix) const
		{
			// The function computes the difference between this matrix and the input matrix
			// The matrix dimensions must be consistent, no checks are made to ensure that. 
			BLAS::ScaleAndAdd(row_count*column_count,1,matrix.entries,-1.0,1,entries);
		}
		Matrix Matrix::operator*(const Matrix& matrix) const
		{
			// The function computes the product of this matrix with the input matrix
			// The matrix dimensions must be consistent, no checks are made to ensure that. 
			Matrix product(row_count,matrix.column_count);
			for(unsigned int i = 0 ; i < matrix.column_count ; i++)
			{
				BLAS::MatrixVectorProduct(row_count,column_count,1.0,entries,matrix.column_count,matrix.entries + i,0.0,matrix.column_count,product.entries + i,1);
			}
			return product;
		}
		Matrix Matrix::operator*(const double& factor) const
		{
			Matrix product(*this);
			BLAS::Scale(row_count*column_count,1,product.entries,factor);
			return product;
		}
		double Matrix::operator^(const Matrix& matrix) const
		{
			// this is a contraction operator, both matrices should have 
			// the same number of rows and columns
			if(matrix.row_count != row_count)			return 0.0;
			if(matrix.column_count != column_count)		return 0.0;
			double contraction = 0.0;
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				for(unsigned int j = 0 ; j < column_count ; j++)
				{
					contraction += (entries[i*column_count + j]*matrix.entries[i*column_count + j]);
				}
			}
			return contraction;
		}
		void Matrix::Blend(const Matrix& matrix,const double& factor)
		{
			// This function forms the weighted sum of this matrix weighted by 
			// the factor and the input matrix weighted by (1 - factor)
			BLAS::Scale(row_count*column_count,1,entries,factor);
			BLAS::ScaleAndAdd(row_count*column_count,1,matrix.entries,1.0 - factor,1,entries);
		}
		Matrix Matrix::DiagonalPreMultiply(const Matrix& matrix_diagonal) const
		{
			Matrix product(row_count,column_count);
			for(unsigned int i = 0 ; i < column_count ; i++)
			{
				BLAS::DiagonalMatrixVectorProduct(row_count,matrix_diagonal.entries,column_count,entries + i,column_count,product.entries + i,2);
			}
			return product;
		}
		Matrix Matrix::DiagonalPostMultiply(const Matrix& matrix_diagonal) const
		{
			Matrix product(row_count,column_count);
			for(unsigned int i = 0 ; i < column_count ; i++)
			{
				BLAS::ScaleAndAdd(row_count,column_count,&entries[i],matrix_diagonal(i,0),column_count,&product.entries[i]);
			}
			return product;
		}
		double Matrix::SquaredNorm() const
		{
			return BLAS::SquaredNorm(row_count*column_count,1,entries);
		}
		Matrix Matrix::Solve(const Matrix& rhs) const{return SolveGE(rhs);}
		Matrix Matrix::SolveGE(const Matrix& rhs) const
		{
			// Solve using Gauss elimination with partial
			// pivoting. very robust but very slow
			// Matrix dimensions are assumed to be consistent, no 
			// checks are made to ensure that. 
			unsigned int rhs_count = rhs.column_count;
			double* scales = new double[row_count];
			// get the maximum for each row and place it in the scales vector
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				scales[i] = (*this)(i,BLAS::MaxAbsoluteIndex(column_count,1,entries + i*column_count));
			}
			// set the working matrix with the row wise scaled version
			// of this matrix
			Matrix working_matrix(row_count,column_count);
			Matrix solution(row_count,rhs_count);
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				BLAS::ScaleAndAdd(column_count,1,entries + i*column_count,1.0/scales[i],1,working_matrix.entries + i*column_count);
				BLAS::ScaleAndAdd(rhs_count,1,rhs.entries + i*rhs_count,1.0/scales[i],1,solution.entries + i*rhs_count);
			}
			delete [] scales;
			unsigned int swap_index = 0;
			double temp = 0.0;
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				temp = fabs(working_matrix(i,i));
				swap_index = i;
				for(unsigned int j = i + 1 ; j < row_count ; j++)
				{
					if(fabs(working_matrix(j,i)) > temp)
					{
						temp = fabs(working_matrix(j,i));
						swap_index = j;
					}
				}
				if(temp < 1.0E-50)
				{
					// mark the solution as invalid
					solution.valid = false;
					return solution;
				}
				// swap both rows, physically, in memory
				BLAS::Swap(column_count,1,working_matrix.entries + i*column_count,1,working_matrix.entries + swap_index*column_count);
				BLAS::Swap(rhs_count,1,solution.entries + i*rhs_count,1,solution.entries + swap_index*rhs_count);
				for(unsigned int j = i + 1 ; j < row_count; j++)
				{
					temp = -working_matrix(j,i)/working_matrix(i,i);
					BLAS::ScaleAndAdd(column_count - i,1,working_matrix.entries + i*column_count + i,temp,1,working_matrix.entries + j*column_count + i);
					BLAS::ScaleAndAdd(rhs_count,1,solution.entries + i*rhs_count,temp,1,solution.entries + j*rhs_count);
				}
			}
			// back substitution to get the final solution
			for(unsigned int k = 0 ; k < rhs_count ; k++)
			{	
				solution(row_count - 1,k,solution(row_count - 1,k)/working_matrix(row_count - 1,row_count - 1));
			}
			for(unsigned int k = 0 ; k < rhs_count ; k++)
			{
				for(int i = row_count - 2 ; i >= 0 ; i--)
				{
					temp = BLAS::DotProduct(column_count - i - 1,1,working_matrix.entries + i*column_count + i + 1,rhs_count,solution.entries + (i + 1)*rhs_count + k);
					solution(i,k,(solution(i,k) - temp)/working_matrix(i,i));
				}
			}
			return solution;
		}
		Matrix Matrix::Transpose() const
		{
			Matrix transpose(column_count,row_count);
			for(unsigned int i = 0 ; i < column_count ; i++)
			{
				for(unsigned int j = 0 ; j < row_count ; j++)
				{
					transpose(i,j,entries[j*column_count + i]);
				}
			}
			return transpose;
		}
		void Matrix::Print() const
		{
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				for(unsigned int j = 0 ; j < column_count ; j++)
				{
					printf("%e\t",entries[i*column_count + j]);
				}
				printf("\n");
			}
		}
		void Matrix::SumRows(double* sums) const
		{
			double sum = 0.0;
			unsigned int index = 0;
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				sum = 0.0;
				for(unsigned int j = 0 ; j < column_count ; j++)
				{
					sum += entries[index++];
				}
				sums[i] = sum;
			}
		}
		bool Matrix::WriteRow(const unsigned int& row,FILE* file) const
		{
			return (column_count == fwrite(&entries[row*column_count],sizeof(double),column_count,file));
		}
		bool Matrix::ReadRow(const unsigned int& row,FILE* file) const
		{
			return (column_count == fread(&entries[row*column_count],sizeof(double),column_count,file));
		}
		void Matrix::ZeroEntries(){memset(entries,0,row_count*column_count*sizeof(double));}
		const double* Matrix::Row(const unsigned int& row_index) const{return &entries[row_index*column_count];}
		const double* Matrix::Entries() const{return entries;}
		void Matrix::Initialize(){entries = 0;}

		SparseMatrix::SparseMatrix(){Initialize();}
		SparseMatrix::SparseMatrix(const SparseMatrix& matrix) : BaseMatrix(matrix)
		{
			*this = matrix;
		}
		SparseMatrix::SparseMatrix(const unsigned int& size)
		{
			Initialize();
			Allocate(size,size);
		}
		SparseMatrix::SparseMatrix(const unsigned int& target_row_count,const unsigned int& target_column_count)
		{
			Initialize();
			Allocate(target_row_count,target_column_count);
		}
		SparseMatrix::~SparseMatrix(){Reset();}
		SparseMatrix& SparseMatrix::operator=(const SparseMatrix& matrix)
		{
			Allocate(matrix.row_count,matrix.column_count);
			for(std::map<unsigned int,double>::const_iterator entry = matrix.entries.begin() ; entry != matrix.entries.end() ; entry++)
			{
				entries[entry->first] = entry->second;
			}
			return *this;
		}
		void SparseMatrix::Reset()
		{
			entries.clear();
			BaseMatrix::Reset();
			Initialize();
		}
		void SparseMatrix::Allocate(const unsigned int& target_row_count,const unsigned int& target_column_count)
		{
			Reset();
			row_count = target_row_count;
			column_count = target_column_count;
			entries.clear();
		}
		void SparseMatrix::Randomize()
		{
			entries.clear();
			double density = 0.01;
			unsigned int size = row_count*column_count;
			unsigned int entry_count = (unsigned int)floor(density*size + 0.5);
			for(unsigned int i = 0 ; i < entry_count ; i++)
			{
				entries[Random::UniformInteger(0,size)] = Random::Uniform(-10.0,10.0);
			}
		}
		void SparseMatrix::RandomizeDiagonallyDominant()
		{
			Randomize();
			double* row_sums = new double[row_count];
			double* diagonals = new double[row_count];
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				row_sums[i] = 0.0;
				diagonals[i] = 0.0;
			}
			// sum the entries per row
			unsigned int row = 0;
			for(std::map<unsigned int,double>::const_iterator entry = entries.begin() ; entry != entries.end() ; entry++)
			{
				row = entry->first/column_count;
				row_sums[row] += fabs(entry->second);
				if(row == (entry->first%column_count))		diagonals[row] = entry->second;
			}
			for(unsigned int i = 0 ; i < row_count ; i++)
			{
				// subtract the main diagonal entry from the row sum
				row_sums[i] -= diagonals[i];
				// update the main diagonal term by scaling the sum with a factor that 
				// is greater than one
				if(row_sums[i] < 1.0e-3)			row_sums[i] = 1.0;
				entries[i*column_count + i] = Random::Uniform(1.0,10.0)*row_sums[i];
			}
		}
		double SparseMatrix::operator()(const unsigned int& row_index,const unsigned int& column_index) const
		{
			unsigned int index = row_index*column_count + column_index;
			std::map<unsigned int,double>::const_iterator entry = entries.find(index);
			if(entry == entries.end())		return 0.0;
			return entry->second;
		}
		void SparseMatrix::operator()(const unsigned int& row_index,const unsigned int& column_index,const double& value)
		{
			unsigned int index = row_index*column_count + column_index;
			entries[index] = value;
		}
		Matrix SparseMatrix::operator*(const Matrix& matrix) const
		{
			unsigned int row = 0;
			unsigned int column = 0;
			unsigned int matrix_columns = matrix.ColumnCount();
			Matrix product(row_count,matrix_columns);
			for(std::map<unsigned int,double>::const_iterator entry = entries.begin() ; entry != entries.end() ; entry++)
			{
				row = entry->first/column_count;
				column = entry->first%column_count;
				for(unsigned int j = 0 ; j < matrix_columns ; j++)
				{
					product.Increment(row,j,entry->second*matrix(column,j));
				}
			}
			return product;
		}
		void SparseMatrix::SumRows(double* sums) const
		{
			for(std::map<unsigned int,double>::const_iterator entry = entries.begin() ; entry != entries.end() ; entry++)
			{
				sums[entry->first/column_count] += entry->second;
			}
		}
		void SparseMatrix::Initialize(){entries.clear();}
	}
}

