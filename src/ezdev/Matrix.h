// Matrix.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 06/16/2019

#ifndef MATRIX_H_
#define MATRIX_H_

#include "map"
#include "stdio.h"

namespace EZ
{
	namespace Math
	{
		class Matrix;
		class BaseMatrix
		{
		public:
			BaseMatrix();
			virtual ~BaseMatrix();
			virtual void Reset();
			unsigned int RowCount() const;
			unsigned int ColumnCount() const;
			virtual void Allocate(const unsigned int& target_row_count,const unsigned int& target_column_count) = 0;
			virtual void Randomize() = 0;
			virtual void RandomizeDiagonallyDominant() = 0;
			virtual double operator()(const unsigned int& row_index,const unsigned int& column_index) const = 0;
			virtual void operator()(const unsigned int& row_index,const unsigned int& column_index,const double& value) = 0;
			virtual Matrix operator*(const double& factor) const = 0;
			virtual Matrix operator*(const Matrix& matrix) const = 0;
			virtual Matrix Solve(const Matrix& rhs) const;
			Matrix SolveJacobi(const Matrix& rhs) const;
			Matrix SolveRelaxedJacobi(const Matrix& rhs,const double& relaxation_factor) const;
			Matrix Invert() const;
			bool Valid() const;
			virtual void SumRows(double* sums) const = 0;
			
		private:
			void Initialize();

		protected:
			BaseMatrix(const BaseMatrix& matrix);
			virtual BaseMatrix& operator=(const BaseMatrix& matrix);
			unsigned int row_count;
			unsigned int column_count;
			bool valid;
		};

		class Matrix : public BaseMatrix
		{
		public:
			Matrix();
			Matrix(const Matrix& matrix);
			Matrix(const unsigned int& size);
			Matrix(const unsigned int& target_row_count,const unsigned int& target_column_count);
			~Matrix();
			Matrix& operator=(const Matrix& matrix);
			void Reset();
			void Allocate(const unsigned int& target_row_count,const unsigned int& target_column_count);
			void Randomize();
			void RandomizeDiagonallyDominant();
			void SetRow(const unsigned int& row_index,double* row_entries);
			double operator()(const unsigned int& row_index) const;
			double operator()(const unsigned int& row_index,const unsigned int& column_index) const;
			void operator()(const unsigned int& row_index,const unsigned int& column_index,const double& value);
			void Increment(const unsigned int& row_index,const unsigned int& column_index,const double& value);
			void Decrement(const unsigned int& row_index,const unsigned int& column_index,const double& value);
			Matrix operator+(const Matrix& matrix) const;
			void operator+=(const Matrix& matrix) const;
			Matrix operator-(const Matrix& matrix) const;
			void operator-=(const Matrix& matrix) const;
			Matrix operator*(const Matrix& matrix) const;
			Matrix operator*(const double& factor) const;
			double operator^(const Matrix& matrix) const;
			void Blend(const Matrix& matrix,const double& factor);
			Matrix DiagonalPreMultiply(const Matrix& matrix_diagonal) const;
			Matrix DiagonalPostMultiply(const Matrix& matrix_diagonal) const;
			double SquaredNorm() const;
			Matrix Solve(const Matrix& rhs) const;
			Matrix SolveGE(const Matrix& rhs) const;
			Matrix Transpose() const;
			void Print() const;
			void SumRows(double* sums) const;
			bool WriteRow(const unsigned int& row,FILE* file) const;
			bool ReadRow(const unsigned int& row,FILE* file) const;
			void ZeroEntries();
			const double* Row(const unsigned int& row_index) const;
			const double* Entries() const;

		private:
			void Initialize();
			double* entries;
		};

		class SparseMatrix : public BaseMatrix
		{
		public:
			SparseMatrix();
			SparseMatrix(const SparseMatrix& matrix);
			SparseMatrix(const unsigned int& size);
			SparseMatrix(const unsigned int& target_row_count,const unsigned int& target_column_count);
			~SparseMatrix();
			SparseMatrix& operator=(const SparseMatrix& matrix);
			void Reset();
			void Allocate(const unsigned int& target_row_count,const unsigned int& target_column_count);
			void Randomize();
			void RandomizeDiagonallyDominant();
			double operator()(const unsigned int& row_index,const unsigned int& column_index) const;
			void operator()(const unsigned int& row_index,const unsigned int& column_index,const double& value);
			Matrix operator*(const Matrix& matrix) const;
			void SumRows(double* sums) const;
			
		private:
			void Initialize();
			std::map<unsigned int,double> entries;
		};
	}
}

#endif

