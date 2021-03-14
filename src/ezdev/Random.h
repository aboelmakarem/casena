// Random.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 05/08/2019

#ifndef RANDOM_H_
#define RANDOM_H_

#include "vector"

namespace EZ
{
	class Random
	{
	public:
		~Random();
		static double Uniform();
		static double Uniform(const double& min,const double& max);
		static int UniformInteger(const int& min,const int& max);
		static int UniformSign();
		static double Normal();
		static double Normal(const double& mean,const double& standard_deviation);
		static double Exponential(const double& mean);
		static void Shuffle(const unsigned int& size,unsigned int* shuffled_array,const unsigned int& passes = 1);
		static void Seed(unsigned int seed);

	private:
		Random();
		Random(const Random& random);
		Random& operator=(const Random& random);
		static double MTRNG();
		static unsigned int TemperingShiftU(const unsigned int& y);
		static unsigned int TemperingShiftS(const unsigned int& y);
		static unsigned int TemperingShiftT(const unsigned int& y);
		static unsigned int TemperingShiftL(const unsigned int& y);
		#define MT_N 624
		#define MT_M 397
		static unsigned int MT_i;
		static unsigned int MT_matrix_A;
		static unsigned int MT_magnitude[2];
		static unsigned int MT_numbers[MT_N];
		static unsigned int MT_upper_mask;
		static unsigned int MT_lower_mask;
		static unsigned int MT_tempering_mask_b;
		static unsigned int MT_tempering_mask_c;
	};
}

#endif

