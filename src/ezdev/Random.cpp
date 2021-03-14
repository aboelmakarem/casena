// Random.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 05/08/2019

#include "Random.h"
#include "math.h"
#include "time.h"

namespace EZ
{
	unsigned int Random::MT_i = MT_N + 1;
	unsigned int Random::MT_matrix_A = 2567483615;
	unsigned int Random::MT_magnitude[2] = {0,Random::MT_matrix_A};
	unsigned int Random::MT_numbers[MT_N];
	unsigned int Random::MT_upper_mask = 2147483648;
	unsigned int Random::MT_lower_mask = 2147483647;
	unsigned int Random::MT_tempering_mask_b = 2636928640;
	unsigned int Random::MT_tempering_mask_c = 4022730752;
	Random::Random(){}
	Random::Random(const Random& random){*this = random;}
	Random::~Random(){}
	Random& Random::operator=(const Random& random){return *this;}
	double Random::Uniform(){return MTRNG();}
	double Random::Uniform(const double& min,const double& max)
	{
		if(max < min)		return 0.0;
		return (min + Uniform()*(max - min));
	}
	int Random::UniformInteger(const int& min,const int& max)
	{
		bool accepted = false;
		int result = 0;
		while(!accepted)
		{
			result = (int)floor(Uniform(min - 1,max + 1) + 0.5);
			if(result != (min - 1) && result != (max + 1))		accepted = true;
		}
		return result;
	}
	int Random::UniformSign()
	{
		if(UniformInteger(0,1) == 0)		return -1;
		return 1;
	}
	double Random::Normal()
	{
		// Use Box-Muller transform
		double x1 = 0.0;
		double x2 = 0.0;
		double w = 0.0;
		double y1 = 0.0;
		//double y2 = 0.0;
		do
		{
			x1 = 2.0*Uniform() - 1.0;
			x2 = 2.0*Uniform() - 1.0;
			w = x1*x1 + x2*x2;
		}
		while(w >= 1.0);
		w = sqrt((-2.0*log(w))/w);
		y1 = x1*w;
		//y2 = x2*w;
		return y1;
	}
	double Random::Normal(const double& mean,const double& standard_deviation){return (mean + standard_deviation*Normal());}
	double Random::Exponential(const double& mean){return (-mean*log(1.0 - Uniform()));}
	void Random::Shuffle(const unsigned int& size,unsigned int* shuffled_array,const unsigned int& passes)
	{
		// This function shuffles an array of unsigned integers that runs from 1 to size
		// The function allocates the array, the function user is responsible for deallocating 
		// it. 
		if(shuffled_array != 0)		delete [] shuffled_array;
		shuffled_array = new unsigned int[size];
		for(unsigned int i = 0 ; i < size ; i++)
		{
			shuffled_array[i] = i + 1;
		}
		unsigned int index = 0;
		unsigned int temp = 0;
		for(unsigned int i = 0 ; i < passes ; i++)
		{
			for(unsigned int j = 0 ; j < size ; j++)
			{
				index = UniformInteger(1,size) - 1;
				temp = shuffled_array[j];
				shuffled_array[j] = shuffled_array[index];
				shuffled_array[index] = temp;
			}
		}
	}
	void Random::Seed(unsigned int seed)
	{
		MT_numbers[0] = seed;
		for(MT_i = 1 ; MT_i < MT_N ; MT_i++)
		{
			MT_numbers[MT_i] = (69069*MT_numbers[MT_i - 1]);
		}
	}
	double Random::MTRNG()
	{
		// Use the Mersenne-Twister algorithm to generate random numbers
		unsigned int y = 0;
		if(MT_i >= MT_N)
		{
			unsigned int k = 0;
			if(MT_i == MT_N + 1)
			{
				Seed((unsigned int)time(0));
			}
			for(k = 0; k < MT_N - MT_M ; k++)
			{
				y = (MT_numbers[k] & MT_upper_mask)|(MT_numbers[k + 1] & MT_lower_mask);
				MT_numbers[k] = MT_numbers[k + MT_M]^(y >> 1)^MT_magnitude[y & 0x1];
			}
			for(; k < MT_N - 1; k++)
			{
				y = (MT_numbers[k] & MT_upper_mask)|(MT_numbers[k + 1] & MT_lower_mask);
				MT_numbers[k] = MT_numbers[k + (MT_M - MT_N)]^(y >> 1)^MT_magnitude[y & 0x1];
			}
			y = (MT_numbers[MT_N - 1] & MT_upper_mask)|(MT_numbers[0] & MT_lower_mask);
			MT_numbers[MT_N - 1] = MT_numbers[MT_M - 1]^(y >> 1)^MT_magnitude[y & 0x1];
			MT_i = 0;
		}

		y = MT_numbers[MT_i++];
		y = y^TemperingShiftU(y);
		y = y^(TemperingShiftS(y) & MT_tempering_mask_b);
		y = y^(TemperingShiftT(y) & MT_tempering_mask_c);
		y = y^TemperingShiftL(y);
		return ((double)y)/((unsigned int)0xffffffff);
	}
	unsigned int Random::TemperingShiftU(const unsigned int& y){return (y >> 11);}
	unsigned int Random::TemperingShiftS(const unsigned int& y){return (y << 7);}
	unsigned int Random::TemperingShiftT(const unsigned int& y){return (y << 15);}
	unsigned int Random::TemperingShiftL(const unsigned int& y){return (y >> 18);}
}

