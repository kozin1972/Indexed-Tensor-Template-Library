/*
 * smplmath.h
 *
 *  Created on: 23 мая 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef SMPLMATH_H_
#define SMPLMATH_H_

#include <math.h>
//#include <cstdlib>

namespace tpp
{
	template <typename T>
	inline T sign(T v)
	{
		return v?(v>0.0?1.0:-1.0):0.0;
	}

	template <typename T>
	T sqrt(T v);

	template<>
	inline double sqrt<double>(double v)
	{
		return ::sqrt(v);
	}

	template <typename T>
	T log(T v);

	template<>
	inline double log<double>(double v)
	{
		return ::log(v);
	}

	template <typename T>
	T sin(T v);

	template<>
	inline double sin<double>(double v)
	{
		return ::sin(v);
	}

	template <typename T>
	T cos(T v);

	template<>
	inline double cos<double>(double v)
	{
		return ::cos(v);
	}

	template <typename T>
	T abs(T v);

	template<>
	inline double abs<double>(double v)
	{
		return ::fabs(v);
	}

	template <typename T>
	inline T rand();

	template <>
	inline double rand<double>()
	{
		return (double)(((double)::rand())/(double)RAND_MAX+::rand())/(double)(RAND_MAX);
	}

	template <typename T>
	T norm_rand();

	template<>
	inline double norm_rand<double>()
	{
		return sqrt<double>(-2.0*log<double>(1.0-rand<double>()))*cos<double>(2.0*M_PI*rand<double>());
	}


};



#endif /* SMPLMATH_H_ */
