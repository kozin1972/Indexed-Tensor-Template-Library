/*
 * scal.h
 *
 *  Created on: 24 мая 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef TSCAL_H_
#define TSCAL_H_

#include <tvalence.h>
#include <blas_tmpl.h>
#include <smplmath.h>

namespace iTTL
{

	template <typename T, typename STR>
	struct scal_common_loop
	{
		const STR& str;
		const T alpha;
	public:
		scal_common_loop(const STR& str): str(str), alpha(1)
		{}
		scal_common_loop(const STR& str, T alpha): str(str), alpha(alpha)
		{}
		void run_scal(T *r)
		{
			*r*=alpha;
		}
		void run_shift(T *r)
		{
			*r+=alpha;
		}
	};


	template <typename T, typename STR, typename ITTYPE, typename BASE_LOOP>
	struct scal_spread_loop: public BASE_LOOP
	{
	private:
	public:
		scal_spread_loop(const STR& str): BASE_LOOP(str)
		{}
		scal_spread_loop(const STR& str, T alpha): BASE_LOOP(str, alpha)
		{}
		void run_scal(T *r)
		{
			ITTYPE it(this->str,r);
			for (;it.not_end(r);it.move_one(r))
				BASE_LOOP::run_scal(r);
		}
		void run_shift(T *r)
		{
			ITTYPE it(this->str,r);
			for (;it.not_end(r);it.move_one(r))
				BASE_LOOP::run_shift(r);
		}
	};

	template <typename T, typename STR, typename SPREAD, typename BASE_LOOP>
	struct scal_create_general_loop
	{
		using type=BASE_LOOP;
	};

	template <typename T, typename STR, typename HEAD, typename ... SPREAD, typename BASE_LOOP>
	struct scal_create_general_loop<T, STR, type_sequence<HEAD, SPREAD...>, BASE_LOOP>:
		public scal_create_general_loop<T, STR, type_sequence<SPREAD...>, scal_spread_loop<T, STR, typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T>, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename VD>
	struct scal_loop
	{
		const STR& str;
		const T alpha;
	private:
		BLAS_INTEGER N;
		BLAS_INTEGER incr;
	public:
		scal_loop(const STR& str): str(str), alpha(1),
		N(vd_iterator_getter<STR, VD, 0>::template iterator_type<T>::length(str)),
		incr(vd_iterator_getter<STR, VD, 0>::template iterator_type<T>::step(str))
		{}
		scal_loop(const STR& str, T alpha): str(str), alpha(alpha),
		N(vd_iterator_getter<STR, VD, 0>::template iterator_type<T>::length(str)),
		incr(vd_iterator_getter<STR, VD, 0>::template iterator_type<T>::step(str))
		{}
		void run_scal(T *r)
		{
			scal<T>(&N, &alpha, r, &incr);
		}
		void run_shift(T *r)
		{
			axpy<T>(&N, &alpha, &type_constants<T>::one, &type_constants<BLAS_INTEGER>::zero, r, &incr);
		}
	};

	template <typename T, typename STR, typename VD>
	struct scal_runner;

//	DONE VECTOR
	template <typename T, typename STR, typename CONT_MASK, int S_V_TYPE, int S_C_MASK, int S_NEXT_VALENCE, size_t ... S_POS, typename ... M10, typename ... SPREAD>
	struct scal_runner<T, STR, type_sequence<CONT_MASK, type_sequence<type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, M10...>,SPREAD...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public scal_create_general_loop<T, STR, type_sequence<SPREAD...>, scal_loop<T, STR, type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, M10...> > >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename CONT_MASK, int S_V_TYPE, int S_C_MASK, int S_V_MASK, int S_NEXT_VALENCE, size_t ... S_POS, typename ... M10, typename ... SPREAD>
	struct scal_runner<T, STR, type_sequence<CONT_MASK, type_sequence<type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, S_V_MASK, S_NEXT_VALENCE, S_POS...>, M10...>,SPREAD...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public scal_create_general_loop<T, STR, type_sequence<type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, S_V_MASK, S_NEXT_VALENCE, S_POS...>, M10...>, SPREAD...>, scal_common_loop<T, STR> >
	{
	};

};

#endif /* TSCAL_H_ */
