/*
 * tdiv.h
 *
 *  Created on: 1 июн. 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef TDIV_H_
#define TDIV_H_

#include <tvalence.h>
#include <blas_tmpl.h>
#include <smplmath.h>

namespace iTTL
{

	template <typename T, typename STR, typename STA, typename VD, typename BASE_LOOP>
	struct div_suma_loop
	{
	};

	template <typename T, typename STR, typename STA, typename VD, typename BASE_LOOP>
	struct div_ra_loop: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STR, VD, 0>::template iterator_type<T> ritype;
		typedef typename vd_iterator_getter<STA, VD, 1>::template iterator_type<T> aitype;
		div_ra_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		div_ra_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void run_sign(T *r, const T *a)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a))
				BASE_LOOP::run_sign(r,a);
		}
		void run_div(T *r, const T *a)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a))
				BASE_LOOP::run_div(r,a);
		}
		void run_scal(T *r, const T *a)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a))
				BASE_LOOP::run_scal(r,a);
		}
		void spread(const T *r0, T *ri)
		{
			ritype rit(this->str, ri);
			for (;rit.not_end(ri);rit.move_one(ri))
				BASE_LOOP::spread(r0,ri);
		}
	};

	template <typename T, typename STR, typename STA>
	struct div_common_loop
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	public:
		div_common_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1)
		{}
		div_common_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha)
		{}
		void run_sign(T *r, const T *a)
		{
			*r=sign<T>(*a);
		}
		void run_div(T *r, const T *a)
		{
			*r/=*a;
		}
		void run_scal(T *r, const T *a)
		{
			*r*=*a;
		}
		void spread(const T *r0, T *ri)
		{
			*ri=*r0;
		}
	};

	template <typename T, typename STR, typename STA, typename ITTYPE, bool IS_FIRST, typename BASE_LOOP>
	struct div_spread_loop;

// LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename ITTYPE, typename BASE_LOOP>
	struct div_spread_loop<T, STR, STA, ITTYPE, false, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		div_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		div_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void spread_first(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::spread(r0,ri);
		}
		void spread(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::spread(r0,ri);
		}
	};

// LAST; FIRST
	template <typename T, typename STR, typename STA, typename ITTYPE, typename BASE_LOOP>
	struct div_spread_loop<T, STR, STA, ITTYPE, true, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		div_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		div_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void spread_first(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::spread(r0,ri);
		}
		void spread(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::spread(r0,ri);
		}
		void run_sign(T *r, const T *a)
		{
			T *r0=r;
			T *ri=r;
			BASE_LOOP::run_sign(r,a);
			ri=r0;
			spread_first(r0, ri);
		}
		void run_div(T *r, const T *a)
		{
			ITTYPE it(this->str,r);
			for (;it.not_end(r);it.move_one(r))
				BASE_LOOP::run_div(r,a);
		}
		void run_scal(T *r, const T *a)
		{
			ITTYPE it(this->str,r);
			for (;it.not_end(r);it.move_one(r))
				BASE_LOOP::run_scal(r,a);
		}
	};




// NOT LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename ITTYPE, typename HEAD, typename BASE_BASE_LOOP>
	struct div_spread_loop<T, STR, STA, ITTYPE, false, div_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> >: public div_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef div_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
	public:
		div_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		div_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void spread_first(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			BASE_LOOP::spread_first(r0,ri);
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::spread(r0,ri);
		}
		void spread(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::spread(r0,ri);
		}
	};

// NOT LAST; FIRST
	template <typename T, typename STR, typename STA, typename ITTYPE, typename HEAD, typename BASE_BASE_LOOP>
	struct div_spread_loop<T, STR, STA, ITTYPE, true, div_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> >: public div_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef div_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
	public:
		div_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		div_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void spread_first(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			BASE_LOOP::spread_first(r0,ri);
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::spread(r0,ri);
		}
		void spread(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::spread(r0,ri);
		}
		void run_sign(T *r, const T *a)
		{
			T *r0=r;
			T *ri=r;
			BASE_LOOP::run_sign(r,a);
			ri=r0;
			spread_first(r0, ri);
		}
		void run_div(T *r, const T *a)
		{
			ITTYPE it(this->str,r);
			for (;it.not_end(r);it.move_one(r))
				BASE_LOOP::run_div(r,a);
		}
		void run_scal(T *r, const T *a)
		{
			ITTYPE it(this->str,r);
			for (;it.not_end(r);it.move_one(r))
				BASE_LOOP::run_scal(r,a);
		}
	};

	template <typename T, typename STR, typename STA, typename SPREAD, typename SUMA, typename RA, typename BASE_LOOP>
	struct div_create_general_loop
	{
		using type=BASE_LOOP;
	};

	template <typename T, typename STR, typename STA, typename SPREAD, typename HEAD, typename ... SUMA, typename RA, typename BASE_LOOP>
	struct div_create_general_loop<T, STR, STA, SPREAD, type_sequence<HEAD, SUMA...>, RA, BASE_LOOP>:
		public div_create_general_loop<T, STR, STA, SPREAD, type_sequence<SUMA...>, RA, div_suma_loop<T, STR, STA, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SPREAD, typename HEAD, typename ... RA, typename BASE_LOOP>
	struct div_create_general_loop<T, STR, STA, SPREAD, type_sequence<>, type_sequence<HEAD, RA...>, BASE_LOOP>:
		public div_create_general_loop<T, STR, STA, SPREAD, type_sequence<>, type_sequence<RA...>, div_ra_loop<T, STR, STA, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename HEAD, typename ... SPREAD, typename BASE_LOOP>
	struct div_create_general_loop<T, STR, STA, type_sequence<HEAD, SPREAD...>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public div_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<>, type_sequence<>, div_spread_loop<T, STR, STA, typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T>, false, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename HEAD, typename BASE_LOOP>
	struct div_create_general_loop<T, STR, STA, type_sequence<HEAD>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public div_create_general_loop<T, STR, STA, type_sequence<>, type_sequence<>, type_sequence<>, div_spread_loop<T, STR, STA, typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T>, true, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename VD>
	struct div_just_spread_loop
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	private:
		BLAS_INTEGER l;
		BLAS_INTEGER s;
	public:
		using ittype=typename vd_iterator_getter<STR, VD, 0>::template iterator_type<T>;
	public:
		div_just_spread_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1),
			l(ittype::length(str)), s(ittype::step(str))
		{}
		div_just_spread_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha),
			l(ittype::length(str)), s(ittype::step(str))
		{}
		void run_sign(T *r, const T *a)
		{
			T v=sign<T>(*a);
			copy<T>(&l, &v, &type_constants<BLAS_INTEGER>::zero, r, &s);
		}
		void run_div(T *r, const T *a)
		{
			T alpha=1.0/(*a);
			scal<T>(&l,&alpha,r,&s);
		}
		void run_scal(T *r, const T *a)
		{
			scal<T>(&l,a,r,&s);
		}
		void spread(const T *r0, T *ri)
		{
			copy<T>(&l, r0, &s, ri, &s);
		}
	};


	template <typename T, typename STR, typename STA, typename VD>
	struct div_runner;

//	DONE just spread
	template <typename T, typename STR, typename STA, typename CONT_MASK, int S_V_TYPE, int S_C_MASK, int S_NEXT_VALENCE, size_t ... S_POS, typename ... M10, typename ... M1>
	struct div_runner<T, STR, STA, type_sequence<CONT_MASK, type_sequence<type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, M10...>, M1...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public div_create_general_loop<T, STR, STA, type_sequence<M1...>, type_sequence<>, type_sequence<>, div_just_spread_loop<T, STR, STA, type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, M10...> > >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, typename ... RA>
	struct div_runner<T, STR, STA, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public div_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, div_common_loop<T, STR, STA> >
	{
	};

};




#endif /* TDIV_H_ */
