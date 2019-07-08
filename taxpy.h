/*
 * taxpy.h
 *
 *  Created on: 1 июн. 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef TAXPY_H_
#define TAXPY_H_

#include <tvalence.h>
#include <blas_tmpl.h>
#include <smplmath.h>

namespace tpp
{

	template <typename T, typename STR, typename STA, typename HEAD, typename BASE_LOOP>
	struct axpy_suma_loop: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STA, HEAD, 1>::template iterator_type<T> itype;
		axpy_suma_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		axpy_suma_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void run_axpy(T *r, const T *a)
		{
			itype it(this->sta, a);
			for (;it.not_end(a);it.move_one(a))
				BASE_LOOP::run_axpy(r,a);
		}
	};

	template <typename T, typename STR, typename STA, typename HEAD, typename BASE_LOOP>
	struct axpy_ra_loop: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T> ritype;
		typedef typename vd_iterator_getter<STA, HEAD, 1>::template iterator_type<T> aitype;
		axpy_ra_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		axpy_ra_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void run_axpy(T *r, const T *a)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a))
				BASE_LOOP::run_axpy(r,a);
		}
		void add0(const T *r0, T *ri)
		{
			ritype rit0(this->str, r0);
			ritype riti(this->str, ri);
			for (;rit0.not_end(r0);rit0.move_one(r0),riti.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void sub0(const T *r0, T *ri)
		{
			ritype rit0(this->str, r0);
			ritype riti(this->str, ri);
			for (;rit0.not_end(r0);rit0.move_one(r0),riti.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
	};

	template <typename T, typename STR, typename STA>
	struct axpy_common_loop
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	public:
		axpy_common_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1)
		{}
		axpy_common_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha)
		{}
		void run_axpy(T *r, const T *a)
		{
			*r+=*a*alpha;
		}
		void add0(const T *r0, T *ri)
		{
			*ri+=*r0;
		}
		void sub0(const T *r0, T *ri)
		{
			*ri-=*r0;
		}
	};

	template <typename T, typename STR, typename STA, typename ITTYPE, bool IS_FIRST, typename BASE_LOOP>
	struct axpy_spread_loop;

// LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename ITTYPE, typename BASE_LOOP>
	struct axpy_spread_loop<T, STR, STA, ITTYPE, false, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		axpy_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		axpy_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void sub0_first(T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void sub0(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void add0_first(T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void add0(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
	};

// LAST; FIRST
	template <typename T, typename STR, typename STA, typename ITTYPE, typename BASE_LOOP>
	struct axpy_spread_loop<T, STR, STA, ITTYPE, true, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		axpy_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		axpy_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void sub0_first(T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void sub0(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void add0_first(T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void add0(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void run_axpy(T *r, const T *a)
		{
			T *r0=r;
			T *ri=r;
			sub0_first(r0, ri);
			BASE_LOOP::run_axpy(r,a);
			ri=r0;
			add0_first(r0, ri);
		}
	};

// NOT LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename ITTYPE, typename HEAD, typename BASE_BASE_LOOP>
	struct axpy_spread_loop<T, STR, STA, ITTYPE, false, axpy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> >: public axpy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef axpy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
	public:
		axpy_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		axpy_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void sub0_first(T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			BASE_LOOP::sub0_first(r0,ri);
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void sub0(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void add0_first(T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			BASE_LOOP::add0_first(r0,ri);
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void add0(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
	};

// NOT LAST; FIRST
	template <typename T, typename STR, typename STA, typename ITTYPE, typename HEAD, typename BASE_BASE_LOOP>
	struct axpy_spread_loop<T, STR, STA, ITTYPE, true, axpy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> >: public axpy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef axpy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
	public:
		axpy_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		axpy_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void sub0_first(T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			BASE_LOOP::sub0_first(r0,ri);
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void sub0(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void add0_first(T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			if (!it.not_end(ri))
				return;
			BASE_LOOP::add0_first(r0,ri);
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void add0(const T *r0, T *ri)
		{
			ITTYPE it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void run_axpy(T *r, const T *a)
		{
			T *r0=r;
			T *ri=r;
			sub0_first(r0, ri);
			BASE_LOOP::run_axpy(r,a);
			ri=r0;
			add0_first(r0, ri);
		}
	};

	template <typename T, typename STR, typename STA, typename SPREAD, typename SUMA, typename RA, typename BASE_LOOP>
	struct axpy_create_general_loop
	{
		using type=BASE_LOOP;
	};

	template <typename T, typename STR, typename STA, typename SPREAD, typename HEAD, typename ... SUMA, typename RA, typename BASE_LOOP>
	struct axpy_create_general_loop<T, STR, STA, SPREAD, type_sequence<HEAD, SUMA...>, RA, BASE_LOOP>:
		public axpy_create_general_loop<T, STR, STA, SPREAD, type_sequence<SUMA...>, RA, axpy_suma_loop<T, STR, STA, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SPREAD, typename HEAD, typename ... RA, typename BASE_LOOP>
	struct axpy_create_general_loop<T, STR, STA, SPREAD, type_sequence<>, type_sequence<HEAD, RA...>, BASE_LOOP>:
		public axpy_create_general_loop<T, STR, STA, SPREAD, type_sequence<>, type_sequence<RA...>, axpy_ra_loop<T, STR, STA, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename HEAD, typename ... SPREAD, typename BASE_LOOP>
	struct axpy_create_general_loop<T, STR, STA, type_sequence<HEAD, SPREAD...>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public axpy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<>, type_sequence<>, axpy_spread_loop<T, STR, STA, typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T>, false, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename HEAD, typename BASE_LOOP>
	struct axpy_create_general_loop<T, STR, STA, type_sequence<HEAD>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public axpy_create_general_loop<T, STR, STA, type_sequence<>, type_sequence<>, type_sequence<>, axpy_spread_loop<T, STR, STA, typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T>, true, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename RAV>
	struct axpy_loop;

	template <typename T, typename STR, typename STA, int RAV_V_TYPE, int RAV_C_MASK, int RAV_NEXT_VALENCE, size_t ... RAV_POS, typename ... RA>
	struct axpy_loop<T, STR, STA, type_sequence<valence_data<RAV_V_TYPE, 3, RAV_C_MASK, 3, RAV_NEXT_VALENCE, RAV_POS...>, RA...> >
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	private:
		typedef type_sequence<valence_data<RAV_V_TYPE, 3, RAV_C_MASK, 3, RAV_NEXT_VALENCE, RAV_POS...>, RA...> VD;
		BLAS_INTEGER N;
		BLAS_INTEGER incr;
		BLAS_INTEGER inca;
	public:
		axpy_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1),
		N(vd_iterator_getter<STR, VD, 0>::template iterator_type<T>::length(str)),
		incr(vd_iterator_getter<STR, VD, 0>::template iterator_type<T>::step(str)),
		inca(vd_iterator_getter<STA, VD, 1>::template iterator_type<T>::step(sta))
		{}
		axpy_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha),
		N(vd_iterator_getter<STR, VD, 0>::template iterator_type<T>::length(str)),
		incr(vd_iterator_getter<STR, VD, 0>::template iterator_type<T>::step(str)),
		inca(vd_iterator_getter<STA, VD, 1>::template iterator_type<T>::step(sta))
		{}
		void run_axpy(T *r, const T *a)
		{
			axpy<T>(&N, &alpha, a, &inca, r, &incr);
		}
		void add0(const T *r0, T *ri)
		{
			axpy<T>(&N, &type_constants<T>::one, r0, &incr, ri, &incr);
		}
		void sub0(const T *r0, T *ri)
		{
			axpy<T>(&N, &type_constants<T>::m_one, r0, &incr, ri, &incr);
		}
	};

	template <typename T, typename STR, typename STA, typename S>
	struct axpy_just_spread_loop;

	template <typename T, typename STR, typename STA, int S_V_TYPE, int S_C_MASK, int S_NEXT_VALENCE, size_t ... S_POS, typename ... SPREAD>
	struct axpy_just_spread_loop<T, STR, STA, type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, SPREAD...> >
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	private:
		typedef type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, SPREAD...> VD;
		BLAS_INTEGER l;
		BLAS_INTEGER s;
	public:
		using ittype=typename vd_iterator_getter<STR, VD, 0>::template iterator_type<T>;
	public:
		axpy_just_spread_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1),
			l(ittype::length(str)), s(ittype::step(str))
		{}
		axpy_just_spread_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha),
			l(ittype::length(str)), s(ittype::step(str))
		{}
		void run_axpy(T *r, const T *a)
		{
			T res=*a*alpha;
			axpy(&l, &alpha, a, &type_constants<BLAS_INTEGER>::zero, r, &s);
		}
		void add0(const T *r0, T *ri)
		{
			axpy<T>(&l, &type_constants<T>::one, r0, &s, ri, &s);
		}
		void sub0(const T *r0, T *ri)
		{
			axpy<T>(&l, &type_constants<T>::m_one, r0, &s, ri, &s);
		}
	};

	template <typename T, typename STR, typename STA, typename S>
	struct axpy_sum_loop;

	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int NEXT_VALENCE, size_t ... POS, typename ... SUMA>
	struct axpy_sum_loop<T, STR, STA, type_sequence<valence_data<V_TYPE, 2, C_MASK, 2, NEXT_VALENCE, POS...>, SUMA...> >
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	private:
		typedef type_sequence<valence_data<V_TYPE, 2, C_MASK, 2, NEXT_VALENCE, POS...>, SUMA...> VD;
		BLAS_INTEGER l;
		BLAS_INTEGER s;
	public:
		using ittype=typename vd_iterator_getter<STA, VD, 1>::template iterator_type<T>;
	public:
		axpy_sum_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1),
			l(ittype::length(sta)), s(ittype::step(sta))
		{}
		axpy_sum_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha),
			l(ittype::length(sta)), s(ittype::step(sta))
		{}
		void run_axpy(T *r, const T *a)
		{
			*r+=dot<T>(&l, a, &s, &alpha, &type_constants<BLAS_INTEGER>::zero);
		}
		void add0(const T *r0, T *ri)
		{
			*ri+=*r0;
		}
		void sub0(const T *r0, T *ri)
		{
			*ri-=*r0;
		}
	};

	template <typename T, typename STR, typename STA, typename SPREAD, typename SUMA, typename RA>
	struct axpy_test_sum;

//	DONE sum
	template <typename T, typename STR, typename STA, typename ... SPREAD, int V_TYPE, int C_MASK, int NEXT_VALENCE, size_t ... POS, typename ... SUMA0, typename ... SUMA, typename ... RA>
	struct axpy_test_sum<T, STR, STA, type_sequence<SPREAD...>, type_sequence<type_sequence<valence_data<V_TYPE, 2, C_MASK, 2, NEXT_VALENCE, POS...>, SUMA0...>, SUMA...>, type_sequence<RA...> >:
		public axpy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, axpy_sum_loop<T, STR, STA, type_sequence<valence_data<V_TYPE, 2, C_MASK, 2, NEXT_VALENCE, POS...>, SUMA0...> > >
//		public gemvA<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<valence_data<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...> >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename ... SPREAD, typename ... SUMA, typename ... RA>
	struct axpy_test_sum<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...> >:
		public axpy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, axpy_common_loop<T, STR, STA> >
//		public common_final<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...> >
	{
	};

	template <typename T, typename STR, typename STA, typename VD>
	struct axpy_runner;

//	DONE axpy
	template <typename T, typename STR, typename STA, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, int RA_V_TYPE, int RA_C_MASK, int RA_NEXT_VALENCE, size_t ... RA_POS, typename ... RA0, typename ... RA>
	struct axpy_runner<T, STR, STA, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>, RA...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public axpy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, axpy_loop<T, STR, STA, type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...> > >
	{
	};

//	DONE just spread
	template <typename T, typename STR, typename STA, typename CONT_MASK, int S_V_TYPE, int S_C_MASK, int S_NEXT_VALENCE, size_t ... S_POS, typename ... SPREAD0, typename ... SPREAD>
	struct axpy_runner<T, STR, STA, type_sequence<CONT_MASK, type_sequence<type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, SPREAD0...>, SPREAD...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public axpy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<>, type_sequence<>, axpy_just_spread_loop<T, STR, STA, type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, SPREAD0...> > >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, typename ... RA>
	struct axpy_runner<T, STR, STA, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public axpy_test_sum<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...> >
	{
	};

};




#endif /* TAXPY_H_ */
