/*
 * tcopy.h
 *
 *  Created on: 21 мая 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef TCOPY_H_
#define TCOPY_H_


#include <tvalence.h>
#include <blas_tmpl.h>
#include <smplmath.h>

namespace iTTL
{

	template <typename T, typename STR, typename STA, typename HEAD, typename BASE_LOOP>
	struct copy_suma_loop: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STA, HEAD, 1>::template iterator_type<T> itype;
		copy_suma_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		copy_suma_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void run_copy(T *r, const T *a)
		{
			itype it(this->sta, a);
			if (!it.not_end(a))
			{
				*r=0;
				return;
			}
			BASE_LOOP::run_copy(r,a);
			for (it.move_one(a);it.not_end(a);it.move_one(a))
				BASE_LOOP::add(r,a);
		}
		void run_asum(T *r, const T *a)
		{
			itype it(this->sta, a);
			if (!it.not_end(a))
			{
				*r=0;
				return;
			}
			BASE_LOOP::run_asum(r,a);
			for (it.move_one(a);it.not_end(a);it.move_one(a))
				BASE_LOOP::add_asum(r,a);
		}
		void add(T *r, const T *a)
		{
			itype it(this->sta, a);
			for (;it.not_end(a);it.move_one(a))
				BASE_LOOP::add(r,a);
		}
		void add_asum(T *r, const T *a)
		{
			itype it(this->sta, a);
			for (;it.not_end(a);it.move_one(a))
				BASE_LOOP::add_asum(r,a);
		}
	};

	template <typename T, typename STR, typename STA, typename HEAD, typename BASE_LOOP>
	struct copy_ra_loop: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T> ritype;
		typedef typename vd_iterator_getter<STA, HEAD, 1>::template iterator_type<T> aitype;
		copy_ra_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		copy_ra_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
		{}
		void run_copy(T *r, const T *a)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a))
				BASE_LOOP::run_copy(r,a);
		}
		void run_asum(T *r, const T *a)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a))
				BASE_LOOP::run_asum(r,a);
		}
		void add(T *r, const T *a)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a))
				BASE_LOOP::add(r,a);
		}
		void add_asum(T *r, const T *a)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a))
				BASE_LOOP::add_asum(r,a);
		}
		void spread(const T *r0, T *ri)
		{
			ritype r0t(this->str, r0);
			ritype rit(this->str, ri);
			for (;rit.not_end(ri);r0t.move_one(r0),rit.move_one(ri))
				BASE_LOOP::spread(r0,ri);
		}
	};

	template <typename T, typename STR, typename STA>
	struct copy_common_loop
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	public:
		copy_common_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1)
		{}
		copy_common_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha)
		{}
		void run_copy(T *r, const T *a)
		{
			*r=*a;
		}
		void run_asum(T *r, const T *a)
		{
			*r=abs<T>(*a);
		}
		void add(T *r, const T *a)
		{
			*r+=*a;
		}
		void add_asum(T *r, const T *a)
		{
			*r+=abs<T>(*a);
		}
		void spread(const T *r0, T *ri)
		{
			*ri=*r0;
		}
	};

	template <typename T, typename STR, typename STA, typename ITTYPE, bool IS_FIRST, typename BASE_LOOP>
	struct copy_spread_loop;

// LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename ITTYPE, typename BASE_LOOP>
	struct copy_spread_loop<T, STR, STA, ITTYPE, false, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		copy_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		copy_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
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
	struct copy_spread_loop<T, STR, STA, ITTYPE, true, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		copy_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		copy_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
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
		void run_copy(T *r, const T *a)
		{
			T *r0=r;
			T *ri=r;
			BASE_LOOP::run_copy(r,a);
			ri=r0;
			spread_first(r0, ri);
		}
		void run_asum(T *r, const T *a)
		{
			T *r0=r;
			T *ri=r;
			BASE_LOOP::run_asum(r,a);
			ri=r0;
			spread_first(r0, ri);
		}
		void add(T *r, const T *a)
		{
			ITTYPE it(this->str, r);
			for (;it.not_end(r);it.move_one(r))
				BASE_LOOP::add(r,a);
		}
		void add_asum(T *r, const T *a)
		{
			ITTYPE it(this->str, r);
			for (;it.not_end(r);it.move_one(r))
				BASE_LOOP::add_asum(r,a);
		}
	};




// NOT LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename ITTYPE, typename HEAD, typename BASE_BASE_LOOP>
	struct copy_spread_loop<T, STR, STA, ITTYPE, false, copy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> >: public copy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef copy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
	public:
		copy_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		copy_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
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
	struct copy_spread_loop<T, STR, STA, ITTYPE, true, copy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> >: public copy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef copy_spread_loop<T, STR, STA, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
	public:
		copy_spread_loop(const STR& str, const STA& sta): BASE_LOOP(str, sta)
		{}
		copy_spread_loop(const STR& str, const STA& sta, T alpha): BASE_LOOP(str, sta, alpha)
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
		void run_copy(T *r, const T *a)
		{
			T *r0=r;
			T *ri=r;
			BASE_LOOP::run_copy(r,a);
			ri=r0;
			spread_first(r0, ri);
		}
		void run_asum(T *r, const T *a)
		{
			T *r0=r;
			T *ri=r;
			BASE_LOOP::run_asum(r,a);
			ri=r0;
			spread_first(r0, ri);
		}
		void add(T *r, const T *a)
		{
			ITTYPE it(this->str, r);
			for (;it.not_end(r);it.move_one(r))
				BASE_LOOP::add(r,a);
		}
		void add_asum(T *r, const T *a)
		{
			ITTYPE it(this->str, r);
			for (;it.not_end(r);it.move_one(r))
				BASE_LOOP::add_asum(r,a);
		}
	};

	template <typename T, typename STR, typename STA, typename SPREAD, typename SUMA, typename RA, typename BASE_LOOP>
	struct copy_create_general_loop
	{
		using type=BASE_LOOP;
	};

	template <typename T, typename STR, typename STA, typename SPREAD, typename HEAD, typename ... SUMA, typename RA, typename BASE_LOOP>
	struct copy_create_general_loop<T, STR, STA, SPREAD, type_sequence<HEAD, SUMA...>, RA, BASE_LOOP>:
		public copy_create_general_loop<T, STR, STA, SPREAD, type_sequence<SUMA...>, RA, copy_suma_loop<T, STR, STA, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SPREAD, typename HEAD, typename ... RA, typename BASE_LOOP>
	struct copy_create_general_loop<T, STR, STA, SPREAD, type_sequence<>, type_sequence<HEAD, RA...>, BASE_LOOP>:
		public copy_create_general_loop<T, STR, STA, SPREAD, type_sequence<>, type_sequence<RA...>, copy_ra_loop<T, STR, STA, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename HEAD, typename ... SPREAD, typename BASE_LOOP>
	struct copy_create_general_loop<T, STR, STA, type_sequence<HEAD, SPREAD...>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public copy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<>, type_sequence<>, copy_spread_loop<T, STR, STA, typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T>, false, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename HEAD, typename BASE_LOOP>
	struct copy_create_general_loop<T, STR, STA, type_sequence<HEAD>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public copy_create_general_loop<T, STR, STA, type_sequence<>, type_sequence<>, type_sequence<>, copy_spread_loop<T, STR, STA, typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T>, true, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename RAC, typename RAV>
	struct lacpy_loop
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER lda;
		BLAS_INTEGER ldb;
	public:
		lacpy_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1),
		M(vd_iterator_getter<STR, RAC, 0>::template iterator_type<T>::length(str)),
		N(vd_iterator_getter<STR, RAV, 0>::template iterator_type<T>::length(str)),
		lda(vd_iterator_getter<STR, RAV, 0>::template iterator_type<T>::step(str)),
		ldb(vd_iterator_getter<STA, RAV, 1>::template iterator_type<T>::step(sta))
		{}
		lacpy_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha),
		M(vd_iterator_getter<STR, RAC, 0>::template iterator_type<T>::length(str)),
		N(vd_iterator_getter<STR, RAV, 0>::template iterator_type<T>::length(str)),
		lda(vd_iterator_getter<STR, RAV, 0>::template iterator_type<T>::step(str)),
		ldb(vd_iterator_getter<STA, RAV, 1>::template iterator_type<T>::step(sta))
		{}
		void run_copy(T *r, const T *a)
		{
			lacpy<T>("O", &M, &N, a, &ldb, r, &lda);
		}
		void add(T *r, const T *a)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
				axpy<T>(&M, &type_constants<T>::one, a+i*ldb, &type_constants<BLAS_INTEGER>::one, r+i*lda, &type_constants<BLAS_INTEGER>::one);
		}
		void spread(const T *r0, T *ri)
		{
			lacpy<T>("O", &M, &N, r0, &lda, ri, &lda);
		}
	};

	template <typename T, typename STR, typename STA, typename RAV>
	struct copy_loop
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	private:
//		segment<RAV_V_TYPE, USAGE_ENUM, 0, 0> srrav;
		BLAS_INTEGER N;
		BLAS_INTEGER incr;
		BLAS_INTEGER inca;
	public:
		copy_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1),
		N(vd_iterator_getter<STR, RAV, 0>::template iterator_type<T>::length(str)),
		incr(vd_iterator_getter<STR, RAV, 0>::template iterator_type<T>::step(str)),
		inca(vd_iterator_getter<STA, RAV, 1>::template iterator_type<T>::step(sta))
		{}
		copy_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha),
		N(vd_iterator_getter<STR, RAV, 0>::template iterator_type<T>::length(str)),
		incr(vd_iterator_getter<STR, RAV, 0>::template iterator_type<T>::step(str)),
		inca(vd_iterator_getter<STA, RAV, 1>::template iterator_type<T>::step(sta))
		{}
		void run_copy(T *r, const T *a)
		{
			copy<T>(&N, a, &inca, r, &incr);
		}
		void add(T *r, const T *a)
		{
			axpy<T>(&N, &type_constants<T>::one, a, &inca, r, &incr);
		}
		void spread(const T *r0, T *ri)
		{
			copy<T>(&N, r0, &incr, ri, &incr);
		}
	};

	template <typename T, typename STR, typename STA, typename S>
	struct just_spread_loop
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	private:
		BLAS_INTEGER l;
		BLAS_INTEGER s;
	public:
		using ittype=typename vd_iterator_getter<STR, S, 0>::template iterator_type<T>;
	public:
		just_spread_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1),
			l(ittype::length(str)), s(ittype::step(str))
		{}
		just_spread_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha),
			l(ittype::length(str)), s(ittype::step(str))
		{}
		void run_copy(T *r, const T *a)
		{
			copy<T>(&l, a, &type_constants<BLAS_INTEGER>::zero, r, &s);
		}
		void run_asum(T *r, const T *a)
		{
			T absa=abs<T>(*a);
			copy<T>(&l, &absa, &type_constants<BLAS_INTEGER>::zero, r, &s);
		}
		void add(T *r, const T *a)
		{
			axpy<T>(&l, &type_constants<T>::one, a, &type_constants<BLAS_INTEGER>::zero, r, &s);
		}
		void add_asum(T *r, const T *a)
		{
			T absa=abs<T>(*a);
			axpy<T>(&l, &type_constants<T>::one, &absa, &type_constants<BLAS_INTEGER>::zero, r, &s);
		}
		void spread(const T *r0, T *ri)
		{
			copy<T>(&l, r0, &s, ri, &s);
		}
	};

	template <typename T, typename STR, typename STA, typename S>
	struct asum_loop
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	private:
		BLAS_INTEGER l;
		BLAS_INTEGER s;
	public:
		using ittype=typename vd_iterator_getter<STA, S, 1>::template iterator_type<T>;
	public:
		asum_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1),
			l(ittype::length(sta)), s(ittype::step(sta))
		{}
		asum_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha),
			l(ittype::length(sta)), s(ittype::step(sta))
		{}
		void run_copy(T *r, const T *a)
		{
			*r=dot<T>(&l, a, &s, &type_constants<T>::one, &type_constants<BLAS_INTEGER>::zero);
		}
		void run_asum(T *r, const T *a)
		{
			*r=asum<T>(&l, a, &s);
		}
		void add(T *r, const T *a)
		{
			*r+=dot<T>(&l, a, &s, &type_constants<T>::one, &type_constants<BLAS_INTEGER>::zero);
		}
		void add_asum(T *r, const T *a)
		{
			*r+=asum<T>(&l, a, &s);
		}
		void spread(const T *r0, T *ri)
		{
			*ri=*r0;
		}
	};


	template <typename T, typename STR, typename STA, typename SPREAD, typename SUMA, typename RA>
	struct copy_test_sum;

//	DONE sum
	template <typename T, typename STR, typename STA, typename ... SPREAD, int V_TYPE, int C_MASK, int NEXT_VALENCE, size_t ... POS, typename ... SUMA0, typename ... SUMA, typename ... RA>
	struct copy_test_sum<T, STR, STA, type_sequence<SPREAD...>, type_sequence<type_sequence<valence_data<V_TYPE, 2, C_MASK, 2, NEXT_VALENCE, POS...>, SUMA0...>, SUMA...>, type_sequence<RA...> >:
		public copy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, asum_loop<T, STR, STA, type_sequence<valence_data<V_TYPE, 2, C_MASK, 2, NEXT_VALENCE, POS...>, SUMA0...> > >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename ... SPREAD, typename ... RA, typename ... SUMA>
	struct copy_test_sum<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...> >:
		public copy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, copy_common_loop<T, STR, STA> >
	{
	};

	template <typename T, typename STR, typename STA, typename VD>
	struct copy_runner;

//	DONE lacpy
	template <typename T, typename STR, typename STA, typename ... SPREAD, typename ... SUMA, int RAC_V_TYPE, int RAC_NEXT_VALENCE, size_t ... RAC_POS, int RAV_V_TYPE, int RAV_C_MASK, int RAV_NEXT_VALENCE, size_t ... RAV_POS, typename ... RA0, typename ... RA1, typename ... RA>
	struct copy_runner<T, STR, STA, type_sequence<std::integral_constant<int,3>, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<type_sequence<valence_data<RAC_V_TYPE, 3, 3, 3, RAC_NEXT_VALENCE, RAC_POS...>, RA0...>, type_sequence<valence_data<RAV_V_TYPE, 3, RAV_C_MASK, 3, RAV_NEXT_VALENCE, RAV_POS...>, RA1...>, RA...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public copy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, lacpy_loop<T, STR, STA, type_sequence<valence_data<RAC_V_TYPE, 3, 3, 3, RAC_NEXT_VALENCE, RAC_POS...>, RA0...>, type_sequence<valence_data<RAV_V_TYPE, 3, RAV_C_MASK, 3, RAV_NEXT_VALENCE, RAV_POS...>, RA1...> > >
	{
	};
//	DONE copy
	template <typename T, typename STR, typename STA, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, int RA_V_TYPE, int RA_C_MASK, int RA_NEXT_VALENCE, size_t ... RA_POS, typename ... RA0, typename ... RA>
	struct copy_runner<T, STR, STA, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>, RA...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public copy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, copy_loop<T, STR, STA, type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...> > >
	{
	};

//	DONE just spread
	template <typename T, typename STR, typename STA, typename CONT_MASK, int S_V_TYPE, int S_C_MASK, int S_NEXT_VALENCE, size_t ... S_POS, typename ... SPREAD0, typename ... SPREAD>
	struct copy_runner<T, STR, STA, type_sequence<CONT_MASK, type_sequence<type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, SPREAD0...>, SPREAD...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public copy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<>, type_sequence<>, just_spread_loop<T, STR, STA, type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, SPREAD0...> > >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, typename ... RA>
	struct copy_runner<T, STR, STA, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public copy_test_sum<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...> >
	{
	};


	template <typename T, typename STR, typename STA, typename VD>
	struct asum_runner;

// DONE ASUM
	template <typename T, typename STR, typename STA, typename CONT_MASK, typename ... SPREAD, int V_TYPE, int C_MASK, int NEXT_VALENCE, size_t ... POS, typename ... SUMA0, typename ... SUMA, typename ... RA>
	struct asum_runner<T, STR, STA, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<type_sequence<valence_data<V_TYPE, 2, C_MASK, 2, NEXT_VALENCE, POS...>, SUMA0...>, SUMA...>, type_sequence<RA...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public copy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, asum_loop<T, STR, STA, type_sequence<valence_data<V_TYPE, 2, C_MASK, 2, NEXT_VALENCE, POS...>, SUMA0...> > >
	{
	};

//	DONE just spread
	template <typename T, typename STR, typename STA, typename CONT_MASK, int S_V_TYPE, int S_C_MASK, int S_NEXT_VALENCE, size_t ... S_POS, typename ... SPREAD0, typename ... SPREAD>
	struct asum_runner<T, STR, STA, type_sequence<CONT_MASK, type_sequence<type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, SPREAD0...>, SPREAD...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public copy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<>, type_sequence<>, just_spread_loop<T, STR, STA, type_sequence<valence_data<S_V_TYPE, 1, S_C_MASK, 1, S_NEXT_VALENCE, S_POS...>, SPREAD0...> > >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, typename ... RA>
	struct asum_runner<T, STR, STA, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public copy_create_general_loop<T, STR, STA, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, copy_common_loop<T, STR, STA> >
	{
	};

};

#endif /* TCOPY_H_ */
