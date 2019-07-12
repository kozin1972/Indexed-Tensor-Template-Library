/*
 * gem.h
 *
 * General Tensor Multiplication
 * May call gemm, gemv, dot, ger, axpy, scal
 *
 *  Created on: 13 мая 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef TGEM_H_
#define TGEM_H_

#include <tvalence.h>
#include <blas_tmpl.h>

namespace iTTL
{

	template <typename T, typename STR, typename STA, typename STB, typename HEAD, typename BASE_LOOP>
	struct gem_sumb_loop: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STB, HEAD, 2>::template iterator_type<T> itype;
		gem_sumb_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): BASE_LOOP(str, sta, stb, alpha, beta)
		{}
		void run(T *r, const T *a, const T *b)
		{
			itype it(this->stb, b);
			if (!it.not_end(b))
				return;
			BASE_LOOP::run(r,a,b);
			for (it.move_one(b);it.not_end(b);it.move_one(b))
				BASE_LOOP::add(r,a,b);
		}
		void add(T *r, const T *a, const T *b)
		{
			itype it(this->stb, b);
			for (;it.not_end(b);it.move_one(b))
				BASE_LOOP::add(r,a,b);
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename HEAD, typename BASE_LOOP>
	struct gem_ab_loop: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STA, HEAD, 1>::template iterator_type<T> aitype;
		typedef typename vd_iterator_getter<STB, HEAD, 2>::template iterator_type<T> bitype;
		gem_ab_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): BASE_LOOP(str, sta, stb, alpha, beta)
		{}
		void run(T *r, const T *a, const T *b)
		{
			aitype ait(this->sta, a);
			bitype bit(this->stb, b);
			if (!ait.not_end(a))
				return;
			BASE_LOOP::run(r,a,b);
			for (ait.move_one(a),bit.move_one(b);ait.not_end(a);ait.move_one(a),bit.move_one(b))
				BASE_LOOP::add(r,a,b);
		}
		void add(T *r, const T *a, const T *b)
		{
			aitype ait(this->sta, a);
			bitype bit(this->stb, b);
			for (;ait.not_end(a);ait.move_one(a),bit.move_one(b))
				BASE_LOOP::add(r,a,b);
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename HEAD, typename BASE_LOOP>
	struct gem_rb_loop: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T> ritype;
		typedef typename vd_iterator_getter<STB, HEAD, 2>::template iterator_type<T> bitype;
		gem_rb_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): BASE_LOOP(str, sta, stb, alpha, beta)
		{}
		void run(T *r, const T *a, const T *b)
		{
			ritype rit(this->str, r);
			bitype bit(this->stb, b);
			for (;rit.not_end(r);rit.move_one(r),bit.move_one(b))
				BASE_LOOP::run(r,a,b);
		}
		void add(T *r, const T *a, const T *b)
		{
			ritype rit(this->str, r);
			bitype bit(this->stb, b);
			for (;rit.not_end(r);rit.move_one(r),bit.move_one(b))
				BASE_LOOP::add(r,a,b);
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

	template <typename T, typename STR, typename STA, typename STB, typename HEAD, typename BASE_LOOP>
	struct gem_suma_loop: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STA, HEAD, 1>::template iterator_type<T> itype;
		gem_suma_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): BASE_LOOP(str, sta, stb, alpha, beta)
		{}
		void run(T *r, const T *a, const T *b)
		{
			itype it(this->sta, a);
			if (!it.not_end(a))
				return;
			BASE_LOOP::run(r,a,b);
			for (it.move_one(a);it.not_end(a);it.move_one(a))
				BASE_LOOP::add(r,a,b);
		}
		void add(T *r, const T *a, const T *b)
		{
			itype it(this->sta, a);
			for (;it.not_end(a);it.move_one(a))
				BASE_LOOP::add(r,a,b);
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename HEAD, typename BASE_LOOP>
	struct gem_ra_loop: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T> ritype;
		typedef typename vd_iterator_getter<STA, HEAD, 1>::template iterator_type<T> aitype;
		gem_ra_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): BASE_LOOP(str, sta, stb, alpha, beta)
		{}
		void run(T *r, const T *a, const T *b)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a))
				BASE_LOOP::run(r,a,b);
		}
		void add(T *r, const T *a, const T *b)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a))
				BASE_LOOP::add(r,a,b);
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

	template <typename T, typename STR, typename STA, typename STB, typename HEAD, typename BASE_LOOP>
	struct gem_com_loop: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T> ritype;
		typedef typename vd_iterator_getter<STA, HEAD, 1>::template iterator_type<T> aitype;
		typedef typename vd_iterator_getter<STB, HEAD, 2>::template iterator_type<T> bitype;
		gem_com_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): BASE_LOOP(str, sta, stb, alpha, beta)
		{}
		void run(T *r, const T *a, const T *b)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			bitype bit(this->stb, b);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a),bit.move_one(b))
				BASE_LOOP::run(r,a,b);
		}
		void add(T *r, const T *a, const T *b)
		{
			ritype rit(this->str, r);
			aitype ait(this->sta, a);
			bitype bit(this->stb, b);
			for (;rit.not_end(r);rit.move_one(r),ait.move_one(a),bit.move_one(b))
				BASE_LOOP::add(r,a,b);
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

	template <typename T, typename STR, typename STA, typename STB, typename HEAD, bool IS_FIRST, typename BASE_LOOP>
	struct gem_spread_loop;

// LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename STB, typename HEAD, typename BASE_LOOP>
	struct gem_spread_loop<T, STR, STA, STB, HEAD, false, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T> itype;
		gem_spread_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta):
			BASE_LOOP(str, sta, stb, alpha, beta)
		{}
		void sub0_first(T *r0, T *ri)
		{
			itype it(this->str,ri);
			if (!it.not_end(ri))
				return;
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void sub0(const T *r0, T *ri)
		{
			itype it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void add0_first(T *r0, T *ri)
		{
			itype it(this->str,ri);
			if (!it.not_end(ri))
				return;
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void add0(const T *r0, T *ri)
		{
			itype it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
	};

// LAST; FIRST
	template <typename T, typename STR, typename STA, typename STB, typename HEAD, typename BASE_LOOP>
	struct gem_spread_loop<T, STR, STA, STB, HEAD, true, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename vd_iterator_getter<STR, HEAD, 0>::template iterator_type<T> itype;
		gem_spread_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta):
			BASE_LOOP(str, sta, stb, alpha, beta)
		{}
		void sub0_first(T *r0, T *ri)
		{
			itype it(this->str,ri);
			if (!it.not_end(ri))
				return;
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void sub0(const T *r0, T *ri)
		{
			itype it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void add0_first(T *r0, T *ri)
		{
			itype it(this->str,ri);
			if (!it.not_end(ri))
				return;
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void add0(const T *r0, T *ri)
		{
			itype it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void run(T *r, const T *a, const T *b)
		{
			T *r0=r;
			T *ri=r;
			sub0_first(r0, ri);
			BASE_LOOP::run(r,a,b);
			ri=r0;
			add0_first(r0, ri);
		}
	};

// NOT LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename STB, typename H0, typename HEAD, typename BASE_BASE_LOOP>
	struct gem_spread_loop<T, STR, STA, STB, H0, false, gem_spread_loop<T, STR, STA, STB, HEAD, false, BASE_BASE_LOOP> >: public gem_spread_loop<T, STR, STA, STB, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef gem_spread_loop<T, STR, STA, STB, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
	public:
		typedef typename vd_iterator_getter<STR, H0, 0>::template iterator_type<T> itype;
		gem_spread_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta):
			BASE_LOOP(str, sta, stb, alpha, beta)
		{}
		void sub0_first(T *r0, T *ri)
		{
			itype it(this->str,ri);
			if (!it.not_end(ri))
				return;
			BASE_LOOP::sub0_first(r0,ri);
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void sub0(const T *r0, T *ri)
		{
			itype it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void add0_first(T *r0, T *ri)
		{
			itype it(this->str,ri);
			if (!it.not_end(ri))
				return;
			BASE_LOOP::add0_first(r0,ri);
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void add0(const T *r0, T *ri)
		{
			itype it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
	};

// NOT LAST; FIRST
	template <typename T, typename STR, typename STA, typename STB, typename H0, typename HEAD, typename BASE_BASE_LOOP>
	struct gem_spread_loop<T, STR, STA, STB, H0, true, gem_spread_loop<T, STR, STA, STB, HEAD, false, BASE_BASE_LOOP> >: public gem_spread_loop<T, STR, STA, STB, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef gem_spread_loop<T, STR, STA, STB, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
	public:
		typedef typename vd_iterator_getter<STR, H0, 0>::template iterator_type<T> itype;
		gem_spread_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta):
			BASE_LOOP(str, sta, stb, alpha, beta)
		{}
		void sub0_first(T *r0, T *ri)
		{
			itype it(this->str,ri);
			if (!it.not_end(ri))
				return;
			BASE_LOOP::sub0_first(r0,ri);
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void sub0(const T *r0, T *ri)
		{
			itype it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::sub0(r0,ri);
		}
		void add0_first(T *r0, T *ri)
		{
			itype it(this->str,ri);
			if (!it.not_end(ri))
				return;
			BASE_LOOP::add0_first(r0,ri);
			it.move_one(ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void add0(const T *r0, T *ri)
		{
			itype it(this->str,ri);
			for (;it.not_end(ri);it.move_one(ri))
				BASE_LOOP::add0(r0,ri);
		}
		void run(T *r, const T *a, const T *b)
		{
			T *r0=r;
			T *ri=r;
			sub0_first(r0, ri);
			BASE_LOOP::run(r,a,b);
			ri=r0;
			add0_first(r0, ri);
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename SPREAD, typename SUMA, typename RA, typename SUMB, typename RB, typename AB, typename COMMON, typename BASE_LOOP>
	struct gem_create_general_loop
	{
		using type=BASE_LOOP;
	};

	template <typename T, typename STR, typename STA, typename STB, typename SPREAD, typename SUMA, typename RA, typename HEAD, typename ... SUMB, typename RB, typename AB, typename COMMON, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SPREAD, SUMA, RA, type_sequence<HEAD, SUMB...>, RB, AB, COMMON, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SPREAD, SUMA, RA, type_sequence<SUMB...>, RB, AB, COMMON, gem_sumb_loop<T, STR, STA, STB, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SPREAD, typename SUMA, typename RA, typename RB, typename HEAD, typename ... AB, typename COMMON, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SPREAD, SUMA, RA, type_sequence<>, RB, type_sequence<HEAD, AB...>, COMMON, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SPREAD, SUMA, RA, type_sequence<>, RB, type_sequence<AB...>, COMMON, gem_ab_loop<T, STR, STA, STB, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SPREAD, typename SUMA, typename RA, typename HEAD, typename ... RB, typename COMMON, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SPREAD, SUMA, RA, type_sequence<>, type_sequence<HEAD, RB...>, type_sequence<>, COMMON, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SPREAD, SUMA, RA, type_sequence<>, type_sequence<RB...>, type_sequence<>, COMMON, gem_rb_loop<T, STR, STA, STB, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SPREAD, typename HEAD, typename ... SUMA, typename RA, typename COMMON, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SPREAD, type_sequence<HEAD, SUMA...>, RA, type_sequence<>, type_sequence<>, type_sequence<>, COMMON, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SPREAD, type_sequence<SUMA...>, RA, type_sequence<>, type_sequence<>, type_sequence<>, COMMON, gem_suma_loop<T, STR, STA, STB, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SPREAD, typename HEAD, typename ... RA, typename COMMON, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SPREAD, type_sequence<>, type_sequence<HEAD, RA...>, type_sequence<>, type_sequence<>, type_sequence<>, COMMON, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SPREAD, type_sequence<>, type_sequence<RA...>, type_sequence<>, type_sequence<>, type_sequence<>, COMMON, gem_ra_loop<T, STR, STA, STB, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SPREAD, typename HEAD, typename ... COMMON, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SPREAD, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<HEAD, COMMON...>, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SPREAD, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<COMMON...>, gem_com_loop<T, STR, STA, STB, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename HEAD, typename ... SPREAD, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, type_sequence<HEAD, SPREAD...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, gem_spread_loop<T, STR, STA, STB, HEAD, false, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename HEAD, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, type_sequence<HEAD>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, gem_spread_loop<T, STR, STA, STB, HEAD, true, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB, typename AB, int VARIANT>
	struct gemm_loop;

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB, typename AB>
	struct gemm_loop<T, STR, STA, STB, RA, RB, AB, 6>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER common;
		BLAS_INTEGER lda;
		BLAS_INTEGER ldb;
		BLAS_INTEGER ldc;
	public:
		gemm_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::length(str)),
		N(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::length(str)),
		common(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::length(stb)),
		lda(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::step(sta)),
		ldb(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::step(stb)),
		ldc(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemm<T>("N", "N", &M, &N, &common, &alpha, b, &ldb, a, &lda, &beta, r, &ldc);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemm<T>("N", "N", &M, &N, &common, &alpha, b, &ldb, a, &lda, &type_constants<T>::one, r, &ldc);
		}
		void add0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
				axpy<T>(&M,&type_constants<T>::one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
		}
		void sub0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
			{
				axpy<T>(&M,&type_constants<T>::m_one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
				scal<T>(&M,&beta,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
			}
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB, typename AB>
	struct gemm_loop<T, STR, STA, STB, RA, RB, AB, 7>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER common;
		BLAS_INTEGER lda;
		BLAS_INTEGER ldb;
		BLAS_INTEGER ldc;
	public:
		gemm_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::length(sta)),
		N(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::length(str)),
		common(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::length(stb)),
		lda(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::step(sta)),
		ldb(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::step(stb)),
		ldc(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemm<T>("T", "T", &M, &N, &common, &alpha, a, &lda, b, &ldb, &beta, r, &ldc);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemm<T>("T", "T", &M, &N, &common, &alpha, a, &lda, b, &ldb, &type_constants<T>::one, r, &ldc);
		}
		void add0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
				axpy<T>(&M,&type_constants<T>::one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
		}
		void sub0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
			{
				axpy<T>(&M,&type_constants<T>::m_one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
				scal<T>(&M,&beta,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
			}
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB, typename AB>
	struct gemm_loop<T, STR, STA, STB, RA, RB, AB, 4>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER common;
		BLAS_INTEGER lda;
		BLAS_INTEGER ldb;
		BLAS_INTEGER ldc;
	public:
		gemm_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::length(str)),
		N(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::length(str)),
		common(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::length(stb)),
		lda(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::step(sta)),
		ldb(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::step(stb)),
		ldc(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemm<T>("N", "T", &M, &N, &common, &alpha, b, &ldb, a, &lda, &beta, r, &ldc);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemm<T>("N", "T", &M, &N, &common, &alpha, b, &ldb, a, &lda, &type_constants<T>::one, r, &ldc);
		}
		void add0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
				axpy<T>(&M,&type_constants<T>::one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
		}
		void sub0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
			{
				axpy<T>(&M,&type_constants<T>::m_one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
				scal<T>(&M,&beta,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
			}
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB, typename AB>
	struct gemm_loop<T, STR, STA, STB, RA, RB, AB, 5>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER common;
		BLAS_INTEGER lda;
		BLAS_INTEGER ldb;
		BLAS_INTEGER ldc;
	public:
		gemm_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::length(str)),
		N(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::length(str)),
		common(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::length(stb)),
		lda(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::step(sta)),
		ldb(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::step(stb)),
		ldc(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemm<T>("N", "T", &M, &N, &common, &alpha, a, &lda, b, &ldb, &beta, r, &ldc);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemm<T>("N", "T", &M, &N, &common, &alpha, a, &lda, b, &ldb, &type_constants<T>::one, r, &ldc);
		}
		void add0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
				axpy<T>(&M,&type_constants<T>::one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
		}
		void sub0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
			{
				axpy<T>(&M,&type_constants<T>::m_one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
				scal<T>(&M,&beta,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
			}
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB, typename AB>
	struct gemm_loop<T, STR, STA, STB, RA, RB, AB, 2>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER common;
		BLAS_INTEGER lda;
		BLAS_INTEGER ldb;
		BLAS_INTEGER ldc;
	public:
		gemm_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::length(stb)),
		N(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::length(str)),
		common(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::length(sta)),
		lda(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::step(sta)),
		ldb(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::step(stb)),
		ldc(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemm<T>("T", "N", &M, &N, &common, &alpha, b, &ldb, a, &lda, &beta, r, &ldc);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemm<T>("T", "N", &M, &N, &common, &alpha, b, &ldb, a, &lda, &type_constants<T>::one, r, &ldc);
		}
		void add0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
				axpy<T>(&M,&type_constants<T>::one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
		}
		void sub0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
			{
				axpy<T>(&M,&type_constants<T>::m_one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
				scal<T>(&M,&beta,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
			}
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB, typename AB>
	struct gemm_loop<T, STR, STA, STB, RA, RB, AB, 3>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER common;
		BLAS_INTEGER lda;
		BLAS_INTEGER ldb;
		BLAS_INTEGER ldc;
	public:
		gemm_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::length(sta)),
		N(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::length(stb)),
		common(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::length(sta)),
		lda(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::step(sta)),
		ldb(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::step(stb)),
		ldc(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemm<T>("T", "N", &M, &N, &common, &alpha, a, &lda, b, &ldb, &beta, r, &ldc);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemm<T>("T", "N", &M, &N, &common, &alpha, a, &lda, b, &ldb, &type_constants<T>::one, r, &ldc);
		}
		void add0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
				axpy<T>(&M,&type_constants<T>::one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
		}
		void sub0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
			{
				axpy<T>(&M,&type_constants<T>::m_one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
				scal<T>(&M,&beta,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
			}
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB, typename AB>
	struct gemm_loop<T, STR, STA, STB, RA, RB, AB, 0>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER common;
		BLAS_INTEGER lda;
		BLAS_INTEGER ldb;
		BLAS_INTEGER ldc;
	public:
		gemm_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::length(stb)),
		N(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::length(str)),
		common(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::length(sta)),
		lda(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::step(sta)),
		ldb(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::step(stb)),
		ldc(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemm<T>("T", "T", &M, &N, &common, &alpha, b, &ldb, a, &lda, &beta, r, &ldc);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemm<T>("T", "T", &M, &N, &common, &alpha, b, &ldb, a, &lda, &type_constants<T>::one, r, &ldc);
		}
		void add0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
				axpy<T>(&M,&type_constants<T>::one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
		}
		void sub0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
			{
				axpy<T>(&M,&type_constants<T>::m_one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
				scal<T>(&M,&beta,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
			}
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB, typename AB>
	struct gemm_loop<T, STR, STA, STB, RA, RB, AB, 1>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER common;
		BLAS_INTEGER lda;
		BLAS_INTEGER ldb;
		BLAS_INTEGER ldc;
	public:
		gemm_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::length(str)),
		N(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::length(stb)),
		common(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::length(sta)),
		lda(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::step(sta)),
		ldb(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::step(stb)),
		ldc(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemm<T>("N", "N", &M, &N, &common, &alpha, a, &lda, b, &ldb, &beta, r, &ldc);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemm<T>("N", "N", &M, &N, &common, &alpha, a, &lda, b, &ldb, &type_constants<T>::one, r, &ldc);
		}
		void add0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
				axpy<T>(&M,&type_constants<T>::one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
		}
		void sub0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
			{
				axpy<T>(&M,&type_constants<T>::m_one,r0+i*ldc,&type_constants<BLAS_INTEGER>::one,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
				scal<T>(&M,&beta,ri+i*ldc,&type_constants<BLAS_INTEGER>::one);
			}
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename AB, int VARIANT>
	struct gemvA_loop;

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename AB>
	struct gemvA_loop<T, STR, STA, STB, RA, AB, 0>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER LDA;
		BLAS_INTEGER INCX;
		BLAS_INTEGER INCY;
	public:
		gemvA_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::length(sta)),
		N(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::length(sta)),
		LDA(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::step(sta)),
		INCX(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::step(stb)),
		INCY(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemv<T>("N", &M, &N, &alpha, a, &LDA, b, &INCX, &beta, r, &INCY);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemv<T>("N", &M, &N, &alpha, a, &LDA, b, &INCX, &type_constants<T>::one, r, &INCY);
		}
		void add0(const T *r0, T *ri)
		{
			axpy<T>(&M,&type_constants<T>::one,r0,&INCY,ri,&INCY);
		}
		void sub0(const T *r0, T *ri)
		{
			axpy<T>(&M,&type_constants<T>::m_one,r0,&INCY,ri,&INCY);
			scal<T>(&M,&beta,ri,&INCY);
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename AB>
	struct gemvA_loop<T, STR, STA, STB, RA, AB, 2>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER LDA;
		BLAS_INTEGER INCX;
		BLAS_INTEGER INCY;
	public:
		gemvA_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::length(sta)),
		N(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::length(sta)),
		LDA(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::step(sta)),
		INCX(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::step(stb)),
		INCY(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemv<T>("T", &M, &N, &alpha, a, &LDA, b, &INCX, &beta, r, &INCY);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemv<T>("T", &M, &N, &alpha, a, &LDA, b, &INCX, &type_constants<T>::one, r, &INCY);
		}
		void add0(const T *r0, T *ri)
		{
			axpy<T>(&M,&type_constants<T>::one,r0,&INCY,ri,&INCY);
		}
		void sub0(const T *r0, T *ri)
		{
			axpy<T>(&M,&type_constants<T>::m_one,r0,&INCY,ri,&INCY);
			scal<T>(&M,&beta,ri,&INCY);
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RB, typename AB, int VARIANT>
	struct gemvB_loop;

	template <typename T, typename STR, typename STA, typename STB, typename RB, typename AB>
	struct gemvB_loop<T, STR, STA, STB, RB, AB, 0>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER LDB;
		BLAS_INTEGER INCX;
		BLAS_INTEGER INCY;
	public:
		gemvB_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::length(stb)),
		N(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::length(stb)),
		LDB(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::step(stb)),
		INCX(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::step(sta)),
		INCY(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemv<T>("T", &M, &N, &alpha, b, &LDB, a, &INCX, &beta, r, &INCY);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemv<T>("T", &M, &N, &alpha, b, &LDB, a, &INCX, &type_constants<T>::one, r, &INCY);
		}
		void add0(const T *r0, T *ri)
		{
			axpy<T>(&M,&type_constants<T>::one,r0,&INCY,ri,&INCY);
		}
		void sub0(const T *r0, T *ri)
		{
			axpy<T>(&M,&type_constants<T>::m_one,r0,&INCY,ri,&INCY);
			scal<T>(&M,&beta,ri,&INCY);
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RB, typename AB>
	struct gemvB_loop<T, STR, STA, STB, RB, AB, 4>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER M;
		BLAS_INTEGER N;
		BLAS_INTEGER LDB;
		BLAS_INTEGER INCX;
		BLAS_INTEGER INCY;
	public:
		gemvB_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::length(stb)),
		N(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::length(stb)),
		LDB(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::step(stb)),
		INCX(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::step(sta)),
		INCY(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			gemv<T>("N", &M, &N, &alpha, b, &LDB, a, &INCX, &beta, r, &INCY);
		}
		void add(T *r, const T *a, const T *b)
		{
			gemv<T>("N", &M, &N, &alpha, b, &LDB, a, &INCX, &type_constants<T>::one, r, &INCY);
		}
		void add0(const T *r0, T *ri)
		{
			axpy<T>(&M,&type_constants<T>::one,r0,&INCY,ri,&INCY);
		}
		void sub0(const T *r0, T *ri)
		{
			axpy<T>(&M,&type_constants<T>::m_one,r0,&INCY,ri,&INCY);
			scal<T>(&M,&beta,ri,&INCY);
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename AB>
	struct gem_dot_loop
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER qty;
		BLAS_INTEGER incA;
		BLAS_INTEGER incB;
	public:
		gem_dot_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		qty(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::length(sta)),
		incA(vd_iterator_getter<STA, AB, 1>::template iterator_type<T>::step(sta)),
		incB(vd_iterator_getter<STB, AB, 2>::template iterator_type<T>::step(stb))
		{}
		void run(T *r, const T *a, const T *b)
		{
			*r*=beta;
			*r+=alpha*dot<T>(&qty,a,&incA,b,&incB);
		}
		void add(T *r, const T *a, const T *b)
		{
			*r+=alpha*dot<T>(&qty,a,&incA,b,&incB);
		}
		void add0(const T *r0, T *ri)
		{
			*ri+=*r0;
		}
		void sub0(const T *r0, T *ri)
		{
			*ri-=*r0;
			*ri*=beta;
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB, int VARIANT>
	struct ger_loop;

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB>
	struct ger_loop<T, STR, STA, STB, RA, RB, 0>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER	M;
		BLAS_INTEGER	N;
		BLAS_INTEGER	INCA;
		BLAS_INTEGER INCB;
		BLAS_INTEGER	LDR;
	public:
		ger_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::length(stb)),
		N(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::length(str)),
		INCA(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::step(sta)),
		INCB(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::step(stb)),
		LDR(vd_iterator_getter<STR, RA, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			if (beta==1.0)
				ger(&M, &N, &alpha, b, &INCB, a, &INCA, r, &LDR);
			else
			{
				BLAS_INTEGER one=1;
				for (BLAS_INTEGER i=0;i<N;i++)
					scal<T>(&M, &beta, r+i*LDR, &one);
				ger(&M, &N, &alpha, b, &INCB, a, &INCA, r, &LDR);
			}
//			gemm<T>(&transB, &transA, &M, &N, &common, &alpha, b, &ldb, a, &lda, &beta, r, &ldc);
		}
		void add(T *r, const T *a, const T *b)
		{
			ger(&M, &N, &alpha, b, &INCB, a, &INCA, r, &LDR);
		}
		void add0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
				axpy<T>(&M,&type_constants<T>::one,r0+i*LDR,&type_constants<BLAS_INTEGER>::one,ri+i*LDR,&type_constants<BLAS_INTEGER>::one);
		}
		void sub0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
			{
				axpy<T>(&M,&type_constants<T>::m_one,r0+i*LDR,&type_constants<BLAS_INTEGER>::one,ri+i*LDR,&type_constants<BLAS_INTEGER>::one);
				scal<T>(&M,&beta,ri+i*LDR,&type_constants<BLAS_INTEGER>::one);
			}
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename RA, typename RB>
	struct ger_loop<T, STR, STA, STB, RA, RB, 1>
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER	M;
		BLAS_INTEGER	N;
		BLAS_INTEGER	INCA;
		BLAS_INTEGER INCB;
		BLAS_INTEGER	LDR;
	public:
		ger_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		M(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::length(sta)),
		N(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::length(str)),
		INCA(vd_iterator_getter<STA, RA, 1>::template iterator_type<T>::step(sta)),
		INCB(vd_iterator_getter<STB, RB, 2>::template iterator_type<T>::step(stb)),
		LDR(vd_iterator_getter<STR, RB, 0>::template iterator_type<T>::step(str))
		{}
		void run(T *r, const T *a, const T *b)
		{
			if (beta==1.0)
				ger(&M, &N, &alpha, a, &INCA, b, &INCB, r, &LDR);
			else
			{
				BLAS_INTEGER one=1;
				for (BLAS_INTEGER i=0;i<N;i++)
					scal<T>(&M, &beta, r+i*LDR, &one);
				ger(&M, &N, &alpha, a, &INCA, b, &INCB, r, &LDR);
			}
		}
		void add(T *r, const T *a, const T *b)
		{
			ger(&M, &N, &alpha, a, &INCA, b, &INCB, r, &LDR);
		}
		void add0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
				axpy<T>(&M,&type_constants<T>::one,r0+i*LDR,&type_constants<BLAS_INTEGER>::one,ri+i*LDR,&type_constants<BLAS_INTEGER>::one);
		}
		void sub0(const T *r0, T *ri)
		{
			for (BLAS_INTEGER i=0;i<N;i++)
			{
				axpy<T>(&M,&type_constants<T>::m_one,r0+i*LDR,&type_constants<BLAS_INTEGER>::one,ri+i*LDR,&type_constants<BLAS_INTEGER>::one);
				scal<T>(&M,&beta,ri+i*LDR,&type_constants<BLAS_INTEGER>::one);
			}
		}
	};

	template <typename T, typename STR, typename STA, typename STB>
	struct gem_common_loop
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	public:
		gem_common_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta)
		{}
		void run(T *r, const T *a, const T *b)
		{
			*r*=beta;
			*r+=*a*(*b)*alpha;
		}
		void add(T *r, const T *a, const T *b)
		{
			*r+=*a*(*b)*alpha;
		}
		void add0(const T *r0, T *ri)
		{
			*ri+=*r0;
		}
		void sub0(const T *r0, T *ri)
		{
			*ri-=*r0;
			*ri*=beta;
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename COMMON>
	struct gem_sbmv_loop
	{
		const STR& str;
		const STA& sta;
		const STB& stb;
		T alpha;
		T beta;
	private:
		BLAS_INTEGER N;
		BLAS_INTEGER INCR;
		BLAS_INTEGER INCA;
		BLAS_INTEGER INCB;
		using rittype=typename vd_iterator_getter<STR, COMMON, 0>::template iterator_type<T>;
		using aittype=typename vd_iterator_getter<STA, COMMON, 1>::template iterator_type<T>;
		using bittype=typename vd_iterator_getter<STB, COMMON, 2>::template iterator_type<T>;
	public:
		gem_sbmv_loop(const STR& str, const STA& sta, const STB& stb, T alpha, T beta): str(str), sta(sta), stb(stb),
		alpha(alpha),beta(beta),
		N(rittype::length(str)),
		INCR(rittype::step(str)),
		INCA(aittype::step(sta)),
		INCB(bittype::step(stb))
		{}
		void run(T *r, const T *a, const T *b)
		{
			sbmv<T>("U", &N, &type_constants<BLAS_INTEGER>::zero, &alpha, a, &INCA, b, &INCB, &beta, r, &INCR);
		}
		void add(T *r, const T *a, const T *b)
		{
			sbmv<T>("U", &N, &type_constants<BLAS_INTEGER>::zero, &alpha, a, &INCA, b, &INCB, &type_constants<T>::one, r, &INCR);
		}
		void add0(const T *r0, T *ri)
		{
			axpy<T>(&N, &type_constants<T>::one, r0, &INCR, ri, &INCR);
		}
		void sub0(const T *r0, T *ri)
		{
			axpy<T>(&N, &type_constants<T>::m_one, r0, &INCR, ri, &INCR);
			scal<T>(&N, &beta, ri, &INCR);
		}
	};

	template <typename T, typename STR, typename STA, typename STB, typename SPREAD, typename SUMA, typename RA, typename SUMB, typename RB, typename AB, typename COMMON, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_sbmv_test;
	template <typename T, typename STR, typename STA, typename STB, typename ... SPREAD, typename ... SUMA, typename ... RA, typename ... SUMB, typename ... RB, typename ... AB, typename ... COMMON, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_sbmv_test<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<COMMON...>, R_CONT, A_CONT, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<COMMON...>, gem_common_loop<T, STR, STA, STB> >
	{
	};
	template <typename T, typename STR, typename STA, typename STB, typename ... SPREAD, typename ... SUMA, typename ... RA, typename ... SUMB, typename ... RB, typename ... AB, int V_TYPE, int C_MASK, int NEXT_VALENCE, size_t ... POS, typename ... COMMON0, typename ... COMMON, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_sbmv_test<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<type_sequence<valence_data<V_TYPE, 7, C_MASK, 7, NEXT_VALENCE, POS...>, COMMON0...>, COMMON...>, R_CONT, A_CONT, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<COMMON...>, gem_sbmv_loop<T, STR, STA, STB, type_sequence<valence_data<V_TYPE, 7, C_MASK, 7, NEXT_VALENCE, POS...>, COMMON0...> > >
	{

	};

	template <typename T, typename STR, typename STA, typename STB, typename VD, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_mask;

//	DONE GEMM
	template <typename T, typename STR, typename STA, typename STB, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, int RA_V_TYPE, int RA_C_MASK, int RA_NEXT_VALENCE, size_t ... RA_POS, typename ... RA0, typename ... RA, typename ... SUMB, int RB_V_TYPE, int RB_C_MASK, int RB_NEXT_VALENCE, size_t ... RB_POS, typename ... RB0, typename ... RB, int AB_V_TYPE, int AB_C_MASK, int AB_NEXT_VALENCE, size_t ... AB_POS, typename ... AB0, typename ... AB, typename ... COMMON>
	struct gem_mask<T, STR, STA, STB, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>,RA...>, type_sequence<SUMB...>, type_sequence<type_sequence<valence_data<RB_V_TYPE, 5, RB_C_MASK, 5, RB_NEXT_VALENCE, RB_POS...>, RB0...>, RB...>, type_sequence<type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...>, AB...>, type_sequence<COMMON...> >, true, true, true>:
		public gem_create_general_loop<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<COMMON...>, gemm_loop<T, STR, STA, STB, type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>, type_sequence<valence_data<RB_V_TYPE, 5, RB_C_MASK, 5, RB_NEXT_VALENCE, RB_POS...>, RB0...>, type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...>, (RA_C_MASK & 1) + (AB_C_MASK & 2) + (RB_C_MASK & 4) > >
	{
	};
//	DONE GEMV A(B)
	template <typename T, typename STR, typename STA, typename STB, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, int RA_V_TYPE, int RA_C_MASK, int RA_NEXT_VALENCE, size_t ... RA_POS, typename ... RA0, typename ... RA, typename ... SUMB, int RB_V_TYPE, int RB_C_MASK, int RB_NEXT_VALENCE, size_t ... RB_POS, typename ... RB0, typename ... RB, int AB_V_TYPE, int AB_C_MASK, int AB_NEXT_VALENCE, size_t ... AB_POS, typename ... AB0, typename ... AB, typename ... COMMON, bool R_CONT, bool B_CONT>
	struct gem_mask<T, STR, STA, STB, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>,RA...>, type_sequence<SUMB...>, type_sequence<type_sequence<valence_data<RB_V_TYPE, 5, RB_C_MASK, 5, RB_NEXT_VALENCE, RB_POS...>, RB0...>, RB...>, type_sequence<type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...>, AB...>, type_sequence<COMMON...> >, R_CONT, true, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<type_sequence<valence_data<RB_V_TYPE, 5, RB_C_MASK, 5, RB_NEXT_VALENCE, RB_POS...>, RB0...>, RB...>, type_sequence<AB...>, type_sequence<COMMON...>, gemvA_loop<T, STR, STA, STB, type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>, type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...>, (AB_C_MASK & 2) > >
	{
	};
//  DONE GEMV A
	template <typename T, typename STR, typename STA, typename STB, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, int RA_V_TYPE, int RA_C_MASK, int RA_NEXT_VALENCE, size_t ... RA_POS, typename ... RA0, typename ... RA, typename ... SUMB, typename ... RB, int AB_V_TYPE, int AB_C_MASK, int AB_NEXT_VALENCE, size_t ... AB_POS, typename ... AB0, typename ... AB, typename ... COMMON, bool R_CONT, bool B_CONT>
	struct gem_mask<T, STR, STA, STB, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>,RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...>, AB...>, type_sequence<COMMON...> >, R_CONT, true, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<COMMON...>, gemvA_loop<T, STR, STA, STB, type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>, type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...>, (AB_C_MASK & 2) > >
	{
	};

//	DONE GEMV B
	template <typename T, typename STR, typename STA, typename STB, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, typename ... RA, typename ... SUMB, int RB_V_TYPE, int RB_C_MASK, int RB_NEXT_VALENCE, size_t ... RB_POS, typename ... RB0, typename ... RB, int AB_V_TYPE, int AB_C_MASK, int AB_NEXT_VALENCE, size_t ... AB_POS, typename ... AB0, typename ... AB, typename ... COMMON, bool R_CONT, bool A_CONT>
	struct gem_mask<T, STR, STA, STB, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<type_sequence<valence_data<RB_V_TYPE, 5, RB_C_MASK, 5, RB_NEXT_VALENCE, RB_POS...>, RB0...>, RB...>, type_sequence<type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...>, AB...>, type_sequence<COMMON...> >, R_CONT, A_CONT, true>:
		public gem_create_general_loop<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<COMMON...>, gemvB_loop<T, STR, STA, STB, type_sequence<valence_data<RB_V_TYPE, 5, RB_C_MASK, 5, RB_NEXT_VALENCE, RB_POS...>, RB0...>, type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...>, (RB_C_MASK & 4) > >
	{
	};

//	DONE DOT
	template <typename T, typename STR, typename STA, typename STB, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, int RA_V_TYPE, int RA_C_MASK, int RA_NEXT_VALENCE, size_t ... RA_POS, typename ... RA0, typename ... RA, typename ... SUMB, int RB_V_TYPE, int RB_C_MASK, int RB_NEXT_VALENCE, size_t ... RB_POS, typename ... RB0, typename ... RB, int AB_V_TYPE, int AB_C_MASK, int AB_NEXT_VALENCE, size_t ... AB_POS, typename ... AB0, typename ... AB, typename ... COMMON, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_mask<T, STR, STA, STB, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>,RA...>, type_sequence<SUMB...>, type_sequence<type_sequence<valence_data<RB_V_TYPE, 5, RB_C_MASK, 5, RB_NEXT_VALENCE, RB_POS...>, RB0...>, RB...>, type_sequence<type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...>, AB...>, type_sequence<COMMON...> >, R_CONT, A_CONT, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>,RA...>, type_sequence<SUMB...>, type_sequence<type_sequence<valence_data<RB_V_TYPE, 5, RB_C_MASK, 5, RB_NEXT_VALENCE, RB_POS...>, RB0...>, RB...>, type_sequence<AB...>, type_sequence<COMMON...>, gem_dot_loop<T, STR, STA, STB, type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...> > >
	{
	};

//	DONE DOT
	template <typename T, typename STR, typename STA, typename STB, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, typename ... RA, typename ... SUMB, typename ... RB, int AB_V_TYPE, int AB_C_MASK, int AB_NEXT_VALENCE, size_t ... AB_POS, typename ... AB0, typename ... AB, typename ... COMMON, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_mask<T, STR, STA, STB, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...>, AB...>, type_sequence<COMMON...> >, R_CONT, A_CONT, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<COMMON...>, gem_dot_loop<T, STR, STA, STB, type_sequence<valence_data<AB_V_TYPE, 6, AB_C_MASK, 6, AB_NEXT_VALENCE, AB_POS...>, AB0...> > >
	{
	};

//	DONE GER
	template <typename T, typename STR, typename STA, typename STB, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, int RA_V_TYPE, int RA_C_MASK, int RA_NEXT_VALENCE, size_t ... RA_POS, typename ... RA0, typename ... RA, typename ... SUMB, int RB_V_TYPE, int RB_C_MASK, int RB_NEXT_VALENCE, size_t ... RB_POS, typename ... RB0, typename ... RB, typename ... AB, typename ... COMMON, bool A_CONT, bool B_CONT>
	struct gem_mask<T, STR, STA, STB, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>,RA...>, type_sequence<SUMB...>, type_sequence<type_sequence<valence_data<RB_V_TYPE, 5, RB_C_MASK, 5, RB_NEXT_VALENCE, RB_POS...>, RB0...>, RB...>, type_sequence<AB...>, type_sequence<COMMON...> >, true, A_CONT, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<COMMON...>, ger_loop<T, STR, STA, STB, type_sequence<valence_data<RA_V_TYPE, 3, RA_C_MASK, 3, RA_NEXT_VALENCE, RA_POS...>, RA0...>, type_sequence<valence_data<RB_V_TYPE, 5, RB_C_MASK, 5, RB_NEXT_VALENCE, RB_POS...>, RB0...>, (RA_C_MASK & 1) > >
	{
	};

//	DONE COMMON OR SBMV
	template <typename T, typename STR, typename STA, typename STB, typename CONT_MASK, typename ... SPREAD, typename ... SUMA, typename ... RA, typename ... SUMB, typename ... RB, typename ... AB, typename ... COMMON, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_mask<T, STR, STA, STB, type_sequence<CONT_MASK, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<COMMON...> >, R_CONT, A_CONT, B_CONT>:
		public gem_sbmv_test<T, STR, STA, STB, type_sequence<SPREAD...>, type_sequence<SUMA...>, type_sequence<RA...>, type_sequence<SUMB...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<COMMON...>, R_CONT, A_CONT, B_CONT>
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename VD>
	struct gem_runner;

	template <typename T, typename STR, typename STA, typename STB, int CONT_MASK, typename ... VD>
	struct gem_runner<T, STR, STA, STB, type_sequence<std::integral_constant<int, CONT_MASK>, VD...> >:
		public gem_mask<T, STR, STA, STB, type_sequence<std::integral_constant<int, CONT_MASK>, VD...>, (CONT_MASK & 1)==1, (CONT_MASK & 2)==2, (CONT_MASK & 4)==4>
	{
	};

};
#endif /* TGEM_H_ */
