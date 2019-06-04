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

namespace tpp
{

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename HEAD, typename BASE_LOOP>
	struct gem_sumb_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct gem_sumb_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STB, SI, 2, V_TYPE>::template iterator_type<T> itype;
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename HEAD, typename BASE_LOOP>
	struct gem_ab_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct gem_ab_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STA, SI, 1, V_TYPE>::template iterator_type<T> aitype;
		typedef typename iterator_getter_by_src_v_type<STB, SI, 2, V_TYPE>::template iterator_type<T> bitype;
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename HEAD, typename BASE_LOOP>
	struct gem_rb_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct gem_rb_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T> ritype;
		typedef typename iterator_getter_by_src_v_type<STB, SI, 2, V_TYPE>::template iterator_type<T> bitype;
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename HEAD, typename BASE_LOOP>
	struct gem_suma_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct gem_suma_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STA, SI, 1, V_TYPE>::template iterator_type<T> itype;
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename HEAD, typename BASE_LOOP>
	struct gem_ra_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct gem_ra_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T> ritype;
		typedef typename iterator_getter_by_src_v_type<STA, SI, 1, V_TYPE>::template iterator_type<T> aitype;
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename HEAD, typename BASE_LOOP>
	struct gem_com_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct gem_com_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T> ritype;
		typedef typename iterator_getter_by_src_v_type<STA, SI, 1, V_TYPE>::template iterator_type<T> aitype;
		typedef typename iterator_getter_by_src_v_type<STB, SI, 2, V_TYPE>::template iterator_type<T> bitype;
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename HEAD, bool IS_FIRST, typename BASE_LOOP>
	struct gem_spread_loop;

// LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename STB, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct gem_spread_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, false, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T> itype;
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
	template <typename T, typename STR, typename STA, typename STB, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct gem_spread_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, true, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T> itype;
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
//		void add(T *r, const T *a, const T *b)
//		{
//			itype it(this->str, r);
//			for (;it.not_end(r);it.move_one(r))
//				BASE_LOOP::add(r,a,b);
//		}
	};

// NOT LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename STB, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename HEAD, typename BASE_BASE_LOOP>
	struct gem_spread_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, false, gem_spread_loop<T, STR, STA, STB, SI, HEAD, false, BASE_BASE_LOOP> >: public gem_spread_loop<T, STR, STA, STB, SI, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef gem_spread_loop<T, STR, STA, STB, SI, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
	public:
		typedef typename iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T> itype;
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
	template <typename T, typename STR, typename STA, typename STB, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename HEAD, typename BASE_BASE_LOOP>
	struct gem_spread_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, true, gem_spread_loop<T, STR, STA, STB, SI, HEAD, false, BASE_BASE_LOOP> >: public gem_spread_loop<T, STR, STA, STB, SI, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef gem_spread_loop<T, STR, STA, STB, SI, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
	public:
		typedef typename iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T> itype;
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename SPREAD, typename COMMON, typename RA, typename RB, typename AB, typename SUMA, typename SUMB, typename BASE_LOOP>
	struct gem_create_general_loop
	{
		using type=BASE_LOOP;
	};

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename SPREAD, typename COMMON, typename RA, typename RB, typename AB, typename SUMA, typename HEAD, typename ... SUMB, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, COMMON, RA, RB, AB, SUMA, type_sequence<HEAD, SUMB...>, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, COMMON, RA, RB, AB, SUMA, type_sequence<SUMB...>, gem_sumb_loop<T, STR, STA, STB, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename SPREAD, typename COMMON, typename RA, typename RB, typename HEAD, typename ... AB, typename SUMA, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, COMMON, RA, RB, type_sequence<HEAD, AB...>, SUMA, type_sequence<>, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, COMMON, RA, RB, type_sequence<AB...>, SUMA, type_sequence<>, gem_ab_loop<T, STR, STA, STB, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename SPREAD, typename COMMON, typename RA, typename HEAD, typename ... RB, typename SUMA, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, COMMON, RA, type_sequence<HEAD, RB...>, type_sequence<>, SUMA, type_sequence<>, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, COMMON, RA, type_sequence<RB...>, type_sequence<>, SUMA, type_sequence<>, gem_rb_loop<T, STR, STA, STB, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename SPREAD, typename COMMON, typename RA, typename HEAD, typename ... SUMA, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, COMMON, RA, type_sequence<>, type_sequence<>, type_sequence<HEAD, SUMA...>, type_sequence<>, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, COMMON, RA, type_sequence<>, type_sequence<>, type_sequence<SUMA...>, type_sequence<>, gem_suma_loop<T, STR, STA, STB, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename SPREAD, typename COMMON, typename HEAD, typename ... RA, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, COMMON, type_sequence<HEAD, RA...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, COMMON, type_sequence<RA...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, gem_ra_loop<T, STR, STA, STB, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename SPREAD, typename HEAD, typename ... COMMON, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, type_sequence<HEAD, COMMON...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SI, SPREAD, type_sequence<COMMON...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, gem_com_loop<T, STR, STA, STB, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename HEAD, typename ... SPREAD, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<HEAD, SPREAD...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, gem_spread_loop<T, STR, STA, STB, SI, HEAD, false, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename HEAD, typename BASE_LOOP>
	struct gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<HEAD>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, gem_spread_loop<T, STR, STA, STB, SI, HEAD, true, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename RA, typename RB, typename AB, int VARIANT>
	struct gemm_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemm_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 6>
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
		M(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::length(str)),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::length(str)),
		common(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::length(stb)),
		lda(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::step(sta)),
		ldb(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::step(stb)),
		ldc(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemm_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 7>
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
		M(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::length(sta)),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::length(str)),
		common(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::length(stb)),
		lda(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::step(sta)),
		ldb(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::step(stb)),
		ldc(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemm_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 4>
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
		M(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::length(str)),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::length(str)),
		common(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::length(stb)),
		lda(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::step(sta)),
		ldb(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::step(stb)),
		ldc(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemm_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 5>
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
		M(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::length(str)),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::length(str)),
		common(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::length(stb)),
		lda(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::step(sta)),
		ldb(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::step(stb)),
		ldc(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemm_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 2>
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
		M(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::length(stb)),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::length(str)),
		common(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::length(sta)),
		lda(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::step(sta)),
		ldb(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::step(stb)),
		ldc(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemm_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 3>
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
		M(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::length(sta)),
		N(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::length(stb)),
		common(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::length(sta)),
		lda(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::step(sta)),
		ldb(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::step(stb)),
		ldc(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemm_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 0>
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
		M(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::length(stb)),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::length(str)),
		common(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::length(sta)),
		lda(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::step(sta)),
		ldb(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::step(stb)),
		ldc(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemm_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 1>
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
		M(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::length(str)),
		N(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::length(stb)),
		common(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::length(sta)),
		lda(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::step(sta)),
		ldb(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::step(stb)),
		ldc(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename RA, typename AB, int VARIANT>
	struct gemvA_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemvA_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 0>
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
		M(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::length(sta)),
		N(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::length(sta)),
		LDA(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::step(sta)),
		INCX(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::step(stb)),
		INCY(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemvA_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 2>
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
		M(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::length(sta)),
		N(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::length(sta)),
		LDA(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::step(sta)),
		INCX(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::step(stb)),
		INCY(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename RB, typename AB, int VARIANT>
	struct gemvB_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemvB_loop<T, STR, STA, STB, SI, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 0>
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
		M(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::length(stb)),
		N(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::length(stb)),
		LDB(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::step(stb)),
		INCX(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::step(sta)),
		INCY(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gemvB_loop<T, STR, STA, STB, SI, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, 4>
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
		M(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::length(stb)),
		N(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::length(stb)),
		LDB(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::step(stb)),
		INCX(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::step(sta)),
		INCY(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename AB>
	struct gem_dot_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER>
	struct gem_dot_loop<T, STR, STA, STB, SI, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER> >
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
		qty(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::length(sta)),
		incA(iterator_getter_by_src_v_type<STA, SI, 1, AB_V_TYPE>::template iterator_type<T>::step(sta)),
		incB(iterator_getter_by_src_v_type<STB, SI, 2, AB_V_TYPE>::template iterator_type<T>::step(stb))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename RA, typename RB, int VARIANT>
	struct ger_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER>
	struct ger_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, 0>
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
		M(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::length(stb)),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::length(str)),
		INCA(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::step(sta)),
		INCB(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::step(stb)),
		LDR(iterator_getter_by_src_v_type<STR, SI, 0, RA_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER>
	struct ger_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, 1>
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
		M(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::length(sta)),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::length(str)),
		INCA(iterator_getter_by_src_v_type<STA, SI, 1, RA_V_TYPE>::template iterator_type<T>::step(sta)),
		INCB(iterator_getter_by_src_v_type<STB, SI, 2, RB_V_TYPE>::template iterator_type<T>::step(stb)),
		LDR(iterator_getter_by_src_v_type<STR, SI, 0, RB_V_TYPE>::template iterator_type<T>::step(str))
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

	template <typename T, typename STR, typename STA, typename STB, typename SI>
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

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename COMMON>
	struct gem_sbmv_loop;

	template <typename T, typename STR, typename STA, typename STB, typename SI, int V_TYPE, int C_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER>
	struct gem_sbmv_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, 7, C_MASK, 7, IS_JOINABLE, JOIN_WITH, ORDER> >
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
		using rittype=typename iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T>;
		using aittype=typename iterator_getter_by_src_v_type<STA, SI, 1, V_TYPE>::template iterator_type<T>;
		using bittype=typename iterator_getter_by_src_v_type<STB, SI, 2, V_TYPE>::template iterator_type<T>;
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



	template <typename ELEMENT, typename VI, typename RES>
	struct insert_mask;

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, int H_V_TYPE, int H_MASK, int H_C_MASK, int H_V_MASK, bool H_IS_JOINABLE, int H_JOIN_WITH, size_t H_ORDER, typename ... VI, typename ... RES>
	struct insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>, type_sequence<RES...> >:
		public std::conditional<
			(MASK==V_MASK && H_MASK!=H_V_MASK)?true:
				(MASK!=V_MASK && H_MASK==H_V_MASK)?false:
					(C_MASK > H_C_MASK)?true:
						(C_MASK < H_C_MASK)?false:
							(ORDER<H_ORDER),
			type_sequence<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>,
			insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<VI...>, type_sequence<RES..., valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER> > >
		>::type
	{
	};

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... RES>
	struct insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<>, type_sequence<RES...> >:
		public type_sequence<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB, typename VISI, typename SPREAD, typename COMMON, typename RA, typename RB, typename AB, typename SUMA, typename SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_sbmv_test;

	template <typename T, typename STR, typename STA, typename STB, typename SI, typename ... SPREAD, typename ... COMMON, typename ... RA, typename ... RB, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_sbmv_test<T, STR, STA, STB, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, gem_common_loop<T, STR, STA, STB, SI> >
	{
	};
	template <typename T, typename STR, typename STA, typename STB, typename SI, typename ... SPREAD, int V_TYPE, int C_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... COMMON, typename ... RA, typename ... RB, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_sbmv_test<T, STR, STA, STB, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<valence_info<V_TYPE, 7, C_MASK, 7, IS_JOINABLE, JOIN_WITH, ORDER>, COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, gem_sbmv_loop<T, STR, STA, STB, SI, valence_info<V_TYPE, 7, C_MASK, 7, IS_JOINABLE, JOIN_WITH, ORDER> > >
	{

	};

	template <typename T, typename STR, typename STA, typename STB, typename VISI, typename SPREAD, typename COMMON, typename RA, typename RB, typename AB, typename SUMA, typename SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_by_mask;

//  SPREAD
	template <typename T, typename STR, typename STA, typename STB, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... COMMON, typename ... RA, typename ... RB, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<VI...>, SI>, typename insert_mask<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<SPREAD...>, type_sequence<> >::type, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>
	{
	};
//	COMMON
	template <typename T, typename STR, typename STA, typename STB, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... COMMON, typename ... RA, typename ... RB, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<valence_info<V_TYPE, 7, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<VI...>, SI>, type_sequence<SPREAD...>, typename insert_mask<valence_info<V_TYPE, 7, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<COMMON...>, type_sequence<> >::type, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>
	{
	};
//  RA
	template <typename T, typename STR, typename STA, typename STB, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... COMMON, typename ... RA, typename ... RB, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, typename insert_mask<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<RA...>, type_sequence<> >::type, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT || ((V_MASK==3 && (C_MASK & 1))), A_CONT || ((V_MASK==3 && (C_MASK & 2))), B_CONT>
	{
	};
//  RB
	template <typename T, typename STR, typename STA, typename STB, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... COMMON, typename ... RA, typename ... RB, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<valence_info<V_TYPE, 5, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, typename insert_mask<valence_info<V_TYPE, 5, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<RB...>, type_sequence<> >::type, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT || ((V_MASK==5 && (C_MASK & 1))), A_CONT, B_CONT || ((V_MASK==5 && (C_MASK & 4)))>
	{
	};
//  AB
	template <typename T, typename STR, typename STA, typename STB, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... COMMON, typename ... RA, typename ... RB, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<valence_info<V_TYPE, 6, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, typename insert_mask<valence_info<V_TYPE, 6, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<AB...>, type_sequence<> >::type, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT || ((V_MASK==6 && (C_MASK & 2))), B_CONT || ((V_MASK==6 && (C_MASK & 4)))>
	{
	};
//	SUMA
	template <typename T, typename STR, typename STA, typename STB, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... COMMON, typename ... RA, typename ... RB, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, typename insert_mask<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<SUMA...>, type_sequence<> >::type, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>
	{
	};
//	SUMB
	template <typename T, typename STR, typename STA, typename STB, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... COMMON, typename ... RA, typename ... RB, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<valence_info<V_TYPE, 4, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<VI...>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, typename insert_mask<valence_info<V_TYPE, 4, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<SUMB...>, type_sequence<> >::type, R_CONT, A_CONT, B_CONT>
	{
	};
//	DONE GEMM
	template <typename T, typename STR, typename STA, typename STB, typename SI, typename ... SPREAD, typename ... COMMON, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, typename ... RA, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, typename ... RB, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER, typename ... AB, typename ... SUMA, typename ... SUMB>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, true, true, true>:
		public gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, gemm_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, (RA_C_MASK & 1) + (AB_C_MASK & 2) + (RB_C_MASK & 4) > >
//		public gemm_final<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, (RA_C_MASK & 1) + (AB_C_MASK & 2) + (RB_C_MASK & 4) >
	{
	};
//	DONE GEMV A(B)
	template <typename T, typename STR, typename STA, typename STB, typename SI, typename ... SPREAD, typename ... COMMON, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, typename ... RA, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, typename ... RB, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, true, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, gemvA_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, (AB_C_MASK & 2) > >
//		public gemvA<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...> >
	{
	};
//  DONE GEMV A
	template <typename T, typename STR, typename STA, typename STB, typename SI, typename ... SPREAD, typename ... COMMON, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, typename ... RA, typename ... RB, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, true, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, gemvA_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, (AB_C_MASK & 2) > >
//		public gemvA<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...> >
	{
	};

//	DONE GEMV B
	template <typename T, typename STR, typename STA, typename STB, typename SI, typename ... SPREAD, typename ... COMMON, typename ... RA, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, typename ... RB, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, true>:
		public gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, gemvB_loop<T, STR, STA, STB, SI, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, (RB_C_MASK & 4) > >
//		public gemvB<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...> >
	{
	};

//	DONE DOT
	template <typename T, typename STR, typename STA, typename STB, typename SI, typename ... SPREAD, typename ... COMMON, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, typename ... RA, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, typename ... RB, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, gem_dot_loop<T, STR, STA, STB, SI, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER> > >
//		public dot_final<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...> >
	{
	};

//	DONE DOT
	template <typename T, typename STR, typename STA, typename STB, typename SI, typename ... SPREAD, typename ... COMMON, typename ... RA, typename ... RB, int AB_V_TYPE, int AB_C_MASK, bool AB_IS_JOINABLE, int AB_JOIN_WITH, size_t AB_ORDER, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, gem_dot_loop<T, STR, STA, STB, SI, valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER> > >
//		public dot_final<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_sequence<SUMA...>, type_sequence<SUMB...> >
	{
	};

//	DONE GER
	template <typename T, typename STR, typename STA, typename STB, typename SI, typename ... SPREAD, typename ... COMMON, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, typename ... RA, int RB_V_TYPE, int RB_C_MASK, bool RB_IS_JOINABLE, int RB_JOIN_WITH, size_t RB_ORDER, typename ... RB, typename ... AB, typename ... SUMA, typename ... SUMB, bool A_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, true, A_CONT, B_CONT>:
		public gem_create_general_loop<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, ger_loop<T, STR, STA, STB, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, (RA_C_MASK & 1) > >
//		public ger_final<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...> >
	{
	};

//	DONE COMMON OR SBMV
	template <typename T, typename STR, typename STA, typename STB, typename SI, typename ... SPREAD, typename ... COMMON, typename ... RA, typename ... RB, typename ... AB, typename ... SUMA, typename ... SUMB, bool R_CONT, bool A_CONT, bool B_CONT>
	struct gem_by_mask<T, STR, STA, STB, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>:
		public gem_sbmv_test<T, STR, STA, STB, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...>, R_CONT, A_CONT, B_CONT>
//		public common_final<T, STR, STA, STB, SI, type_sequence<SPREAD...>, type_sequence<COMMON...>, type_sequence<RA...>, type_sequence<RB...>, type_sequence<AB...>, type_sequence<SUMA...>, type_sequence<SUMB...> >
	{
	};

	template <typename T, typename STR, typename STA, typename STB>
	struct gem_runner: public gem_by_mask<T, STR, STA, STB, typename make_valence_info_and_join<STR, STA, STB>::type, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, false, false, false>
	{
	};

};
#endif /* TGEM_H_ */
