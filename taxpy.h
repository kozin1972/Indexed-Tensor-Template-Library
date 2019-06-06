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

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename BASE_LOOP>
	struct axpy_suma_loop;

	template <typename T, typename STR, typename STA, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct axpy_suma_loop<T, STR, STA, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STA, SI, 1, V_TYPE>::template iterator_type<T> itype;
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

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename BASE_LOOP>
	struct axpy_ra_loop;

	template <typename T, typename STR, typename STA, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct axpy_ra_loop<T, STR, STA, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T> ritype;
		typedef typename iterator_getter_by_src_v_type<STA, SI, 1, V_TYPE>::template iterator_type<T> aitype;
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

	template <typename T, typename STR, typename STA, typename SI>
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

	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, bool IS_FIRST, typename BASE_LOOP>
	struct axpy_spread_loop;

// LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename BASE_LOOP>
	struct axpy_spread_loop<T, STR, STA, SI, ITTYPE, false, BASE_LOOP>: public BASE_LOOP
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
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename BASE_LOOP>
	struct axpy_spread_loop<T, STR, STA, SI, ITTYPE, true, BASE_LOOP>: public BASE_LOOP
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
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename HEAD, typename BASE_BASE_LOOP>
	struct axpy_spread_loop<T, STR, STA, SI, ITTYPE, false, axpy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> >: public axpy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef axpy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
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
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename HEAD, typename BASE_BASE_LOOP>
	struct axpy_spread_loop<T, STR, STA, SI, ITTYPE, true, axpy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> >: public axpy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef axpy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
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

	template <typename T, typename STR, typename STA, typename SI, typename SPREAD, typename RA, typename SUMA, typename BASE_LOOP>
	struct axpy_create_general_loop
	{
		using type=BASE_LOOP;
	};

	template <typename T, typename STR, typename STA, typename SI, typename SPREAD, typename RA, typename HEAD, typename ... SUMA, typename BASE_LOOP>
	struct axpy_create_general_loop<T, STR, STA, SI, SPREAD, RA, type_pack<HEAD, SUMA...>, BASE_LOOP>:
		public axpy_create_general_loop<T, STR, STA, SI, SPREAD, RA, type_pack<SUMA...>, axpy_suma_loop<T, STR, STA, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename SPREAD, typename HEAD, typename ... RA, typename BASE_LOOP>
	struct axpy_create_general_loop<T, STR, STA, SI, SPREAD, type_pack<HEAD, RA...>, type_pack<>, BASE_LOOP>:
		public axpy_create_general_loop<T, STR, STA, SI, SPREAD, type_pack<RA...>, type_pack<>, axpy_ra_loop<T, STR, STA, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename ... SPREAD, typename BASE_LOOP>
	struct axpy_create_general_loop<T, STR, STA, SI, type_pack<HEAD, SPREAD...>, type_pack<>, type_pack<>, BASE_LOOP>:
		public axpy_create_general_loop<T, STR, STA, SI, type_pack<SPREAD...>, type_pack<>, type_pack<>, axpy_spread_loop<T, STR, STA, SI, typename iterator_getter_by_src_v_type<STR, SI, 0, HEAD::v_type>::template iterator_type<T>, false, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename BASE_LOOP>
	struct axpy_create_general_loop<T, STR, STA, SI, type_pack<HEAD>, type_pack<>, type_pack<>, BASE_LOOP>:
		public axpy_create_general_loop<T, STR, STA, SI, type_pack<>, type_pack<>, type_pack<>, axpy_spread_loop<T, STR, STA, SI, typename iterator_getter_by_src_v_type<STR, SI, 0, HEAD::v_type>::template iterator_type<T>, true, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename RAV>
	struct axpy_loop;

	template <typename T, typename STR, typename STA, typename SI, int RAV_V_TYPE, int RAV_C_MASK, bool RAV_IS_JOINABLE, int RAV_JOIN_WITH, size_t RAV_ORDER>
	struct axpy_loop<T, STR, STA, SI, valence_info<RAV_V_TYPE, 3, RAV_C_MASK, 3, RAV_IS_JOINABLE, RAV_JOIN_WITH, RAV_ORDER> >
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
		axpy_loop(const STR& str, const STA& sta): str(str), sta(sta), alpha(1),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::length(str)),
		incr(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::step(str)),
		inca(iterator_getter_by_src_v_type<STA, SI, 1, RAV_V_TYPE>::template iterator_type<T>::step(sta))
		{}
		axpy_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::length(str)),
		incr(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::step(str)),
		inca(iterator_getter_by_src_v_type<STA, SI, 1, RAV_V_TYPE>::template iterator_type<T>::step(sta))
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

	template <typename T, typename STR, typename STA, typename SI, typename S>
	struct axpy_just_spread_loop;

	template <typename T, typename STR, typename STA, typename SI, int S_V_TYPE, int S_C_MASK, bool S_IS_JOINABLE, int S_JOIN_WITH, size_t S_ORDER>
	struct axpy_just_spread_loop<T, STR, STA, SI, valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER> >
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	private:
		BLAS_INTEGER l;
		BLAS_INTEGER s;
	public:
		using ittype=typename iterator_getter_by_src_v_type<STR, SI, 0, S_V_TYPE>::template iterator_type<T>;
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

	template <typename T, typename STR, typename STA, typename SI, typename S>
	struct axpy_sum_loop;

	template <typename T, typename STR, typename STA, typename SI, int V_TYPE, int C_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER>
	struct axpy_sum_loop<T, STR, STA, SI, valence_info<V_TYPE, 2, C_MASK, 2, IS_JOINABLE, JOIN_WITH, ORDER> >
	{
		const STR& str;
		const STA& sta;
		const T alpha;
	private:
		BLAS_INTEGER l;
		BLAS_INTEGER s;
	public:
		using ittype=typename iterator_getter_by_src_v_type<STA, SI, 1, V_TYPE>::template iterator_type<T>;
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
//	template <typename SI, typename RA, typename RB, typename AB, int VARIANT>
//	struct create_gemm_loop;

//	template <typename T, typename STR, typename STA, typename STB, typename SI, typename SPREAD, typename COMMON, typename RA, typename RB, typename AB, typename SUMA, typename SUMB, int VARIANT>
//	struct gemm_final;
//
//	template <typename T, typename STR, typename STA, typename STB, typename SI, typename SPREAD, typename COMMON, typename RAH, typename ... RA, typename RBH, typename ... RB, typename ABH, typename ... AB, typename SUMA, typename SUMB, int VARIANT>
//	struct gemm_final<T, STR, STA, STB, SI, SPREAD, COMMON, type_pack<RAH, RA...>, type_pack<RBH, RB...>, type_pack<ABH, AB...>, SUMA, SUMB, VARIANT>:
//		public create_general_loop<T, STR, STA, STB, SI, SPREAD, COMMON, type_pack<RA...>, type_pack<RB...>, type_pack<AB...>, SUMA, SUMB, gemm_loop<T, STR, STA, STB, SI, RAH, RBH, ABH, VARIANT> >
//	{
//
//	};



	template <typename ELEMENT, typename VI, typename RES>
	struct axpy_insert_mask;

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, int H_V_TYPE, int H_MASK, int H_C_MASK, int H_V_MASK, bool H_IS_JOINABLE, int H_JOIN_WITH, size_t H_ORDER, typename ... VI, typename ... RES>
	struct axpy_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>, type_pack<RES...> >:
		public std::conditional<
			(MASK==V_MASK && H_MASK!=H_V_MASK)?true:
				(MASK!=V_MASK && H_MASK==H_V_MASK)?false:
					(C_MASK > H_C_MASK)?true:
						(C_MASK < H_C_MASK)?false:
							(ORDER<H_ORDER),
			type_pack<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>,
			axpy_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<VI...>, type_pack<RES..., valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER> > >
		>::type
	{
	};

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... RES>
	struct axpy_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<>, type_pack<RES...> >:
		public type_pack<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER> >
	{
	};

	template <typename T, typename STR, typename STA, typename VISI, typename SPREAD, typename RA, typename SUMA, bool R_CONT, bool A_CONT>
	struct axpy_test_sum;

//	DONE sum
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, typename ... RA, int V_TYPE, int C_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct axpy_test_sum<T, STR, STA, type_pack<type_pack<>, SI>, type_pack<SPREAD...>, type_pack<RA...>, type_pack<valence_info<V_TYPE, 2, C_MASK, 2, IS_JOINABLE, JOIN_WITH, ORDER>, SUMA...>, R_CONT, A_CONT>:
		public axpy_create_general_loop<T, STR, STA, SI, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, axpy_sum_loop<T, STR, STA, SI, valence_info<V_TYPE, 2, C_MASK, 2, IS_JOINABLE, JOIN_WITH, ORDER> > >
//		public gemvA<T, STR, STA, STB, SI, type_pack<SPREAD...>, type_pack<COMMON...>, type_pack<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_pack<valence_info<RB_V_TYPE, 5, RB_C_MASK, 5, RB_IS_JOINABLE, RB_JOIN_WITH, RB_ORDER>, RB...>, type_pack<valence_info<AB_V_TYPE, 6, AB_C_MASK, 6, AB_IS_JOINABLE, AB_JOIN_WITH, AB_ORDER>, AB...>, type_pack<SUMA...>, type_pack<SUMB...> >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct axpy_test_sum<T, STR, STA, type_pack<type_pack<>, SI>, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>:
		public axpy_create_general_loop<T, STR, STA, SI, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, axpy_common_loop<T, STR, STA, SI> >
//		public common_final<T, STR, STA, STB, SI, type_pack<SPREAD...>, type_pack<COMMON...>, type_pack<RA...>, type_pack<RB...>, type_pack<AB...>, type_pack<SUMA...>, type_pack<SUMB...> >
	{
	};

	template <typename T, typename STR, typename STA, typename VISI, typename SPREAD, typename RA, typename SUMA, bool R_CONT, bool A_CONT>
	struct axpy_by_mask;

//  SPREAD
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct axpy_by_mask<T, STR, STA, type_pack<type_pack<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>:
		public axpy_by_mask<T, STR, STA, type_pack<type_pack<VI...>, SI>, typename axpy_insert_mask<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<SPREAD...>, type_pack<> >::type, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>
	{
	};
//  RA
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct axpy_by_mask<T, STR, STA, type_pack<type_pack<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>:
		public axpy_by_mask<T, STR, STA, type_pack<type_pack<VI...>, SI>, type_pack<SPREAD...>, typename axpy_insert_mask<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<RA...>, type_pack<> >::type, type_pack<SUMA...>, R_CONT || ((V_MASK==3 && (C_MASK & 1))), A_CONT || ((V_MASK==3 && (C_MASK & 2)))>
	{
	};
//	SUMA
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct axpy_by_mask<T, STR, STA, type_pack<type_pack<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>:
		public axpy_by_mask<T, STR, STA, type_pack<type_pack<VI...>, SI>, type_pack<SPREAD...>, type_pack<RA...>, typename axpy_insert_mask<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<SUMA...>, type_pack<> >::type, R_CONT, A_CONT>
	{
	};
//	DONE axpy
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct axpy_by_mask<T, STR, STA, type_pack<type_pack<>, SI>, type_pack<SPREAD...>, type_pack<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_pack<SUMA...>, R_CONT, A_CONT>:
		public axpy_create_general_loop<T, STR, STA, SI, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, axpy_loop<T, STR, STA, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER> > >
	{
	};

//	DONE just spread
	template <typename T, typename STR, typename STA, typename SI, int S_V_TYPE, int S_C_MASK, bool S_IS_JOINABLE, int S_JOIN_WITH, size_t S_ORDER, typename ... SPREAD, bool R_CONT, bool A_CONT>
	struct axpy_by_mask<T, STR, STA, type_pack<type_pack<>, SI>, type_pack<valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER>, SPREAD...>, type_pack<>, type_pack<>, R_CONT, A_CONT>:
		public axpy_create_general_loop<T, STR, STA, SI, type_pack<SPREAD...>, type_pack<>, type_pack<>, axpy_just_spread_loop<T, STR, STA, SI, valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER> > >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct axpy_by_mask<T, STR, STA, type_pack<type_pack<>, SI>, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>:
		public axpy_test_sum<T, STR, STA, type_pack<type_pack<>, SI>, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>
	{
	};

	template <typename T, typename STR, typename STA>
	struct axpy_runner: public axpy_by_mask<T, STR, STA, typename make_valence_info_and_join<STR, STA>::type, type_pack<>, type_pack<>, type_pack<>, false, false>
	{
	};

};




#endif /* TAXPY_H_ */