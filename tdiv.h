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

namespace tpp
{

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename BASE_LOOP>
	struct div_suma_loop
	{
//		static_assert(false,"Free argument index is not allowed");
	};

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename BASE_LOOP>
	struct div_ra_loop;

	template <typename T, typename STR, typename STA, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct div_ra_loop<T, STR, STA, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T> ritype;
		typedef typename iterator_getter_by_src_v_type<STA, SI, 1, V_TYPE>::template iterator_type<T> aitype;
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

	template <typename T, typename STR, typename STA, typename SI>
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

	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, bool IS_FIRST, typename BASE_LOOP>
	struct div_spread_loop;

// LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename BASE_LOOP>
	struct div_spread_loop<T, STR, STA, SI, ITTYPE, false, BASE_LOOP>: public BASE_LOOP
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
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename BASE_LOOP>
	struct div_spread_loop<T, STR, STA, SI, ITTYPE, true, BASE_LOOP>: public BASE_LOOP
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
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename HEAD, typename BASE_BASE_LOOP>
	struct div_spread_loop<T, STR, STA, SI, ITTYPE, false, div_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> >: public div_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef div_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
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
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename HEAD, typename BASE_BASE_LOOP>
	struct div_spread_loop<T, STR, STA, SI, ITTYPE, true, div_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> >: public div_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef div_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
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

	template <typename T, typename STR, typename STA, typename SI, typename SPREAD, typename RA, typename SUMA, typename BASE_LOOP>
	struct div_create_general_loop
	{
		using type=BASE_LOOP;
	};

	template <typename T, typename STR, typename STA, typename SI, typename SPREAD, typename RA, typename HEAD, typename ... SUMA, typename BASE_LOOP>
	struct div_create_general_loop<T, STR, STA, SI, SPREAD, RA, type_pack<HEAD, SUMA...>, BASE_LOOP>:
		public div_create_general_loop<T, STR, STA, SI, SPREAD, RA, type_pack<SUMA...>, div_suma_loop<T, STR, STA, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename SPREAD, typename HEAD, typename ... RA, typename BASE_LOOP>
	struct div_create_general_loop<T, STR, STA, SI, SPREAD, type_pack<HEAD, RA...>, type_pack<>, BASE_LOOP>:
		public div_create_general_loop<T, STR, STA, SI, SPREAD, type_pack<RA...>, type_pack<>, div_ra_loop<T, STR, STA, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename ... SPREAD, typename BASE_LOOP>
	struct div_create_general_loop<T, STR, STA, SI, type_pack<HEAD, SPREAD...>, type_pack<>, type_pack<>, BASE_LOOP>:
		public div_create_general_loop<T, STR, STA, SI, type_pack<SPREAD...>, type_pack<>, type_pack<>, div_spread_loop<T, STR, STA, SI, typename iterator_getter_by_src_v_type<STR, SI, 0, HEAD::v_type>::template iterator_type<T>, false, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename BASE_LOOP>
	struct div_create_general_loop<T, STR, STA, SI, type_pack<HEAD>, type_pack<>, type_pack<>, BASE_LOOP>:
		public div_create_general_loop<T, STR, STA, SI, type_pack<>, type_pack<>, type_pack<>, div_spread_loop<T, STR, STA, SI, typename iterator_getter_by_src_v_type<STR, SI, 0, HEAD::v_type>::template iterator_type<T>, true, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename S>
	struct div_just_spread_loop;

	template <typename T, typename STR, typename STA, typename SI, int S_V_TYPE, int S_C_MASK, bool S_IS_JOINABLE, int S_JOIN_WITH, size_t S_ORDER>
	struct div_just_spread_loop<T, STR, STA, SI, valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER> >
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


	template <typename ELEMENT, typename VI, typename RES>
	struct div_insert_mask;

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, int H_V_TYPE, int H_MASK, int H_C_MASK, int H_V_MASK, bool H_IS_JOINABLE, int H_JOIN_WITH, size_t H_ORDER, typename ... VI, typename ... RES>
	struct div_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>, type_pack<RES...> >:
		public std::conditional<
			(MASK==V_MASK && H_MASK!=H_V_MASK)?true:
				(MASK!=V_MASK && H_MASK==H_V_MASK)?false:
					(C_MASK > H_C_MASK)?true:
						(C_MASK < H_C_MASK)?false:
							(ORDER<H_ORDER),
			type_pack<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>,
			div_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<VI...>, type_pack<RES..., valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER> > >
		>::type
	{
	};

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... RES>
	struct div_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<>, type_pack<RES...> >:
		public type_pack<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER> >
	{
	};

	template <typename T, typename STR, typename STA, typename VISI, typename SPREAD, typename RA, typename SUMA, bool R_CONT, bool A_CONT>
	struct div_by_mask;

//  SPREAD
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct div_by_mask<T, STR, STA, type_pack<type_pack<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>:
		public div_by_mask<T, STR, STA, type_pack<type_pack<VI...>, SI>, typename div_insert_mask<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<SPREAD...>, type_pack<> >::type, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>
	{
	};
//  RA
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct div_by_mask<T, STR, STA, type_pack<type_pack<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>:
		public div_by_mask<T, STR, STA, type_pack<type_pack<VI...>, SI>, type_pack<SPREAD...>, typename div_insert_mask<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<RA...>, type_pack<> >::type, type_pack<SUMA...>, R_CONT || ((V_MASK==3 && (C_MASK & 1))), A_CONT || ((V_MASK==3 && (C_MASK & 2)))>
	{
	};
//	SUMA
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct div_by_mask<T, STR, STA, type_pack<type_pack<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>:
		public div_by_mask<T, STR, STA, type_pack<type_pack<VI...>, SI>, type_pack<SPREAD...>, type_pack<RA...>, typename div_insert_mask<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_pack<SUMA...>, type_pack<> >::type, R_CONT, A_CONT>
	{
	};

//	DONE just spread
	template <typename T, typename STR, typename STA, typename SI, int S_V_TYPE, int S_C_MASK, bool S_IS_JOINABLE, int S_JOIN_WITH, size_t S_ORDER, typename ... SPREAD, bool R_CONT, bool A_CONT>
	struct div_by_mask<T, STR, STA, type_pack<type_pack<>, SI>, type_pack<valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER>, SPREAD...>, type_pack<>, type_pack<>, R_CONT, A_CONT>:
		public div_create_general_loop<T, STR, STA, SI, type_pack<SPREAD...>, type_pack<>, type_pack<>, div_just_spread_loop<T, STR, STA, SI, valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER> > >
	{
		using vi_by_mask=type_pack<
				type_pack<>,
				type_pack<valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER>, SPREAD...>,
				type_pack<>,
				type_pack<> >;
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct div_by_mask<T, STR, STA, type_pack<type_pack<>, SI>, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, R_CONT, A_CONT>:
		public div_create_general_loop<T, STR, STA, SI, type_pack<SPREAD...>, type_pack<RA...>, type_pack<SUMA...>, div_common_loop<T, STR, STA, SI> >
	{
		using vi_by_mask=type_pack<
				type_pack<>,
				type_pack<SPREAD...>,
				type_pack<SUMA...>,
				type_pack<RA...> >;
	};

	template <typename T, typename STR, typename STA>
	struct div_runner: public div_by_mask<T, STR, STA, typename make_valence_info_and_join<STR, STA>::type, type_pack<>, type_pack<>, type_pack<>, false, false>
	{
	};

};




#endif /* TDIV_H_ */
