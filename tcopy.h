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

namespace tpp
{

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename BASE_LOOP>
	struct copy_suma_loop;

	template <typename T, typename STR, typename STA, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct copy_suma_loop<T, STR, STA, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STA, SI, 1, V_TYPE>::template iterator_type<T> itype;
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

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename BASE_LOOP>
	struct copy_ra_loop;

	template <typename T, typename STR, typename STA, typename SI, int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename BASE_LOOP>
	struct copy_ra_loop<T, STR, STA, SI, valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, BASE_LOOP>: public BASE_LOOP
	{
	private:
	public:
		typedef typename iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T> ritype;
		typedef typename iterator_getter_by_src_v_type<STA, SI, 1, V_TYPE>::template iterator_type<T> aitype;
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
			ritype rit(this->str, ri);
			for (;rit.not_end(ri);rit.move_one(ri))
				BASE_LOOP::spread(r0,ri);
		}
	};

	template <typename T, typename STR, typename STA, typename SI>
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

	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, bool IS_FIRST, typename BASE_LOOP>
	struct copy_spread_loop;

// LAST; NOT FIRST
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename BASE_LOOP>
	struct copy_spread_loop<T, STR, STA, SI, ITTYPE, false, BASE_LOOP>: public BASE_LOOP
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
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename BASE_LOOP>
	struct copy_spread_loop<T, STR, STA, SI, ITTYPE, true, BASE_LOOP>: public BASE_LOOP
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
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename HEAD, typename BASE_BASE_LOOP>
	struct copy_spread_loop<T, STR, STA, SI, ITTYPE, false, copy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> >: public copy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef copy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
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
	template <typename T, typename STR, typename STA, typename SI, typename ITTYPE, typename HEAD, typename BASE_BASE_LOOP>
	struct copy_spread_loop<T, STR, STA, SI, ITTYPE, true, copy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> >: public copy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP>
	{
	private:
		typedef copy_spread_loop<T, STR, STA, SI, HEAD, false, BASE_BASE_LOOP> BASE_LOOP;
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

	template <typename T, typename STR, typename STA, typename SI, typename SPREAD, typename RA, typename SUMA, typename BASE_LOOP>
	struct copy_create_general_loop
	{
		using type=BASE_LOOP;
	};

	template <typename T, typename STR, typename STA, typename SI, typename SPREAD, typename RA, typename HEAD, typename ... SUMA, typename BASE_LOOP>
	struct copy_create_general_loop<T, STR, STA, SI, SPREAD, RA, type_sequence<HEAD, SUMA...>, BASE_LOOP>:
		public copy_create_general_loop<T, STR, STA, SI, SPREAD, RA, type_sequence<SUMA...>, copy_suma_loop<T, STR, STA, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename SPREAD, typename HEAD, typename ... RA, typename BASE_LOOP>
	struct copy_create_general_loop<T, STR, STA, SI, SPREAD, type_sequence<HEAD, RA...>, type_sequence<>, BASE_LOOP>:
		public copy_create_general_loop<T, STR, STA, SI, SPREAD, type_sequence<RA...>, type_sequence<>, copy_ra_loop<T, STR, STA, SI, HEAD, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename ... SPREAD, typename BASE_LOOP>
	struct copy_create_general_loop<T, STR, STA, SI, type_sequence<HEAD, SPREAD...>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public copy_create_general_loop<T, STR, STA, SI, type_sequence<SPREAD...>, type_sequence<>, type_sequence<>, copy_spread_loop<T, STR, STA, SI, typename iterator_getter_by_src_v_type<STR, SI, 0, HEAD::v_type>::template iterator_type<T>, false, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename HEAD, typename BASE_LOOP>
	struct copy_create_general_loop<T, STR, STA, SI, type_sequence<HEAD>, type_sequence<>, type_sequence<>, BASE_LOOP>:
		public copy_create_general_loop<T, STR, STA, SI, type_sequence<>, type_sequence<>, type_sequence<>, copy_spread_loop<T, STR, STA, SI, typename iterator_getter_by_src_v_type<STR, SI, 0, HEAD::v_type>::template iterator_type<T>, true, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename STA, typename SI, typename RAC, typename RAV>
	struct lacpy_loop;

	template <typename T, typename STR, typename STA, typename SI, int RAC_V_TYPE, int RAC_C_MASK, bool RAC_IS_JOINABLE, int RAC_JOIN_WITH, size_t RAC_ORDER, int RAV_V_TYPE, int RAV_C_MASK, bool RAV_IS_JOINABLE, int RAV_JOIN_WITH, size_t RAV_ORDER>
	struct lacpy_loop<T, STR, STA, SI, valence_info<RAC_V_TYPE, 3, RAC_C_MASK, 3, RAC_IS_JOINABLE, RAC_JOIN_WITH, RAC_ORDER>, valence_info<RAV_V_TYPE, 3, RAV_C_MASK, 3, RAV_IS_JOINABLE, RAV_JOIN_WITH, RAV_ORDER> >
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
		M(iterator_getter_by_src_v_type<STR, SI, 0, RAC_V_TYPE>::template iterator_type<T>::length(str)),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::length(str)),
		lda(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::step(str)),
		ldb(iterator_getter_by_src_v_type<STA, SI, 1, RAV_V_TYPE>::template iterator_type<T>::step(sta))
		{}
		lacpy_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha),
		M(iterator_getter_by_src_v_type<STR, SI, 0, RAC_V_TYPE>::template iterator_type<T>::length(str)),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::length(str)),
		lda(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::step(str)),
		ldb(iterator_getter_by_src_v_type<STA, SI, 1, RAV_V_TYPE>::template iterator_type<T>::step(sta))
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

	template <typename T, typename STR, typename STA, typename SI, typename RAV>
	struct copy_loop;

	template <typename T, typename STR, typename STA, typename SI, int RAV_V_TYPE, int RAV_C_MASK, bool RAV_IS_JOINABLE, int RAV_JOIN_WITH, size_t RAV_ORDER>
	struct copy_loop<T, STR, STA, SI, valence_info<RAV_V_TYPE, 3, RAV_C_MASK, 3, RAV_IS_JOINABLE, RAV_JOIN_WITH, RAV_ORDER> >
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
		N(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::length(str)),
		incr(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::step(str)),
		inca(iterator_getter_by_src_v_type<STA, SI, 1, RAV_V_TYPE>::template iterator_type<T>::step(sta))
		{}
		copy_loop(const STR& str, const STA& sta, T alpha): str(str), sta(sta), alpha(alpha),
		N(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::length(str)),
		incr(iterator_getter_by_src_v_type<STR, SI, 0, RAV_V_TYPE>::template iterator_type<T>::step(str)),
		inca(iterator_getter_by_src_v_type<STA, SI, 1, RAV_V_TYPE>::template iterator_type<T>::step(sta))
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

	template <typename T, typename STR, typename STA, typename SI, typename S>
	struct just_spread_loop;

	template <typename T, typename STR, typename STA, typename SI, int S_V_TYPE, int S_C_MASK, bool S_IS_JOINABLE, int S_JOIN_WITH, size_t S_ORDER>
	struct just_spread_loop<T, STR, STA, SI, valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER> >
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

	template <typename T, typename STR, typename STA, typename SI, typename S>
	struct asum_loop;

	template <typename T, typename STR, typename STA, typename SI, int V_TYPE, int C_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER>
	struct asum_loop<T, STR, STA, SI, valence_info<V_TYPE, 2, C_MASK, 2, IS_JOINABLE, JOIN_WITH, ORDER> >
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



	template <typename ELEMENT, typename VI, typename RES>
	struct copy_insert_mask;

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, int H_V_TYPE, int H_MASK, int H_C_MASK, int H_V_MASK, bool H_IS_JOINABLE, int H_JOIN_WITH, size_t H_ORDER, typename ... VI, typename ... RES>
	struct copy_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>, type_sequence<RES...> >:
		public std::conditional<
			(MASK==V_MASK && H_MASK!=H_V_MASK)?true:
				(MASK!=V_MASK && H_MASK==H_V_MASK)?false:
					(C_MASK > H_C_MASK)?true:
						(C_MASK < H_C_MASK)?false:
							(ORDER<H_ORDER),
			type_sequence<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>,
			copy_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<VI...>, type_sequence<RES..., valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER> > >
		>::type
	{
	};

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... RES>
	struct copy_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<>, type_sequence<RES...> >:
		public type_sequence<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER> >
	{
	};

	template <typename T, typename STR, typename STA, typename VISI, typename SPREAD, typename RA, typename SUMA, bool R_CONT, bool A_CONT>
	struct copy_test_sum;

//	DONE sum
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, typename ... RA, int V_TYPE, int C_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct copy_test_sum<T, STR, STA, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<valence_info<V_TYPE, 2, C_MASK, 2, IS_JOINABLE, JOIN_WITH, ORDER>, SUMA...>, R_CONT, A_CONT>:
		public copy_create_general_loop<T, STR, STA, SI, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, asum_loop<T, STR, STA, SI, valence_info<V_TYPE, 2, C_MASK, 2, IS_JOINABLE, JOIN_WITH, ORDER> > >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct copy_test_sum<T, STR, STA, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>:
		public copy_create_general_loop<T, STR, STA, SI, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, copy_common_loop<T, STR, STA, SI> >
	{
	};

	template <typename T, typename STR, typename STA, typename VISI, typename SPREAD, typename RA, typename SUMA, bool R_CONT, bool A_CONT>
	struct copy_by_mask;

//  SPREAD
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct copy_by_mask<T, STR, STA, type_sequence<type_sequence<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>:
		public copy_by_mask<T, STR, STA, type_sequence<type_sequence<VI...>, SI>, typename copy_insert_mask<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<SPREAD...>, type_sequence<> >::type, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>
	{
	};
//  RA
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct copy_by_mask<T, STR, STA, type_sequence<type_sequence<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>:
		public copy_by_mask<T, STR, STA, type_sequence<type_sequence<VI...>, SI>, type_sequence<SPREAD...>, typename copy_insert_mask<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<RA...>, type_sequence<> >::type, type_sequence<SUMA...>, R_CONT || ((V_MASK==3 && (C_MASK & 1))), A_CONT || ((V_MASK==3 && (C_MASK & 2)))>
	{
	};
//	SUMA
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct copy_by_mask<T, STR, STA, type_sequence<type_sequence<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>:
		public copy_by_mask<T, STR, STA, type_sequence<type_sequence<VI...>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, typename copy_insert_mask<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<SUMA...>, type_sequence<> >::type, R_CONT, A_CONT>
	{
	};
//	DONE lacpy
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, int RAC_V_TYPE, bool RAC_IS_JOINABLE, int RAC_JOIN_WITH, size_t RAC_ORDER, int RAV_V_TYPE, int RAV_C_MASK, bool RAV_IS_JOINABLE, int RAV_JOIN_WITH, size_t RAV_ORDER, typename ... RA, typename ... SUMA>
	struct copy_by_mask<T, STR, STA, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<valence_info<RAC_V_TYPE, 3, 3, 3, RAC_IS_JOINABLE, RAC_JOIN_WITH, RAC_ORDER>, valence_info<RAV_V_TYPE, 3, RAV_C_MASK, 3, RAV_IS_JOINABLE, RAV_JOIN_WITH, RAV_ORDER>, RA...>, type_sequence<SUMA...>, true, true>:
		public copy_create_general_loop<T, STR, STA, SI, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, lacpy_loop<T, STR, STA, SI, valence_info<RAC_V_TYPE, 3, 3, 3, RAC_IS_JOINABLE, RAC_JOIN_WITH, RAC_ORDER>, valence_info<RAV_V_TYPE, 3, RAV_C_MASK, 3, RAV_IS_JOINABLE, RAV_JOIN_WITH, RAV_ORDER> > >
	{
	};
//	DONE copy
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, int RA_V_TYPE, int RA_C_MASK, bool RA_IS_JOINABLE, int RA_JOIN_WITH, size_t RA_ORDER, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct copy_by_mask<T, STR, STA, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER>, RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>:
		public copy_create_general_loop<T, STR, STA, SI, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, copy_loop<T, STR, STA, SI, valence_info<RA_V_TYPE, 3, RA_C_MASK, 3, RA_IS_JOINABLE, RA_JOIN_WITH, RA_ORDER> > >
	{
	};

//	DONE just spread
	template <typename T, typename STR, typename STA, typename SI, int S_V_TYPE, int S_C_MASK, bool S_IS_JOINABLE, int S_JOIN_WITH, size_t S_ORDER, typename ... SPREAD, bool R_CONT, bool A_CONT>
	struct copy_by_mask<T, STR, STA, type_sequence<type_sequence<>, SI>, type_sequence<valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER>, SPREAD...>, type_sequence<>, type_sequence<>, R_CONT, A_CONT>:
		public copy_create_general_loop<T, STR, STA, SI, type_sequence<SPREAD...>, type_sequence<>, type_sequence<>, just_spread_loop<T, STR, STA, SI, valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER> > >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct copy_by_mask<T, STR, STA, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>:
		public copy_test_sum<T, STR, STA, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>
	{
	};

	template <typename T, typename STR, typename STA>
	struct copy_runner: public copy_by_mask<T, STR, STA, typename make_valence_info_and_join<STR, STA>::type, type_sequence<>, type_sequence<>, type_sequence<>, false, false>
	{
	};

	template <typename T, typename STR, typename STA, typename VISI, typename SPREAD, typename RA, typename SUMA, bool R_CONT, bool A_CONT>
	struct asum_by_mask;

//  SPREAD
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct asum_by_mask<T, STR, STA, type_sequence<type_sequence<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>:
		public asum_by_mask<T, STR, STA, type_sequence<type_sequence<VI...>, SI>, typename copy_insert_mask<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<SPREAD...>, type_sequence<> >::type, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>
	{
	};
//  RA
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct asum_by_mask<T, STR, STA, type_sequence<type_sequence<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>:
		public asum_by_mask<T, STR, STA, type_sequence<type_sequence<VI...>, SI>, type_sequence<SPREAD...>, typename copy_insert_mask<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<RA...>, type_sequence<> >::type, type_sequence<SUMA...>, R_CONT || ((V_MASK==3 && (C_MASK & 1))), A_CONT || ((V_MASK==3 && (C_MASK & 2)))>
	{
	};
//	SUMA
	template <typename T, typename STR, typename STA, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct asum_by_mask<T, STR, STA, type_sequence<type_sequence<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>:
		public asum_by_mask<T, STR, STA, type_sequence<type_sequence<VI...>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, typename copy_insert_mask<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<SUMA...>, type_sequence<> >::type, R_CONT, A_CONT>
	{
	};

// DONE ASUM
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, typename ... RA, int V_TYPE, int C_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct asum_by_mask<T, STR, STA, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<valence_info<V_TYPE, 2, C_MASK, 2, IS_JOINABLE, JOIN_WITH, ORDER>, SUMA...>, R_CONT, A_CONT>:
		public copy_create_general_loop<T, STR, STA, SI, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, asum_loop<T, STR, STA, SI, valence_info<V_TYPE, 2, C_MASK, 2, IS_JOINABLE, JOIN_WITH, ORDER> > >
	{
		using vi_by_mask=type_sequence<
				type_sequence<>,
				type_sequence<SPREAD...>,
				type_sequence<valence_info<V_TYPE, 2, C_MASK, 2, IS_JOINABLE, JOIN_WITH, ORDER>, SUMA...>,
				type_sequence<RA...> >;
	};

//	DONE just spread
	template <typename T, typename STR, typename STA, typename SI, int S_V_TYPE, int S_C_MASK, bool S_IS_JOINABLE, int S_JOIN_WITH, size_t S_ORDER, typename ... SPREAD, bool R_CONT, bool A_CONT>
	struct asum_by_mask<T, STR, STA, type_sequence<type_sequence<>, SI>, type_sequence<valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER>, SPREAD...>, type_sequence<>, type_sequence<>, R_CONT, A_CONT>:
		public copy_create_general_loop<T, STR, STA, SI, type_sequence<SPREAD...>, type_sequence<>, type_sequence<>, just_spread_loop<T, STR, STA, SI, valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER> > >
	{
		using vi_by_mask=type_sequence<
				type_sequence<>,
				type_sequence<valence_info<S_V_TYPE, 1, S_C_MASK, 1, S_IS_JOINABLE, S_JOIN_WITH, S_ORDER>, SPREAD...>,
				type_sequence<>,
				type_sequence<> >;
	};

//	DONE COMMON
	template <typename T, typename STR, typename STA, typename SI, typename ... SPREAD, typename ... RA, typename ... SUMA, bool R_CONT, bool A_CONT>
	struct asum_by_mask<T, STR, STA, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, R_CONT, A_CONT>:
		public copy_create_general_loop<T, STR, STA, SI, type_sequence<SPREAD...>, type_sequence<RA...>, type_sequence<SUMA...>, copy_common_loop<T, STR, STA, SI> >
	{
		using vi_by_mask=type_sequence<
				type_sequence<>,
				type_sequence<SPREAD...>,
				type_sequence<SUMA...>,
				type_sequence<RA...> >;
	};

	template <typename T, typename STR, typename STA>
	struct asum_runner: public asum_by_mask<T, STR, STA, typename make_valence_info_and_join<STR, STA>::type, type_sequence<>, type_sequence<>, type_sequence<>, false, false>
	{
	};


};

#endif /* TCOPY_H_ */
