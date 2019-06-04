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

namespace tpp
{

	template <typename T, typename STR, typename SI>
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


	template <typename T, typename STR, typename SI, typename ITTYPE, typename BASE_LOOP>
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

	template <typename T, typename STR, typename SI, typename SPREAD, typename BASE_LOOP>
	struct scal_create_general_loop
	{
		using type=BASE_LOOP;
	};

	template <typename T, typename STR, typename SI, typename HEAD, typename ... SPREAD, typename BASE_LOOP>
	struct scal_create_general_loop<T, STR, SI, type_sequence<HEAD, SPREAD...>, BASE_LOOP>:
		public scal_create_general_loop<T, STR, SI, type_sequence<SPREAD...>, scal_spread_loop<T, STR, SI, typename iterator_getter_by_src_v_type<STR, SI, 0, HEAD::v_type>::template iterator_type<T>, BASE_LOOP> >
	{
	};

	template <typename T, typename STR, typename SI, typename RAV>
	struct scal_loop;

	template <typename T, typename STR, typename SI, int V_TYPE, int C_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER>
	struct scal_loop<T, STR, SI, valence_info<V_TYPE, 1, C_MASK, 1, IS_JOINABLE, JOIN_WITH, ORDER> >
	{
		const STR& str;
		const T alpha;
	private:
		BLAS_INTEGER N;
		BLAS_INTEGER incr;
	public:
		scal_loop(const STR& str): str(str), alpha(1),
		N(iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T>::length(str)),
		incr(iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T>::step(str))
		{}
		scal_loop(const STR& str, T alpha): str(str), alpha(alpha),
		N(iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T>::length(str)),
		incr(iterator_getter_by_src_v_type<STR, SI, 0, V_TYPE>::template iterator_type<T>::step(str))
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


	template <typename ELEMENT, typename VI, typename RES>
	struct scal_insert_mask;

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, int H_V_TYPE, int H_MASK, int H_C_MASK, int H_V_MASK, bool H_IS_JOINABLE, int H_JOIN_WITH, size_t H_ORDER, typename ... VI, typename ... RES>
	struct scal_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>, type_sequence<RES...> >:
		public std::conditional<
			(MASK==V_MASK && H_MASK!=H_V_MASK)?true:
				(MASK!=V_MASK && H_MASK==H_V_MASK)?false:
					(C_MASK > H_C_MASK)?true:
						(C_MASK < H_C_MASK)?false:
							(ORDER<H_ORDER),
			type_sequence<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>,
			scal_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<VI...>, type_sequence<RES..., valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER> > >
		>::type
	{
	};

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... RES>
	struct scal_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<>, type_sequence<RES...> >:
		public type_sequence<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER> >
	{
	};

	template <typename T, typename STR, typename VISI, typename SPREAD, bool R_CONT, bool A_CONT>
	struct scal_by_mask;

//  SPREAD
	template <typename T, typename STR, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... SPREAD, bool R_CONT, bool A_CONT>
	struct scal_by_mask<T, STR, type_sequence<type_sequence<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<SPREAD...>, R_CONT, A_CONT>:
		public scal_by_mask<T, STR, type_sequence<type_sequence<VI...>, SI>, typename scal_insert_mask<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<SPREAD...>, type_sequence<> >::type, R_CONT, A_CONT>
	{
	};
//	DONE VECTOR
	template <typename T, typename STR, typename SI, typename ... SPREAD, int V_TYPE, int C_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, bool R_CONT, bool A_CONT>
	struct scal_by_mask<T, STR, type_sequence<type_sequence<>, SI>, type_sequence<valence_info<V_TYPE, 1, C_MASK, 1, IS_JOINABLE, JOIN_WITH, ORDER>, SPREAD...>, R_CONT, A_CONT>:
		public scal_create_general_loop<T, STR, SI, type_sequence<SPREAD...>, scal_loop<T, STR, SI, valence_info<V_TYPE, 1, C_MASK, 1, IS_JOINABLE, JOIN_WITH, ORDER> > >
	{
	};

//	DONE COMMON
	template <typename T, typename STR, typename SI, typename ... SPREAD, bool R_CONT, bool A_CONT>
	struct scal_by_mask<T, STR, type_sequence<type_sequence<>, SI>, type_sequence<SPREAD...>, R_CONT, A_CONT>:
		public scal_create_general_loop<T, STR, SI, type_sequence<SPREAD...>, scal_common_loop<T, STR, SI> >
	{
	};

	template <typename T, typename STR>
	struct scal_runner: public scal_by_mask<T, STR, typename make_valence_info_and_join<STR>::type, type_sequence<>, false, false>
	{
	};

};

#endif /* TSCAL_H_ */
