/*
 * gesv.h
 *
 *  Created on: 24 мая 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef GESV_H_
#define GESV_H_

#include <tvalence.h>
#include <blas_tmpl.h>
#include <smplmath.h>

namespace tpp
{

	template <typename T, typename STM, typename STV, typename SI, typename M_FREE, typename COMMON, typename V_FREE, bool R_CONT>
	struct gesv_matrix
	{
		static_assert(true, "Internal error during instantiation of gesv_matrix");
	};

	template <typename T, typename STM, typename STV, typename SI,
		int M_V_TYPE, int M_C_MASK, int M_JOIN_WITH, size_t M_ORDER,
		int C_V_TYPE, int C_C_MASK, int C_JOIN_WITH, size_t C_ORDER,
		int V_V_TYPE, int V_C_MASK, int V_JOIN_WITH, size_t V_ORDER,
		bool R_CONT>
	struct gesv_matrix<T, STM, STV, SI,
		type_sequence<valence_info<M_V_TYPE, 1, M_C_MASK, 1, false, M_JOIN_WITH, M_ORDER> >,
		type_sequence<valence_info<C_V_TYPE, 3, C_C_MASK, 3, false, C_JOIN_WITH, C_ORDER> >,
		type_sequence<valence_info<V_V_TYPE, 2, V_C_MASK, 2, false, V_JOIN_WITH, V_ORDER> >, R_CONT>
	{
		static_assert(C_C_MASK & 2,"gesv: Linked dimension of arguments should be continuous");
		static BLAS_INTEGER run_gesv(const STM& stm, const STV& stv, T *M, T *V, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::length(stm);
			if ((size_t)N!=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::step(stm);
			BLAS_INTEGER LDB=iterator_getter_by_src_v_type<STV, SI, 1, V_V_TYPE>::template iterator_type<T>::step(stv);
			BLAS_INTEGER NRHS=iterator_getter_by_src_v_type<STV, SI, 1, V_V_TYPE>::template iterator_type<T>::length(stv);
			if (!IPIV)
				IPIV=(BLAS_INTEGER *)alloca((N>1?N:1)*sizeof(BLAS_INTEGER));
			BLAS_INTEGER INFO;
			getrf<T>(&N, &N, M, &LDA, IPIV, &INFO);
			if (INFO)
				return INFO;
			getrs<T>("N", &N, &NRHS, M, &LDA, IPIV, V, &LDB, &INFO);
			return INFO;
		}
		static BLAS_INTEGER run_solve(const STM& stm, const STV& stv, T *M, T *V, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::length(stm);
			if ((size_t)N!=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::step(stm);
			BLAS_INTEGER LDB=iterator_getter_by_src_v_type<STV, SI, 1, V_V_TYPE>::template iterator_type<T>::step(stv);
			BLAS_INTEGER NRHS=iterator_getter_by_src_v_type<STV, SI, 1, V_V_TYPE>::template iterator_type<T>::length(stv);
			BLAS_INTEGER INFO;
			getrs<T>("N", &N, &NRHS, M, &LDA, IPIV, V, &LDB, &INFO);
			return INFO;
		}
	};

	template <typename T, typename STM, typename STV, typename SI,
	int M_V_TYPE, int M_C_MASK, int M_JOIN_WITH, size_t M_ORDER,
	int C_V_TYPE, int C_C_MASK, int C_JOIN_WITH, size_t C_ORDER,
	int V_V_TYPE, int V_C_MASK, int V_JOIN_WITH, size_t V_ORDER>
	struct gesv_matrix<T, STM, STV, SI,
	type_sequence<valence_info<M_V_TYPE, 1, M_C_MASK, 1, false, M_JOIN_WITH, M_ORDER> >,
	type_sequence<valence_info<C_V_TYPE, 3, C_C_MASK, 3, false, C_JOIN_WITH, C_ORDER> >,
	type_sequence<valence_info<V_V_TYPE, 2, V_C_MASK, 2, false, V_JOIN_WITH, V_ORDER> >, 0>
	{
		static_assert(C_C_MASK & 2,"gesv: Linked dimension of arguments should be continuous");
		static BLAS_INTEGER run_gesv(const STM& stm, const STV& stv, T *M, T *V, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::length(stm);
			if ((size_t)N!=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::step(stm);
			BLAS_INTEGER LDB=iterator_getter_by_src_v_type<STV, SI, 1, V_V_TYPE>::template iterator_type<T>::step(stv);
			BLAS_INTEGER NRHS=iterator_getter_by_src_v_type<STV, SI, 1, V_V_TYPE>::template iterator_type<T>::length(stv);
			if (!IPIV)
				IPIV=(BLAS_INTEGER *)alloca((N>1?N:1)*sizeof(BLAS_INTEGER));
			BLAS_INTEGER INFO;
			getrf<T>(&N, &N, M, &LDA, IPIV, &INFO);
			if (INFO)
				return INFO;
			getrs<T>("T", &N, &NRHS, M, &LDA, IPIV, V, &LDB, &INFO);
			return INFO;
		}
		static BLAS_INTEGER run_solve(const STM& stm, const STV& stv, T *M, T *V, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::length(stm);
			if ((size_t)N!=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::step(stm);
			BLAS_INTEGER LDB=iterator_getter_by_src_v_type<STV, SI, 1, V_V_TYPE>::template iterator_type<T>::step(stv);
			BLAS_INTEGER NRHS=iterator_getter_by_src_v_type<STV, SI, 1, V_V_TYPE>::template iterator_type<T>::length(stv);
			BLAS_INTEGER INFO;
			getrs<T>("T", &N, &NRHS, M, &LDA, IPIV, V, &LDB, &INFO);
			return INFO;
		}
	};


	template <typename T, typename STM, typename STV, typename SI, typename M_FREE, typename COMMON, bool R_CONT>
	struct gesv_vector
	{
		static_assert(true, "Internal error during instantiation of gesv_vector");
	};

	template <typename T, typename STM, typename STV, typename SI,
		int M_V_TYPE, int M_C_MASK, int M_JOIN_WITH, size_t M_ORDER,
		int C_V_TYPE, int C_C_MASK, int C_JOIN_WITH, size_t C_ORDER,
		bool R_CONT>
	struct gesv_vector<T, STM, STV, SI,
	type_sequence<valence_info<M_V_TYPE, 1, M_C_MASK, 1, false, M_JOIN_WITH, M_ORDER> >,
	type_sequence<valence_info<C_V_TYPE, 3, C_C_MASK, 3, false, C_JOIN_WITH, C_ORDER> >, R_CONT>
	{
		static_assert(C_C_MASK & 2,"gesv: Linked dimension of arguments should be continuous");
		static BLAS_INTEGER run_gesv(const STM& stm, const STV& stv, T *M, T *V, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::length(stm);
			if ((size_t)N!=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::step(stm);
			if (!IPIV)
				IPIV=(BLAS_INTEGER *)alloca((N>1?N:1)*sizeof(BLAS_INTEGER));
			BLAS_INTEGER INFO;
			getrf<T>(&N, &N, M, &LDA, IPIV, &INFO);
			if (INFO)
				return INFO;
			getrs<T>("N", &N, &type_constants<BLAS_INTEGER>::one, M, &LDA, IPIV, V, &N, &INFO);
			return INFO;
		}
		static BLAS_INTEGER run_solve(const STM& stm, const STV& stv, T *M, T *V, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::length(stm);
			if ((size_t)N!=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::step(stm);
			BLAS_INTEGER INFO;
			getrs<T>("N", &N, &type_constants<BLAS_INTEGER>::one, M, &LDA, IPIV, V, &N, &INFO);
			return INFO;
		}
	};

	template <typename T, typename STM, typename STV, typename SI,
		int M_V_TYPE, int M_C_MASK, int M_JOIN_WITH, size_t M_ORDER,
		int C_V_TYPE, int C_C_MASK, int C_JOIN_WITH, size_t C_ORDER>
	struct gesv_vector<T, STM, STV, SI,
		type_sequence<valence_info<M_V_TYPE, 1, M_C_MASK, 1, false, M_JOIN_WITH, M_ORDER> >,
		type_sequence<valence_info<C_V_TYPE, 3, C_C_MASK, 3, false, C_JOIN_WITH, C_ORDER> >, 0>
	{
		static_assert(C_C_MASK & 2,"gesv: Linked dimension of arguments should be continuous");
		static BLAS_INTEGER run_gesv(const STM& stm, const STV& stv, T *M, T *V, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::length(stm);
			if ((size_t)N!=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::step(stm);
			if (!IPIV)
				IPIV=(BLAS_INTEGER *)alloca((N>1?N:1)*sizeof(BLAS_INTEGER));
			BLAS_INTEGER INFO;
			getrf<T>(&N, &N, M, &LDA, IPIV, &INFO);
			if (INFO)
				return INFO;
			getrs<T>("T", &N, &type_constants<BLAS_INTEGER>::one, M, &LDA, IPIV, V, &N, &INFO);
			return INFO;
		}
		static BLAS_INTEGER run_solve(const STM& stm, const STV& stv, T *M, T *V, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=iterator_getter_by_src_v_type<STM, SI, 0, M_V_TYPE>::template iterator_type<T>::length(stm);
			if ((size_t)N!=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=iterator_getter_by_src_v_type<STM, SI, 0, C_V_TYPE>::template iterator_type<T>::step(stm);
			BLAS_INTEGER INFO;
			getrs<T>("T", &N, &type_constants<BLAS_INTEGER>::one, M, &LDA, IPIV, V, &N, &INFO);
			return INFO;
		}
	};

	template <typename ELEMENT, typename VI, typename RES>
	struct gesv_insert_mask;

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, int H_V_TYPE, int H_MASK, int H_C_MASK, int H_V_MASK, bool H_IS_JOINABLE, int H_JOIN_WITH, size_t H_ORDER, typename ... VI, typename ... RES>
	struct gesv_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>, type_sequence<RES...> >:
		public std::conditional<
			(MASK==V_MASK && H_MASK!=H_V_MASK)?true:
				(MASK!=V_MASK && H_MASK==H_V_MASK)?false:
					(C_MASK > H_C_MASK)?true:
						(C_MASK < H_C_MASK)?false:
							(ORDER<H_ORDER),
			type_sequence<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER>, VI...>,
			gesv_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<VI...>, type_sequence<RES..., valence_info<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_IS_JOINABLE, H_JOIN_WITH, H_ORDER> > >
		>::type
	{
	};

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... RES>
	struct gesv_insert_mask<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<>, type_sequence<RES...> >:
		public type_sequence<RES..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER> >
	{
	};


	template <typename T, typename STM, typename STV, typename VISI, typename M_FREE, typename COMMON, typename V_FREE, bool R_CONT, bool M_CONT, bool V_CONT>
	struct gesv_by_mask;

//  M_FREE
	template <typename T, typename STM, typename STV, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... M_FREE, typename ... COMMON, typename ... V_FREE, bool R_CONT, bool M_CONT, bool V_CONT>
	struct gesv_by_mask<T, STM, STV, type_sequence<type_sequence<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<M_FREE...>, type_sequence<COMMON...>, type_sequence<V_FREE...>, R_CONT, M_CONT, V_CONT>:
		public gesv_by_mask<T, STM, STV, type_sequence<type_sequence<VI...>, SI>, typename gesv_insert_mask<valence_info<V_TYPE, 1, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<M_FREE...>, type_sequence<> >::type, type_sequence<COMMON...>, type_sequence<V_FREE...>, R_CONT, M_CONT || (C_MASK & 1), V_CONT>
	{
		static_assert(V_MASK==1,"gesv: Free matrix dimension cannot be optimized. Use defaultIndex or segmentIndex only");
	};
//  COMMON
	template <typename T, typename STM, typename STV, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... M_FREE, typename ... COMMON, typename ... V_FREE, bool R_CONT, bool M_CONT, bool V_CONT>
	struct gesv_by_mask<T, STM, STV, type_sequence<type_sequence<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<M_FREE...>, type_sequence<COMMON...>, type_sequence<V_FREE...>, R_CONT, M_CONT, V_CONT>:
		public gesv_by_mask<T, STM, STV, type_sequence<type_sequence<VI...>, SI>, type_sequence<M_FREE...>, typename gesv_insert_mask<valence_info<V_TYPE, 3, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<COMMON...>, type_sequence<> >::type, type_sequence<V_FREE...>, R_CONT || ((V_MASK==3 && (C_MASK & 1))), M_CONT || (C_MASK & 1), V_CONT || (C_MASK & 2)>
	{
		static_assert(V_MASK & 1,"gesv: Common matrix dimension on matrix side cannot be optimized. Use defaultIndex or segmentIndex only");
		static_assert(V_MASK & 2,"gesv: Common matrix dimension on argument side cannot be optimized. Use defaultIndex or segmentIndex only");
	};
//	V_FREE
	template <typename T, typename STM, typename STV, int V_TYPE, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VI, typename SI, typename ... M_FREE, typename ... COMMON, typename ... V_FREE, bool R_CONT, bool M_CONT, bool V_CONT>
	struct gesv_by_mask<T, STM, STV, type_sequence<type_sequence<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VI...>, SI>, type_sequence<M_FREE...>, type_sequence<COMMON...>, type_sequence<V_FREE...>, R_CONT, M_CONT, V_CONT>:
		public gesv_by_mask<T, STM, STV, type_sequence<type_sequence<VI...>, SI>, type_sequence<M_FREE...>, type_sequence<COMMON...>, typename gesv_insert_mask<valence_info<V_TYPE, 2, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, type_sequence<V_FREE...>, type_sequence<> >::type, R_CONT, M_CONT, V_CONT || (C_MASK & 2)>
	{
		static_assert(V_MASK==2,"gesv: Free argument dimension cannot be optimized. Use defaultIndex or segmentIndex only");
	};

//  DONE MATRIX
	template <typename T, typename STM, typename STV, typename SI, typename ... M_FREE, typename ... COMMON, typename ... V_FREE, bool R_CONT, bool M_CONT, bool V_CONT>
	struct gesv_by_mask<T, STM, STV, type_sequence<type_sequence<>, SI>, type_sequence<M_FREE...>, type_sequence<COMMON...>, type_sequence<V_FREE...>, R_CONT, M_CONT, V_CONT>:
		public gesv_matrix<T, STM, STV, SI, type_sequence<M_FREE...>, type_sequence<COMMON...>, type_sequence<V_FREE...>, R_CONT>
	{
		static_assert(M_CONT,"gesv: One of dimensions of Matrix of coefficients should be continuous");
		static_assert(sizeof...(M_FREE)==1,"gesv: Number of free indices on Matrix coefficient side should be 1");
		static_assert(sizeof...(V_FREE)==1,"gesv: Number of free indices on argument side cannot be more than 1");
		static_assert(sizeof...(COMMON)==1,"gesv: Number of common indices should be 1");
	};

//  DONE VECTOR
	template <typename T, typename STM, typename STV, typename SI, typename ... M_FREE, typename ... COMMON, bool R_CONT, bool M_CONT, bool V_CONT>
	struct gesv_by_mask<T, STM, STV, type_sequence<type_sequence<>, SI>, type_sequence<M_FREE...>, type_sequence<COMMON...>, type_sequence<>, R_CONT, M_CONT, V_CONT>:
		public gesv_vector<T, STM, STV, SI, type_sequence<M_FREE...>, type_sequence<COMMON...>, R_CONT>
	{
		static_assert(M_CONT,"gesv: One of dimensions of Matrix of coefficients should be continuous");
		static_assert(V_CONT,"gesv: One of dimensions of arguments should be continuous");
		static_assert(sizeof...(M_FREE)==1,"gesv: Number of free indices on Matrix coefficient side should be 1");
		static_assert(sizeof...(COMMON)==1,"gesv: Number of common indices should be 1");
	};


	template <typename T, typename STM, typename STV>
	struct gesv_runner;

	template <typename T, size_t M_SNUM, size_t M_CONT, bool M_IS_INDEXED, typename M_OST, typename ... M_SHAPES, size_t V_SNUM, size_t V_CONT, bool V_IS_INDEXED, typename V_OST, typename ... V_SHAPES>
	struct gesv_runner<T, stuple<M_SNUM, M_CONT, M_IS_INDEXED, M_OST, M_SHAPES...>, stuple<V_SNUM, V_CONT, V_IS_INDEXED, V_OST, V_SHAPES...> >:
		public gesv_by_mask<T, stuple<M_SNUM, M_CONT, M_IS_INDEXED, M_OST, M_SHAPES...>, stuple<V_SNUM, V_CONT, V_IS_INDEXED, V_OST, V_SHAPES...>,
			typename make_valence_info<stuple<M_SNUM, M_CONT, M_IS_INDEXED, M_OST, M_SHAPES...>, stuple<V_SNUM, V_CONT, V_IS_INDEXED, V_OST, V_SHAPES...> >::type,
			type_sequence<>, type_sequence<>, type_sequence<>, false, false, false>
	{
		static_assert(M_SNUM==2, "Trying to call gesv with non-matrix tensor (dimension≠2) as a linear system coefficients");
		static_assert(V_SNUM==2 || V_SNUM==1, "Trying to call gesv with non-vector and non-matrix parameter (dimension≠1 and dimension≠2)");
		static_assert(M_IS_INDEXED,"Matrix should be indexed to call gesv");
		static_assert(V_IS_INDEXED,"Argument should be indexed to call gesv");
	};
};

#endif /* GESV_H_ */
