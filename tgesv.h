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

	template <typename T, typename STM, typename STV, typename M_FREE, typename V_FREE, typename COMMON, int COMARG, int VVM2, bool NORMAL>
	struct gesv_matrix
	{
		static_assert(true, "Internal error during instantiation of gesv_matrix");
	};

	template <typename T, typename STM, typename STV,
		int M_V_TYPE, int M_C_MASK, int M_NEXT_VALENCE, size_t M_POS,
		int V_V_TYPE, int V_C_MASK, int V_V_MASK, int V_NEXT_VALENCE, size_t V_POS,
		int C_V_TYPE, int C_C_MASK, int C_V_MASK, int C_NEXT_VALENCE, size_t CM_POS, size_t CV_POS,
		bool NORMAL>
	struct gesv_matrix<T, STM, STV,
		valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>,
		valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>,
		valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 3, 2,
		NORMAL>
	{
		static BLAS_INTEGER get_LDA(const STM& stm)
		{
			return NORMAL?
					vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::step(stm)
					:vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::step(stm);
		}
		static const char *TRANS() { return NORMAL?"N":"T"; }

		static BLAS_INTEGER run_solve(const STM& stm, const STV& stv, T *M, T *V, T *tmp, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::length(stm);
			if ((size_t)N!=vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=get_LDA(stm);
			BLAS_INTEGER LDB=vd_iterator_getter<STV, valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>, 1>::template iterator_type<T>::step(stv);
			BLAS_INTEGER NRHS=vd_iterator_getter<STV, valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>, 1>::template iterator_type<T>::length(stv);
			BLAS_INTEGER INFO;
			getrs<T>(TRANS(), &N, &NRHS, M, &LDA, IPIV, V, &LDB, &INFO);
			return INFO;
		}
	};

	template <typename T, typename STM, typename STV,
		int M_V_TYPE, int M_C_MASK, int M_NEXT_VALENCE, size_t M_POS,
		int V_V_TYPE, int V_C_MASK, int V_V_MASK, int V_NEXT_VALENCE, size_t V_POS,
		int C_V_TYPE, int C_C_MASK, int C_V_MASK, int C_NEXT_VALENCE, size_t CM_POS, size_t CV_POS,
		bool NORMAL>
	struct gesv_matrix<T, STM, STV,
		valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>,
		valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>,
		valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 3, 0,
		NORMAL>
	{
		static BLAS_INTEGER get_LDA(const STM& stm)
		{
			return NORMAL?
					vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::step(stm)
					:vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::step(stm);
		}
		static const char *TRANS() { return NORMAL?"N":"T"; }

		static BLAS_INTEGER run_solve(const STM& stm, const STV& stv, T *M, T *V, T *tmp, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::length(stm);
			if ((size_t)N!=vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=get_LDA(stm);
			BLAS_INTEGER INFO;
			typedef typename vd_iterator_getter<STV, valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>, 1>::template iterator_type<T> RHS_it;
			for (RHS_it it(stv, V);it.not_end(V);it.move_one(V))
			{
				getrs<T>(TRANS(), &N, &type_constants<BLAS_INTEGER>::one, M, &LDA, IPIV, V, &N, &INFO);
			}
			return INFO;
		}
	};

	template <typename T, typename STM, typename STV,
		int M_V_TYPE, int M_C_MASK, int M_NEXT_VALENCE, size_t M_POS,
		int V_V_TYPE, int V_C_MASK, int V_V_MASK, int V_NEXT_VALENCE, size_t V_POS,
		int C_V_TYPE, int C_C_MASK, int C_V_MASK, int C_NEXT_VALENCE, size_t CM_POS, size_t CV_POS, int VVM2,
		bool NORMAL>
	struct gesv_matrix<T, STM, STV,
		valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>,
		valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>,
		valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 1, VVM2,
		NORMAL>
	{
		static BLAS_INTEGER get_LDA(const STM& stm)
		{
			return NORMAL?
					vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::step(stm)
					:vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::step(stm);
		}
		static const char *TRANS() { return NORMAL?"N":"T"; }
		static BLAS_INTEGER run_solve(const STM& stm, const STV& stv, T *M, T *V, T *tmp, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::length(stm);
			if ((size_t)N!=vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=get_LDA(stm);
//			BLAS_INTEGER LDB=vd_iterator_getter<STV, valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>, 1>::template iterator_type<T>::step(stv);
			BLAS_INTEGER LDB=vd_iterator_getter<STV, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 1>::template iterator_type<T>::step(stv);
//			BLAS_INTEGER NRHS=vd_iterator_getter<STV, valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>, 1>::template iterator_type<T>::length(stv);
			BLAS_INTEGER INFO;
			typedef typename vd_iterator_getter<STV, valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>, 1>::template iterator_type<T> RHS_it;
			for (RHS_it it(stv, V);it.not_end(V);it.move_one(V))
//			for (BLAS_INTEGER i=0;i<NRHS;i++)
			{
				copy<T>(&N, V, &LDB, tmp, &type_constants<BLAS_INTEGER>::one);
				getrs<T>(TRANS(), &N, &type_constants<BLAS_INTEGER>::one, M, &LDA, IPIV, tmp, &N, &INFO);
				if (INFO)
					return INFO;
				copy<T>(&N, tmp, &type_constants<BLAS_INTEGER>::one, V, &LDB);
			}
			return INFO;
		}
	};

	template <typename T, typename STM, typename STV,
		int M_V_TYPE, int M_C_MASK, int M_NEXT_VALENCE, size_t M_POS,
		int V_V_TYPE, int V_C_MASK, int V_V_MASK, int V_NEXT_VALENCE, size_t V_POS,
		int C_V_TYPE, int C_C_MASK, int C_V_MASK, int C_NEXT_VALENCE, size_t CM_POS, size_t CV_POS, int VVM2,
		bool NORMAL>
	struct gesv_matrix<T, STM, STV,
		valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>,
		valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>,
		valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0, VVM2,
		NORMAL>
	{
		static BLAS_INTEGER get_LDA(const STM& stm)
		{
			return NORMAL?
					vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::step(stm)
					:vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::step(stm);
		}
		static const char *TRANS() { return NORMAL?"N":"T"; }
		static BLAS_INTEGER run_solve(const STM& stm, const STV& stv, T *M, T *V, T *tmp, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::length(stm);
			if ((size_t)N!=vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=get_LDA(stm);
			BLAS_INTEGER INFO;
			typedef typename vd_iterator_getter<STV, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 1>::template iterator_type<T> V_it;
			typedef typename vd_iterator_getter<STV, valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>, 1>::template iterator_type<T> RHS_it;
			for (RHS_it it(stv, V);it.not_end(V);it.move_one(V))
			{
				T *VC=V;
				V_it v_it(stv, VC);
				for (size_t i=0;v_it.not_end(VC);i++,v_it.move_one(VC))
				{
					tmp[i]=*VC;
				}
				getrs<T>(TRANS(), &N, &type_constants<BLAS_INTEGER>::one, M, &LDA, IPIV, tmp, &N, &INFO);
				if (INFO)
					return INFO;
				VC=V;
				V_it v_itb(stv, VC);
				for (size_t i=0;v_itb.not_end(VC);i++,v_itb.move_one(VC))
				{
					*VC=tmp[i];
				}
			}
			return INFO;
		}
	};

	template <typename T, typename STM, typename STV, typename M_FREE, typename COMMON, int COMARG, bool NORMAL>
	struct gesv_vector
	{
		static_assert(true, "Internal error during instantiation of gesv_vector");
	};

	template <typename T, typename STM, typename STV,
		int M_V_TYPE, int M_C_MASK, int M_NEXT_VALENCE, size_t M_POS,
		int C_V_TYPE, int C_C_MASK, int C_V_MASK, int C_NEXT_VALENCE, size_t CM_POS, size_t CV_POS,
		bool NORMAL>
	struct gesv_vector<T, STM, STV,
	valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>,
	valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 3, NORMAL>
	{
		static BLAS_INTEGER get_LDA(const STM& stm)
		{
			return NORMAL?
					vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::step(stm)
					:vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::step(stm);
		}
		static const char *TRANS() { return NORMAL?"N":"T"; }
		static BLAS_INTEGER run_solve(const STM& stm, const STV& stv, T *M, T *V, T *tmp, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::length(stm);
			if ((size_t)N!=vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=get_LDA(stm);
			BLAS_INTEGER INFO;
			getrs<T>(TRANS(), &N, &type_constants<BLAS_INTEGER>::one, M, &LDA, IPIV, V, &N, &INFO);
			return INFO;
		}
	};

	template <typename T, typename STM, typename STV,
		int M_V_TYPE, int M_C_MASK, int M_NEXT_VALENCE, size_t M_POS,
		int C_V_TYPE, int C_C_MASK, int C_V_MASK, int C_NEXT_VALENCE, size_t CM_POS, size_t CV_POS,
		bool NORMAL>
	struct gesv_vector<T, STM, STV,
	valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>,
	valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 1, NORMAL>
	{
		static BLAS_INTEGER get_LDA(const STM& stm)
		{
			return NORMAL?
					vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::step(stm)
					:vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::step(stm);
		}
		static const char *TRANS() { return NORMAL?"N":"T"; }
		static BLAS_INTEGER run_solve(const STM& stm, const STV& stv, T *M, T *V, T *tmp, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::length(stm);
			if ((size_t)N!=vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=get_LDA(stm);
			BLAS_INTEGER LDB=vd_iterator_getter<STV, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 1>::template iterator_type<T>::step(stv);
			BLAS_INTEGER INFO;
			copy<T>(&N, V, &LDB, tmp, &type_constants<BLAS_INTEGER>::one);
			getrs<T>(TRANS(), &N, &type_constants<BLAS_INTEGER>::one, M, &LDA, IPIV, tmp, &N, &INFO);
			copy<T>(&N, tmp, &type_constants<BLAS_INTEGER>::one, V, &LDB);
			return INFO;
		}
	};

	template <typename T, typename STM, typename STV,
		int M_V_TYPE, int M_C_MASK, int M_NEXT_VALENCE, size_t M_POS,
		int C_V_TYPE, int C_C_MASK, int C_V_MASK, int C_NEXT_VALENCE, size_t CM_POS, size_t CV_POS,
		bool NORMAL>
	struct gesv_vector<T, STM, STV,
	valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>,
	valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0, NORMAL>
	{
		static BLAS_INTEGER get_LDA(const STM& stm)
		{
			return NORMAL?
					vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::step(stm)
					:vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::step(stm);
		}
		static const char *TRANS() { return NORMAL?"N":"T"; }
		static BLAS_INTEGER run_solve(const STM& stm, const STV& stv, T *M, T *V, T *tmp, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=vd_iterator_getter<STM, valence_data<M_V_TYPE, 1, M_C_MASK, 1, M_NEXT_VALENCE, M_POS>, 0>::template iterator_type<T>::length(stm);
			if ((size_t)N!=vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=get_LDA(stm);
			BLAS_INTEGER INFO;
			typedef typename vd_iterator_getter<STV, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 1>::template iterator_type<T> V_it;
			T *VC=V;
			V_it v_it(stv, VC);
			for (size_t i=0;v_it.not_end(VC);i++,v_it.move_one(VC))
			{
				tmp[i]=*VC;
			}
			getrs<T>(TRANS(), &N, &type_constants<BLAS_INTEGER>::one, M, &LDA, IPIV, tmp, &N, &INFO);
			VC=V;
			V_it v_itb(stv, VC);
			for (size_t i=0;v_itb.not_end(VC);i++,v_itb.move_one(VC))
			{
				*VC=tmp[i];
			}
			return INFO;
		}
	};


	template <typename T, typename STM, typename STV, typename VD>
	struct gesv_runner
	{
		static_assert(1,"Internal error");
	};

//  DONE MATRIX
	template <typename T, typename STM, typename STV, int CONT_MASK, typename M_FREE, int C_V_TYPE, int C_C_MASK, int C_V_MASK, int C_NEXT_VALENCE, size_t CM_POS, size_t CV_POS, int V_V_TYPE, int V_C_MASK, int V_V_MASK, int V_NEXT_VALENCE, size_t V_POS>
	struct gesv_runner<T, STM, STV, type_sequence<std::integral_constant<int, CONT_MASK>, type_sequence<M_FREE>, type_sequence<valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS> >, type_sequence<valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS> >, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public gesv_matrix<T, STM, STV, M_FREE, valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, (C_C_MASK & 2) + ((C_V_MASK & 2)>>1), (V_V_MASK & 2), CM_POS==1>
	{
		static BLAS_INTEGER get_LDA(const STM& stm)
		{
			return CM_POS==1?
					vd_iterator_getter<STM, M_FREE, 0>::template iterator_type<T>::step(stm)
					:vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::step(stm);
		}
		static const char *TRANS() { return CM_POS==1?"N":"T"; }
		static BLAS_INTEGER run_gesv(const STM& stm, const STV& stv, T *M, T *V, BLAS_INTEGER *IPIV)
		{
			static_assert(V_V_MASK==2, "Cannot optimize free argument dimension. Use defaultIndex or segmentIndex only");
			BLAS_INTEGER N=vd_iterator_getter<STM, M_FREE, 0>::template iterator_type<T>::length(stm);
			if ((size_t)N!=vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=get_LDA(stm);
			BLAS_INTEGER LDB=vd_iterator_getter<STV, valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>, 1>::template iterator_type<T>::step(stv);
			BLAS_INTEGER NRHS=vd_iterator_getter<STV, valence_data<V_V_TYPE, 2, V_C_MASK, V_V_MASK, V_NEXT_VALENCE, V_POS>, 1>::template iterator_type<T>::length(stv);
			if (!IPIV)
				IPIV=(BLAS_INTEGER *)alloca((N>1?N:1)*sizeof(BLAS_INTEGER));
			BLAS_INTEGER INFO;
			getrf<T>(&N, &N, M, &LDA, IPIV, &INFO);
			if (INFO)
				return INFO;
			getrs<T>(TRANS(), &N, &NRHS, M, &LDA, IPIV, V, &LDB, &INFO);
			return INFO;
		}
	};

//  DONE VECTOR
	template <typename T, typename STM, typename STV, int CONT_MASK, typename M_FREE, int C_V_TYPE, int C_C_MASK, int C_V_MASK, int C_NEXT_VALENCE, size_t CM_POS, size_t CV_POS>
	struct gesv_runner<T, STM, STV, type_sequence<std::integral_constant<int, CONT_MASK>, type_sequence<M_FREE>, type_sequence<>, type_sequence<valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS> >, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<> > >:
		public gesv_vector<T, STM, STV, M_FREE, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, (C_C_MASK & 2) + ((C_V_MASK & 2)>>1), CM_POS==1>
	{
		static BLAS_INTEGER get_LDA(const STM& stm)
		{
			return CM_POS==1?
					vd_iterator_getter<STM, M_FREE, 0>::template iterator_type<T>::step(stm)
					:vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::step(stm);
		}
		static const char *TRANS() { return CM_POS==1?"N":"T"; }
		static BLAS_INTEGER run_gesv(const STM& stm, const STV& stv, T *M, T *V, BLAS_INTEGER *IPIV)
		{
			BLAS_INTEGER N=vd_iterator_getter<STM, M_FREE, 0>::template iterator_type<T>::length(stm);
			if ((size_t)N!=vd_iterator_getter<STM, valence_data<C_V_TYPE, 3, C_C_MASK, C_V_MASK, C_NEXT_VALENCE, CM_POS, CV_POS>, 0>::template iterator_type<T>::length(stm))
			{
				size_t shape[2];
				stm.get_shape(shape);
				throw MatrixIsNotSquareException(shape);
			}
			BLAS_INTEGER LDA=get_LDA(stm);
			if (!IPIV)
				IPIV=(BLAS_INTEGER *)alloca((N>1?N:1)*sizeof(BLAS_INTEGER));
			BLAS_INTEGER INFO;
			getrf<T>(&N, &N, M, &LDA, IPIV, &INFO);
			if (INFO)
				return INFO;
			getrs<T>(TRANS(), &N, &type_constants<BLAS_INTEGER>::one, M, &LDA, IPIV, V, &N, &INFO);
			return INFO;
		}
	};

};

#endif /* GESV_H_ */
