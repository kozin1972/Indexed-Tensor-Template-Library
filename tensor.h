/*
 * tensor.h
 *
 *  Created on: 29 апр. 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef TENSOR_H_
#define TENSOR_H_

#include <tvalence.h>
#include <tgem.h>
#include <tcopy.h>
#include <taxpy.h>
#include <tdiv.h>
#include <tgesv.h>
#include <tscal.h>
#include <texpress.h>
#include <dim_init.h>

namespace tpp
{

	template <typename T, size_t SNUM>
	struct initializer_list_type
	{
		typedef std::initializer_list<typename initializer_list_type<T, SNUM-1>::type> type;
	};

	template <typename T>
	struct initializer_list_type<T,0>
	{
		typedef T type;
	};

	template <typename T, typename ST>
	class tensor;
	template <typename T, typename ST, typename INDEX_ORDER>
	class ind_tensor;
	template <typename T, typename ST>
	class LU;

	template <size_t INUM, typename BPP_TYPE>
	struct shapes_iterator;

	template <size_t INUM, typename T, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_iterator<INUM, tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> > >
//	:
//		public shapes_iterator_base<T, position_number<INUM,stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >::pos, typename std::tuple_element<position_number<INUM,stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >::pos, std::tuple<SHAPES...> >::type, check_shape<typename std::tuple_element<position_number<INUM,stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >::pos, std::tuple<SHAPES...> >::type>::need_parent, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >
	{
		static const size_t POS=position_number<INUM,stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >::pos;
//		static const size_t POS=position_number<INUM,stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >();
		using type=shapes_iterator_base<T, POS, (CONT>0 && POS==sizeof...(SHAPES)-1), typename std::tuple_element<POS, std::tuple<SHAPES...> >::type, check_shape<typename std::tuple_element<POS, std::tuple<SHAPES...> >::type>::need_parent, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >;
//		inline shapes_iterator<INUM, tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> > >	(const stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>& st, const T *p):sib_type(st,p){}
//		inline shapes_iterator<INUM, tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> > >	(const tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >& tensor_v, T *&p):sib_type(tensor_v,p){}
	};

	template <typename T, typename SHAPES>
	class base_tensor;
	template <typename T, typename SHAPES>
	class non_trivial;

	template <typename T, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	class base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >:
		public stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>
	{
		template <typename,typename> friend class base_tensor;
		template <typename,typename> friend class tensor;
		template <typename,typename> friend class LU;
		template <typename, typename, typename, typename, size_t, size_t, typename, typename ...> friend struct index_shift_applyer;
		template <size_t, typename> friend struct shapes_iterator;
		T *data_handle;
		T *data;
		T *allocate_new(const size_t (&aw)[SNUM])
		{
			return shared_array<T>::alloc(segment_shapes_initializer<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >::init(*this,aw));
		}
		base_tensor():ST(),data_handle(NULL),data(NULL) {}
		template <typename OLD_ST, typename ... Ts>
		base_tensor(const tensor<T, OLD_ST>& src, const std::tuple<Ts...>& ind):ST(static_cast<OLD_ST>(src),ind,typename std::integral_constant<bool, defaultIndexVariant<Ts...>::defaultIndices>::type()),data_handle(src.data_handle),data(src.data) { shared_array<T>::share(data_handle); }
		base_tensor(const base_tensor& src):ST(src),data_handle(src.data_handle), data(src.data) { shared_array<T>::share(data_handle); }
		base_tensor(const size_t (&aw)[SNUM]):ST(),data_handle(allocate_new(aw)),data(data_handle) {}
		template <typename PREVST, int ... V_TYPE>
		base_tensor(const PREVST& ost, tpp::integer_sequence<int, V_TYPE...> vn): ST(ost, vn),data_handle(ost.data_handle), data(ost.data) { shared_array<T>::share(data_handle); }
		void move_ptrs(T *data_handle, T *data)
		{
			shared_array<T>::free(this->data_handle);
			this->data_handle=data_handle;
			this->data=data;
			shared_array<T>::share(data_handle);
		}
		~base_tensor() { shared_array<T>::free(data_handle); }
	public:
		using ST=stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>;
		bool is_allocated() const noexcept { return this->data_handle!=0; }
		void free() { shared_array<T>::free(this->data_handle); this->data_handle=NULL; }
		size_t size() const noexcept { return ST::size(); }
//		void get_shape(size_t (&aw)[SNUM]) const noexcept { ST::shape(aw); }
		T *get_step(size_t (&ast)[SNUM]) const { step_getter<SNUM, ST>::step(*this, ast); return data; }
	};

	template <typename T, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	class tensor<T, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >
	{
		template <typename,typename> friend class tensor;
		template <typename,typename> friend class base_tensor;
		template <typename OLD_ST, typename ... Ts>
		tensor(const tensor<T, OLD_ST>& src, const std::tuple<Ts...>& ind):base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >(src,ind){}
		template <typename PREVST, int ... V_TYPE>
		tensor(const PREVST& ost, tpp::integer_sequence<int, V_TYPE...> vn): base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >(ost, vn) {}
//		template <typename ... Ts>
//		inline std::array<size_t,SNUM> get_std_array(Ts...args)
//		{
//			std::array<size_t,SNUM> ret;
//			aw_getter::get_aw(*((size_t (*)[SNUM])ret.data()), args...);
//			return ret;
//		}

		template <bool TrueValue=true>
		typename std::enable_if<(SNUM==1 && TrueValue),tensor&>::type assign(int POS, std::initializer_list<T> list)
		{
			size_t nrows=get<0>(*this).length();
			size_t row=0;
			if (list.size()>nrows)
				throw outOfBounds(list.size()-1,POS,nrows);
			for (auto it=list.begin();it!=list.end() && row<nrows;it++,row++)
				operator()(row)=*it;
			for (;row<nrows;row++)
				operator()(row)=0;
			return *this;
		}

		template <bool TrueValue=true>
		typename std::enable_if<(SNUM>=2 && TrueValue),tensor&>::type assign(size_t POS, typename initializer_list_type<T, SNUM>::type list)
		{
			size_t nrows=get<0>(*this).length();
			size_t row[1]={0};
			if (list.size()>nrows)
				throw outOfBounds(list.size()-1,POS,nrows);
			for (auto it=list.begin();it!=list.end() && row[0]<nrows;it++,row[0]++)
			{
				tensor<T, typename valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, check_shape<typename ST::template element_type<0>>::v_type>::type> res(*this, tpp::integer_sequence<int, check_shape<typename ST::template element_type<0>>::v_type>());
				valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, check_shape<typename ST::template element_type<0>>::v_type>::apply_numeric_index(*this,res.data,row);
				res.assign(POS+1,*it);
			}
			for (;row[0]<nrows;row[0]++)
			{
				tensor<T, typename valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, check_shape<typename ST::template element_type<0>>::v_type>::type> res(*this, tpp::integer_sequence<int, check_shape<typename ST::template element_type<0>>::v_type>());
				valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, check_shape<typename ST::template element_type<0>>::v_type>::apply_numeric_index(*this,res.data,row);
				res=0;
			}
			return *this;
		}
public:
		static const bool is_continuous=(CONT>=sizeof...(SHAPES));
		static const bool is_indexed=IS_INDEXED;
		typedef LU<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> > LU_type;
		tensor() = default;
		tensor(const size_t (&aw)[SNUM]):base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >(aw){}
//		template <size_t ... SIZE>
//		tensor(const m_array<T, SIZE...>&& idata) :base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >({SIZE...})
//		{
//			static_assert(SNUM==sizeof...(SIZE),"Incorrect initializers");
//			BLAS_INTEGER size=this->size();
//			copy<T>(&size, &idata[0], &type_constants<BLAS_INTEGER>::one, this->data, &type_constants<BLAS_INTEGER>::one);
//		}
		template <size_t S0>
		tensor(const T (&idata)[S0]):base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >({S0})
		{
			static_assert(SNUM==1,"Cannot initialize non-vector with 1-d array");
			BLAS_INTEGER size=this->size();
			copy<T>(&size, &idata[0], &type_constants<BLAS_INTEGER>::one, this->data, &type_constants<BLAS_INTEGER>::one);
		}
		template <size_t S0, size_t S1>
		tensor(const T (&idata)[S0][S1]):base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >({S0,S1})
		{
			static_assert(SNUM==2,"Cannot initialize non-matrix with 2-d array");
			BLAS_INTEGER size=this->size();
			copy<T>(&size, &idata[0][0], &type_constants<BLAS_INTEGER>::one, this->data, &type_constants<BLAS_INTEGER>::one);
		}
		template <size_t S0, size_t S1, size_t S2>
		tensor(const T (&idata)[S0][S1][S2]):base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >({S0,S1,S2})
		{
			static_assert(SNUM==3,"Cannot initialize non-3d tensor with 3-d array");
			BLAS_INTEGER size=this->size();
			copy<T>(&size, &idata[0][0][0], &type_constants<BLAS_INTEGER>::one, this->data, &type_constants<BLAS_INTEGER>::one);
		}
		template <size_t S0, size_t S1, size_t S2, size_t S3>
		tensor(const T (&idata)[S0][S1][S2][S3]):base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >({S0,S1,S2,S3})
		{
			static_assert(SNUM==4,"Cannot initialize non-4d tensor with 4-d array");
			BLAS_INTEGER size=this->size();
			copy<T>(&size, &idata[0][0][0][0], &type_constants<BLAS_INTEGER>::one, this->data, &type_constants<BLAS_INTEGER>::one);
		}
//		template <size_t S0, size_t S1, size_t S2, size_t S3, size_t S4>
//		tensor(const T (&idata)[S0][S1][S2][S3][S4]):base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >({S0,S1,S2,S3,S4})
//		{
//			static_assert(SNUM==5,"Cannot initialize non-5d tensor with 5-d array");
//			BLAS_INTEGER size=this->size();
//			copy<T>(&size, &idata[0][0][0][0][0], &type_constants<BLAS_INTEGER>::one, this->data, &type_constants<BLAS_INTEGER>::one);
//		}
//		template <size_t S0, size_t S1, size_t S2, size_t S3, size_t S4, size_t S5>
//		tensor(const T (&idata)[S0][S1][S2][S3][S4][S5]):base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >({S0,S1,S2,S3,S4,S5})
//		{
//			static_assert(SNUM==6,"Cannot initialize non-6d tensor with 6-d array");
//			BLAS_INTEGER size=this->size();
//			copy<T>(&size, &idata[0][0][0][0][0][0], &type_constants<BLAS_INTEGER>::one, this->data, &type_constants<BLAS_INTEGER>::one);
//		}
//		template <size_t S0, size_t S1, size_t S2, size_t S3, size_t S4, size_t S5, size_t S6>
//		tensor(const T (&idata)[S0][S1][S2][S3][S4][S5][S6]):base_tensor<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> >({S0,S1,S2,S3,S4,S5})
//		{
//			static_assert(SNUM==7,"Cannot initialize non-7d tensor with 7-d array");
//			BLAS_INTEGER size=this->size();
//			copy<T>(&size, &idata[0][0][0][0][0][0][0], &type_constants<BLAS_INTEGER>::one, this->data, &type_constants<BLAS_INTEGER>::one);
//		}

		template <bool TrueValue=true>
		typename std::enable_if<(SNUM>=2 && TrueValue),tensor&>::type operator=(typename initializer_list_type<T, SNUM>::type list)
		{
			size_t nrows=get<0>(*this).length();
			size_t row[1]={0};
			if (list.size()>nrows)
				throw outOfBounds(list.size()-1,0,nrows);
			for (auto it=list.begin();it!=list.end() && row[0]<nrows;it++,row[0]++)
			{
				tensor<T, typename valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, check_shape<typename ST::template element_type<0>>::v_type>::type> res(*this, tpp::integer_sequence<int, check_shape<typename ST::template element_type<0>>::v_type>());
				valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, check_shape<typename ST::template element_type<0>>::v_type>::apply_numeric_index(*this,res.data,row);
				res.assign(1,*it);
			}
			for (;row[0]<nrows;row[0]++)
			{
				tensor<T, typename valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, check_shape<typename ST::template element_type<0>>::v_type>::type> res(*this, tpp::integer_sequence<int, check_shape<typename ST::template element_type<0>>::v_type>());
				valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, check_shape<typename ST::template element_type<0>>::v_type>::apply_numeric_index(*this,res.data,row);
				res=0;
			}
			return *this;
		}

		template <bool TrueValue=true>
		typename std::enable_if<(SNUM==1 && TrueValue),tensor&>::type operator=(std::initializer_list<T> list)
		{
			static_assert(SNUM==1,"Trying to initialize non-vector with one-dimansional initilizer list");
			size_t aw[SNUM];
			this->get_shape(aw);
			size_t row=0;
			if (list.size()>aw[0])
				throw outOfBounds(list.size()-1,0,aw[0]);
			for (auto it=list.begin();it!=list.end() && row<aw[0];it++,row++)
			{
				operator()(row)=*it;
			}
			for (;row<aw[0];row++)
				operator()(row)=0;
			return *this;
		}


//		template <typename ... Ts>
//		void init_data(Ts...args)
//		{
//			size_t aw[SNUM];
//			bool data_exists=aw_getter::get_aw(aw, args...);
////			this->data_handle=this->allocate_new(aw);
////			this->data=this->data_handle;
//			if (data_exists)
//			{
//				aw_getter::init_data(aw, this->data, args...);
//			}
//		}
		using ST=stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>;
		T* data_ptr() const noexcept { static_assert((CONT>=sizeof...(SHAPES)),"You cannot access data if is_continuous static data member of tensor is not true"); return this->data; }
		template <typename ... Ts>
		inline tensor<T, typename parenthesis<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, Ts...>::type> operator[](const std::tuple<Ts...>& ind) const
		{
			static_assert(SNUM>0,"Indexing is available only for non-trivial tensors");
			if (!this->is_allocated())
				throw outOfBounds();
			using ret_type=tensor<T, typename parenthesis<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, Ts...>::type>;
			ret_type ret(*this,ind);
			index_shift_applyer<T, ST, typename ret_type::ost_type, std::tuple<Ts...>, 0, 0, type_sequence<SHAPES...>, Ts...>::shift_data(*this,ret.post(),ret.data,ind);
			return ret;
		}
		template <typename ... Ts>
		inline tensor<T, typename parenthesis<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, Ts...>::type> operator()(Ts...inds) const
		{
			static_assert(SNUM>0,"Indexing is available only for non-trivial tensors");
			std::tuple<Ts...> ind(inds...);
			return operator[](ind);
		}
		template <size_t O_NUM, int ... V_TYPE>
		tensor<T, typename valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, V_TYPE...>::type> remove_valences(const size_t (&o)[O_NUM], tpp::integer_sequence<int, V_TYPE...>)
		{
			static_assert(sizeof...(V_TYPE)==O_NUM,"Invalid number of offsets");
			static_assert(IS_INDEXED,"Valence could be removed only from indexed tensors");
			tensor<T, typename valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, V_TYPE...>::type> res(*this, tpp::integer_sequence<int, V_TYPE...>());
//			tensor<T, typename valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, V_TYPE...>::type> res;
			valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, V_TYPE...>::apply_numeric_index(*this,res.data,o);
			return res;
		}
//		template <size_t O_NUM, int ... V_TYPE>
//		tensor<T, typename valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, V_TYPE...>::type> remove_valences(const size_t (&o)[O_NUM])
//		{
//			static_assert(sizeof...(V_TYPE)==O_NUM,"Invalid number of offsets");
//			static_assert(IS_INDEXED,"Valence could be removed only from indexed tensors");
//			tensor<T, typename valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, V_TYPE...>::type> res(*this, tpp::integer_sequence<int, V_TYPE...>());
////			tensor<T, typename valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, V_TYPE...>::type> res;
//			valence_remover<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, V_TYPE...>::apply_numeric_index(*this,res.data,o);
//			return res;
//		}
		tensor& operator=(const tensor& A)
		{
			typedef typename valence_parser<true, ST, ST>::type vd_type;
			check_shape_length<vd_type>(*this, A);
			typedef typename join_valence_data<vd_type>::type jvd_type;
//			test_shape_length<ST, ST>::test(std::tuple<ST, ST>(*this,A));
			typename copy_runner<T, ST, ST, jvd_type>::type runner(*this, A);
			runner.run_copy(this->data, A.data);
			return *this;
		}
		template <typename STA>
		tensor& operator=(const tensor<T, STA>& A)
		{
			static_assert(IS_INDEXED,"Only non-indexed tensor of the same dimension can be copied to non-indexed tensor");
			typedef typename valence_parser<true, ST, STA>::type vd_type;
			check_shape_length<vd_type>(*this, A);
			typedef typename join_valence_data<vd_type>::type jvd_type;
//			typename tseq_element<3,jvd_type>::head jvdv;
//			test_shape_length<ST, STA>::test(std::tuple<ST, STA>(*this,A));
			typename copy_runner<T, ST, STA, jvd_type>::type runner(*this, A);
			runner.run_copy(this->data, A.data);
			return *this;
		}
		template <typename STA>
		tensor& axpy(const tensor<T, STA>& A, T alpha=1.0)
		{
			typedef typename valence_parser<true, ST, STA>::type vd_type;
			check_shape_length<vd_type>(*this, A);
			typedef typename join_valence_data<vd_type>::type jvd_type;
//			test_shape_length<ST, STA>::test(std::tuple<ST, STA>(*this,A));
			typename axpy_runner<T, ST, STA, jvd_type>::type runner(*this, A, alpha);
			runner.run_axpy(this->data, A.data);
			return *this;
		}
		template <typename STA, typename STB>
		tensor& gem(const tensor<T, STA>& A, const tensor<T, STB>& B, T alpha=1.0, T beta=0.0)
		{
			static_assert(IS_INDEXED,"Tensor is not indexed yet. Multiplication is prohibited.");
			typedef typename valence_parser<true, ST, STA, STB>::type vd_type;
			check_shape_length<vd_type>(*this, A, B);
			typedef typename join_valence_data<vd_type>::type jvd_type;
//			test_shape_length<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA, STB>::test(std::tuple<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA, STB>(*this,A,B));
			typename gem_runner<T, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA, STB, jvd_type>::type runner(*this, A, B, alpha, beta);
//			jvd_type jvd;
//			int cm=jvd_type::head::value;
//			size_t s2=tpp::tseq_element<2,jvd_type>::size;
//			size_t s4=tpp::tseq_element<4,jvd_type>::size;
//			int sumat=tpp::tseq_element<2,jvd_type>::head::head::v_type;
//			int sumbt=tpp::tseq_element<4,jvd_type>::head::head::v_type;
			runner.run(this->data, A.data, B.data);
			return *this;
		}
		template <typename STA>
		T dot(const tensor<T, STA>& A) const
		{
			static_assert(IS_INDEXED,"Tensor is not indexed yet. Multiplication is prohibited.");
//			test_shape_length<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>::test(std::tuple<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>(*this,A));
			typedef stuple<0, 0, true, void> VST;
			static VST vst;
			typedef typename valence_parser<true, VST, ST, STA>::type vd_type;
			check_shape_length<vd_type>(vst, *this, A);
			typedef typename join_valence_data<vd_type>::type jvd_type;
			typename gem_runner<T, VST, stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA, jvd_type>::type runner(vst, *this, A, 1.0, 0.0);
			T res=0.0;
			runner.run(&res, this->data, A.data);
			return res;
		}
		T sum() const
		{
			typedef stuple<0, 0, true, void> VST;
			static VST vst;
			typedef typename valence_parser<true, VST, ST>::type vd_type;
			typedef typename join_valence_data<vd_type>::type jvd_type;
			typename copy_runner<T, VST, ST, jvd_type>::type runner(vst, *this);
			T res;
			runner.run_copy(&res, this->data);
			return res;
		}
		T asum() const
		{
			typedef stuple<0, 0, true, void> VST;
			static VST vst;
			typedef typename valence_parser<true, VST, ST>::type vd_type;
			typedef typename join_valence_data<vd_type>::type jvd_type;
			typename asum_runner<T, VST, ST, jvd_type>::type runner(vst, *this);
			T res;
			runner.run_asum(&res, this->data);
			return res;
		}
		template <typename STA>
		tensor& asum(const tensor<T, STA>& A)
		{
			static_assert(IS_INDEXED,"Tensor is not indexed yet. Multiplication is prohibited.");
			typedef typename valence_parser<true, ST, STA>::type vd_type;
			check_shape_length<vd_type>(*this, A);
			typedef typename join_valence_data<vd_type>::type jvd_type;
//			test_shape_length<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>::test(std::tuple<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>(*this,A));
			typename asum_runner<T, ST, STA, jvd_type>::type runner(*this, A);
			runner.run_asum(this->data, A.data);
			return *this;
		}
		template <typename STA>
		tensor& sign(const tensor<T, STA>& A)
		{
			static_assert(IS_INDEXED,"Tensor is not indexed yet. Multiplication is prohibited.");
			typedef typename valence_parser<true, ST, STA>::type vd_type;
			static_assert(tpp::tseq_element<2,vd_type>::size==0, "Free argument index in not allowed");
			check_shape_length<vd_type>(*this, A);
			typedef typename join_valence_data<vd_type>::type jvd_type;
//			using vd_type=typename valence_parser_join<true, ST, STA>::type;
//			test_shape_length<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>::test(std::tuple<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>(*this,A));
			typename div_runner<T, ST, STA, jvd_type>::type runner(*this, A);
//			static_assert(tpp::tseq_element<2,typename div_runner<T, ST, STA>::vd_type>::size==0, "Free argument index in not allowed");
			runner.run_sign(this->data, A.data);
			return *this;
		}
		template <typename STA>
		tensor& div(const tensor<T, STA>& A)
		{
//			static_assert(IS_INDEXED,"Tensor is not indexed yet. Multiplication is prohibited.");
//			test_shape_length<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>::test(std::tuple<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>(*this,A));
//			typename div_runner<T, ST, STA>::type runner(*this, A);
//			static_assert(tpp::tseq_element<2,typename div_runner<T, ST, STA>::vd_type>::size==0, "Free argument index in not allowed");
			static_assert(IS_INDEXED,"Tensor is not indexed yet. Multiplication is prohibited.");
			typedef typename valence_parser<true, ST, STA>::type vd_type;
			static_assert(tpp::tseq_element<2,vd_type>::size==0, "Free argument index in not allowed");
			check_shape_length<vd_type>(*this, A);
			typedef typename join_valence_data<vd_type>::type jvd_type;
//			using vd_type=typename valence_parser_join<true, ST, STA>::type;
//			test_shape_length<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>::test(std::tuple<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>(*this,A));
			typename div_runner<T, ST, STA, jvd_type>::type runner(*this, A);
//			static_assert(tpp::tseq_element<2,typename div_runner<T, ST, STA>::vd_type>::size==0, "Free argument index in not allowed");
			runner.run_div(this->data, A.data);
			return *this;
		}
		template <typename STA>
		tensor& scal(const tensor<T, STA>& A)
		{
//			static_assert(IS_INDEXED,"Tensor is not indexed yet. Multiplication is prohibited.");
//			test_shape_length<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>::test(std::tuple<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>(*this,A));
//			typename div_runner<T, ST, STA>::type runner(*this, A);
//			static_assert(size(typename decltype(get<2>(typename div_runner<T, ST, STA>::vi_by_mask()))::type())==0, "Free argument index in not allowed");
			static_assert(IS_INDEXED,"Tensor is not indexed yet. Multiplication is prohibited.");
			typedef typename valence_parser<true, ST, STA>::type vd_type;
			static_assert(tpp::tseq_element<2,vd_type>::size==0, "Free argument index in not allowed");
			check_shape_length<vd_type>(*this, A);
			typedef typename join_valence_data<vd_type>::type jvd_type;
//			using vd_type=typename valence_parser_join<true, ST, STA>::type;
//			test_shape_length<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>::test(std::tuple<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STA>(*this,A));
			typename div_runner<T, ST, STA, jvd_type>::type runner(*this, A);
//			static_assert(tpp::tseq_element<2,typename div_runner<T, ST, STA>::vd_type>::size==0, "Free argument index in not allowed");
			runner.run_scal(this->data, A.data);
			return *this;
		}
		tensor& scal(T scale)
		{
			typedef typename valence_parser<true, ST>::type vd_type;
			typedef typename join_valence_data<vd_type>::type jvd_type;
			typename scal_runner<T, ST, jvd_type>::type runner(*this, scale);
			runner.run_scal(this->data);
			return *this;
		}
		tensor& shift(T alpha)
		{
			typedef typename valence_parser<true, ST>::type vd_type;
			typedef typename join_valence_data<vd_type>::type jvd_type;
			typename scal_runner<T, ST, jvd_type>::type runner(*this, alpha);
			runner.run_shift(this->data);
			return *this;
		}

// gesv - solve linear equation A(I,J)*X(J)=V(I) or A(I,J)*X(I)=V(J)
// A is square invertible matrix
// A is overwritten after execution
// Initial value of V is overwritten by the result value of X after execution
		template <typename STV>
		BLAS_INTEGER gesv(tensor<T, STV>&& V, BLAS_INTEGER *IPIV=NULL)
		{
			static_assert(ST::snum==2, "Trying to call gesv with non-matrix tensor (number of valences ≠ 2) as a linear system coefficients");
			static_assert(STV::snum==2 || STV::snum==1, "Trying to call gesv with non-vector and non-matrix parameter (number of valences ≠ 1 and number of valences ≠ 2)");
			static_assert(IS_INDEXED,"Matrix is not indexed yet. Call of gesv is prohibited.");
			static_assert(STV::is_indexed, "Argument should be indexed to call gesv");
			typedef typename valence_parser<true, ST, STV>::type vd_type;
			static_assert(tseq_element<1,vd_type>::size==1, "gesv: The matrix should have exactly 1 free valence");
			static_assert(tseq_element<2,vd_type>::size<=1, "gesv: The argument should not have more than 1 free valences");
			static_assert(tseq_element<3,vd_type>::size==1, "gesv: There should be exactly 1 common valence");
			static_assert(tseq_element<1,vd_type>::head::v_mask==1, "gesv: Free matrix dimension cannot be optimized. Use defaultIndex or segmentIndex only");
			static_assert(tseq_element<3,vd_type>::head::v_mask & 1,"gesv: Common matrix dimension on matrix side cannot be optimized. Use defaultIndex or segmentIndex only");
			static_assert(tseq_element<3,vd_type>::head::c_mask & 2,"gesv: Common matrix dimension on argument side should be continuous. Use defaultIndex only");
			static_assert(tseq_element<1,vd_type>::head::c_mask==1 || (tseq_element<3,vd_type>::head::c_mask & 1), "gesv: One of valences of Matrix of coefficients should be continuous");
			check_shape_length<vd_type>(*this, V);
//			test_shape_length<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STV>::test(std::tuple<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...>, STV>(*this,V));
			return gesv_runner<T, ST, STV, vd_type>::run_gesv(*this, V, this->data, V.data, IPIV);
		}
		template <size_t newTotq>
		using reshaped_type=tensor<T, typename initial_shapes<newTotq>::type>;
		template <size_t newTotq>
		reshaped_type<newTotq> reshape(const size_t (&aw)[newTotq]) const
		{
			static_assert((CONT>=sizeof...(SHAPES)),"You cannot reshape a tensor if is_continuous static data member of tensor is not true");
			reshaped_type<newTotq> ret;
			ret.move_ptrs(this->data_handle, this->data);
			size_t rsize=initial_shapes<newTotq>::init(ret,aw);
			if (this->size()<rsize)
			{
				size_t oaw[SNUM];
				this->get_shape(oaw);
				throw reshapeException(oaw,aw);
			}
			return ret;
		}
		template <size_t newTotq>
		reshaped_type<newTotq> reshape(const int (&awi)[newTotq]) const
		{
			static_assert((CONT>=sizeof...(SHAPES)),"You cannot reshape a tensor if is_continuous static data member is not true");
			size_t aw[newTotq];
			for (size_t i=0;i<newTotq;i++) aw[i]=awi[i];
			return reshape(aw);
		}
		template <typename ... Ts>
		using reshaped_type2=tensor<T, typename initial_shapes<sizeof...(Ts)+1>::type>;
		template <typename ... Ts>
		reshaped_type2<Ts...> reshape(size_t w0, Ts... wr) const
		{
			static_assert((CONT>=sizeof...(SHAPES)),"You cannot reshape a tensor if is_continuous static data member is not true");
			static const size_t dim=sizeof...(Ts)+1;
			size_t aw[dim]={w0,wr...};
			return reshape<dim>(aw);
		}
		template <typename ... Ts>
		reshaped_type2<Ts...> reshape(int w0, Ts... wr) const
		{
			static_assert((CONT>=sizeof...(SHAPES)),"You cannot reshape a tensor if is_continuous static data member is not true");
			static const size_t dim=sizeof...(Ts)+1;
			int awi[dim]={w0,wr...};
			size_t aw[dim];
			for (size_t i=0;i<dim;i++) aw[i]=awi[i];
			return reshape<dim>(aw);
		}
		tensor& redim(const size_t (&aw)[SNUM])
		{
			static_assert(!IS_INDEXED,"redim cannot be applied to indexed tensors");
			shared_array<T>::free(this->data_handle);
			this->data_handle=shared_array<T>::alloc(initial_shapes<SNUM>::init(*this,aw));
			this->data=this->data_handle;
			return *this;
		}
		tensor& redim(const int (&awi)[SNUM])
		{
			static_assert(!IS_INDEXED,"redim cannot be applied to indexed tensors");
			const size_t (&aw)[SNUM];
			for (size_t i=0;i<SNUM;i++) aw[i]=awi[i];
			shared_array<T>::free(this->data_handle);
			this->data_handle=shared_array<T>::alloc(initial_shapes<SNUM>::init(*this,aw));
			this->data=this->data_handle;
			return *this;
		}
		template <typename ... Ts>
		tensor& redim(size_t w0, Ts... wr)
		{
			static_assert(!IS_INDEXED,"redim cannot be applied to indexed tensors");
			static const size_t dim=sizeof...(wr)+1;
			static_assert(dim==SNUM,"Number of initializers of dimensions is not equal to number of dimensions");
			size_t aw[SNUM]={w0,wr...};
			return redim(aw);
		}
		template <typename ... Ts>
		tensor& redim(int w0, Ts... wr)
		{
			static_assert(!IS_INDEXED,"redim cannot be applied to indexed tensors");
			static const size_t dim=sizeof...(wr)+1;
			static_assert(dim==SNUM,"Number of initializers of dimensions is not equal to number of dimensions");
			int awi[SNUM]={w0,wr...};
			size_t aw[SNUM];
			for (size_t i=0;i<SNUM;i++) aw[i]=awi[i];
			return redim(aw);
		}
		const T& value() const { static_assert(SNUM==0,"Data element is not specified completely"); return *this->data; }
		T& value() { static_assert(SNUM==0,"Data element is not specified completely"); return *this->data; }
		operator const T&() const { static_assert(SNUM==0,"Data element is not specified completely"); return *this->data; }
		operator T&() { static_assert(SNUM==0,"Data element is not specified completely"); return *this->data; }
		inline tensor& operator=(const T& v)
		{
			if (SNUM==0)
				*this->data=v;
			else
			{
				typedef stuple<0, 0, true, void> VST;
				static VST vst;
				typedef typename valence_parser<true, ST, VST>::type vd_type;
				typedef typename join_valence_data<vd_type>::type jvd_type;
				typename copy_runner<T, ST, VST, jvd_type>::type runner(*this, vst);
				runner.run_copy(this->data, &v);
			}
			return *this;
		}
		LU_type lu() { return LU_type(*this); }
		tensor<T, typename stuple_like<ST>::type> empty_like()
		{
			size_t shapes[SNUM];
			this->get_shape(shapes);
			tensor<T, typename stuple_like<ST>::type> tmp(shapes);
			typename defaultIndices<ST>::type di{};
			return tmp[di];
		}
		template <typename STA, typename STB>
		tensor& operator=(const number_tt_c<T, STA, STB>& nt) { return gem(nt.t0, nt.t1, nt.number, 0.0); }
		template <typename STA, typename STB>
		tensor& operator+=(const number_tt_c<T, STA, STB>& nt) { return gem(nt.t0, nt.t1, nt.number, 1.0); }
		template <typename STA, typename STB>
		tensor& operator-=(const number_tt_c<T, STA, STB>& nt) { return gem(nt.t0, nt.t1, -nt.number, 1.0); }
		template <typename STA>
		tensor& operator=(const number_tensor_c<T,STA>& nt) { operator=(nt.t); return scal(nt.number); }
//		template <typename STA>
//		tensor& operator=(number_tensor<T,STA>&& nt) { operator=(nt.t); return scal(nt.number); }
		template <typename STA>
		tensor& operator+=(const number_tensor_c<T,STA>& nt) { return axpy(nt.t,nt.number); }
//		template <typename STA>
//		tensor& operator+=(number_tensor<T,STA>&& nt) { return axpy(nt.t,nt.number); }
		template <typename STA>
		tensor& operator+=(const tensor<T,STA>& t) { return axpy(t,1.0); }
		template <typename STA>
		tensor& operator+=(tensor<T,STA>&& t) { return axpy(t,1.0); }
		template <typename STA>
		tensor& operator-=(const number_tensor_c<T,STA>& nt) { return axpy(nt.t,-nt.number); }
//		template <typename STA>
//		tensor& operator-=(number_tensor<T,STA>&& nt) { return axpy(nt.t,-nt.number); }
		template <typename STA>
		tensor& operator-=(const tensor<T,STA>& t) { return axpy(t,-1.0); }
		template <typename STA>
		tensor& operator-=(tensor<T,STA>&& t) { return axpy(t,-1.0); }
		tensor& operator*=(T number) { return scal(number); }
		tensor& operator*=(int number) { return scal(number); }
		tensor& operator/=(T number) { return scal(1.0/number); }
		tensor& operator/=(int number) { return scal(1.0/(T)number); }
		tensor& operator+=(T number) { return shift(number); }
		tensor& operator-=(T number) { return shift(-number); }
		tensor& operator+=(int number) { return shift(number); }
		tensor& operator-=(int number) { return shift(-number); }
		template <typename T2>
		tensor& operator=(const sum_num_struct<T, T2>& sn) { operator=(sn.t2); if (sn.multiplier!=1.0) scal(sn.multiplier); return operator+=(sn.summand); }
		template <typename T2>
		tensor& axpy(const number_tensor_c<T, T2>& nt, T number) { return axpy(nt.t, nt.number*number); }
		template <typename T2>
		tensor& axpy(const sum_num_struct<T, T2>& sn, T number) { axpy(sn.t2, sn.multiplier*number); return operator+=(sn.summand*number); }
		template <typename T1, typename T2>
		tensor& axpy(const sum_struct<T, T1, T2>& ss, T number) { axpy(ss.t1, ss.m1*number); return axpy(ss.t2, ss.m2*number); }
		template <typename T2>
		tensor& operator+=(const sum_num_struct<T, T2>& sn) { axpy(sn.t2, sn.multiplier); return operator+=(sn.summand); }
		template <typename T2>
		tensor& operator-=(const sum_num_struct<T, T2>& sn) { axpy(sn.t2, -sn.multiplier); return operator-=(sn.summand); }
		template <typename T1, typename T2>
		tensor& operator=(const sum_struct<T, T1, T2>& ss) { operator=(ss.t1); if (ss.m1!=1.0) scal(ss.m1); return axpy(ss.t2, ss.m2); }
		template <typename T1, typename T2>
		tensor& operator+=(const sum_struct<T, T1, T2>& ss) { axpy(ss.t1, ss.m1); return axpy(ss.t2, ss.m2); }
		template <typename T1, typename T2>
		tensor& operator-=(const sum_struct<T, T1, T2>& ss) { axpy(ss.t1, -ss.m1); return axpy(ss.t2, -ss.m2); }
		template <int ... V_ORDER>
		ind_tensor<T, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, typename get_index_order<tpp::integer_sequence<int, V_ORDER...>, type_sequence<SHAPES...> >::type> order_indices() const
		{
			static_assert(IS_INDEXED,"Order indices nakes sense only if tensor is indexed");
			static_assert(sizeof...(V_ORDER)==SNUM,"Count of valences is not equal to number of dimensions");
			static_assert(is_unique(tpp::integer_sequence<int, V_ORDER...>()),"Valences must be unique");
			return static_cast<ind_tensor<T, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, typename get_index_order<tpp::integer_sequence<int, V_ORDER...>, type_sequence<SHAPES...> >::type> >(*this);
		}
	};

	template <typename T, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, typename INDEX_ORDER>
	class ind_tensor<T, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, INDEX_ORDER>: public tensor<T, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >
	{
	public:
		inline ind_tensor(const tensor<T, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >& t): tensor<T, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >(t) {}
		inline ind_tensor(tensor<T, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >&& t): tensor<T, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >(t) {}
		template <typename ... Ts>
		inline tensor<T, typename get_ordered_tuple<INDEX_ORDER, std::tuple<Ts...> >::template parenthesis<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> > > operator()(const std::tuple<Ts...>& ind) const
		{
			static_assert(SNUM>0,"Indexing is available only for non-trivial tensors");
			return tensor<T, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >::operator[](get_ordered_tuple<INDEX_ORDER, std::tuple<Ts...> >::template create_tuple<>(ind));
		}
		template <typename ... Ts>
		inline tensor<T, typename get_ordered_tuple<INDEX_ORDER, std::tuple<Ts...> >::template parenthesis<stuple<SNUM,CONT, IS_INDEXED, OST, SHAPES...> > > operator()(Ts...inds) const
		{
			static_assert(SNUM>0,"Indexing is available only for non-trivial tensors");
			std::tuple<Ts...> ind(inds...);
			return tensor<T, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >::operator[](get_ordered_tuple<INDEX_ORDER, std::tuple<Ts...> >::template create_tuple<>(ind));
		}
	};

	template <typename T, typename ST>
	class LU
	{
		typedef typename stuple_like<ST>::type NST;
		tensor<T, NST> lu;
		T *TMP_VECTOR;
		BLAS_INTEGER *IPIV;
		BLAS_INTEGER INFO;
//		static size_t BLAS_int_qty(size_t sz) {	return (sizeof(T)*sz+sizeof(BLAS_INTEGER)-1)/sizeof(BLAS_INTEGER); }
//		static size_t new_sh0(size_t sz, size_t sh1) { return (sizeof(BLAS_INTEGER)*(BLAS_int_qty(sz)+sh1)+sizeof(T)*sh1-1)/(sizeof(T)*sh1); }
		static size_t T_qty(size_t blas_qty) { return (sizeof(BLAS_INTEGER)*blas_qty+sizeof(T)-1)/sizeof(T); }
		static size_t T_rows(size_t blas_qty, size_t row_size) { return (T_qty(blas_qty)+row_size-1)/row_size; }
		static tensor<T, NST> get_lu(const tensor<T, ST>& M)
		{
			size_t shapes[2];
			M.get_shape(shapes);
			size_t sh0=shapes[0];
			shapes[0]=sh0+1+T_rows(shapes[1], shapes[1]);
			tensor<T, typename initial_shapes<2>::type> tmp(shapes);
			get<0>(tmp).init(sh0,shapes[1]);
			typename defaultIndices<ST>::type di{};
			tmp[di]=M;
			return tmp[di];
		}
//		static BLAS_INTEGER *get_IPIV(const tensor<T, ST>& M, const tensor<T, NST>& lu)
//		{
//
//			size_t shapes[2];
//			M.get_shape(shapes);
//			return ((BLAS_INTEGER *)(lu.data_ptr()+(shapes[0]+1)*shapes[1])); //+BLAS_int_qty((shapes[0]+1)*shapes[1]);
//		}
	public:
		LU(const tensor<T, ST>& M): lu(get_lu(M)), TMP_VECTOR(lu.data_ptr()+M.size()), IPIV((BLAS_INTEGER *)(TMP_VECTOR+tpp::get<1>(M).length())), INFO(0)
		{
			static_assert(ST::snum==2, "LU factorization is available only for matrices");
			static_assert(NST::is_indexed,"Matrix is not indexed yet. Call of gesv is prohibited.");
//			typedef typename valence_parser<true, NST>::type vd_type;
//			static_assert(tseq_element<1,vd_type>::head::v_mask==1, "gesv: Free matrix dimension cannot be optimized. Use defaultIndex or segmentIndex only");
			size_t shapes[2];
			M.get_shape(shapes);
			BLAS_INTEGER s1=shapes[1];
			BLAS_INTEGER s0=shapes[0];
			getrf<T>(&s1, &s0, lu.data_ptr(), &s1, IPIV, &INFO);
		}
		template <typename STV>
		BLAS_INTEGER solve(tensor<T, STV>&& V)
		{
			static_assert(NST::snum==2, "Trying to call gesv with non-matrix tensor (number of valences ≠ 2) as a linear system coefficients");
			static_assert(STV::snum==2 || STV::snum==1, "Trying to call gesv with non-vector and non-matrix parameter (number of valences ≠ 1 and number of valences ≠ 2)");
			static_assert(STV::is_indexed, "Argument should be indexed to call gesv");
			typedef typename valence_parser<true, NST, STV>::type vd_type;
			static_assert(tseq_element<1,vd_type>::size==1, "gesv: The matrix should have exactly 1 free valence");
			static_assert(tseq_element<2,vd_type>::size<=1, "gesv: The argument should not have more than 1 free valences");
			static_assert(tseq_element<3,vd_type>::size==1, "gesv: There should be exactly 1 common valence");
			static_assert(tseq_element<3,vd_type>::head::v_mask & 1,"gesv: Common valence on matrix side cannot be optimized. Use defaultIndex or segmentIndex only");
//			static_assert(tseq_element<3,vd_type>::head::c_mask & 2,"gesv: Common matrix dimension on argument side should be continuous. Use defaultIndex only");
//			static_assert(tseq_element<1,vd_type>::head::c_mask==1 || (tseq_element<3,vd_type>::head::c_mask & 1), "gesv: One of valences of Matrix of coefficients should be continuous");
//			static_assert(tseq_element<3,vd_type>::head::v_mask & 2,"gesv: Common valence on argument side cannot be optimized. Use defaultIndex or segmentIndex only. defaultIndex is faster");
			check_shape_length<vd_type>(lu, V);


//			test_shape_length<NST, STV>::test(std::tuple<NST, STV>(lu,V));
//			if (tseq_element<3,vd_type>::head::c_mask & 2)
//				return gesv_runner<T, NST, STV, vd_type>::run_solve(lu, V, lu.data_ptr(), V.data, IPIV);
//			else
			return gesv_runner<T, NST, STV, vd_type>::run_solve(lu, V, lu.data_ptr(), V.data, TMP_VECTOR, IPIV);
		}
		BLAS_INTEGER info() { return INFO; }
	};

	template <size_t DIM, typename T=double>
	using TENSOR=tensor<T, typename initial_shapes<DIM>::type>;
	template <typename T=double>
	using VECTOR=tensor<T, typename initial_shapes<1>::type>;
	template <typename T=double>
	using MATRIX=tensor<T, typename initial_shapes<2>::type>;
};


#endif /* TENSOR_H_ */
