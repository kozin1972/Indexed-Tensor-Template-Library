/*
 * stuple.h
 *
 *  Created on: 25 апр. 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef STUPLE_H_
#define STUPLE_H_

#include <tppindex.h>
#include <tsequence.h>

namespace tpp
{

//	template <typename T, size_t ... SIZE>
//	struct m_array;
//	template <typename T>
//	struct m_array<T>
//	{
//		typedef T type;
//		type data;
//	};
//	template <typename T, size_t S0, size_t ... SIZE>
//	struct m_array<T, S0, SIZE...>
//	{
//		typedef typename m_array<T, SIZE...>::type type[S0];
//		type data;
//	};

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct stuple;

	template <size_t I, typename ... SHAPES>
	struct stuple_base;

	template <std::size_t __i, typename ST>
	struct stuple_element;
	/**
	* Recursive case for tuple_element: strip off the first element in
	* the tuple and retrieve the (i-1)th element of the remaining tuple.
	*/
	template <std::size_t __i, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename _Head, typename... _Tail>
	struct stuple_element<__i, stuple<SNUM, CONT, IS_INDEXED, OST, _Head, _Tail...> >
	: stuple_element<__i - 1, stuple<SNUM, CONT, IS_INDEXED, OST,_Tail...> > { };

	/**
	* Basis case for tuple_element: The first element is the one we're seeking.
	*/
	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename _Head, typename... _Tail>
	struct stuple_element<0, stuple<SNUM, CONT, IS_INDEXED, OST,_Head, _Tail...> >
	{
		typedef _Head type;
	};
    // Duplicate of C++14's tuple_element_t for internal use in C++11 mode
    template<std::size_t __i, typename _Tp>
	    using __stuple_element_t = typename stuple_element<__i, _Tp>::type;

    template<std::size_t __i, typename _Tp>
	struct stuple_element<__i, const _Tp>
	{
	    typedef typename std::add_const<__stuple_element_t<__i, _Tp>>::type type;
	};

    template<std::size_t __i, typename _Tp>
	struct stuple_element<__i, volatile _Tp>
	{
	    typedef typename std::add_volatile<__stuple_element_t<__i, _Tp>>::type type;
	};

    template<std::size_t __i, typename _Tp>
	struct stuple_element<__i, const volatile _Tp>
	{
	    typedef typename std::add_cv<__stuple_element_t<__i, _Tp>>::type type;
	};

#if __cplusplus > 201103L
#define __cpp_lib_tuple_element_t 201402

    template<std::size_t __i, typename _Tp>
	  using stuple_element_t = typename stuple_element<__i, _Tp>::type;
#endif

    template<std::size_t __i, typename _Head, typename... _Tail>
	    constexpr const _Head&
	    __get_helper(const stuple_base<__i, _Head, _Tail...>& __t) noexcept
	    { return stuple_base<__i, _Head, _Tail...>::head(__t); }
    template<std::size_t __i, typename _Head, typename... _Tail>
	    constexpr _Head&
	    __get_helper(stuple_base<__i, _Head, _Tail...>& __t) noexcept
	    { return stuple_base<__i, _Head, _Tail...>::head(__t); }

 /// Return a reference to the ith element of a tuple.
    template<std::size_t __i, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename... _Elements>
	constexpr const __stuple_element_t<__i, stuple<SNUM, CONT, IS_INDEXED, OST, _Elements...>>&
	get(const stuple<SNUM, CONT, IS_INDEXED, OST, _Elements...>& __t) noexcept
	{ return tpp::__get_helper<__i>(__t); }

    template<std::size_t __i, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename... _Elements>
	constexpr __stuple_element_t<__i, stuple<SNUM, CONT, IS_INDEXED, OST, _Elements...>>&
	get(stuple<SNUM, CONT, IS_INDEXED, OST, _Elements...>& __t) noexcept
	{ return tpp::__get_helper<__i>(__t); }


	template <typename ST, typename ... IND>
	struct parenthesis;

//	template <size_t INUM, typename ST>
//	constexpr size_t position_number() { return tpp::get<INUM>(typename ST::active_shape_indices{}); }

	template <size_t INUM, typename SHAPES>
	struct position_number;

	template <size_t INUM, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename HEAD, typename ... SHAPES>
	struct position_number<INUM, stuple<SNUM, CONT, IS_INDEXED, OST, HEAD, SHAPES...> >
	{
		static const size_t pos=position_number<INUM-1, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >::pos+1;
	};

	template <size_t INUM, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, template <int, enum dSU, size_t, size_t> class BSC, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct position_number<INUM, stuple<SNUM, CONT, IS_INDEXED, OST, BSC<V_TYPE, USAGE_SLAVE, USE_ALSO, PARENT>, SHAPES...> >
	{
		static const size_t pos=position_number<INUM, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >::pos+1;
	};

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, template <int, enum dSU, size_t, size_t> class BSC, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct position_number<0, stuple<SNUM, CONT, IS_INDEXED, OST, BSC<V_TYPE, USAGE_SLAVE, USE_ALSO, PARENT>, SHAPES...> >
	{
		static const size_t pos=position_number<0, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >::pos+1;
	};

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename HEAD, typename ... SHAPES>
	struct position_number<0, stuple<SNUM, CONT, IS_INDEXED, OST, HEAD, SHAPES...> >
	{
		static const size_t pos=0;
	};

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST>
	struct position_number<0, stuple<SNUM, CONT, IS_INDEXED, OST> >
	{
		static const size_t pos=0;
	};

	template <size_t POS, typename ... SHAPES>
	struct index_number;
	template <>
	struct index_number<0>
	{
		static const size_t inum=0;
	};
	template <template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct index_number<0, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>
	{
		static const size_t inum=0;
	};
	template <template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct index_number<0, sClass<V_TYPE, USAGE_FULL, USE_ALSO, PARENT>, SHAPES...>
	{
		static const size_t inum=0;
	};
	template <size_t POS, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct index_number<POS, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>
	{
		static const size_t inum=index_number<POS-1, SHAPES...>::inum+1;
	};
	template <size_t POS, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct index_number<POS, sClass<V_TYPE, USAGE_SLAVE, USE_ALSO, PARENT>, SHAPES...>
	{
		static const size_t inum=index_number<POS-1, SHAPES...>::inum;
	};


	template <typename SHAPE> struct length_getter;

	template <template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	struct length_getter<shapeClass<V_TYPE,USAGE,USE_ALSO,PARENT> >
	{
		inline static size_t length(const shapeClass<V_TYPE,USAGE,USE_ALSO,PARENT>& s) { return s.length(); }
	};

	template <template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, size_t USE_ALSO, size_t PARENT>
	struct length_getter<shapeClass<V_TYPE,USAGE_SLAVE,USE_ALSO,PARENT> >
	{
		inline static size_t length(const shapeClass<V_TYPE,USAGE_SLAVE,USE_ALSO,PARENT>& s) noexcept { return 1; }
	};

	template <size_t CUR, typename ST>
	struct size_getter
	{
		inline static size_t size(const ST& st) { return length_getter<typename ST::template element_type<CUR-1> >::length(get<CUR-1>(st))*size_getter<CUR-1,ST>::size(st); }
	};

	template <typename ST>
	struct size_getter<0, ST>
	{
		inline static size_t size(const ST& st) noexcept { return 1; }
	};


	template <size_t CUR, typename ST>
	struct shape_getter
	{
		static const size_t pos=position_number<CUR-1,ST>::pos;
//		static const size_t pos=position_number<CUR-1,ST>();
		inline static void shape(const ST& st, size_t (&aw)[ST::snum]) { aw[CUR-1]=length_getter<typename ST::template element_type<pos> >::length(get<CUR-1>(st)); shape_getter<CUR-1,ST>::shape(st,aw); }
	};

	template <typename ST>
	struct shape_getter<0, ST>
	{
		inline static void shape(const ST& st, size_t (&ast)[ST::snum]) noexcept {}
	};

	template <size_t CUR, typename ST>
	struct step_getter
	{
		static const size_t pos=position_number<CUR-1,ST>();
		inline static void step(const ST& st, size_t (&ast)[ST::snum]) { ast[CUR-1]=get<CUR-1>(st).step(); step_getter<CUR-1,ST>::step(st,ast); }
	};

	template <typename ST>
	struct step_getter<0, ST>
	{
		inline static void step(const ST& st, size_t (&aw)[ST::snum]) noexcept {}
	};

	template <typename ST, int ... V_TYPE>
	struct valence_remover;


	template <size_t I, typename HEAD, typename ... SHAPES>
	struct stuple_base<I, HEAD, SHAPES...>: private HEAD, public stuple_base<I+1, SHAPES...>
	{
		template <size_t, size_t, typename> friend struct stuple_getter;
		typedef HEAD element_type;
		inline static constexpr const element_type& head(const stuple_base& st) noexcept { return st; }
		inline static constexpr element_type& head(stuple_base& st) noexcept { return st; }
		typedef stuple_base<I+1, SHAPES...> base_type;
		inline static constexpr const base_type& base(const stuple_base& st) noexcept { return st; }
		inline static constexpr base_type& base(stuple_base& st) noexcept { return st; }
		typedef typename std::conditional<check_shape<HEAD>::usage==USAGE_SLAVE,typename base_type::active_shape_indices,decltype(push_front<size_t, I>(typename base_type::active_shape_indices{})) >::type active_shape_indices;
		template <typename OST, typename ... Ts>
		inline stuple_base(const OST& ost, const std::tuple<Ts...>& inds):
			element_type(parenthesis<OST, Ts...>::template index<I>(ost, inds),
					parenthesis<OST, Ts...>::template parent<I>(ost, inds)),
			base_type(ost, inds) {}
		template <typename OST, typename ... Ts>
		inline stuple_base(const std::tuple<Ts...>& inds, const OST& ost):
					element_type(get<I>(ost)),
					base_type(inds, ost) {}
		template <typename OST, int ... V_TYPE>
		inline stuple_base(const OST& ost, tpp::integer_sequence<int, V_TYPE...> vn):
					element_type(valence_remover<OST, V_TYPE...>::template shape<I>(ost)),
					base_type(ost, vn) {}
		stuple_base() = default;
	};

	template <size_t I>
	struct stuple_base<I>
	{
		typedef index_sequence<> active_shape_indices;
		template <typename OST, typename ... Ts>
		inline stuple_base(const OST& ost, const std::tuple<Ts...>& inds) noexcept {}
		template <typename OST, typename ... Ts>
		inline stuple_base(const std::tuple<Ts...>& inds, const OST& ost) noexcept {}
		template <typename OST, int ... V_TYPE>
		inline stuple_base(const OST& ost, tpp::integer_sequence<int, V_TYPE...> vn) noexcept {}
		stuple_base() = default;
	};


	template <typename CHK, typename ... IND>
	struct check_same_variance_or_numeric
	{
		static const bool value=true;
	};
	template <int V_TYPE, int V_TYPE2, typename ... IND>
	struct check_same_variance_or_numeric<defaultIndex<V_TYPE>, defaultIndex<V_TYPE2>, IND...>:
		public check_same_variance_or_numeric<defaultIndex<V_TYPE>, IND...>
	{
	};
	template <int V_TYPE, typename ... IND>
	struct check_same_variance_or_numeric<defaultIndex<V_TYPE>, defaultIndex<V_TYPE>, IND...>
	{
		static const bool value=true;
	};
	template <int V_TYPE>
	struct check_same_variance_or_numeric<defaultIndex<V_TYPE> >
	{
		static const bool value=false;
	};

	template <typename ... IND>
	struct defaultIndexVariant;

	template <int V_TYPE, typename ... IND>
	struct defaultIndexVariant<defaultIndex<V_TYPE>, IND...>
	{
		static const bool defaultIndices=(!check_same_variance_or_numeric<defaultIndex<V_TYPE>, IND...>::value) && defaultIndexVariant<IND...>::defaultIndices;
	};

	template <typename HEAD, typename ... IND>
	struct defaultIndexVariant<HEAD, IND...>
	{
		static const bool defaultIndices=false;
	};

	template<>
	struct defaultIndexVariant<>
	{
		static const bool defaultIndices=true;
	};

	template <typename ST, typename ... IND>
	struct parenthesis;

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct stuple: public stuple_base<0, SHAPES...>
	{
		typedef OST ost_type;
		OST ost;
#if __cplusplus > 201103L
		inline constexpr ost_type *post() noexcept { return &ost; }
#else
		inline ost_type *post() noexcept { return &ost; }
#endif
		using base_type=stuple_base<0, SHAPES...>;
		inline constexpr const base_type& base() const noexcept { return *this; }
		typedef stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> tpl_type;
		template <typename ... Ts>
		inline stuple(const OST& ost, const std::tuple<Ts...>& inds, std::integral_constant<bool, false>): base_type(ost, inds), ost(ost) { parenthesis<OST, Ts...>::finalize(*this, ost, inds); }
		template <typename OLD_ST, typename ... Ts>
		inline stuple(const OLD_ST& ost, const std::tuple<Ts...>& inds, std::integral_constant<bool, true>): base_type(inds,ost), ost(ost.ost) {}
		template <typename VN, size_t O_SNUM, size_t O_CONT, bool O_IS_INDEXED, typename ... O_SHAPES>
		inline stuple(const stuple<O_SNUM, O_CONT, O_IS_INDEXED, OST, O_SHAPES...>& st, VN vn): base_type(st, vn), ost(ost) { }
		stuple() = default;
		static const size_t snum=SNUM;
		static const bool is_indexed=IS_INDEXED;
		template <size_t I>
		using element_type=typename stuple_element<I, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >::type;
		inline void get_shape(size_t (&aw)[SNUM]) const noexcept { shape_getter<SNUM, tpl_type >::shape(*this,aw); }
		inline size_t size() const noexcept { return size_getter<sizeof...(SHAPES), tpl_type>::size(*this); }
	};

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename ... SHAPES>
	struct stuple<SNUM, CONT, IS_INDEXED, void, SHAPES...>: public stuple_base<0, SHAPES...>
	{
		typedef void ost_type;
		inline constexpr ost_type *post() noexcept { return NULL; }
		using base_type=stuple_base<0, SHAPES...>;
		inline constexpr const base_type& base() const noexcept { return *this; }
		typedef stuple<SNUM, CONT, IS_INDEXED, void, SHAPES...> tpl_type;
		template <typename OST, typename ... Ts>
		inline stuple(const OST& ost, const std::tuple<Ts...>& inds, std::integral_constant<bool, false>): base_type(ost, inds) { parenthesis<OST, Ts...>::finalize(*this, ost, inds); }
		template <typename OST, typename ... Ts>
		inline stuple(const OST& ost, const std::tuple<Ts...>& inds, std::integral_constant<bool, true>): base_type(inds, ost) {}
		template <typename VN, size_t O_SNUM, size_t O_CONT, bool O_IS_INDEXED, typename ... O_SHAPES>
		inline stuple(const stuple<O_SNUM, O_CONT, O_IS_INDEXED, void, O_SHAPES...>& st, VN vn): base_type(st, vn) { }
		stuple() = default;
		static const size_t snum=SNUM;
		static const bool is_indexed=IS_INDEXED;
		template <size_t I>
		using element_type=typename stuple_element<I, stuple<SNUM, CONT, IS_INDEXED, void, SHAPES...> >::type;
		inline void get_shape(size_t (&aw)[SNUM]) const noexcept { shape_getter<SNUM, tpl_type >::shape(*this,aw); }
		inline size_t size() const noexcept { return size_getter<sizeof...(SHAPES), tpl_type>::size(*this); }
	};


	template <size_t DIM, typename START_TUPLE>
	struct add_initial_shapes;

	template <size_t DIM, size_t SNUM, typename OST, typename ... SHAPES>
	struct add_initial_shapes<DIM, stuple<SNUM, SNUM, false, OST, SHAPES...> >
	{
		using type=typename add_initial_shapes<DIM-1, stuple<SNUM, SNUM, false, OST, segment<-DIM,USAGE_FULL,0,0>, SHAPES...> >::type;
	};

	template <size_t SNUM, typename OST, typename ... SHAPES>
	struct add_initial_shapes<0, stuple<SNUM, SNUM, false, OST, SHAPES...> >
	{
		using type=stuple<SNUM, SNUM, false, OST, SHAPES...>;
	};

	template <size_t DIM, size_t CUR, typename SHAPES>
	struct shapes_initializer;

	template <size_t DIM, size_t CUR, size_t SNUM, typename HEAD, typename ... SHAPES>
	struct shapes_initializer<DIM, CUR, stuple<SNUM, SNUM, false, void, HEAD, SHAPES...> >
	{
		inline static size_t init(stuple<SNUM, SNUM, false, void, HEAD, SHAPES...>& sh, const size_t *aw)
		{
			size_t s=shapes_initializer<DIM, CUR+1, stuple<SNUM, SNUM, false, void, HEAD, SHAPES...> >::init(sh,aw);
			get<CUR>(sh).init(aw[CUR],s);
			return s*aw[CUR];
		}
	};

	template <size_t DIM, size_t SNUM, typename ... SHAPES>
	struct shapes_initializer<DIM, DIM, stuple<SNUM, SNUM, false, void, SHAPES...> >
	{
		inline static size_t init(stuple<SNUM, SNUM, false, void, SHAPES...>& sh, const size_t *aw) noexcept
		{
			return 1;
		}
	};

	template <size_t DIM, size_t SNUM, typename HEAD, typename ... SHAPES>
	struct shapes_initializer<DIM, DIM, stuple<SNUM, SNUM, false, void, HEAD, SHAPES...> >
	{
		inline static size_t init(stuple<SNUM, SNUM, false, void, HEAD, SHAPES...>& sh, const size_t *aw) noexcept
		{
			return 1;
		}
	};

	template <size_t DIM>
	struct initial_shapes
	{
		using type=typename add_initial_shapes<DIM, stuple<DIM, DIM, false, void> >::type;
		inline static size_t init(type& sh, const size_t (&aw)[DIM])
		{
			return shapes_initializer<DIM,0,type>::init(sh,aw);
		}
	};

	template <size_t CUR, typename ST>
	struct segment_shapes_initializer_base;

	template <size_t CUR, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct segment_shapes_initializer_base<CUR, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >
	{
		inline static size_t init(stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, const size_t (&aw)[SNUM])
		{
			size_t s=segment_shapes_initializer_base<CUR+1, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >::init(st,aw);
			get<CUR>(st).init(aw[CUR],s);
			return s*aw[CUR];
		}
	};

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct segment_shapes_initializer_base<SNUM, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >
	{
		inline static size_t init(stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, const size_t (&aw)[SNUM])
		{
			return 1;
		}
	};

	template <typename ST>
	struct segment_shapes_initializer: public segment_shapes_initializer_base<0, ST>
	{
	};


	template <size_t SNUM, typename TS, typename ... SHAPES>
	struct create_shapes_like;
	template <size_t SNUM, typename ... TS, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct create_shapes_like<SNUM, type_sequence<TS...>, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>:
		public create_shapes_like<SNUM+1, type_sequence<TS..., segment<V_TYPE, USAGE_FULL, 0,0> >, SHAPES...>
	{
	};
	template <size_t SNUM, typename ... TS>
	struct create_shapes_like<SNUM, type_sequence<TS...> >
	{
		using type=stuple<SNUM, SNUM, true, void, TS...>;
	};

	template <typename ST>
	struct stuple_like;

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct stuple_like<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public create_shapes_like<0, type_sequence<>, SHAPES...>
	{
	};

	template <size_t MAX_CONT, typename SC, typename ST>
	struct finalize_parenthesis;

	template <size_t MAX_CONT, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SRCSHAPES, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct finalize_parenthesis<MAX_CONT, type_sequence<sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SRCSHAPES...>, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public finalize_parenthesis<(CONT<MAX_CONT?CONT:MAX_CONT), type_sequence<SRCSHAPES...>, stuple<SNUM, CONT, IS_INDEXED, OST, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...> >
	{
	};

	template <size_t MAX_CONT, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SRCSHAPES, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct finalize_parenthesis<MAX_CONT, type_sequence<segment<V_TYPE, USAGE_FULL, USE_ALSO, PARENT>, SRCSHAPES...>, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public finalize_parenthesis<MAX_CONT, type_sequence<SRCSHAPES...>, stuple<SNUM, (CONT+1<MAX_CONT?CONT+1:MAX_CONT), IS_INDEXED, OST, segment<V_TYPE, USAGE_FULL, USE_ALSO, PARENT>, SHAPES...> >
	{
	};

	template <size_t MAX_CONT, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SRCSHAPES, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES >
	struct finalize_parenthesis<MAX_CONT, type_sequence<segment<V_TYPE, USAGE_CONT, USE_ALSO, PARENT>, SRCSHAPES...>, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public finalize_parenthesis<(CONT+1<MAX_CONT?CONT+1:MAX_CONT), type_sequence<SRCSHAPES...>, stuple<SNUM, (CONT+1<MAX_CONT?CONT+1:MAX_CONT), IS_INDEXED, OST, segment<V_TYPE, USAGE_CONT, USE_ALSO, PARENT>, SHAPES...> >
	{
	};

	template <size_t MAX_CONT, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct finalize_parenthesis<MAX_CONT, type_sequence<>, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >
	{
		using type=stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>;
	};

	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t SNUM, size_t INUM, bool NEED_PARENT, typename NSC, typename INUMS, typename OPOSS, typename OST, typename SC, typename INDT, typename ... IND>
	struct parenthesis_ind_skip_slave;

	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t GP, typename OSC, typename OINUMS, typename NSC, typename INUMS, typename OPOSS, typename OSHAPE, typename ITYPE, typename CREATED_SHAPE, typename NINUM, size_t SNUM, size_t INUM, bool NEED_PARENT, typename OST, typename HEAD, typename SC, typename INDT, typename ... IND>
	struct parenthesis_slave;

// common
	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t GP, typename NEW_SHAPE, typename INUM0, typename ... TNSC, typename ... TINUMS, typename ... NS, typename ... NIN, typename OPOSS, template <int, enum dSU, size_t, size_t> class sClass, int O_V_TYPE, enum dSU O_USAGE, size_t O_USE_ALSO, size_t O_PARENT, template <int> class iClass, int I_V_TYPE, typename CREATED_SHAPE, typename NINUM, size_t SNUM, size_t INUM, bool NEED_PARENT, typename OST, typename HEAD, typename SC, typename INDT, typename ... IND>
	struct parenthesis_slave<CONT, MIN_USAGE, OPOS, POS, GP, type_sequence<NEW_SHAPE, TNSC...>, std::tuple<INUM0,TINUMS...>, type_sequence<NS...>, std::tuple<NIN...>, OPOSS, sClass<O_V_TYPE, O_USAGE, O_USE_ALSO, O_PARENT>, iClass<I_V_TYPE>, CREATED_SHAPE, NINUM, SNUM, INUM, NEED_PARENT, OST, HEAD, SC, INDT, IND...>:
		public parenthesis_slave<CONT, MIN_USAGE, OPOS, POS, GP-1, type_sequence<TNSC...>, std::tuple<TINUMS...>, type_sequence<NS...,NEW_SHAPE>, std::tuple<NIN..., INUM0>, OPOSS, sClass<O_V_TYPE, O_USAGE, O_USE_ALSO, O_PARENT>, iClass<I_V_TYPE>, CREATED_SHAPE, NINUM, SNUM, INUM, NEED_PARENT, OST, HEAD, SC, INDT, IND...>
	{
//		using update=typename parenthesis_slave<OPOS, POS, GP-1, type_sequence<TNSC...>, type_sequence<TNIC...>, type_sequence<TPSC...>, std::tuple<TINUMS...>, type_sequence<NS...,NEW_SHAPE>, type_sequence<ITYPE, NI...>, type_sequence<PARENT_SHAPE, NP...>, std::tuple<INUM0, NIN...>, sClass<O_V_TYPE, O_USAGE, O_USE_ALSO, O_PARENT>, iClass<I_V_TYPE>, CREATED_SHAPE, NINUM>::update;
	};

// GLUE
	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t GP, enum dSU N_USAGE, size_t N_USE_ALSO, typename INUM0, typename ... TNSC, typename ... TINUMS, typename ... NS, typename ... NIN, typename ... IOPOSS, template <int, enum dSU, size_t, size_t> class sClass, int O_V_TYPE, enum dSU O_USAGE, size_t O_USE_ALSO, size_t O_PARENT, template <int> class iClass, int I_V_TYPE, enum dSU C_USAGE, typename NINUM, size_t SNUM, size_t INUM, bool NEED_PARENT, typename OST, typename HEAD, typename ... SC, typename INDT, typename ... IND>
	struct parenthesis_slave<CONT, MIN_USAGE, OPOS, POS, GP, type_sequence<segment<I_V_TYPE, N_USAGE, N_USE_ALSO, 0>, TNSC...>, std::tuple<INUM0,TINUMS...>, type_sequence<NS...>, std::tuple<NIN...>, std::tuple<IOPOSS...>, sClass<O_V_TYPE, O_USAGE, O_USE_ALSO, O_PARENT>, iClass<I_V_TYPE>, segment<I_V_TYPE, C_USAGE, 0, 0>, NINUM, SNUM, INUM, NEED_PARENT, OST, HEAD, type_sequence<SC...>, INDT, IND...>:
		public parenthesis_ind_skip_slave<(sizeof...(SC)<CONT?sizeof...(SC):CONT), USAGE_CONT, OPOS+1, POS, SNUM, INUM+1, NEED_PARENT, type_sequence<NS..., segment<I_V_TYPE, (N_USAGE<USAGE_ENUM)?USAGE_ENUM:N_USAGE, N_USE_ALSO, 0>, TNSC...>, std::tuple<NIN..., INUM0, TINUMS...>, std::tuple<IOPOSS...>, OST, type_sequence<SC...>, INDT, IND...>
	{
		using base=parenthesis_ind_skip_slave<CONT, USAGE_CONT, OPOS+1, POS, SNUM, INUM+1, NEED_PARENT, type_sequence<NS..., segment<I_V_TYPE, (N_USAGE<USAGE_ENUM)?USAGE_ENUM:N_USAGE, N_USE_ALSO, 0>, TNSC...>, std::tuple<NIN..., INUM0, TINUMS...>, std::tuple<IOPOSS...>, OST, type_sequence<SC...>, INDT, IND...>;
		inline static void finalize(typename base::type& st, const OST& ost, const INDT& ind)
		{
			using created_shape=segment<I_V_TYPE, C_USAGE, 0, 0>;
			const HEAD& h=get<OPOS>(ost);
			std::get<INUM>(ind).check_bounds(INUM, h.length());
			created_shape nsh(std::get<INUM>(ind),get<OPOS>(ost));
			using GPSH=typename base::type::template element_type<GP>;
			GPSH& msh=get<GP>(st);
			segment_gluer<GPSH, created_shape>::update(msh,nsh);
			if (nsh.length()!=msh.length())
				throw incompatibleIndices(INUM0::value, INUM, msh.length(), nsh.length());
			base::finalize(st, ost, ind);
		}
	};

// GLUE && SLAVE
	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t GP, enum dSU N_USAGE, typename INUM0, typename ... TNSC, typename ... TINUMS, typename ... NS, typename ... NIN, typename ... IOPOSS, template <int, enum dSU, size_t, size_t> class sClass, int O_V_TYPE, enum dSU O_USAGE, size_t O_USE_ALSO, size_t O_PARENT, template <int> class iClass, int I_V_TYPE, enum dSU C_USAGE, typename NINUM, size_t SNUM, size_t INUM, bool NEED_PARENT, typename OST, typename HEAD, typename ... SC, typename INDT, typename ... IND>
	struct parenthesis_slave<CONT, MIN_USAGE, OPOS, POS, GP, type_sequence<segment<I_V_TYPE, N_USAGE, 0, 0>, TNSC...>, std::tuple<INUM0,TINUMS...>, type_sequence<NS...>, std::tuple<NIN...>, std::tuple<IOPOSS...>, sClass<O_V_TYPE, O_USAGE, O_USE_ALSO, O_PARENT>, iClass<I_V_TYPE>, segment<I_V_TYPE, C_USAGE, 0, 0>, NINUM, SNUM, INUM, NEED_PARENT, OST, HEAD, type_sequence<SC...>, INDT, IND...>:
		public parenthesis_ind_skip_slave<(sizeof...(SC)<CONT?sizeof...(SC):CONT), USAGE_CONT, OPOS+1, POS, SNUM, INUM+1, NEED_PARENT, type_sequence<NS..., segment<I_V_TYPE, (N_USAGE<USAGE_ENUM)?USAGE_ENUM:N_USAGE, 0, 0>, TNSC...>, std::tuple<NIN..., INUM0, TINUMS...>, std::tuple<IOPOSS...>, OST, type_sequence<SC...>, INDT, IND...>
	{
		using base=parenthesis_ind_skip_slave<CONT, USAGE_CONT, OPOS+1, POS, SNUM, INUM+1, NEED_PARENT, type_sequence<NS..., segment<I_V_TYPE, (N_USAGE<USAGE_ENUM)?USAGE_ENUM:N_USAGE, 0, 0>, TNSC...>, std::tuple<NIN..., INUM0, TINUMS...>, std::tuple<IOPOSS...>, OST, type_sequence<SC...>, INDT, IND...>;
		inline static void finalize(typename base::type& st, const OST& ost, const INDT& ind)
		{
			using created_shape=segment<I_V_TYPE, C_USAGE, 0, 0>;
			const HEAD& h=get<OPOS>(ost);
			std::get<INUM>(ind).check_bounds(INUM, h.length());
			created_shape nsh(std::get<INUM>(ind),get<OPOS>(ost));
			using GPSH=typename base::type::template element_type<GP>;
			GPSH& msh=get<GP>(st);
			segment_gluer<GPSH, created_shape>::update(msh,nsh);
			if (nsh.length()!=msh.length())
				throw incompatibleIndices(INUM0::value, INUM, I_V_TYPE, msh.length(), nsh.length());
			base::finalize(st, ost, ind);
		}
	};


// SLAVE
	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t GP, template <int, enum dSU, size_t, size_t> class csClass, enum dSU N_USAGE, size_t N_PARENT, typename INUM0, typename ... TNSC, typename ... TINUMS, typename ... NS, typename ... NIN, typename ... IOPOSS, template <int, enum dSU, size_t, size_t> class sClass, int O_V_TYPE, enum dSU O_USAGE, size_t O_USE_ALSO, size_t O_PARENT, template <int> class iClass, int I_V_TYPE, template <int, enum dSU, size_t, size_t> class crsClass, enum dSU C_USAGE, size_t C_PARENT, typename NINUM, size_t SNUM, size_t INUM, bool NEED_PARENT, typename OST, typename HEAD, typename ... SC, typename INDT, typename ... IND>
	struct parenthesis_slave<CONT, MIN_USAGE, OPOS, POS, GP, type_sequence<csClass<I_V_TYPE, N_USAGE, 0, N_PARENT>, TNSC...>, std::tuple<INUM0,TINUMS...>, type_sequence<NS...>, std::tuple<NIN...>, std::tuple<IOPOSS...>, sClass<O_V_TYPE, O_USAGE, O_USE_ALSO, O_PARENT>, iClass<I_V_TYPE>, crsClass<I_V_TYPE, C_USAGE, 0, C_PARENT>, NINUM, SNUM, INUM, NEED_PARENT, OST, HEAD, type_sequence<SC...>, INDT, IND...>:
		public parenthesis_ind_skip_slave<(sizeof...(SC)<CONT?sizeof...(SC):CONT), USAGE_CONT, OPOS+1, POS+1, SNUM, INUM+1, NEED_PARENT || C_PARENT?true:false, type_sequence<crsClass<I_V_TYPE, USAGE_SLAVE, 0, C_PARENT>, NS..., csClass<I_V_TYPE, (N_USAGE<USAGE_ENUM)?USAGE_ENUM:N_USAGE, POS, N_PARENT>, TNSC...>, std::tuple<NIN..., INUM0, TINUMS..., NINUM>, std::tuple<IOPOSS...,typename std::integral_constant<size_t,OPOS>::type>, OST, type_sequence<SC...>, INDT, IND...>
	{
		using base=parenthesis_ind_skip_slave<CONT, USAGE_CONT, OPOS+1, POS+1, SNUM, INUM+1, NEED_PARENT || C_PARENT?true:false, type_sequence<crsClass<I_V_TYPE, USAGE_SLAVE, 0, C_PARENT>, NS..., csClass<I_V_TYPE, (N_USAGE<USAGE_ENUM)?USAGE_ENUM:N_USAGE, POS, N_PARENT>, TNSC...>, std::tuple<NIN..., INUM0, TINUMS..., NINUM>, std::tuple<IOPOSS...,typename std::integral_constant<size_t,OPOS>::type>, OST, type_sequence<SC...>, INDT, IND...>;
		inline static void finalize(typename base::type& st, const OST& ost, const INDT& ind)
		{
			const HEAD& h=get<OPOS>(ost);
			std::get<INUM>(ind).check_bounds(INUM, h.length());
			using GPSH=typename base::type::template element_type<GP>;
			GPSH& msh=get<GP>(st);
			using NSH=typename base::type::template element_type<POS>;
			NSH& nsh=get<POS>(st);
			if (nsh.length()!=msh.length())
				throw incompatibleIndices(INUM0::value, INUM, I_V_TYPE, msh.length(), nsh.length());
			base::finalize(st, ost, ind);
		}
	};

// last
	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t GP, typename ... NS, typename ... NIN, typename ... IOPOSS, template <int, enum dSU, size_t, size_t> class sClass, int O_V_TYPE, enum dSU O_USAGE, size_t O_USE_ALSO, size_t O_PARENT, template <int> class iClass, int I_V_TYPE, template <int, enum dSU, size_t, size_t> class crsClass, enum dSU C_USAGE, size_t C_PARENT, typename NINUM, size_t SNUM, size_t INUM, bool NEED_PARENT, typename OST, typename HEAD, typename SC, typename INDT, typename ... IND>
	struct parenthesis_slave<CONT, MIN_USAGE, OPOS, POS, GP, type_sequence<>, std::tuple<>, type_sequence<NS...>, std::tuple<NIN...>, std::tuple<IOPOSS...>, sClass<O_V_TYPE, O_USAGE, O_USE_ALSO, O_PARENT>, iClass<I_V_TYPE>, crsClass<I_V_TYPE, C_USAGE, 0, C_PARENT>, NINUM, SNUM, INUM, NEED_PARENT, OST, HEAD, SC, INDT, IND...>:
		public parenthesis_ind_skip_slave<CONT, USAGE_FULL, OPOS+1, POS+1, SNUM+1, INUM+1, (NEED_PARENT || C_PARENT>0), type_sequence<crsClass<I_V_TYPE, (MIN_USAGE<C_USAGE)?C_USAGE:MIN_USAGE, 0, C_PARENT>, NS...>, std::tuple<NIN..., NINUM>, std::tuple<IOPOSS...,typename std::integral_constant<size_t,OPOS>::type>, OST, SC, INDT, IND...>
	{
		using base=parenthesis_ind_skip_slave<CONT, USAGE_FULL, OPOS+1, POS+1, SNUM+1, INUM+1, (NEED_PARENT || C_PARENT>0), type_sequence<crsClass<I_V_TYPE, (MIN_USAGE<C_USAGE)?C_USAGE:MIN_USAGE, 0, C_PARENT>, NS...>, std::tuple<NIN..., NINUM>, std::tuple<IOPOSS...,typename std::integral_constant<size_t,OPOS>::type>, OST, SC, INDT, IND...>;
		inline static void finalize(typename base::type& st, const OST& ost, const INDT& ind)
		{
			const HEAD& h=get<OPOS>(ost);
			std::get<INUM>(ind).check_bounds(INUM, h.length());
			base::finalize(st, ost, ind);
		}
	};



	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t SNUM, size_t INUM, bool NEED_PARENT, typename NSC, typename INUMS, typename OPOSS, typename OST, typename SC, typename INDT, typename ... IND>
	struct parenthesis_ind;

	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t SNUM, size_t INUM, bool NEED_PARENT, typename NSC, typename INUMS, typename OPOSS, typename OST, typename ... SHAPES, typename INDT, typename ... IND>
	struct parenthesis_ind_skip_slave<CONT, MIN_USAGE, OPOS, POS, SNUM, INUM, NEED_PARENT, NSC, INUMS, OPOSS, OST, type_sequence<SHAPES...>, INDT, IND...>:
		public parenthesis_ind<CONT, MIN_USAGE, OPOS, POS, SNUM, INUM, NEED_PARENT, NSC, INUMS, OPOSS, OST, type_sequence<SHAPES...>, INDT, IND...>
	{
	};

// USAGE_SLAVE
	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t SNUM, size_t INUM, bool NEED_PARENT, typename NSC, typename INUMS, typename OPOSS, typename OST, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES, typename INDT, typename ... IND>
	struct parenthesis_ind_skip_slave<CONT, MIN_USAGE, OPOS, POS, SNUM, INUM, NEED_PARENT, NSC, INUMS, OPOSS, OST, type_sequence<sClass<V_TYPE, USAGE_SLAVE, USE_ALSO, PARENT>, SHAPES...>, INDT, IND...>:
		public parenthesis_ind_skip_slave<CONT, USAGE_CONT, OPOS+1, POS, SNUM, INUM, NEED_PARENT, NSC, INUMS, OPOSS, OST, type_sequence<SHAPES...>, INDT, IND...>
	{
	};


// Index is not size_t and old shape is not USAGE_SLAVE
	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t SNUM, size_t INUM, bool NEED_PARENT, typename ... NS, typename INUMS, typename OPOSS, typename OST, typename HEAD, typename ... SHAPES, typename INDT, template <int> class iClass, int I_V_TYPE, typename ... IND>
	struct parenthesis_ind<CONT, MIN_USAGE, OPOS, POS, SNUM, INUM, NEED_PARENT, type_sequence<NS...>, INUMS, OPOSS, OST, type_sequence<HEAD, SHAPES...>, INDT, iClass<I_V_TYPE>, IND...>:
		public parenthesis_slave<CONT, MIN_USAGE, OPOS, POS, POS-1, type_sequence<NS...>, INUMS, type_sequence<>, std::tuple<>, OPOSS, HEAD, iClass<I_V_TYPE>, typename iClass<I_V_TYPE>::template applyer<HEAD, OPOS>::new_shape, typename std::integral_constant<size_t, INUM>::type, SNUM, INUM, NEED_PARENT, OST, HEAD, type_sequence<SHAPES...>, INDT, IND...>
	{
	};

// numeric index
	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t SNUM, size_t INUM, bool NEED_PARENT, typename ... NS, typename INUMS, typename OPOSS, typename OST, typename HEAD, typename ... SHAPES, typename INDT, typename NITYPE, typename ... IND>
	struct parenthesis_ind<CONT, MIN_USAGE, OPOS, POS, SNUM, INUM, NEED_PARENT, type_sequence<NS...>, INUMS, OPOSS, OST, type_sequence<HEAD, SHAPES...>, INDT, NITYPE, IND...>:
		public parenthesis_ind_skip_slave<(sizeof...(SHAPES)<CONT?sizeof...(SHAPES):CONT), USAGE_CONT, OPOS+1, POS, SNUM, INUM+1, NEED_PARENT, type_sequence<NS...>, INUMS, OPOSS, OST, type_sequence<SHAPES...>, INDT, IND...>
	{
	};

// No more indices need_parent, rest of old shapes
	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t SNUM, size_t INUM, typename NSC, typename INUMS, typename OPOSS, typename OST, typename ... SHAPES, typename INDT>
	struct parenthesis_ind<CONT, MIN_USAGE, OPOS, POS, SNUM, INUM, true, NSC, INUMS, OPOSS, OST, type_sequence<SHAPES...>, INDT>:
		public finalize_parenthesis<CONT, NSC, stuple<SNUM, 0, true, OST> >
	{
		inline static void finalize(typename finalize_parenthesis<CONT, NSC, stuple<SNUM, 0, true, OST> >::type& st, const OST& ost, const INDT& ind) {}
		template <size_t I>
		inline static constexpr const typename OST::template element_type<std::tuple_element<I,OPOSS>::type::type::value>& parent(const OST& ost, const INDT& ind) noexcept
		{
			return get<std::tuple_element<I,OPOSS>::type::type::value>(ost);
		}
		template <size_t I>
		inline static constexpr const typename std::tuple_element<std::tuple_element<I,INUMS>::type::type::value,INDT>::type& index(const OST& ost, const INDT& ind) noexcept
		{
			return std::get<std::tuple_element<I,INUMS>::type::type::value>(ind);
		}
	};

// No more indices need_parent==false, rest of old shapes
	template <size_t CONT, enum dSU MIN_USAGE, size_t OPOS, size_t POS, size_t SNUM, size_t INUM, typename NSC, typename INUMS, typename OPOSS, typename OST, typename ... SHAPES, typename INDT>
	struct parenthesis_ind<CONT, MIN_USAGE, OPOS, POS, SNUM, INUM, false, NSC, INUMS, OPOSS, OST, type_sequence<SHAPES...>, INDT>:
		public finalize_parenthesis<CONT, NSC, stuple<SNUM, 0, true, void> >
	{
		inline static void finalize(typename finalize_parenthesis<CONT, NSC, stuple<SNUM, 0, true, void> >::type& st, const OST& ost, const INDT& ind) {}
		template <size_t I>
		inline static constexpr const typename OST::template element_type<std::tuple_element<I,OPOSS>::type::type::value>& parent(const OST& ost, const INDT& ind) noexcept
		{
			return get<std::tuple_element<I,OPOSS>::type::type::value>(ost);
		}
		template <size_t I>
		inline static constexpr const typename std::tuple_element<std::tuple_element<I,INUMS>::type::type::value,INDT>::type& index(const OST& ost, const INDT& ind) noexcept
		{
			return std::get<std::tuple_element<I,INUMS>::type::type::value>(ind);
		}
	};


	template <size_t POS, size_t OPOS, typename INDT, typename ST, typename OSC, typename SC, typename ... IND>
	struct parenthesis_special;
// default index
	template <size_t POS, size_t OPOS, typename INDT, typename ST, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES, typename ... NSHAPES, int I_V_TYPE, typename ... IND>
	struct parenthesis_special<POS, OPOS, INDT, ST, type_sequence<sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>, type_sequence<NSHAPES...>, defaultIndex<I_V_TYPE>, IND...>:
		public parenthesis_special<POS+1, OPOS+1, INDT, ST, type_sequence<SHAPES...>, type_sequence<NSHAPES..., sClass<I_V_TYPE, USAGE, USE_ALSO, PARENT> >, IND...>
	{
		inline static constexpr const defaultIndex<I_V_TYPE>& index(const ST& ost, const INDT& ind, typename std::integral_constant<size_t, POS>::type) noexcept
		{
			return std::get<POS>(ind);
		}
		inline static constexpr const sClass<V_TYPE, USAGE, USE_ALSO, PARENT>& parent(const ST& ost, const INDT& ind, typename std::integral_constant<size_t, POS>::type) noexcept
		{
			return get<POS>(ost);
		}
	};
// no more indices
	template <size_t POS, size_t OPOS, typename INDT, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, typename ... NSHAPES>
	struct parenthesis_special<POS, OPOS, INDT, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, type_sequence<>, type_sequence<NSHAPES...> >
	{
		using type=stuple<SNUM, CONT, true, OST, NSHAPES...>;
		inline static void finalize(type& st, const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& ost, const INDT& ind) {}
	};

	template <typename ST, bool NORMALVARIANT, typename ... IND>
	struct parenthesis_chk_normal;

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, typename ... IND>
	struct parenthesis_chk_normal<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, false, IND...>:
		public parenthesis_ind_skip_slave<CONT, USAGE_FULL, 0, 0, 0, 0, false, type_sequence<>, std::tuple<>, std::tuple<>, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, type_sequence<SHAPES...>, std::tuple<IND...>, IND...>
	{
		static const bool defaultIndices=false;
	};
// special case
	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, typename ... IND>
	struct parenthesis_chk_normal<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, true, IND...>:
		public parenthesis_special<0, 0, std::tuple<IND...>, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, type_sequence<SHAPES...>, type_sequence<>, IND...>
	{
		static const bool defaultIndices=true;
	};


	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, typename ... IND>
	struct parenthesis<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, IND...>:
		public parenthesis_chk_normal<
		stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, defaultIndexVariant<IND...>::defaultIndices, IND...>
	{
		static_assert(SNUM==sizeof...(IND),"Incorrect number of indices");
	};

	template <typename T, size_t POS, bool IS_LAST_CONT, typename T_SHAPE, typename SHAPES>
	struct shapes_caller_base;

// parent- also+
	template <typename T, size_t POS, bool IS_LAST_CONT, template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_caller_base<T, POS, IS_LAST_CONT, shapeClass<V_TYPE, USAGE, USE_ALSO, 0>, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >
	{
		using ST=stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>;
		using scba=shapes_caller_base<T, USE_ALSO, (CONT>0 && USE_ALSO==sizeof...(SHAPES)-1), typename std::tuple_element<USE_ALSO, std::tuple<SHAPES...> >::type, ST>;
		inline static void apply_numeric_index(stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, T *&p, ptrdiff_t o)
		{
			get<POS>(st).template apply_numeric_index<T>(p,o);
			scba::apply_numeric_index(st,p,o);
		}
		inline static void apply_numeric_index(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, T *&p, ptrdiff_t o)
		{
			get<POS>(st).template apply_numeric_index<T>(p,o);
			scba::apply_numeric_index(st,p,o);
		}
	};

// parent- also-
	template <typename T, size_t POS, bool IS_LAST_CONT, template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, enum dSU USAGE, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_caller_base<T, POS, IS_LAST_CONT, shapeClass<V_TYPE, USAGE, 0, 0>, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >
	{
		inline static void apply_numeric_index(stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, T *&p, ptrdiff_t o) { get<POS>(st).template apply_numeric_index<T>(p,o); }
		inline static void apply_numeric_index(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, T *&p, ptrdiff_t o) { get<POS>(st).template apply_numeric_index<T>(p,o); }
	};


// parent- also- continuous
	template <typename T, size_t POS, int V_TYPE, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_caller_base<T, POS, true, segment<V_TYPE, USAGE_FULL, 0, 0>, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >
	{
		inline static void apply_numeric_index(stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, T *&p, ptrdiff_t o) noexcept { p+=o; }
		inline static void apply_numeric_index(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, T *&p, ptrdiff_t o) noexcept { p+=o; }
	};

// parent+ also-
	template <typename T, size_t POS, bool IS_LAST_CONT, template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, enum dSU USAGE, size_t PARENT, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_caller_base<T, POS, IS_LAST_CONT, shapeClass<V_TYPE, USAGE, 0, PARENT>, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >
	{
		using PARENT_SHAPE=typename OST::template element_type<PARENT-1>;
		using scbp=shapes_caller_base<T, PARENT-1, (CONT>0 && PARENT-1==sizeof...(SHAPES)-1), PARENT_SHAPE, OST>;
		inline static void apply_numeric_index(stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, T *&p, ptrdiff_t o)
		{
			o=get<POS>(st).template apply_numeric_index<T>(p,o); // ,std::get<POS>(st.ost.tpl)
			scbp::apply_numeric_index(st.ost, p, o);
		}
		inline static void apply_numeric_index(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, T *&p, ptrdiff_t o)
		{
			o=get<POS>(st).template apply_numeric_index<T>(p,o); // ,std::get<POS>(st.ost.tpl)
			scbp::apply_numeric_index(st.ost, p, o);
		}
	};

//parent+ also+
	template <typename T, size_t POS, bool IS_LAST_CONT, template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_caller_base<T, POS, IS_LAST_CONT, shapeClass<V_TYPE, USAGE, USE_ALSO, PARENT>, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >
	{
		using PARENT_SHAPE=typename OST::template element_type<PARENT-1>;
		using ST=stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>;
		using scba=shapes_caller_base<T, USE_ALSO, (CONT>0 && USE_ALSO==sizeof...(SHAPES)-1), typename std::tuple_element<USE_ALSO, std::tuple<SHAPES...> >::type, ST>;
		using scbp=shapes_caller_base<T, PARENT-1, (CONT>0 && PARENT-1==sizeof...(SHAPES)-1), PARENT_SHAPE, OST>;
		inline static void apply_numeric_index(stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, T *&p, ptrdiff_t o)
		{
			scba::apply_numeric_index(st,p,o);
			o=get<POS>(st).template apply_numeric_index<T>(p,o); // ,std::get<POS>(st.ost.tpl)
			scbp::apply_numeric_index(st.ost, p, o);
		}
		inline static void apply_numeric_index(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, T *&p, ptrdiff_t o)
		{
			scba::apply_numeric_index(st,p,o);
			o=get<POS>(st).template apply_numeric_index<T>(p,o);
			scbp::apply_numeric_index(st.ost, p, o);
		}
	};


	template <typename T, size_t POS, typename SHAPES>
	struct shapes_caller;

	template <typename T, size_t POS, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_caller<T, POS, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public shapes_caller_base<T, POS, (CONT>0 && POS==sizeof...(SHAPES)-1), typename std::tuple_element<POS, std::tuple<SHAPES...> >::type, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >
	{
	};


	template <typename T, size_t POS, bool IS_LAST_CONT, typename T_SHAPE, bool NEED_PARENT, typename SHAPES>
	struct shapes_iterator_base;

// parent- also+
	template <typename T, size_t POS, bool IS_LAST_CONT, template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_iterator_base<T, POS, IS_LAST_CONT, shapeClass<V_TYPE, USAGE, USE_ALSO, PARENT>, false, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public shapeClass<V_TYPE, USAGE, USE_ALSO, PARENT>::template iterator<T>
	{
		using ST=stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>;
		using siba=shapes_iterator_base<T, USE_ALSO, (CONT>0 && USE_ALSO==sizeof...(SHAPES)-1), typename std::tuple_element<USE_ALSO, std::tuple<SHAPES...> >::type, check_shape<typename std::tuple_element<USE_ALSO, std::tuple<SHAPES...> >::type>::need_parent, ST>;
		using master_iterator=typename shapeClass<V_TYPE, USAGE, USE_ALSO, PARENT>::template iterator<T>;
		siba also_iterator;
		inline shapes_iterator_base
			(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, const T * const p) : // noexcept ????
			master_iterator(get<POS>(st),p),also_iterator(st,p) {}
		inline void move_one(T *&p) { master_iterator::move_one(p); also_iterator.move_one(p); } // noexcept ????
		inline void move_one(const T *&p) { master_iterator::move_one(p); also_iterator.move_one(p); } // noexcept ????
		inline void move(T *&p, size_t offset) { master_iterator::move(p, offset); also_iterator.move(p, offset); } // noexcept ???
		inline void move(const T *&p, size_t offset) { master_iterator::move(p, offset); also_iterator.move(p, offset); } // noexcept ???
	};

// parent- also-
	template <typename T, size_t POS, bool IS_LAST_CONT, template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, enum dSU USAGE, size_t PARENT, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_iterator_base<T, POS, IS_LAST_CONT, shapeClass<V_TYPE, USAGE, 0, PARENT>, false, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public shapeClass<V_TYPE, USAGE, 0, PARENT>::template iterator<T>
	{
		using ST=stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>;
		using master_iterator=typename shapeClass<V_TYPE, USAGE, 0, PARENT>::template iterator<T>;
		inline shapes_iterator_base
			(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, const T * const p): // noexcept ????
			master_iterator(get<POS>(st),p) {}
	};

// parent- also-, segment
	template <typename T, size_t POS, bool IS_LAST_CONT, int V_TYPE, enum dSU USAGE, size_t PARENT, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_iterator_base<T, POS, IS_LAST_CONT, segment<V_TYPE, USAGE, 0, PARENT>, false, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public segment<V_TYPE, USAGE, 0, PARENT>::template iterator<T>
	{
		using ST=stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>;
		using master_iterator=typename segment<V_TYPE, USAGE, 0, PARENT>::template iterator<T>;
		static size_t length(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st) { return get<POS>(st).length(); }
		static size_t step(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st) { return get<POS>(st).step(); }
		inline shapes_iterator_base
			(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, const T * const p): // noexcept ????
			master_iterator(get<POS>(st),p) {}
	};

// parent- also-, continuous segment
	template <typename T, size_t POS, int V_TYPE, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_iterator_base<T, POS, true, segment<V_TYPE, USAGE_FULL, 0, 0>, false, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public segment<V_TYPE, USAGE_FULL, 0, 0>::template cont_iterator<T>
	{
		using ST=stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>;
		using master_iterator=typename segment<V_TYPE, USAGE_FULL, 0, 0>::template cont_iterator<T>;
		static size_t length(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st) { return get<POS>(st).length(); }
		static size_t step(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st) { return get<POS>(st).step(); }
		inline shapes_iterator_base
			(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, const T * const p): // noexcept ????
			master_iterator(get<POS>(st),p) {}
	};

// parent+ also-
	template <typename T, size_t POS, bool IS_LAST_CONT, template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, enum dSU USAGE, size_t PARENT, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_iterator_base<T, POS, IS_LAST_CONT, shapeClass<V_TYPE, USAGE, 0, PARENT>, true, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public shapeClass<V_TYPE, USAGE, 0, PARENT>::template iterator<T>
	{
		using PARENT_SHAPE=typename OST::template element_type<PARENT-1>;
		using ST=stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>;
		using sibp=shapes_iterator_base<T, PARENT-1, (CONT>0 && PARENT-1==sizeof...(SHAPES)-1), PARENT_SHAPE, check_shape<PARENT_SHAPE>::need_parent, OST>;
		using master_iterator=typename shapeClass<V_TYPE, USAGE, 0, PARENT>::template iterator<T>;
		sibp parent_iterator;
		inline shapes_iterator_base
			(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, const T * const p): // noexcept ????
			master_iterator(get<POS>(st),p),parent_iterator(st.ost,p) {}
		inline void move_one(T *&p) { parent_iterator.move(p,master_iterator::move_one(p)); } // noexcept ????
		inline void move_one(const T *&p) { parent_iterator.move(p,master_iterator::move_one(p)); } // noexcept ????
		inline void move(T *&p, size_t offset) { parent_iterator.move(p, master_iterator::move(p, offset)); } // noexcept ????
		inline void move(const T *&p, size_t offset) { parent_iterator.move(p, master_iterator::move(p, offset)); } // noexcept ????
	};

//parent+ also+
	template <typename T, size_t POS, bool IS_LAST_CONT, template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES>
	struct shapes_iterator_base<T, POS, IS_LAST_CONT, shapeClass<V_TYPE, USAGE, USE_ALSO, PARENT>, true, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >:
		public shapeClass<V_TYPE, USAGE, USE_ALSO, PARENT>::template iterator<T>
	{
		using PARENT_SHAPE=typename OST::template element_type<PARENT-1>;
		using ST=stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>;
		using siba=shapes_iterator_base<T, USE_ALSO, (CONT>0 && USE_ALSO==sizeof...(SHAPES)-1), typename std::tuple_element<USE_ALSO, std::tuple<SHAPES...> >::type, check_shape<typename std::tuple_element<USE_ALSO, std::tuple<SHAPES...> >::type>::need_parent, ST>;
		using sibp=shapes_iterator_base<T, PARENT-1, (CONT>0 && PARENT-1==sizeof...(SHAPES)-1), PARENT_SHAPE, check_shape<PARENT_SHAPE>::need_parent, OST>;
		using master_iterator=typename shapeClass<V_TYPE, USAGE, USE_ALSO, PARENT>::template iterator<T>;
		siba also_iterator;
		sibp parent_iterator;
		inline shapes_iterator_base
			(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>& st, const T * const p): // noexcept ????
			master_iterator(get<POS>(st),p),also_iterator(st,p),parent_iterator(st.ost,p) {}
		inline void move_one(T *&p) { parent_iterator.move(p,master_iterator::move_one(p)); also_iterator.move_one(p); } // noexcept ????
		inline void move_one(const T *&p) { parent_iterator.move(p,master_iterator::move_one(p)); also_iterator.move_one(p); } // noexcept ????
		inline void move(T *&p, size_t offset) { parent_iterator.move(p, master_iterator::move(p, offset)); also_iterator.move(p, offset); } // noexcept ????
		inline void move(const T *&p, size_t offset) { parent_iterator.move(p, master_iterator::move(p, offset)); also_iterator.move(p, offset); } // noexcept ????
	};

	template <typename T, typename OT, typename INDT, size_t POS, size_t INUM, typename SC, typename ... IND>
	struct index_shift_applyer_base;

	template <typename T, typename OT, typename INDT, size_t POS, size_t INUM, typename SC, typename ... IND>
	struct index_shift_applyer_base_skip_slave: public index_shift_applyer_base<T, OT, INDT, POS, INUM, SC, IND...>
	{
	};
	template <typename T, typename OT, typename INDT, size_t POS, size_t INUM, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... RSHAPES, typename ... IND>
	struct index_shift_applyer_base_skip_slave<T, OT, INDT, POS, INUM, type_sequence<sClass<V_TYPE, USAGE_SLAVE, USE_ALSO, PARENT>, RSHAPES...>, IND...>:
		public index_shift_applyer_base_skip_slave<T, OT, INDT, POS+1, INUM, type_sequence<RSHAPES...>, IND...>
	{
	};

// NORMAL
	template <typename T, typename OT, typename INDT, size_t POS, size_t INUM, typename HEAD, typename ... RSHAPES, template <int> class iClass, int I_V_TYPE, typename ... IND>
	struct index_shift_applyer_base<T, OT, INDT, POS, INUM, type_sequence<HEAD, RSHAPES...>, iClass<I_V_TYPE>, IND...>
	{
		inline static void shift_data(const OT& ot, T *& p, const INDT& ind)
		{
			shapes_caller<T,POS,OT>::apply_numeric_index(ot,p,std::get<INUM>(ind).first());
			index_shift_applyer_base_skip_slave<T, OT,INDT,POS+1,INUM+1,type_sequence<RSHAPES...>,IND...>::shift_data(ot,p,ind);
		}
		inline static void shift_data(OT& ot, T *& p, const INDT& ind)
		{
			shapes_caller<T,POS,OT>::apply_numeric_index(ot,p,std::get<INUM>(ind).first());
			index_shift_applyer_base_skip_slave<T, OT,INDT,POS+1,INUM+1,type_sequence<RSHAPES...>,IND...>::shift_data(ot,p,ind);
		}
	};
// numeric index
	template <typename T, typename OT, typename INDT, size_t POS, size_t INUM, typename HEAD, typename ... RSHAPES, typename NITYPE, typename ... IND>
	struct index_shift_applyer_base<T, OT, INDT, POS, INUM, type_sequence<HEAD, RSHAPES...>, NITYPE, IND...>
	{
		inline static void shift_data(const OT& ot, T *& p, const INDT& ind)
		{
			size_t v=std::get<INUM>(ind);
			const HEAD& h=get<POS>(ot);
			if (v>=h.length())
				throw outOfBounds(v,INUM,h.length());
			shapes_caller<T,POS,OT>::apply_numeric_index(ot,p,v);
			index_shift_applyer_base_skip_slave<T, OT,INDT,POS+1,INUM+1,type_sequence<RSHAPES...>,IND...>::shift_data(ot,p,ind);
		}
		inline static void shift_data(OT& ot, T *& p, const INDT& ind)
		{
			size_t v=std::get<INUM>(ind);
			const HEAD& h=get<POS>(ot);
			if (v>=h.length())
				throw outOfBounds(v,INUM,h.length());
			shapes_caller<T,POS,OT>::apply_numeric_index(ot,p,v);
			index_shift_applyer_base_skip_slave<T, OT,INDT,POS+1,INUM+1,type_sequence<RSHAPES...>,IND...>::shift_data(ot,p,ind);
		}
	};

	template <typename T, typename OT, typename INDT, size_t POS, size_t INUM, typename SC>
	struct index_shift_applyer_base<T, OT, INDT, POS, INUM, SC>
	{
		inline static void shift_data(const OT& ot, T *& p, const INDT& ind) {}
	};

	template <typename T, typename OT, typename ROST, typename INDT, size_t POS, size_t INUM, typename SC, typename ... IND>
	struct index_shift_applyer
	{
		inline static void shift_data(const OT& ot, ROST *post, T *& p, const INDT& ind)
		{
			index_shift_applyer_base_skip_slave<T, ROST, INDT, POS, INUM, SC, IND...>::shift_data(*post,p,ind);
		}
	};

	template <typename T, typename OT, typename INDT, size_t POS, size_t INUM, typename SC, typename ... IND>
	struct index_shift_applyer<T, OT, void, INDT, POS, INUM, SC, IND...>
	{
		inline static void shift_data(const OT& ot, void *post, T *& p, const INDT& ind)
		{
			index_shift_applyer_base_skip_slave<T, OT, INDT, POS, INUM, SC, IND...>::shift_data(ot,p,ind);
		}
	};

	template <typename ST, typename ... IND>
	struct defaultIndices;

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES, typename ... IND>
	struct defaultIndices <stuple<SNUM, CONT, IS_INDEXED, OST, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>, IND...>:
		public defaultIndices <stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, IND..., defaultIndex<V_TYPE>>
	{
	};
	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... IND>
	struct defaultIndices <stuple<SNUM, CONT, IS_INDEXED, OST>, IND...>
	{
		using type=std::tuple<IND...>;
	};

	template <size_t SNUM, size_t CONT, typename CNSC, typename NSC, typename OST, typename NPOS, size_t POS, typename NPOS_CHK>
	struct CORRECT_USE_ALSO;

	template <size_t SNUM, size_t CONT, typename ... CNSC, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... NSC, typename OST, typename NPOS, size_t POS, size_t N_USE_ALSO, size_t ... NPOS_CHK>
	struct CORRECT_USE_ALSO<SNUM, CONT, type_sequence<CNSC...>, type_sequence<sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, NSC...>, OST, NPOS, POS, index_sequence<N_USE_ALSO, NPOS_CHK...> >:
		public CORRECT_USE_ALSO<SNUM, CONT, type_sequence<CNSC...>, type_sequence<sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, NSC...>, OST, NPOS, POS+1, index_sequence<NPOS_CHK...> >
	{
	};


	template <size_t SNUM, size_t CONT, typename ... CNSC, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... NSC, typename OST, typename NPOS, size_t N_USE_ALSO, size_t ... NPOS_CHK>
	struct CORRECT_USE_ALSO<SNUM, CONT, type_sequence<CNSC...>, type_sequence<sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, NSC...>, OST, NPOS, USE_ALSO, index_sequence<N_USE_ALSO, NPOS_CHK...> >:
		public CORRECT_USE_ALSO<(USAGE==USAGE_SLAVE?SNUM:SNUM+1), (USAGE==USAGE_FULL?CONT:USAGE==USAGE_CONT?(CONT>sizeof...(NSC)+1?sizeof...(NSC)+1:CONT):(CONT>sizeof...(NSC)?sizeof...(NSC):CONT)), type_sequence<CNSC...,sClass<V_TYPE, USAGE, N_USE_ALSO, PARENT> >, type_sequence<NSC...>, OST, NPOS, 0, NPOS>
	{
	};

	template <size_t SNUM, size_t CONT, typename ... CNSC, typename OST, typename NPOS, size_t POS, typename NPOS_CHK>
	struct CORRECT_USE_ALSO<SNUM, CONT, type_sequence<CNSC...>, type_sequence<>, OST, NPOS, POS, NPOS_CHK>
	{
		typedef stuple<SNUM, CONT, true, OST, CNSC...> type;
	};



	template <typename ST, typename NSC, size_t CPOS, size_t INUM, typename NPOS, typename POS, size_t CONT, typename OST, typename VN, typename VN_CHK, bool PREV_REMOVED, typename ... SHAPES>
	struct valence_remover_base;


	template <typename ST, typename ... NSC, size_t CPOS, size_t INUM, size_t ... NPOS, size_t ... POS, size_t CONT, typename OST, int ... VN, int VNC, int ... VN_CHK, bool PREV_REMOVED, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct valence_remover_base<ST, type_sequence<NSC...>, CPOS, INUM, index_sequence<NPOS...>, index_sequence<POS...>, CONT, OST, tpp::integer_sequence<int, VN...>, tpp::integer_sequence<int, VNC, VN_CHK...>, PREV_REMOVED, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>:
		public valence_remover_base<ST, type_sequence<NSC...>, CPOS, INUM, index_sequence<NPOS...>, index_sequence<POS...>, CONT, OST, tpp::integer_sequence<int, VN...>, tpp::integer_sequence<int, VN_CHK...>, PREV_REMOVED, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>
	{
	};

	template <typename ST, typename ... NSC, size_t CPOS, size_t INUM, size_t ... NPOS, size_t ... POS, size_t CONT, typename OST, int ... VN, bool PREV_REMOVED, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct valence_remover_base<ST, type_sequence<NSC...>, CPOS, INUM, index_sequence<NPOS...>, index_sequence<POS...>, CONT, OST, tpp::integer_sequence<int, VN...>, tpp::integer_sequence<int>, PREV_REMOVED, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>:
		public valence_remover_base<ST, type_sequence<NSC..., sClass<V_TYPE, (PREV_REMOVED && USAGE==USAGE_FULL?USAGE_CONT:USAGE), USE_ALSO, PARENT> >, CPOS+1, INUM, tpp::index_sequence<NPOS..., sizeof...(NSC)>, tpp::index_sequence<POS..., CPOS>, CONT, OST, tpp::integer_sequence<int, VN...>, tpp::integer_sequence<int, VN...>, false, SHAPES...>
	{
	};

	template <typename ST, typename ... NSC, size_t CPOS, size_t INUM, size_t... NPOS, size_t ... POS, size_t CONT, typename OST, int ... VN, int ... VN_CHK, bool PREV_REMOVED, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct valence_remover_base<ST, type_sequence<NSC...>, CPOS, INUM, index_sequence<NPOS...>, index_sequence<POS...>, CONT, OST, tpp::integer_sequence<int, VN...>, tpp::integer_sequence<int, V_TYPE, VN_CHK...>, PREV_REMOVED, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>:
		public valence_remover_base<ST, type_sequence<NSC...>, CPOS+1, INUM+1, tpp::index_sequence<NPOS..., 0>, tpp::index_sequence<POS...>, (CONT<sizeof...(SHAPES)?CONT:sizeof...(SHAPES)), OST, tpp::integer_sequence<int, VN...>, tpp::integer_sequence<int, VN...>, true, SHAPES...>
	{
		typedef valence_remover_base<ST, type_sequence<NSC...>, CPOS+1, INUM+1, tpp::index_sequence<NPOS..., 0>, tpp::index_sequence<POS...>, (CONT<sizeof...(SHAPES)?CONT:sizeof...(SHAPES)), OST, tpp::integer_sequence<int, VN...>, tpp::integer_sequence<int, VN...>, true, SHAPES...> base;
		template <typename T>
		static void apply_numeric_index(ST& st, T *&p, const size_t *o)
		{
			static size_t inum=find<int, V_TYPE>(tpp::integer_sequence<int, VN...>{});
			if (o[inum]>=get<CPOS>(st).length())
				throw outOfBounds(o[inum], sizeof...(POS), get<CPOS>(st).length());
			shapes_caller<T, CPOS, ST>::apply_numeric_index(st, p, o[inum]);
			base::apply_numeric_index(st, p, o);
		}
	};

	template <typename ST, typename ... NSC, size_t CPOS, size_t INUM, size_t ... NPOS, size_t ... POS, size_t CONT, typename OST, int ... VN, int ... VN_CHK, bool PREV_REMOVED>
	struct valence_remover_base<ST, type_sequence<NSC...>, CPOS, INUM, tpp::index_sequence<NPOS...>, tpp::index_sequence<POS...>, CONT, OST, tpp::integer_sequence<int, VN...>, tpp::integer_sequence<int, VN_CHK...>, PREV_REMOVED>:
		public CORRECT_USE_ALSO<0, CONT, type_sequence<>, type_sequence<NSC...>, OST, tpp::index_sequence<NPOS...>, 0, tpp::index_sequence<NPOS...> >
	{
		static_assert(INUM==sizeof...(VN), "Cannot remove incorrect valence");
		template <typename T>
		static void apply_numeric_index(ST& st, T *&p, const size_t *o) noexcept {}
		template <size_t I>
		inline static constexpr const typename ST::template element_type<get<I>(tpp::index_sequence<POS...>{})>& shape(const ST& st) noexcept
		{
			return get<get<I>(tpp::index_sequence<POS...>{})>(st);
		}
	};

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, int ... V_TYPE>
	struct valence_remover<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, V_TYPE...>:
		public valence_remover_base<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, type_sequence<>, 0, 0, tpp::index_sequence<>, tpp::index_sequence<>, CONT, OST, tpp::integer_sequence<int, V_TYPE...>, tpp::integer_sequence<int, V_TYPE...>, false, SHAPES...>
	{
		static_assert(is_unique(tpp::integer_sequence<int, V_TYPE...>{}),"Valences must be unique");
		typedef valence_remover_base<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, type_sequence<>, 0, 0, tpp::index_sequence<>, tpp::index_sequence<>, CONT, OST, tpp::integer_sequence<int, V_TYPE...>, tpp::integer_sequence<int, V_TYPE...>, false, SHAPES...> base;
		static typename base::type create_stuple(const stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>&st) { return typename base::type(st, tpp::integer_sequence<int, V_TYPE...>{}); }
	};

};

#endif /* STUPLE_H_ */
