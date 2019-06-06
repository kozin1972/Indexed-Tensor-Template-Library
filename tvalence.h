/*
 * valence.h
 *
 *  Created on: 7 мая 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef VALENCE_H_
#define VALENCE_H_

#include <stuple.h>

namespace tpp
{

	template <int SRC, int V_TYPE, bool IS_CONTINUOUS, bool IS_VECTOR, size_t ... POS>
	struct shape_info
	{
		static const int src=SRC;
		static const int v_type=V_TYPE;
		static const bool is_continuous=IS_CONTINUOUS;
		static const bool is_vector=IS_VECTOR;
	};

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER>
	struct valence_info
	{
		static const int v_type=V_TYPE;
		static const int mask=MASK;
		static const int c_mask=C_MASK;
		static const int v_mask=V_MASK;
		static const bool is_joinable=IS_JOINABLE;
		static const int join_with=JOIN_WITH;
		static const size_t order=ORDER;
	};

	template <int V_TYPE, int SRC, bool IS_CONTINUOUS, bool IS_VECTOR, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename CVI, typename VI>
	struct VI_updater;

	template <int V_TYPE, int SRC, bool IS_CONTINUOUS, bool IS_VECTOR, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... CVTS, typename HEAD, typename ... VTS>
	struct VI_updater<V_TYPE, SRC, IS_CONTINUOUS, IS_VECTOR, IS_JOINABLE, JOIN_WITH, ORDER, type_pack<CVTS...>, type_pack<HEAD, VTS...> >
	{
		using type=
			typename std::conditional<(V_TYPE<HEAD::v_type),
				type_pack<CVTS..., valence_info<V_TYPE,(1<<SRC),(IS_CONTINUOUS?(1<<SRC):0),(IS_VECTOR?(1<<SRC):0), IS_JOINABLE,JOIN_WITH,ORDER>, HEAD, VTS...>,
				typename std::conditional<(V_TYPE==HEAD::v_type),
					typename std::conditional<(IS_JOINABLE && HEAD::is_joinable && JOIN_WITH==HEAD::join_with),
						type_pack<CVTS..., valence_info<V_TYPE,HEAD::mask | (1<<SRC),HEAD::c_mask | (IS_CONTINUOUS?(1<<SRC):0),HEAD::v_mask | (IS_VECTOR?(1<<SRC):0),IS_JOINABLE,JOIN_WITH,(SRC==1?ORDER:((HEAD::mask & 2)?HEAD::order:(SRC==2?ORDER:HEAD::order)))>, VTS...>,
						type_pack<CVTS..., valence_info<V_TYPE,HEAD::mask | (1<<SRC),HEAD::c_mask | (IS_CONTINUOUS?(1<<SRC):0),HEAD::v_mask | (IS_VECTOR?(1<<SRC):0),false,0,(SRC==1?ORDER:((HEAD::mask & 2)?HEAD::order:(SRC==2?ORDER:HEAD::order)))>, VTS...>
					>::type,
					typename VI_updater<V_TYPE, SRC, IS_CONTINUOUS, IS_VECTOR, IS_JOINABLE, JOIN_WITH, ORDER, type_pack<CVTS..., HEAD>, type_pack<VTS...> >::type
				>::type
			>::type;
	};

	template <int V_TYPE, int SRC, bool IS_CONTINUOUS, bool IS_VECTOR, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... CVTS>
	struct VI_updater<V_TYPE, SRC, IS_CONTINUOUS, IS_VECTOR, IS_JOINABLE, JOIN_WITH, ORDER, type_pack<CVTS...>, type_pack<> >
	{
		using type=type_pack<CVTS..., valence_info<V_TYPE,(1<<SRC),(IS_CONTINUOUS?(1<<SRC):0),(IS_VECTOR?(1<<SRC):0),IS_JOINABLE,JOIN_WITH,ORDER> >;
	};

	template <int V_TYPE, typename VI>
	struct VI_getter;

	template <int V_TYPE, int C_V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VTS>
	struct VI_getter<V_TYPE, type_pack<valence_info<C_V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VTS...> >:
		public VI_getter<V_TYPE, type_pack<VTS...> >
	{
	};

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... VTS>
	struct VI_getter<V_TYPE, type_pack<valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, VTS...> >
	{
		using valence_data=valence_info<V_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>;
	};

	template <int V_TYPE>
	struct VI_getter<V_TYPE, type_pack<> >
	{
		static_assert(true,"Internal error.");
	};

	template <typename VI, typename DST>
	struct VI_corrector;

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, int JOIN_WITH, size_t ORDER, typename ... VTS, typename ... DSTV>
	struct VI_corrector<type_pack<valence_info<V_TYPE, MASK, C_MASK, V_MASK, true, JOIN_WITH, ORDER>, VTS...>, type_pack<DSTV...> >:
		public VI_corrector<type_pack<VTS...>, type_pack<DSTV..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, MASK==VI_getter<JOIN_WITH,type_pack<VTS...,DSTV...> >::valence_data::mask, JOIN_WITH, ORDER> > >
	{
	};
	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, int JOIN_WITH, size_t ORDER, typename ... VTS, typename ... DSTV>
	struct VI_corrector<type_pack<valence_info<V_TYPE, MASK, C_MASK, V_MASK, false, JOIN_WITH, ORDER>, VTS...>, type_pack<DSTV...> >:
		public VI_corrector<type_pack<VTS...>, type_pack<DSTV..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, false, JOIN_WITH, ORDER> > >
	{
	};

	template <typename ... DSTV>
	struct VI_corrector<type_pack<>, type_pack<DSTV...> >
	{
		using type=type_pack<DSTV...>;
	};

	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename TS, typename ... SHAPES>
	struct make_shape_info_base;

// SLAVE
	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename ... SEQ, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, size_t USE_ALSO, size_t PARENT, template <int, enum dSU, size_t, size_t> class N_sClass, int N_V_TYPE, enum dSU N_USAGE, size_t N_USE_ALSO, size_t N_PARENT, typename ... SHAPES>
	struct make_shape_info_base<FIRST_CONT, POS, SRC, SNUM, VI, type_pack<SEQ...>, sClass<V_TYPE, USAGE_SLAVE, USE_ALSO, PARENT>, N_sClass<N_V_TYPE, N_USAGE, N_USE_ALSO, N_PARENT>, SHAPES...>:
		public make_shape_info_base<FIRST_CONT, POS+1, SRC, SNUM, VI, type_pack<SEQ...>, N_sClass<N_V_TYPE, N_USAGE, N_USE_ALSO, N_PARENT>, SHAPES...>
	{
	};

// LAST SLAVE
	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename ... SEQ, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, size_t USE_ALSO, size_t PARENT>
	struct make_shape_info_base<FIRST_CONT, POS, SRC, SNUM, VI, type_pack<SEQ...>, sClass<V_TYPE, USAGE_SLAVE, USE_ALSO, PARENT> >:
		public make_shape_info_base<FIRST_CONT, POS+1, SRC, SNUM, VI, type_pack<SEQ...> >
	{
	};

// COMMON
	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename ... SEQ, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, template <int, enum dSU, size_t, size_t> class N_sClass, int N_V_TYPE, enum dSU N_USAGE, size_t N_USE_ALSO, size_t N_PARENT, typename ... SHAPES>
	struct make_shape_info_base<FIRST_CONT, POS, SRC, SNUM, VI, type_pack<SEQ...>, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, N_sClass<N_V_TYPE, N_USAGE, N_USE_ALSO, N_PARENT>, SHAPES...>:
		public make_shape_info_base<FIRST_CONT, POS+1, SRC, SNUM, typename VI_updater<V_TYPE, SRC, false, (USAGE<=USAGE_CONT)?true:false, false, 0, sizeof...(SHAPES)+1, type_pack<>, VI>::type, type_pack<SEQ..., shape_info<SRC, V_TYPE, false, (USAGE<=USAGE_CONT)?true:false, POS> >, N_sClass<N_V_TYPE, N_USAGE, N_USE_ALSO, N_PARENT>, SHAPES...>
	{
	};

// LAST COMMON
	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename ... SEQ, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	struct make_shape_info_base<FIRST_CONT, POS, SRC, SNUM, VI, type_pack<SEQ...>, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>>:
		public make_shape_info_base<FIRST_CONT, POS+1, SRC, SNUM, typename VI_updater<V_TYPE, SRC, (USAGE<=USAGE_CONT && POS>=FIRST_CONT)?true:false,(USAGE<=USAGE_CONT)?true:false, false, 0, 0, type_pack<>, VI>::type, type_pack<SEQ..., shape_info<SRC, V_TYPE, (USAGE<=USAGE_CONT && POS>=FIRST_CONT)?true:false,(USAGE<=USAGE_CONT)?true:false,POS> > >
	{
	};

// FULL
	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename ... SEQ, int V_TYPE, template <int, enum dSU, size_t, size_t> class N_sClass, int N_V_TYPE, size_t N_USE_ALSO, size_t N_PARENT, typename ... SHAPES>
	struct make_shape_info_base<FIRST_CONT, POS, SRC, SNUM, VI, type_pack<SEQ...>, segment<V_TYPE, USAGE_FULL, 0, 0>, N_sClass<N_V_TYPE, USAGE_FULL, N_USE_ALSO, N_PARENT>, SHAPES...>:
		public make_shape_info_base<FIRST_CONT, POS+1, SRC, SNUM, typename VI_updater<V_TYPE, SRC, false, true, true, N_V_TYPE, sizeof...(SHAPES)+1, type_pack<>, VI>::type, type_pack<SEQ..., shape_info<SRC, V_TYPE, false, true, POS> >, N_sClass<N_V_TYPE, USAGE_FULL, N_USE_ALSO, N_PARENT>, SHAPES...>
	{
	};

// LAST FULL
	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename ... SEQ, int V_TYPE>
	struct make_shape_info_base<FIRST_CONT, POS, SRC, SNUM, VI, type_pack<SEQ...>, segment<V_TYPE, USAGE_FULL, 0, 0>>:
		public make_shape_info_base<FIRST_CONT, POS+1, SRC, SNUM, typename VI_updater<V_TYPE, SRC, (POS>=FIRST_CONT), true, false, 0, 0, type_pack<>, VI>::type, type_pack<SEQ..., shape_info<SRC, V_TYPE, (POS>=FIRST_CONT), true, POS> > >
	{
	};

// CONT
	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename ... SEQ, int V_TYPE, template <int, enum dSU, size_t, size_t> class N_sClass, int N_V_TYPE, size_t N_USE_ALSO, size_t N_PARENT, typename ... SHAPES>
	struct make_shape_info_base<FIRST_CONT, POS, SRC, SNUM, VI, type_pack<SEQ...>, segment<V_TYPE, USAGE_CONT, 0, 0>, N_sClass<N_V_TYPE, USAGE_FULL, N_USE_ALSO, N_PARENT>, SHAPES...>:
		public make_shape_info_base<FIRST_CONT, POS+1, SRC, SNUM, typename VI_updater<V_TYPE, SRC, false, true, true, N_V_TYPE, sizeof...(SHAPES)+1, type_pack<>, VI>::type, type_pack<SEQ..., shape_info<SRC, V_TYPE, false, true, POS> >, N_sClass<N_V_TYPE, USAGE_FULL, N_USE_ALSO, N_PARENT>, SHAPES...>
	{
	};

// LAST CONT
	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename ... SEQ, int V_TYPE>
	struct make_shape_info_base<FIRST_CONT, POS, SRC, SNUM, VI, type_pack<SEQ...>, segment<V_TYPE, USAGE_CONT, 0, 0>>:
		public make_shape_info_base<FIRST_CONT, POS+1, SRC, SNUM, typename VI_updater<V_TYPE, SRC, (POS>=FIRST_CONT), true, false, 0, 0, type_pack<>, VI>::type, type_pack<SEQ..., shape_info<SRC, V_TYPE, (POS>=FIRST_CONT), true, POS> > >
	{
	};

// ENUM segment
	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename ... SEQ, int V_TYPE, template <int, enum dSU, size_t, size_t> class N_sClass, int N_V_TYPE, enum dSU N_USAGE, size_t N_USE_ALSO, size_t N_PARENT, typename ... SHAPES>
	struct make_shape_info_base<FIRST_CONT, POS, SRC, SNUM, VI, type_pack<SEQ...>, segment<V_TYPE, USAGE_ENUM, 0, 0>, N_sClass<N_V_TYPE, N_USAGE, N_USE_ALSO, N_PARENT>, SHAPES...>:
		public make_shape_info_base<FIRST_CONT, POS+1, SRC, SNUM, typename VI_updater<V_TYPE, SRC, false, true, false, 0, sizeof...(SHAPES)+1, type_pack<>, VI>::type, type_pack<SEQ..., shape_info<SRC, V_TYPE, false, true, POS> >, N_sClass<N_V_TYPE, N_USAGE, N_USE_ALSO, N_PARENT>, SHAPES...>
	{
	};

// LAST ENUM segment
	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename ... SEQ, int V_TYPE>
	struct make_shape_info_base<FIRST_CONT, POS, SRC, SNUM, VI, type_pack<SEQ...>, segment<V_TYPE, USAGE_ENUM, 0, 0>>:
		public make_shape_info_base<FIRST_CONT, POS+1, SRC, SNUM, typename VI_updater<V_TYPE, SRC, (POS>=FIRST_CONT), true, false, 0, 0, type_pack<>, VI>::type, type_pack<SEQ..., shape_info<SRC, V_TYPE, (POS>=FIRST_CONT), true, POS> > >
	{
	};

	template <size_t FIRST_CONT, int POS, int SRC, size_t SNUM, typename VI, typename ... SEQ>
	struct make_shape_info_base<FIRST_CONT, POS, SRC, SNUM, VI, type_pack<SEQ...> >
	{
		using type=type_pack<typename VI_corrector<VI,type_pack<> >::type, type_pack<SEQ...> >;
	};


	template <int SRC, typename ST, typename VI, typename TS>
	struct make_shape_info;

	template <int SRC, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, typename VI, typename TS>
	struct make_shape_info<SRC, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, VI, TS>:
		public make_shape_info_base<sizeof...(SHAPES)-CONT, 0, SRC, SNUM, VI, TS, SHAPES...>
	{
	};

	template <bool FOUND, int V_TYPE, int JOIN_WITH>
	struct find_result
	{
		static const bool found=FOUND;
		static const int v_type=V_TYPE;
		static const int join_with=JOIN_WITH;
	};

	template <int RV_TYPE, typename VI, typename ... TSVI>
	struct remove_valence_info;

	template <int RV_TYPE, typename ... CVI, int MASK, int C_MASK, int V_MASK, bool IS_JOINABLE, int JOIN_WITH, size_t ORDER, typename ... TVI>
	struct remove_valence_info<RV_TYPE, type_pack<CVI...>, valence_info<RV_TYPE, MASK, C_MASK, V_MASK, IS_JOINABLE, JOIN_WITH, ORDER>, TVI...>:
		public remove_valence_info<RV_TYPE, type_pack<CVI...>, TVI...>
	{
	};

	template <int RV_TYPE, typename ... CVI,  int V_TYPE, int MASK, int C_MASK, int V_MASK, size_t ORDER, typename ... TVI>
	struct remove_valence_info<RV_TYPE, type_pack<CVI...>, valence_info<V_TYPE, MASK, C_MASK, V_MASK, true, RV_TYPE, ORDER>, TVI...>:
		public remove_valence_info<RV_TYPE, type_pack<CVI..., valence_info<V_TYPE, MASK, C_MASK, V_MASK, false, RV_TYPE, ORDER> >, TVI...>
	{
	};

	template <int RV_TYPE, typename ... CVI,  typename HEAD, typename ... TVI>
	struct remove_valence_info<RV_TYPE, type_pack<CVI...>, HEAD, TVI...>:
		public remove_valence_info<RV_TYPE, type_pack<CVI..., HEAD>, TVI...>
	{
	};

	template <int RV_TYPE, typename ... CVI>
	struct remove_valence_info<RV_TYPE, type_pack<CVI...> >
	{
		using type=type_pack<CVI...>;
	};

	template <int V_TYPE, typename TS, typename ... SEQ>
	struct join_shape;

	template <int V_TYPE, typename ... RSEQ, int SRC, bool IS_CONTINUOUS, bool IS_VECTOR, size_t POS, int JOIN_WITH, bool N_IS_CONTINUOUS, bool N_IS_VECTOR, size_t... TPOS, typename ... SEQ>
	struct join_shape<V_TYPE, type_pack<RSEQ...>, shape_info<SRC, V_TYPE, IS_CONTINUOUS, IS_VECTOR, POS>, shape_info<SRC, JOIN_WITH, N_IS_CONTINUOUS, N_IS_VECTOR, TPOS...>, SEQ...>:
		public join_shape<V_TYPE, type_pack<RSEQ...,shape_info<SRC, V_TYPE, N_IS_CONTINUOUS, IS_VECTOR, POS, TPOS...> >, SEQ...>
	{
	};

	template <int V_TYPE, typename ... RSEQ, int SRC_0, int V_TYPE_0, bool IS_CONTINUOUS_0, bool IS_VECTOR_0, size_t ... POS_0, int SRC_1, int V_TYPE_1, bool IS_CONTINUOUS_1, bool IS_VECTOR_1, size_t ... POS_1, typename ... SEQ>
	struct join_shape<V_TYPE, type_pack<RSEQ...>, shape_info<SRC_0, V_TYPE_0, IS_CONTINUOUS_0, IS_VECTOR_0, POS_0...>, shape_info<SRC_1, V_TYPE_1, IS_CONTINUOUS_1, IS_VECTOR_1, POS_1...>, SEQ...>:
		public join_shape<V_TYPE, type_pack<RSEQ..., shape_info<SRC_0, V_TYPE_0, IS_CONTINUOUS_0, IS_VECTOR_0, POS_0...> >, shape_info<SRC_1, V_TYPE_1, IS_CONTINUOUS_1, IS_VECTOR_1, POS_1...>, SEQ...>
	{
	};

	template <int V_TYPE, typename ... RSEQ, int SRC, int CV_TYPE, bool IS_CONTINUOUS, bool IS_VECTOR, size_t ... TPOS>
	struct join_shape<V_TYPE, type_pack<RSEQ...>, shape_info<SRC, CV_TYPE, IS_CONTINUOUS, IS_VECTOR, TPOS...> >
	{
		using type=type_pack<RSEQ...,shape_info<SRC, CV_TYPE, IS_CONTINUOUS, IS_VECTOR, TPOS...> >;
	};

	template <int V_TYPE, typename ... RSEQ>
	struct join_shape<V_TYPE, type_pack<RSEQ...> >
	{
		using type=type_pack<RSEQ...>;
	};

	template <typename VI, typename TS, typename FR>
	struct shape_joiner_base;

	template <typename VI, typename TS>
	struct shape_joiner;

	template <typename VI, typename TS, int V_TYPE, int JOIN_WITH>
	struct join_valence;

	template <typename ... TSVI, typename ... SEQ, int V_TYPE, int JOIN_WITH>
	struct join_valence<type_pack<TSVI...>, type_pack<SEQ...>, V_TYPE, JOIN_WITH>:
		public shape_joiner<typename remove_valence_info<JOIN_WITH, type_pack<>, TSVI...>::type, typename join_shape<V_TYPE, type_pack<>, SEQ...>::type>
	{
	};

	template <typename VI, typename TS, int V_TYPE, int JOIN_WITH>
	struct shape_joiner_base<VI, TS, find_result<true, V_TYPE, JOIN_WITH> >:
		public join_valence<VI, TS, V_TYPE, JOIN_WITH>
	{
	};

	template <typename VI, typename TS, int V_TYPE, int JOIN_WITH>
	struct shape_joiner_base<VI, TS, find_result<false, V_TYPE, JOIN_WITH> >
	{
//		using valence_info_seq=VI;
//		using shape_info_seq=TS;
		using type=type_pack<VI, TS>;
	};

	template <typename FVI, typename VI>
	struct find_last_valence;

	template <typename FVI, int V_TYPE, int JOIN_WITH, bool WITH_JOINABLE, typename ... TVI>
	struct test_joinable;

	template <typename FVI, int V_TYPE, int JOIN_WITH, typename ... TVI>
	struct test_joinable<FVI, V_TYPE, JOIN_WITH, true, TVI...>:
		public find_last_valence<FVI, type_pack<TVI...> >
	{
	};

	template <typename FVI, int V_TYPE, int JOIN_WITH>
	struct test_joinable<FVI, V_TYPE, JOIN_WITH, true>
	{
		using type=find_result<false, V_TYPE, JOIN_WITH>;
	};

	template <typename FVI, int V_TYPE, int JOIN_WITH, typename ... TVI>
	struct test_joinable<FVI, V_TYPE, JOIN_WITH, false, TVI...>
	{
		using type=find_result<true, V_TYPE, JOIN_WITH>;
	};

	template <typename FVI, int V_TYPE, int MASK, int C_MASK, int V_MASK, int JOIN_WITH, size_t ORDER, typename ... TVI>
	struct find_last_valence<FVI, type_pack<valence_info<V_TYPE, MASK, C_MASK, V_MASK, true, JOIN_WITH, ORDER>, TVI...> >:
		public test_joinable<FVI, V_TYPE, JOIN_WITH, VI_getter<JOIN_WITH, FVI>::valence_data::is_joinable, TVI...>
	{
	};

	template <typename FVI, int V_TYPE, int MASK, int C_MASK, int V_MASK, int JOIN_WITH, size_t ORDER, typename ... TVI>
	struct find_last_valence<FVI, type_pack<valence_info<V_TYPE, MASK, C_MASK, V_MASK, false, JOIN_WITH, ORDER>, TVI...> >:
		public find_last_valence<FVI, type_pack<TVI...> >
	{
	};

	template <typename FVI>
	struct find_last_valence<FVI, type_pack<> >
	{
		using type=find_result<false, 0, 0>;
	};

	template <typename VI, typename TS>
	struct shape_joiner: public shape_joiner_base<VI, TS, typename find_last_valence<VI, VI>::type>
	{
	};

	template <int SRC, typename PTS, typename ... TST>
	struct make_valence_info_base;

	template <int SRC, typename VI, typename TS, typename ST, typename ... TST>
	struct make_valence_info_base<SRC, type_pack<VI, TS>, ST, TST...>:
		public make_valence_info_base<SRC+1, typename make_shape_info<SRC, ST, VI, TS>::type, TST...>
	{
	};

	template <int SRC, typename VI, typename TS, typename ST>
	struct make_valence_info_base<SRC, type_pack<VI, TS>, ST>:
		public make_shape_info<SRC, ST, VI, TS>
	{
	};

	template <int SRC, typename VI, typename TS>
	struct make_valence_info_base<SRC, type_pack<VI, TS> >
	{
		using type=type_pack<type_pack<>, type_pack<> >;
	};

	template <typename ... ST>
	struct make_valence_info:
		public make_valence_info_base<0, type_pack<type_pack<>, type_pack<> >, ST...>
	{
	};

	template <int SRC, typename PTS, typename ... TST>
	struct make_valence_info_and_join_base;

	template <int SRC, typename VI, typename TS, typename ST, typename ... TST>
	struct make_valence_info_and_join_base<SRC, type_pack<VI, TS>, ST, TST...>:
		public make_valence_info_and_join_base<SRC+1, typename make_shape_info<SRC, ST, VI, TS>::type, TST...>
	{
	};

	template <int SRC, typename VI, typename TS>
	struct make_valence_info_and_join_base<SRC, type_pack<VI, TS> >:
		public shape_joiner<VI, TS>
	{
	};

	template <typename ... ST>
	struct make_valence_info_and_join:
		public make_valence_info_and_join_base<0, type_pack<type_pack<>, type_pack<>>, ST...>
	{
	};

	template <int V_TYPE, int SRC, size_t POS>
	struct length_info
	{
	};

	template <int SRC, size_t POS, typename ST, typename TS, typename SH, typename ... RST>
	struct test_length;

	template <int SRC, size_t POS, int V_TYPE, typename ST, typename TS, typename DTS, typename SH, typename ... RST>
	struct test_v_type;

	template <int SRC, typename ST, typename TS, typename ... RST>
	struct test_shape_length_base;

// All length_info are searched
	template <int SRC, size_t POS, int V_TYPE, typename ST, typename ... DTS, typename SH, typename ... RST>
	struct test_v_type<SRC, POS, V_TYPE, ST, type_pack<>, type_pack<DTS...>, SH, RST...>:
		public test_length<SRC, POS+1, ST, type_pack<DTS..., length_info<V_TYPE, SRC, POS> >, SH, RST...>
	{
	};

// V_TYPE is different
	template <int SRC, size_t POS, int V_TYPE, typename ST, int S_V_TYPE, int S_SRC, size_t S_POS, typename ... TS, typename ... DTS, typename SH, typename ... RST>
	struct test_v_type<SRC, POS, V_TYPE, ST, type_pack<length_info<S_V_TYPE, S_SRC, S_POS>, TS...>, type_pack<DTS...>, SH, RST...>:
		public test_v_type<SRC, POS, V_TYPE, ST, type_pack<TS...>, type_pack<DTS...,length_info<S_V_TYPE, S_SRC, S_POS> >, SH, RST...>
	{
	};

// We found the same VALENCE
	template <int SRC, size_t POS, int V_TYPE, typename ST, int S_SRC, size_t S_POS, typename ... TS, typename ... DTS, typename SH, typename ... RST>
	struct test_v_type<SRC, POS, V_TYPE, ST, type_pack<length_info<V_TYPE, S_SRC, S_POS>, TS...>, type_pack<DTS...>, SH, RST...>:
		public test_v_type<SRC, POS, V_TYPE, ST, type_pack<TS...>, type_pack<DTS..., length_info<V_TYPE, S_SRC, S_POS> >, SH, RST...>
	{
		static void test(const ST& st)
		{
			if (get<POS>(std::get<SRC>(st)).length()!=get<S_POS>(std::get<S_SRC>(st)).length())
				throw incompatibleIndices(S_POS, S_SRC, POS, SRC, V_TYPE, get<S_POS>(std::get<S_SRC>(st)).length(), get<POS>(std::get<SRC>(st)).length());
			test_v_type<SRC, POS, V_TYPE, ST, type_pack<TS...>, type_pack<DTS..., length_info<V_TYPE, S_SRC, S_POS> >, SH, RST...>::test(st);
		}
	};

// Next shape
	template <int SRC, size_t POS, typename ST, typename TS, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES, typename ... RST>
	struct test_length<SRC, POS, ST, TS, type_pack<sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>, RST...>:
		public test_v_type<SRC, POS, V_TYPE, ST, TS, type_pack<>, type_pack<SHAPES...>, RST...>
	{

	};

// Next tensor
	template <int SRC, size_t POS, typename ST, typename TS, typename ... RST>
	struct test_length<SRC, POS, ST, TS, type_pack<>, RST...>:
		public test_shape_length_base<SRC+1, ST, TS, RST...>
	{
	};

	template <int SRC, typename ST, typename TS, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, typename ... RST>
	struct test_shape_length_base<SRC, ST, TS, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, RST...>:
		public test_length<SRC, 0, ST, TS, type_pack<SHAPES...>, RST...>
	{

	};

	template <int SRC, typename ST, typename TS>
	struct test_shape_length_base<SRC, ST, TS>
	{
		static void test(const ST& st) {}
	};

	template <typename ... ST>
	struct test_shape_length: public test_shape_length_base<0, std::tuple<ST...>, type_pack<>, ST...>
	{
	};

	template <typename T, typename ST, size_t ... POS>
	class joined_iterator;
	template <typename T, typename ST, size_t POS0, size_t ... POS>
	class joined_iterator<T, ST, POS0, POS...>
	{
		size_t s;
		const T *end;
	public:
		static size_t length(const ST& st) { return joined_iterator<T, ST, POS...>::length(st)*get<POS0>(st).length(); }
		static size_t step(const ST& st) { return joined_iterator<T, ST, POS...>::step(st); } //return get<POS0>(st).step(); }
		joined_iterator(const ST& st, const T * const p): s(step(st)), end(p+s*length(st)) {}
		void move_one(T *&p) { p+=s; }
		void move_one(const T *&p) { p+=s; }
		bool not_end(const T *p) { return p<end; }
	};
	template <typename T, typename ST, size_t POS>
	class joined_iterator<T, ST, POS>
	{
		size_t s;
		T *end;
	public:
		static size_t length(const ST& st) { return get<POS>(st).length(); }
		static size_t step(const ST& st) { return get<POS>(st).step(); }
	};

	template <typename ST, typename SI, int SRC, int V_TYPE>
	struct iterator_getter_by_src_v_type;

	template <typename ST, typename HEAD, typename ... SI, int SRC, int V_TYPE>
	struct iterator_getter_by_src_v_type<ST, type_pack<HEAD, SI...>, SRC, V_TYPE>:
		public iterator_getter_by_src_v_type<ST, type_pack<SI...>, SRC, V_TYPE>
	{
	};

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, bool IS_CONTINUOUS, bool IS_VECTOR, size_t POS, typename ... SI, int SRC, int V_TYPE>
	struct iterator_getter_by_src_v_type<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, type_pack<shape_info<SRC, V_TYPE, IS_CONTINUOUS, IS_VECTOR, POS>, SI...>, SRC, V_TYPE>
	{
		template <typename T>
		using iterator_type=shapes_iterator_base<T, POS, (CONT>0 && POS==sizeof...(SHAPES)-1), typename std::tuple_element<POS, std::tuple<SHAPES...> >::type, check_shape<typename std::tuple_element<POS, std::tuple<SHAPES...> >::type>::need_parent, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >;
	};

	template <typename ST, bool IS_CONTINUOUS, size_t POS0, size_t ... POS, typename ... SI, int SRC, int V_TYPE>
	struct iterator_getter_by_src_v_type<ST, type_pack<shape_info<SRC, V_TYPE, IS_CONTINUOUS, true, POS0, POS...>, SI...>, SRC, V_TYPE>
	{
		template <typename T>
		using iterator_type=joined_iterator<T, ST, POS0, POS...>;
	};

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, bool IS_CONTINUOUS, size_t POS, typename ... SI, int SRC, int V_TYPE>
	struct iterator_getter_by_src_v_type<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, type_pack<shape_info<SRC, V_TYPE, IS_CONTINUOUS, true, POS>, SI...>, SRC, V_TYPE>
	{
		template <typename T>
		using iterator_type=shapes_iterator_base<T, POS, (CONT>0 && POS==sizeof...(SHAPES)-1), typename std::tuple_element<POS, std::tuple<SHAPES...> >::type, check_shape<typename std::tuple_element<POS, std::tuple<SHAPES...> >::type>::need_parent, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >;
	};


	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t POS0, size_t ... POS>
	struct valence_data
	{
		static const int v_type=V_TYPE;
		static const int mask=MASK;
		static const int c_mask=C_MASK;
		static const int v_mask=V_MASK;
		static const int next_valence=NEXT_VALENCE;
		static const int pos0=POS0;
	};

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t POS0, size_t POS1, size_t ... POS>
	struct valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS0, POS1, POS...>
	{
		static const int v_type=V_TYPE;
		static const int mask=MASK;
		static const int c_mask=C_MASK;
		static const int v_mask=V_MASK;
		static const int next_valence=NEXT_VALENCE;
		static const int pos0=POS0;
		static const int pos1=POS1;
	};

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t POS0, size_t POS1, size_t POS2, size_t ... POS>
	struct valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS0, POS1, POS2, POS...>
	{
		static const int v_type=V_TYPE;
		static const int mask=MASK;
		static const int c_mask=C_MASK;
		static const int v_mask=V_MASK;
		static const int next_valence=NEXT_VALENCE;
		static const int pos0=POS0;
		static const int pos1=POS1;
		static const int pos2=POS2;
	};


	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t POS, typename CVD, typename VD>
	struct add_valence_data;

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t POS, typename ... CVD, int H_V_TYPE, int H_MASK, int H_C_MASK, int H_V_MASK, int H_NEXT_VALENCE, size_t ... H_POS, typename ... VD>
	struct add_valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS, type_pack<CVD...>, type_pack<valence_data<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_NEXT_VALENCE, H_POS...>, VD...> >:
		public
				std::conditional<(V_TYPE==H_V_TYPE),
					type_pack<CVD..., valence_data<V_TYPE, MASK | H_MASK, C_MASK | H_C_MASK, V_MASK | H_V_MASK, (NEXT_VALENCE==H_NEXT_VALENCE?NEXT_VALENCE:V_TYPE), H_POS..., POS>, VD...>,
					typename add_valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS, type_pack<CVD..., valence_data<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_NEXT_VALENCE, H_POS...> >, type_pack<VD...> >::type
				>::type
	{
	};

	template <int V_TYPE, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t POS, typename ... CVD>
	struct add_valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS, type_pack<CVD...>, type_pack<> >
	{
		typedef type_pack<CVD..., valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS> > type;
	};

//	template <typename CVD, typename EVD, typename VD>
//	struct add_valence_data_by_pos;
//
//	template <typename ... CVD, int H_V_TYPE, int H_MASK, int H_C_MASK, int H_V_MASK, size_t H_HEAD_POS, size_t ... H_POS, typename ... EVD, int V_TYPE, int MASK, int C_MASK, int V_MASK, size_t HEAD_POS, size_t ... POS>
//	struct add_valence_data_by_pos<type_pack<CVD...>, type_pack<valence_data<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_HEAD_POS, H_POS...>, EVD...>, valence_data<V_TYPE, MASK, C_MASK, V_MASK, HEAD_POS, POS...> >:
//		public
//			std::conditional<(HEAD_POS<H_HEAD_POS),
//				type_pack<CVD..., valence_data<V_TYPE, MASK, C_MASK, V_MASK, HEAD_POS, POS...>, valence_data<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_HEAD_POS, H_POS...>, EVD...>,
//				typename add_valence_data_by_pos<type_pack<CVD..., valence_data<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_POS...> >, type_pack<EVD...>, valence_data<V_TYPE, MASK, C_MASK, V_MASK, HEAD_POS, POS...> >::type
//			>::type
//	{
//	};
//
//	template <typename ... CVD, int V_TYPE, int MASK, int C_MASK, int V_MASK, size_t ... POS>
//	struct add_valence_data_by_pos<type_pack<CVD...>, type_pack<>, valence_data<V_TYPE, MASK, C_MASK, V_MASK, POS...> >
//	{
//		typedef type_pack<CVD..., valence_data<V_TYPE, MASK, C_MASK, V_MASK, POS...> > type;
//	};

	template <int FIRST_CONT, int SRC, int POS, typename VD, typename ... SHAPES>
	struct add_shapes_to_valence_data;

	template <int FIRST_CONT, int SRC, int POS, typename VD, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, POS, VD, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, POS+1, typename add_valence_data<V_TYPE, (1<<SRC), 0, 0, V_TYPE, POS, type_pack<>, VD>::type, SHAPES...>
	{
	};

	template <int FIRST_CONT, int SRC, int POS, typename VD, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, POS, VD, sClass<V_TYPE, USAGE_SLAVE, USE_ALSO, PARENT>, SHAPES...>:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, POS+1, VD, SHAPES...>
	{
	};

	template <int FIRST_CONT, int SRC, int POS, typename VD, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, POS, VD, segment<V_TYPE, USAGE, USE_ALSO, PARENT> >:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, POS+1, typename add_valence_data<V_TYPE, (1<<SRC), ((USAGE<=USAGE_CONT && POS>=FIRST_CONT)?(1<<SRC):0), ((USE_ALSO==0 && PARENT==0)?(1<<SRC):0), V_TYPE, POS, type_pack<>, VD>::type>
	{
	};

	template <int FIRST_CONT, int SRC, int POS, typename VD, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, template <int, enum dSU, size_t, size_t> class N_sClass, int N_V_TYPE, enum dSU N_USAGE, size_t N_USE_ALSO, size_t N_PARENT, typename ... SHAPES>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, POS, VD, segment<V_TYPE, USAGE, USE_ALSO, PARENT>, N_sClass<N_V_TYPE, N_USAGE, N_USE_ALSO, N_PARENT>, SHAPES...>:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, POS+1, typename add_valence_data<V_TYPE, (1<<SRC), 0, ((USE_ALSO==0 && PARENT==0)?(1<<SRC):0), V_TYPE, POS, type_pack<>, VD>::type, N_sClass<N_V_TYPE, N_USAGE, N_USE_ALSO, N_PARENT>, SHAPES...>
	{
	};

	template <int FIRST_CONT, int SRC, int POS, typename VD, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, int N_V_TYPE, typename ... SHAPES>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, POS, VD, segment<V_TYPE, USAGE, USE_ALSO, PARENT>, segment<N_V_TYPE, USAGE_FULL, 0, 0>, SHAPES...>:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, POS+1, typename add_valence_data<V_TYPE, (1<<SRC), 0, ((USE_ALSO==0 && PARENT==0)?(1<<SRC):0), (USE_ALSO==0 && PARENT==0 && USAGE<=USAGE_CONT?N_V_TYPE:V_TYPE), POS, type_pack<>, VD>::type, segment<N_V_TYPE, USAGE_FULL, 0, 0>, SHAPES...>
	{
	};

	template <int FIRST_CONT, int SRC, int POS, typename VD, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, POS, VD, segment<V_TYPE, USAGE_SLAVE, USE_ALSO, PARENT>, SHAPES...>:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, POS+1, VD, SHAPES...>
	{
	};

	template <int FIRST_CONT, int SRC, int POS, typename VD>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, POS, VD>
	{
		typedef VD type;
	};


	template <typename VD, typename M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename M7>
	struct sep_valence_by_mask;

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename ... M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename M7>
	struct sep_valence_by_mask<type_pack<valence_data<V_TYPE, 1, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, type_pack<M1...>, M2, M3, M4, M5, M6, M7>:
		public sep_valence_by_mask<type_pack<VD...>, type_pack<valence_data<V_TYPE, 1, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M1...>, M2, M3, M4, M5, M6, M7>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename ... M2, typename M3, typename M4, typename M5, typename M6, typename M7>
	struct sep_valence_by_mask<type_pack<valence_data<V_TYPE, 2, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, type_pack<M2...>, M3, M4, M5, M6, M7>:
		public sep_valence_by_mask<type_pack<VD...>, M1, type_pack<valence_data<V_TYPE, 2, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M2...>, M3, M4, M5, M6, M7>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename M2, typename ... M3, typename M4, typename M5, typename M6, typename M7>
	struct sep_valence_by_mask<type_pack<valence_data<V_TYPE, 3, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, M2, type_pack<M3...>, M4, M5, M6, M7>:
		public sep_valence_by_mask<type_pack<VD...>, M1, M2, type_pack<valence_data<V_TYPE, 3, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M3...>, M4, M5, M6, M7>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename M2, typename M3, typename ... M4, typename M5, typename M6, typename M7>
	struct sep_valence_by_mask<type_pack<valence_data<V_TYPE, 4, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, M2, M3, type_pack<M4...>, M5, M6, M7>:
		public sep_valence_by_mask<type_pack<VD...>, M1, M2, M3, type_pack<valence_data<V_TYPE, 4, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M4...>, M5, M6, M7>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename M2, typename M3, typename M4, typename ... M5, typename M6, typename M7>
	struct sep_valence_by_mask<type_pack<valence_data<V_TYPE, 5, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, M2, M3, M4, type_pack<M5...>, M6, M7>:
		public sep_valence_by_mask<type_pack<VD...>, M1, M2, M3, M4, type_pack<valence_data<V_TYPE, 5, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M5...>, M6, M7>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename M2, typename M3, typename M4, typename M5, typename ... M6, typename M7>
	struct sep_valence_by_mask<type_pack<valence_data<V_TYPE, 6, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, M2, M3, M4, M5, type_pack<M6...>, M7>:
		public sep_valence_by_mask<type_pack<VD...>, M1, M2, M3, M4, M5, type_pack<valence_data<V_TYPE, 6, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M6...>, M7>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename ... M7>
	struct sep_valence_by_mask<type_pack<valence_data<V_TYPE, 7, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, M2, M3, M4, M5, M6, type_pack<M7...> >:
		public sep_valence_by_mask<type_pack<VD...>, M1, M2, M3, M4, M5, M6, type_pack<valence_data<V_TYPE, 7, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M7...> >
	{
	};

	template <typename M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename M7>
	struct sep_valence_by_mask<type_pack<>, M1, M2, M3, M4, M5, M6, M7>
	{
		typedef type_pack<type_pack<>, M1, M2, M3, M4, M5, M6, M7> vd_by_mask;
	};


	template <int SRC, typename VD, typename ... ST>
	struct valence_parser_base;

	template <int SRC, typename VD, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, typename ... ST>
	struct valence_parser_base<SRC, VD, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, ST...>:
		public valence_parser_base<SRC+1, typename add_shapes_to_valence_data<sizeof...(SHAPES)-CONT, SRC, 0, VD, SHAPES...>::type, ST...>
	{
	};

	template <int SRC, typename VD>
	struct valence_parser_base<SRC, VD>:
		public sep_valence_by_mask<VD, type_pack<>, type_pack<>, type_pack<>, type_pack<>, type_pack<>, type_pack<>, type_pack<> >
	{
	};

	template <typename ... ST>
	struct valence_parser: public valence_parser_base<0, type_pack<>, ST...>
	{
	};

	template <typename IO, typename VO, typename SHAPES>
	struct get_index_order_base;

	template <int ... IO, int ... VO, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct get_index_order_base<std::integer_sequence<int, IO...>, std::integer_sequence<int, VO...>, type_pack<sClass<V_TYPE, USAGE, USE_ALSO,PARENT>, SHAPES...> >:
		public get_index_order_base<std::integer_sequence<int, IO..., find<int, V_TYPE>(std::integer_sequence<int, VO...>())>, std::integer_sequence<int, VO...>, type_pack<SHAPES...> >
	{
		static_assert(find<int, V_TYPE>(std::integer_sequence<int, VO...>())<sizeof...(VO),"Tensor valence is not found in the list of valences");
	};

	template <int ... IO, int ... VO, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct get_index_order_base<std::integer_sequence<int, IO...>, std::integer_sequence<int, VO...>, type_pack<sClass<V_TYPE, USAGE_SLAVE, USE_ALSO,PARENT>, SHAPES...> >:
		public get_index_order_base<std::integer_sequence<int, IO...>, std::integer_sequence<int, VO...>, type_pack<SHAPES...> >
	{
	};

	template <int ... IO, typename VO>
	struct get_index_order_base<std::integer_sequence<int, IO...>, VO, type_pack<> >
	{
		typedef std::integer_sequence<int, IO...> type;
	};

	template <typename VO, typename SHAPES>
	struct get_index_order: public get_index_order_base<std::integer_sequence<int>, VO, SHAPES>
	{
	};

	template <typename OT, typename INDEX_ORDER, typename INDT>
	struct get_ordered_tuple_base;

	template <typename ... OT, int FIRST, int ... INDEX_ORDER, typename INDT>
	struct get_ordered_tuple_base<std::tuple<OT...>, std::integer_sequence<int, FIRST, INDEX_ORDER...>, INDT>
	{
		typedef typename std::tuple_element<FIRST,INDT>::type element_type;
		typedef get_ordered_tuple_base<std::tuple<OT..., element_type>, std::integer_sequence<int, INDEX_ORDER...>, INDT> base;
		typedef typename base::type type;
		template <typename ... Ts>
		static type create_tuple(const INDT& src, Ts ... args) { return base::create_tuple(src, args..., std::get<FIRST>(src)); }
		template <typename ST>
		using parenthesis=typename base::template parenthesis<ST>;
	};

	template <typename ... OT, typename INDT>
	struct get_ordered_tuple_base<std::tuple<OT...>, std::integer_sequence<int>, INDT>
	{
		typedef std::tuple<OT...> type;
		template <typename ... Ts>
		static type create_tuple(const INDT& src, OT... args) { return std::tuple<OT...>(args...); }
		template <typename ST>
		using parenthesis=typename parenthesis<ST, OT...>::type;
	};

	template <typename INDEX_ORDER, typename INDT>
	struct get_ordered_tuple: public get_ordered_tuple_base<std::tuple<>, INDEX_ORDER, INDT>
	{
	};

};

#endif /* VALENCE_H_ */
