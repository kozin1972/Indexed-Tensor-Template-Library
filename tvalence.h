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
	public:
		static size_t length(const ST& st) { return get<POS>(st).length(); }
		static size_t step(const ST& st) { return get<POS>(st).step(); }
	};

	template <typename T, typename ST, size_t ... POS>
	class joined_iterator_cont;
	template <typename T, typename ST, size_t POS0, size_t ... POS>
	class joined_iterator_cont<T, ST, POS0, POS...>
	{
		const T *end;
	public:
		static size_t length(const ST& st) { return joined_iterator_cont<T, ST, POS...>::length(st)*get<POS0>(st).length(); }
		static constexpr size_t step(const ST& st) { return 1; }
		joined_iterator_cont(const ST& st, const T * const p): end(p+length(st)) {}
		void move_one(T *&p) { p++; }
		void move_one(const T *&p) { p++; }
		bool not_end(const T *p) { return p<end; }
	};
	template <typename T, typename ST, size_t POS>
	class joined_iterator_cont<T, ST, POS>
	{
	public:
		static size_t length(const ST& st) { return get<POS>(st).length(); }
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


	template <typename CVDV, typename CVD, typename VD, int V_TYPE, bool VECTOR_FIRST, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS>
	struct add_valence_data;

	template <typename ... CVDV, typename ... CVD, int H_V_TYPE, int H_MASK, int H_C_MASK, int H_V_MASK, int H_NEXT_VALENCE, size_t ... H_POS, typename ... VD, int V_TYPE, bool VECTOR_FIRST, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS>
	struct add_valence_data<type_sequence<CVDV...>, type_sequence<CVD...>, type_sequence<valence_data<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_NEXT_VALENCE, H_POS...>, VD...>, V_TYPE, VECTOR_FIRST, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...>:
		public
				std::conditional<(V_TYPE==H_V_TYPE),
					typename std::conditional<((H_MASK==H_V_MASK) == ((MASK | H_MASK)==(V_MASK | H_V_MASK))),
						type_sequence<CVDV..., CVD..., valence_data<V_TYPE, MASK | H_MASK, C_MASK | H_C_MASK, V_MASK | H_V_MASK, (NEXT_VALENCE==H_NEXT_VALENCE?NEXT_VALENCE:V_TYPE), H_POS..., POS...>, VD...>,
						type_sequence<CVDV..., CVD..., VD..., valence_data<V_TYPE, MASK | H_MASK, C_MASK | H_C_MASK, V_MASK | H_V_MASK, (NEXT_VALENCE==H_NEXT_VALENCE?NEXT_VALENCE:V_TYPE), H_POS..., POS...> >
					>::type,
					typename std::conditional<(H_MASK==H_V_MASK),
						add_valence_data<type_sequence<CVDV..., valence_data<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_NEXT_VALENCE, H_POS...> >, type_sequence<CVD...>, type_sequence<VD...>, V_TYPE, VECTOR_FIRST, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...>,
						add_valence_data<type_sequence<CVDV...>, type_sequence<CVD..., valence_data<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, H_NEXT_VALENCE, H_POS...> >, type_sequence<VD...>, V_TYPE, VECTOR_FIRST, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...>
					>::type
				>::type
	{
	};

	template <typename ... CVDV, typename ... CVD, int V_TYPE, bool VECTOR_FIRST, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS>
	struct add_valence_data<type_sequence<CVDV...>, type_sequence<CVD...>, type_sequence<>, V_TYPE, VECTOR_FIRST, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...>
	{
		using type=typename std::conditional<(MASK==V_MASK) && VECTOR_FIRST,
					type_sequence<CVDV..., valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...>, CVD...>,
					type_sequence<CVDV..., CVD..., valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...> >
				>::type;
	};

	template <int FIRST_CONT, int SRC, bool VECTOR_FIRST, size_t POS, typename VD, typename ... SHAPES>
	struct add_shapes_to_valence_data;

	template <int FIRST_CONT, int SRC, bool VECTOR_FIRST, size_t POS, typename VD, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS, VD, sClass<V_TYPE, USAGE, USE_ALSO, PARENT>, SHAPES...>:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS+1, typename add_valence_data<type_sequence<>, type_sequence<>, VD, V_TYPE, VECTOR_FIRST, (1<<SRC), 0, 0, V_TYPE, POS>::type, SHAPES...>
	{
	};

	template <int FIRST_CONT, int SRC, bool VECTOR_FIRST, size_t POS, typename VD, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS, VD, sClass<V_TYPE, USAGE_SLAVE, USE_ALSO, PARENT>, SHAPES...>:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS+1, VD, SHAPES...>
	{
	};

	template <int FIRST_CONT, int SRC, bool VECTOR_FIRST, size_t POS, typename VD, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS, VD, segment<V_TYPE, USAGE, USE_ALSO, PARENT> >:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS+1, typename add_valence_data<type_sequence<>, type_sequence<>, VD, V_TYPE, VECTOR_FIRST, (1<<SRC), ((USAGE<=USAGE_CONT && POS>=FIRST_CONT)?(1<<SRC):0), ((USE_ALSO==0 && PARENT==0)?(1<<SRC):0), V_TYPE, POS>::type>
	{
	};

	template <int FIRST_CONT, int SRC, bool VECTOR_FIRST, size_t POS, typename VD, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, template <int, enum dSU, size_t, size_t> class N_sClass, int N_V_TYPE, enum dSU N_USAGE, size_t N_USE_ALSO, size_t N_PARENT, typename ... SHAPES>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS, VD, segment<V_TYPE, USAGE, USE_ALSO, PARENT>, N_sClass<N_V_TYPE, N_USAGE, N_USE_ALSO, N_PARENT>, SHAPES...>:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS+1, typename add_valence_data<type_sequence<>, type_sequence<>, VD, V_TYPE, VECTOR_FIRST, (1<<SRC), 0, ((USE_ALSO==0 && PARENT==0)?(1<<SRC):0), V_TYPE, POS>::type, N_sClass<N_V_TYPE, N_USAGE, N_USE_ALSO, N_PARENT>, SHAPES...>
	{
	};

	template <int FIRST_CONT, int SRC, bool VECTOR_FIRST, size_t POS, typename VD, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, int N_V_TYPE, typename ... SHAPES>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS, VD, segment<V_TYPE, USAGE, USE_ALSO, PARENT>, segment<N_V_TYPE, USAGE_FULL, 0, 0>, SHAPES...>:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS+1, typename add_valence_data<type_sequence<>, type_sequence<>, VD, V_TYPE, VECTOR_FIRST, (1<<SRC), 0, ((USE_ALSO==0 && PARENT==0)?(1<<SRC):0), (USE_ALSO==0 && PARENT==0 && USAGE<=USAGE_CONT?N_V_TYPE:V_TYPE), POS>::type, segment<N_V_TYPE, USAGE_FULL, 0, 0>, SHAPES...>
	{
	};

	template <int FIRST_CONT, int SRC, bool VECTOR_FIRST, size_t POS, typename VD, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS, VD, segment<V_TYPE, USAGE_SLAVE, USE_ALSO, PARENT>, SHAPES...>:
		public add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS+1, VD, SHAPES...>
	{
	};

	template <int FIRST_CONT, int SRC, bool VECTOR_FIRST, size_t POS, typename VD>
	struct add_shapes_to_valence_data<FIRST_CONT, SRC, VECTOR_FIRST, POS, VD>
	{
		typedef VD type;
	};


	template <typename VD, typename M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename M7, bool CONT_0, bool CONT_1, bool CONT_2>
	struct sep_valence_by_mask;

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename ... M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename M7, bool CONT_0, bool CONT_1, bool CONT_2>
	struct sep_valence_by_mask<type_sequence<valence_data<V_TYPE, 1, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, type_sequence<M1...>, M2, M3, M4, M5, M6, M7, CONT_0, CONT_1, CONT_2>:
		public sep_valence_by_mask<type_sequence<VD...>, type_sequence<valence_data<V_TYPE, 1, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M1...>, M2, M3, M4, M5, M6, M7, CONT_0, CONT_1, CONT_2>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename ... M2, typename M3, typename M4, typename M5, typename M6, typename M7, bool CONT_0, bool CONT_1, bool CONT_2>
	struct sep_valence_by_mask<type_sequence<valence_data<V_TYPE, 2, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, type_sequence<M2...>, M3, M4, M5, M6, M7, CONT_0, CONT_1, CONT_2>:
		public sep_valence_by_mask<type_sequence<VD...>, M1, type_sequence<valence_data<V_TYPE, 2, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M2...>, M3, M4, M5, M6, M7, CONT_0, CONT_1, CONT_2>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename M2, typename ... M3, typename M4, typename M5, typename M6, typename M7, bool CONT_0, bool CONT_1, bool CONT_2>
	struct sep_valence_by_mask<type_sequence<valence_data<V_TYPE, 3, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, M2, type_sequence<M3...>, M4, M5, M6, M7, CONT_0, CONT_1, CONT_2>:
		public sep_valence_by_mask<type_sequence<VD...>, M1, M2, type_sequence<valence_data<V_TYPE, 3, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M3...>, M4, M5, M6, M7, CONT_0 || ((V_MASK==3 && (C_MASK & 1))), CONT_1 || ((V_MASK==3 && (C_MASK & 2))), CONT_2>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename M2, typename M3, typename ... M4, typename M5, typename M6, typename M7, bool CONT_0, bool CONT_1, bool CONT_2>
	struct sep_valence_by_mask<type_sequence<valence_data<V_TYPE, 4, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, M2, M3, type_sequence<M4...>, M5, M6, M7, CONT_0, CONT_1, CONT_2>:
		public sep_valence_by_mask<type_sequence<VD...>, M1, M2, M3, type_sequence<valence_data<V_TYPE, 4, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M4...>, M5, M6, M7, CONT_0, CONT_1, CONT_2>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename M2, typename M3, typename M4, typename ... M5, typename M6, typename M7, bool CONT_0, bool CONT_1, bool CONT_2>
	struct sep_valence_by_mask<type_sequence<valence_data<V_TYPE, 5, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, M2, M3, M4, type_sequence<M5...>, M6, M7, CONT_0, CONT_1, CONT_2>:
		public sep_valence_by_mask<type_sequence<VD...>, M1, M2, M3, M4, type_sequence<valence_data<V_TYPE, 5, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M5...>, M6, M7, CONT_0 || ((V_MASK==5 && (C_MASK & 1))), CONT_1, CONT_2 || ((V_MASK==5 && (C_MASK & 4)))>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename M2, typename M3, typename M4, typename M5, typename ... M6, typename M7, bool CONT_0, bool CONT_1, bool CONT_2>
	struct sep_valence_by_mask<type_sequence<valence_data<V_TYPE, 6, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, M2, M3, M4, M5, type_sequence<M6...>, M7, CONT_0, CONT_1, CONT_2>:
		public sep_valence_by_mask<type_sequence<VD...>, M1, M2, M3, M4, M5, type_sequence<valence_data<V_TYPE, 6, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M6...>, M7, CONT_0, CONT_1 || ((V_MASK==6 && (C_MASK & 2))), CONT_2 || ((V_MASK==6 && (C_MASK & 4)))>
	{
	};

	template <int V_TYPE, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ... POS, typename ... VD, typename M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename ... M7, bool CONT_0, bool CONT_1, bool CONT_2>
	struct sep_valence_by_mask<type_sequence<valence_data<V_TYPE, 7, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, M1, M2, M3, M4, M5, M6, type_sequence<M7...>, CONT_0, CONT_1, CONT_2>:
		public sep_valence_by_mask<type_sequence<VD...>, M1, M2, M3, M4, M5, M6, type_sequence<valence_data<V_TYPE, 7, C_MASK, V_MASK, NEXT_VALENCE, POS...>, M7...>, CONT_0, CONT_1, CONT_2>
	{
	};

	template <typename M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename M7, bool CONT_0, bool CONT_1, bool CONT_2>
	struct sep_valence_by_mask<type_sequence<>, M1, M2, M3, M4, M5, M6, M7, CONT_0, CONT_1, CONT_2>
	{
		typedef type_sequence<std::integral_constant<int, (CONT_0?1:0)+2*(CONT_1?1:0)+4*(CONT_2?1:0)>, M1, M2, M3, M4, M5, M6, M7> type;
	};


	template <int SRC, bool VECTOR_FIRST, typename VD, typename ... ST>
	struct valence_parser_base;

	template <int SRC, bool VECTOR_FIRST, typename VD, size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, typename ... ST>
	struct valence_parser_base<SRC, VECTOR_FIRST, VD, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, ST...>:
		public valence_parser_base<SRC+1, VECTOR_FIRST, typename add_shapes_to_valence_data<sizeof...(SHAPES)-CONT, SRC, VECTOR_FIRST, 0, VD, SHAPES...>::type, ST...>
	{
	};

	template <int SRC, bool VECTOR_FIRST, typename VD>
	struct valence_parser_base<SRC, VECTOR_FIRST, VD>:
		public sep_valence_by_mask<VD, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, type_sequence<>, false, false, false>
	{
	};

	template <bool VECTOR_FIRST, typename ... ST>
	struct valence_parser: public valence_parser_base<0, VECTOR_FIRST, type_sequence<>, ST...>
	{
	};

	template <int SRC0, int SRC1, size_t POS0, size_t POS1, size_t SNUM_0, size_t CONT_0, bool IS_INDEXED_0, typename OST_0, typename ... SHAPES_0, size_t SNUM_1, size_t CONT_1, bool IS_INDEXED_1, typename OST_1, typename ... SHAPES_1>
	void report_incompatible_shapes(const stuple<SNUM_0, CONT_0, IS_INDEXED_0, OST_0, SHAPES_0...>& st0, const stuple<SNUM_1, CONT_1, IS_INDEXED_1, OST_1, SHAPES_1...>& st1)
	{
//		static const size_t inum0=tpp::index_number(std::integral_constant<size_t, POS0>{}, type_sequence<SHAPES_0...>{});
//		static const size_t inum1=tpp::index_number(std::integral_constant<size_t, POS1>{}, type_sequence<SHAPES_1...>{});
		static const size_t inum0=tpp::index_number<POS0, SHAPES_0...>::inum;
		static const size_t inum1=tpp::index_number<POS1, SHAPES_1...>::inum;
		static const int v_type=check_shape<tpp::tseq_element<POS0, type_sequence<SHAPES_0...> > >::v_type;
		throw incompatibleIndices(inum0, SRC0, inum1, SRC1, v_type, tpp::get<POS0>(st0).length(), tpp::get<POS1>(st1).length());
	}

	template <int SRC0, int SRC1, size_t SNUM_0, size_t CONT_0, bool IS_INDEXED_0, typename OST_0, typename ... SHAPES_0, size_t SNUM_1, size_t CONT_1, bool IS_INDEXED_1, typename OST_1, typename ... SHAPES_1>
	void check_shape_length_01(type_sequence<>, const stuple<SNUM_0, CONT_0, IS_INDEXED_0, OST_0, SHAPES_0...>& st0, const stuple<SNUM_1, CONT_1, IS_INDEXED_1, OST_1, SHAPES_1...>& st1)
	{
	}

	template <int SRC0, int SRC1, typename VDH, typename ... VDT, size_t SNUM_0, size_t CONT_0, bool IS_INDEXED_0, typename OST_0, typename ... SHAPES_0, size_t SNUM_1, size_t CONT_1, bool IS_INDEXED_1, typename OST_1, typename ... SHAPES_1>
	void check_shape_length_01(type_sequence<VDH, VDT...>, const stuple<SNUM_0, CONT_0, IS_INDEXED_0, OST_0, SHAPES_0...>& st0, const stuple<SNUM_1, CONT_1, IS_INDEXED_1, OST_1, SHAPES_1...>& st1)
	{
		if (tpp::get<VDH::pos0>(st0).length()!=tpp::get<VDH::pos1>(st1).length())
			report_incompatible_shapes<SRC0, SRC1, VDH::pos0, VDH::pos1>(st0, st1);
		check_shape_length_01<SRC0, SRC1>(type_sequence<VDT...>{}, st0, st1);
	}

	template <size_t SNUM_0, size_t CONT_0, bool IS_INDEXED_0, typename OST_0, typename ... SHAPES_0, size_t SNUM_1, size_t CONT_1, bool IS_INDEXED_1, typename OST_1, typename ... SHAPES_1>
	void check_shape_length_12(type_sequence<>, const stuple<SNUM_0, CONT_0, IS_INDEXED_0, OST_0, SHAPES_0...>& st0, const stuple<SNUM_1, CONT_1, IS_INDEXED_1, OST_1, SHAPES_1...>& st1)
	{
	}

	template <typename VDH, typename ... VDT, size_t SNUM_0, size_t CONT_0, bool IS_INDEXED_0, typename OST_0, typename ... SHAPES_0, size_t SNUM_1, size_t CONT_1, bool IS_INDEXED_1, typename OST_1, typename ... SHAPES_1>
	void check_shape_length_12(type_sequence<VDH, VDT...>, const stuple<SNUM_0, CONT_0, IS_INDEXED_0, OST_0, SHAPES_0...>& st0, const stuple<SNUM_1, CONT_1, IS_INDEXED_1, OST_1, SHAPES_1...>& st1)
	{
		if (tpp::get<VDH::pos1>(st0).length()!=tpp::get<VDH::pos2>(st1).length())
			report_incompatible_shapes<1, 2, VDH::pos1, VDH::pos2>(st0, st1);
		check_shape_length_12(type_sequence<VDT...>{}, st0, st1);
	}

	template <typename VD, size_t SNUM_0, size_t CONT_0, bool IS_INDEXED_0, typename OST_0, typename ... SHAPES_0, size_t SNUM_1, size_t CONT_1, bool IS_INDEXED_1, typename OST_1, typename ... SHAPES_1>
	void check_shape_length(const stuple<SNUM_0, CONT_0, IS_INDEXED_0, OST_0, SHAPES_0...>& st0, const stuple<SNUM_1, CONT_1, IS_INDEXED_1, OST_1, SHAPES_1...>& st1)
	{
		check_shape_length_01<0, 1>(tseq_element<3, VD>{}, st0, st1);
	}

	template <typename VD, size_t SNUM_0, size_t CONT_0, bool IS_INDEXED_0, typename OST_0, typename ... SHAPES_0, size_t SNUM_1, size_t CONT_1, bool IS_INDEXED_1, typename OST_1, typename ... SHAPES_1, size_t SNUM_2, size_t CONT_2, bool IS_INDEXED_2, typename OST_2, typename ... SHAPES_2>
	void check_shape_length(const stuple<SNUM_0, CONT_0, IS_INDEXED_0, OST_0, SHAPES_0...>& st0, const stuple<SNUM_1, CONT_1, IS_INDEXED_1, OST_1, SHAPES_1...>& st1, const stuple<SNUM_2, CONT_2, IS_INDEXED_2, OST_2, SHAPES_2...>& st2)
	{
		check_shape_length_01<0, 1>(tseq_element<3, VD>{}, st0, st1);
		check_shape_length_01<0, 2>(tseq_element<5, VD>{}, st0, st2);
		check_shape_length_01<1, 2>(tseq_element<6, VD>{}, st1, st2);
		check_shape_length_01<0, 1>(tseq_element<7, VD>{}, st0, st1);
		check_shape_length_12(tseq_element<7, VD>{}, st1, st2);
	}

	template <typename JVD, typename VD, typename ... CVD>
	struct join_one_vd;

	template <typename ... JVD, typename LAST, typename ... CVD>
	struct join_one_vd<type_sequence<JVD...>, type_sequence<>, LAST, CVD...>
	{
		typedef type_sequence<JVD..., type_sequence<CVD..., LAST> > type;
	};

	template <typename ... JVD>
	struct join_one_vd<type_sequence<JVD...>, type_sequence<>>
	{
		typedef type_sequence<JVD...> type;
	};

	template <typename ... JVD, typename HEAD, typename ... VD>
	struct join_one_vd<type_sequence<JVD...>, type_sequence<HEAD, VD...> >:
		public join_one_vd<type_sequence<JVD...>, type_sequence<VD...>, HEAD>
	{
	};

	template <typename ... JVD, int H_V_TYPE, int H_MASK, int H_C_MASK, int H_V_MASK, size_t ... H_POS, typename ... VD, int V_TYPE, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ...POS, typename ... CVD>
	struct join_one_vd<type_sequence<JVD...>, type_sequence<valence_data<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, V_TYPE, H_POS...>, VD...>, valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...>, CVD...>:
		public join_one_vd<type_sequence<JVD...>, type_sequence<VD...>, valence_data<H_V_TYPE, H_MASK, H_C_MASK, H_V_MASK, V_TYPE, H_POS...>, CVD..., valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...> >
	{
	};

	template <typename ... JVD, typename HEAD, typename ... VD, typename LAST, typename ... CVD>
	struct join_one_vd<type_sequence<JVD...>, type_sequence<HEAD, VD...>, LAST, CVD...>:
		public join_one_vd<type_sequence<JVD..., type_sequence<CVD...,LAST> >, type_sequence<VD...>, HEAD>
	{
	};

	template <typename VD8>
	struct join_valence_data;

	template <typename M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename M7, int CONT_MASK>
	struct join_valence_data<type_sequence<std::integral_constant<int, CONT_MASK>, M1, M2, M3, M4, M5, M6, M7> >
	{
	    typedef type_sequence<std::integral_constant<int, CONT_MASK>,
				typename join_one_vd<type_sequence<>, M1>::type,
				typename join_one_vd<type_sequence<>, M2>::type,
				typename join_one_vd<type_sequence<>, M3>::type,
				typename join_one_vd<type_sequence<>, M4>::type,
				typename join_one_vd<type_sequence<>, M5>::type,
				typename join_one_vd<type_sequence<>, M6>::type,
				typename join_one_vd<type_sequence<>, M7>::type
				> type;
	};

	template <bool VECTOR_FIRST, typename ... ST>
	struct valence_parser_join: public join_valence_data<typename valence_parser_base<0, VECTOR_FIRST, type_sequence<>, ST...>::type>
	{
	};


	template <int MASK, int SRC, size_t... POS>
	struct vd_pos;

	template <int MASK, int SRC, size_t POS0, size_t... POS>
	struct vd_pos<MASK, SRC, POS0, POS...>:
	    public std::conditional<(MASK & 1),
			vd_pos<(MASK>>1), SRC-1, POS...>,
			vd_pos<(MASK>>1), SRC-1, POS0, POS...>
			>::type
	{
	};

	template <int MASK, size_t POS0, size_t... POS>
	struct vd_pos<MASK, 0, POS0, POS...>
	{
		static const size_t pos=POS0;
	};

	template <typename ST, bool IS_CONT, typename VD, int SRC, size_t ... SPOS>
	struct vd_iterator_getter_mult;

	template <typename ST, bool IS_CONT, int V_TYPE, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ...POS, typename ... NVD, int SRC, size_t ... SPOS>
	struct vd_iterator_getter_mult<ST, IS_CONT, type_sequence<valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...>, NVD...>, SRC, SPOS...>:
		public vd_iterator_getter_mult<ST, IS_CONT, type_sequence<NVD...>, SRC, vd_pos<MASK, SRC, POS...>::pos, SPOS...>
	{
	};

	template <typename ST, int SRC, size_t ... SPOS>
	struct vd_iterator_getter_mult<ST, false, type_sequence<>, SRC, SPOS...>
	{
		template <typename T>
		using iterator_type=joined_iterator<T, ST, SPOS...>;
	};
	template <typename ST, int SRC, size_t ... SPOS>
	struct vd_iterator_getter_mult<ST, true, type_sequence<>, SRC, SPOS...>
	{
		template <typename T>
		using iterator_type=joined_iterator_cont<T, ST, SPOS...>;
	};

	template <typename ST, typename VD, int SRC>
	struct vd_iterator_getter;

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, int V_TYPE, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ...POS, int SRC>
	struct vd_iterator_getter<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...>, SRC>
	{
		static const size_t pos=vd_pos<MASK, SRC, POS...>::pos;
		template <typename T>
		using iterator_type=shapes_iterator_base<T, pos, (CONT>0 && pos==sizeof...(SHAPES)-1), typename std::tuple_element<pos, std::tuple<SHAPES...> >::type, check_shape<typename std::tuple_element<pos, std::tuple<SHAPES...> >::type>::need_parent, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >;
	};

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, int V_TYPE, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ...POS, typename ... VD, int SRC>
	struct vd_iterator_getter<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, type_sequence<valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, SRC>:
		public vd_iterator_getter_mult<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, ((C_MASK & (1<<SRC))!=0), type_sequence<valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...>, VD...>, SRC>
	{
	};

	template <size_t SNUM, size_t CONT, bool IS_INDEXED, typename OST, typename ... SHAPES, int MASK, int C_MASK, int V_MASK, int NEXT_VALENCE, size_t ...POS, int SRC, int V_TYPE>
	struct vd_iterator_getter<stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...>, type_sequence<valence_data<V_TYPE, MASK, C_MASK, V_MASK, NEXT_VALENCE, POS...>>, SRC>
	{
		static const size_t pos=vd_pos<MASK, SRC, POS...>::pos;
		template <typename T>
		using iterator_type=shapes_iterator_base<T, pos, (CONT>0 && pos==sizeof...(SHAPES)-1), typename std::tuple_element<pos, std::tuple<SHAPES...> >::type, check_shape<typename std::tuple_element<pos, std::tuple<SHAPES...> >::type>::need_parent, stuple<SNUM, CONT, IS_INDEXED, OST, SHAPES...> >;
	};

	template <typename IO, typename VO, typename SHAPES>
	struct get_index_order_base;

	template <int ... IO, int ... VO, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct get_index_order_base<tpp::integer_sequence<int, IO...>, tpp::integer_sequence<int, VO...>, type_sequence<sClass<V_TYPE, USAGE, USE_ALSO,PARENT>, SHAPES...> >:
		public get_index_order_base<tpp::integer_sequence<int, IO..., find<int, V_TYPE>(integer_sequence<int,VO...>{})>, tpp::integer_sequence<int, VO...>, type_sequence<SHAPES...> >
	{
		static_assert(find<int, V_TYPE>(integer_sequence<int, VO...>{})<sizeof...(VO),"Tensor valence is not found in the list of valences");
	};

	template <int ... IO, int ... VO, template <int, enum dSU, size_t, size_t> class sClass, int V_TYPE, size_t USE_ALSO, size_t PARENT, typename ... SHAPES>
	struct get_index_order_base<tpp::integer_sequence<int, IO...>, tpp::integer_sequence<int, VO...>, type_sequence<sClass<V_TYPE, USAGE_SLAVE, USE_ALSO,PARENT>, SHAPES...> >:
		public get_index_order_base<tpp::integer_sequence<int, IO...>, tpp::integer_sequence<int, VO...>, type_sequence<SHAPES...> >
	{
	};

	template <int ... IO, typename VO>
	struct get_index_order_base<tpp::integer_sequence<int, IO...>, VO, type_sequence<> >
	{
		typedef tpp::integer_sequence<int, IO...> type;
	};

	template <typename VO, typename SHAPES>
	struct get_index_order: public get_index_order_base<tpp::integer_sequence<int>, VO, SHAPES>
	{
	};

	template <typename OT, typename INDEX_ORDER, typename INDT>
	struct get_ordered_tuple_base;

	template <typename ... OT, int FIRST, int ... INDEX_ORDER, typename INDT>
	struct get_ordered_tuple_base<std::tuple<OT...>, tpp::integer_sequence<int, FIRST, INDEX_ORDER...>, INDT>
	{
		typedef typename std::tuple_element<FIRST,INDT>::type element_type;
		typedef get_ordered_tuple_base<std::tuple<OT..., element_type>, tpp::integer_sequence<int, INDEX_ORDER...>, INDT> base;
		typedef typename base::type type;
		template <typename ... Ts>
		static type create_tuple(const INDT& src, Ts ... args) { return base::create_tuple(src, args..., std::get<FIRST>(src)); }
		template <typename ST>
		using parenthesis=typename base::template parenthesis<ST>;
	};

	template <typename ... OT, typename INDT>
	struct get_ordered_tuple_base<std::tuple<OT...>, tpp::integer_sequence<int>, INDT>
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
