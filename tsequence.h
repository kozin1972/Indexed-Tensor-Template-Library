/*
 * sequence.h
 *
 *  Created on: 30 мая 2019 г.
 *      Author: Alexey Kozin
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <cstddef>

namespace tpp
{
	template <typename ... SEQ>
	struct type_sequence
	{
		typedef type_sequence type;
		static const size_t size=sizeof...(SEQ);
	};

	template <typename HEAD, typename ... SEQ>
	struct type_sequence<HEAD, SEQ...>
	{
		typedef type_sequence type;
		static const size_t size=sizeof...(SEQ)+1;
		typedef HEAD head;
		typedef type_sequence<SEQ...> tail;
	};

	template <size_t POS, typename TS>
	struct tseq_element;

	template <size_t POS, typename HEAD, typename ... SEQ>
	struct tseq_element<POS, type_sequence<HEAD, SEQ...> >:
		public tseq_element<POS-1, type_sequence<SEQ...> >
	{
	};

	template <typename HEAD, typename ... SEQ>
	struct tseq_element<0, type_sequence<HEAD, SEQ...> >
	{
		using type=HEAD;
	};

	template <typename TS, typename RES>
	struct tseq_add_backward;

	template <typename HEAD, typename ... TS, typename ... RES>
	struct tseq_add_backward<type_sequence<HEAD, TS...>, type_sequence<RES...> >:
		public tseq_add_backward<type_sequence<TS...>, type_sequence<HEAD, RES...> >
	{
	};

	template <typename ... RES>
	struct tseq_add_backward<type_sequence<>, type_sequence<RES...> >
	{
		using type=type_sequence<RES...>;
	};

	template <size_t ... SEQ>
	struct size_t_sequence
	{
		typedef size_t_sequence type;
		static const size_t size=sizeof...(SEQ);
	};

	template <size_t HEAD, size_t ... SEQ>
	struct size_t_sequence<HEAD, SEQ...>
	{
		typedef size_t_sequence type;
		static const size_t size=sizeof...(SEQ)+1;
		static const size_t head=HEAD;
		typedef size_t_sequence<SEQ...> tail;
	};

	template <size_t POS, typename SS>
	struct sseq_element;

	template <size_t POS, size_t HEAD, size_t ... SEQ>
	struct sseq_element<POS, size_t_sequence<HEAD, SEQ...> >:
		public sseq_element<POS-1, size_t_sequence<SEQ...> >
	{
	};

	template <size_t HEAD, size_t ... SEQ>
	struct sseq_element<0, size_t_sequence<HEAD, SEQ...> >
	{
		static const size_t value=HEAD;
	};

	template <size_t POS, size_t ELEMENT, size_t ... SEQ>
	struct find_size_t;

	template <size_t POS, size_t ELEMENT, size_t ... SEQ>
	struct find_size_t<POS, ELEMENT, ELEMENT, SEQ...>
	{
		static const bool found=true;
		static const size_t pos=POS;
	};

	template <size_t POS, size_t ELEMENT, size_t HEAD, size_t ... SEQ>
	struct find_size_t<POS, ELEMENT, HEAD, SEQ...>: public find_size_t<POS+1, ELEMENT, SEQ...>
	{
	};

	template <size_t POS, size_t ELEMENT>
	struct find_size_t<POS, ELEMENT>
	{
		static const bool found=false;
	};

	template <size_t ... SEQ>
	struct is_unique_size_t;

	template <size_t ELEMENT, size_t ... SEQ>
	struct is_unique_size_t<ELEMENT, SEQ...>
	{
		static const bool value=find_size_t<0, ELEMENT, SEQ...>::found?false:is_unique_size_t<SEQ...>::value;
	};

	template <>
	struct is_unique_size_t<>
	{
		static const bool value=true;
	};

	template <int ... SEQ>
	struct int_sequence
	{
		typedef int_sequence type;
		static const size_t size=sizeof...(SEQ);
	};

	template <int HEAD, int ... SEQ>
	struct int_sequence<HEAD, SEQ...>
	{
		typedef int_sequence type;
		static const size_t size=sizeof...(SEQ)+1;
		static const int head=HEAD;
		typedef int_sequence<SEQ...> tail;
	};

	template <size_t POS, typename SS>
	struct iseq_element;

	template <size_t POS, int HEAD, int ... SEQ>
	struct iseq_element<POS, int_sequence<HEAD, SEQ...> >:
		public iseq_element<POS-1, int_sequence<SEQ...> >
	{
	};

	template <int HEAD, int ... SEQ>
	struct iseq_element<0, int_sequence<HEAD, SEQ...> >
	{
		static const int value=HEAD;
	};

	template <size_t POS, int ELEMENT, int ... SEQ>
	struct find_int;

	template <size_t POS, int ELEMENT, int ... SEQ>
	struct find_int<POS, ELEMENT, ELEMENT, SEQ...>
	{
		static const bool found=true;
		static const size_t pos=POS;
	};

	template <size_t POS, int ELEMENT, int HEAD, int ... SEQ>
	struct find_int<POS, ELEMENT, HEAD, SEQ...>: public find_int<POS+1, ELEMENT, SEQ...>
	{
	};

	template <size_t POS, int ELEMENT>
	struct find_int<POS, ELEMENT>
	{
		static const bool found=false;
	};

	template <int ... SEQ>
	struct is_unique_int;

	template <int ELEMENT, int ... SEQ>
	struct is_unique_int<ELEMENT, SEQ...>
	{
		static const bool value=find_int<0, ELEMENT, SEQ...>::found?false:is_unique_int<SEQ...>::value;
	};

	template <>
	struct is_unique_int<>
	{
		static const bool value=true;
	};

}

#endif /* SEQUENCE_H_ */
