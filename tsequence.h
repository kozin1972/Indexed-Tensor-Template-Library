/*
 * sequence.h
 *
 *  Created on: 30 мая 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 *
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <cstddef>
#include <type_traits>
#include <utility>


namespace iTTL
{

#if __cplusplus <= 201103L
	template <typename T, T ... Is>
	struct integer_sequence {};
#else
	template <typename T, T ... Is>
	using integer_sequence = std::integer_sequence<T, Is...>;
#endif
	template <size_t ... Is>
	using index_sequence=integer_sequence<size_t, Is...>;

	template <class T>
	struct just_type { typedef T type; };
//
//	template <class... Ts>
//	struct type_pack { typedef type_pack type; };
//
//	typedef type_pack<> empty_pack;
//
//	template <class... Ts>
//	constexpr size_t size(type_pack<Ts...> ) { return sizeof...(Ts); }
//
//	template <class... Ts>
//	constexpr bool is_empty(type_pack<Ts...> tp) { return sizeof...(Ts) == 0; }
//
//	template <class T, class... Ts>
//	constexpr just_type<T> head(type_pack<T, Ts...>) { return {}; }
//
//	template <class T, class... Ts>
//	constexpr type_pack<Ts...> tail(type_pack<T, Ts...>) { return {}; }
//
//	template <class T, class... Ts>
//	constexpr size_t find(type_pack<Ts...> tp)
//	{
//		bool bs[] = {std::is_same<T, Ts>::value...};
//		for (size_t i = 0; i < sizeof...(Ts); ++i)
//		{
//			if (bs[i])
//				return i;
//		}
//		return sizeof...(Ts);
//	}
//
//	template <class T, T value, T... Is>
//	constexpr size_t find(std::integer_sequence<T, Is...> is)
//	{
//		bool bs[] = {value==Is...};
//		for (size_t i = 0; i < sizeof...(Is); ++i)
//		{
//			if (bs[i])
//				return i;
//		}
//		return sizeof...(Is);
//	}
//
//	template <class T, T... Is>
//	constexpr bool is_unique(std::integer_sequence<T, Is...> is)
//	{
//		size_t pos[]={find<T, Is>(is)...};
//		for (size_t i = 0; i < sizeof...(Is); ++i)
//		{
//			if (pos[i]<i)
//				return false;
//		}
//		return true;
//	}
//
//	template<class T, size_t POS, T... Ints>
//	constexpr T get(std::integer_sequence<T, Ints...>) {
//	    constexpr T arr[] = {Ints...};
//	    return arr[POS];
//	}
//
////	template<size_t POS, class... Ts>
////	constexpr auto get(type_pack<Ts...>) {
////		constexpr std::tuple<just_type<Ts>...> tpl;
////		return std::get<POS>(tpl);
////	}
//
//	template <size_t I, class T>
//	struct indexed_type
//	{
//		static constexpr size_t value = I;
//		using type = T;
//	};
//
//	template <class IS, class... Ts>
//	struct indexed_types;
//
//	template <size_t... Is, class... Ts>
//	struct indexed_types<std::index_sequence<Is...>, Ts...>
//	{
//		struct type : indexed_type<Is, Ts>... {};
//	};
//
//	template <class... Ts>
//	using indexed_types_for = typename indexed_types<std::index_sequence_for<Ts...>, Ts...>::type;
//
//	template <size_t I, class T>
//	constexpr just_type<T> get_indexed_type(indexed_type<I, T>) { return {}; }
//
//	template <size_t I, class... Ts>
//	constexpr auto get(type_pack<Ts...>)
//	{
//		return get_indexed_type<I>(indexed_types_for<Ts...>{});
//	}
//
//	template <size_t I, class... Ts>
//	constexpr auto get(just_type<type_pack<Ts...> >)
//	{
//		return get_indexed_type<I>(indexed_types_for<Ts...>{});
//	}
//
//	template <size_t I, class TP>
//	struct type_pack_element_getter;
//
//	template <size_t I, class ... Ts>
//	struct type_pack_element_getter<I, type_pack<Ts...> >
//	{
//		using type = typename decltype(get<I>(type_pack<Ts...>()))::type;
//	};
//
//	template <size_t POS, typename TP>
//	using type_pack_element = typename type_pack_element_getter<POS, TP>::type;

	template <typename ... Ts>
	struct type_sequence
	{
		static const size_t size=sizeof...(Ts);
		using type=type_sequence<Ts...>;
	};

	template <typename T0, typename ... Ts>
	struct type_sequence<T0, Ts...>
	{
		using type=type_sequence<T0, Ts...>;
		using head=T0;
		using tail=type_sequence<Ts...>;
		static const size_t size=sizeof...(Ts)+1;
	};

	template <size_t I, typename ... Ts>
	struct tseq_getter;

	template <typename T0, typename ... Ts>
	struct tseq_getter<0, type_sequence<T0, Ts...> >
	{
		using type=T0;
	};
	template <size_t I, typename T0, typename ... Ts>
	struct tseq_getter<I, type_sequence<T0, Ts...> >
	{
		using type=typename tseq_getter<I-1, type_sequence<Ts...> >::type;
	};

	template <size_t I, typename TS>
	using tseq_element = typename tseq_getter<I, TS>::type;

//	template <class T, T toFind, T... Is>
//	struct find_helper;
//
//	template <class T, T toFind>
//	struct find_helper<T, toFind>
//	{
//		static const size_t value=0;
//	};
//
//	template <class T, T toFind, T I0, T... Is>
//	struct find_helper<T, toFind, I0, Is...>
//	{
//		static const size_t value=((toFind==I0)?0:(find_helper<T, toFind, Is...>::value+1));
//	};
//
	template <class T, T value>
	constexpr size_t find(integer_sequence<T> is) noexcept { return 0; }

	template <class T, T value, T I0, T... Is>
	constexpr size_t find(integer_sequence<T, I0, Is...> is) noexcept
	{
//		return find_helper<T, value, I0, Is...>::value;
		return value==I0?0:(find<T,value>(integer_sequence<T, Is...>{})+1);
	}
//
//	template <class T, T value>
//	constexpr size_t find() { return 0; }
//
//	template <class T, T value, T I0, T... Is>
//	constexpr size_t find()
//	{
//		return value==I0?0:(find<T,value,Is...>()+1);
//	}
//
//	template <class T, T value>
//	constexpr bool in_sequence() { return false; }
//	template <class T, T value, T I0, T... Is>
//	constexpr bool in_sequence()
//	{
//		return value==I0?true:in_sequence<T,value,Is...>();
//	}

	template <typename ... Ts>
	constexpr size_t size(type_sequence<Ts...>) noexcept { return sizeof...(Ts); }

	template <typename T, T NEW_ELEMENT, T... Is>
	constexpr integer_sequence<T, Is..., NEW_ELEMENT> push_back(integer_sequence<T, Is...>) noexcept { return {}; }

	template <typename T, T NEW_ELEMENT, T... Is>
	constexpr integer_sequence<T, NEW_ELEMENT, Is...> push_front(integer_sequence<T, Is...>) noexcept { return {}; }

#if __cplusplus <= 201103L
	template<size_t POS, class T, T I0, T... Ints, typename std::enable_if<POS==0>::type...>
	constexpr T get(integer_sequence<T, I0, Ints...>)
	{
	    return I0;
	}
	template<size_t POS, class T, T I0, T... Ints, typename std::enable_if<POS!=0>::type...>
	constexpr T get(integer_sequence<T, I0, Ints...>)
	{
	    return get<POS-1>(integer_sequence<T, Ints...>{});
	}
	template <size_t I>
	constexpr bool value_is_not_less_than_index(index_sequence<>)
	{
		return true;
	}
	template <size_t I, size_t V0, size_t ... V>
	constexpr bool value_is_not_less_than_index(index_sequence<V0, V...>)
	{
		return V0<I?false:value_is_not_less_than_index<I+1>(index_sequence<V...>{});
	}
	template <class T, T... Is>
	constexpr bool is_unique(integer_sequence<T, Is...> is)
	{
		return value_is_not_less_than_index<0>(index_sequence<find<T, Is>(is)...>{});
	}
#else
	template<size_t POS, class T, T... Ints>
	constexpr T get(std::integer_sequence<T, Ints...>) noexcept
	{
	    constexpr T arr[] = {Ints...};
	    return arr[POS];
	}
	template <class T, T... Is>
	constexpr bool is_unique(std::integer_sequence<T, Is...> is) noexcept
	{
		size_t pos[]={find<T, Is>(is)...};
		for (size_t i = 0; i < sizeof...(Is); ++i)
		{
			if (pos[i]<i)
				return false;
		}
		return true;
	}
#endif

}


#endif /* SEQUENCE_H_ */
