/*
 * sequence.h
 *
 *  Created on: 30 мая 2019 г.
 *      Author: Alexey Kozin
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <cstddef>
#include <type_traits>
#include <utility>

namespace tpp
{

	template <class T>
	struct just_type { typedef T type; };

	template <class... Ts>
	struct type_pack { typedef type_pack type; };

	typedef type_pack<> empty_pack;

	template <class... Ts>
	constexpr size_t size(type_pack<Ts...> ) { return sizeof...(Ts); }

	template <class... Ts>
	constexpr bool is_empty(type_pack<Ts...> tp) { return sizeof...(Ts) == 0; }

	template <class T, class... Ts>
	constexpr just_type<T> head(type_pack<T, Ts...>) { return {}; }

	template <class T, class... Ts>
	constexpr type_pack<Ts...> tail(type_pack<T, Ts...>) { return {}; }

	template <class T, class... Ts>
	constexpr size_t find(type_pack<Ts...> tp)
	{
		bool bs[] = {std::is_same<T, Ts>::value...};
		for (size_t i = 0; i < sizeof...(Ts); ++i)
		{
			if (bs[i])
				return i;
		}
		return sizeof...(Ts);
	}

	template <class T, T value, T... Is>
	constexpr size_t find(std::integer_sequence<T, Is...> is)
	{
		bool bs[] = {value==Is...};
		for (size_t i = 0; i < sizeof...(Is); ++i)
		{
			if (bs[i])
				return i;
		}
		return sizeof...(Is);
	}

	template <class T, T... Is>
	constexpr bool is_unique(std::integer_sequence<T, Is...> is)
	{
		size_t pos[]={find<T, Is>(is)...};
		for (size_t i = 0; i < sizeof...(Is); ++i)
		{
			if (pos[i]<i)
				return false;
		}
		return true;
	}

	template<class T, size_t POS, T... Ints>
	constexpr T get(std::integer_sequence<T, Ints...>) {
	    constexpr T arr[] = {Ints...};
	    return arr[POS];
	}

//	template<size_t POS, class... Ts>
//	constexpr auto get(type_pack<Ts...>) {
//		constexpr std::tuple<just_type<Ts>...> tpl;
//		return std::get<POS>(tpl);
//	}

	template <size_t I, class T>
	struct indexed_type
	{
		static constexpr size_t value = I;
		using type = T;
	};

	template <class IS, class... Ts>
	struct indexed_types;

	template <size_t... Is, class... Ts>
	struct indexed_types<std::index_sequence<Is...>, Ts...>
	{
		struct type : indexed_type<Is, Ts>... {};
	};

	template <class... Ts>
	using indexed_types_for = typename indexed_types<std::index_sequence_for<Ts...>, Ts...>::type;

	template <size_t I, class T>
	constexpr just_type<T> get_indexed_type(indexed_type<I, T>) { return {}; }

	template <size_t I, class... Ts>
	constexpr auto get(type_pack<Ts...>)
	{
		return get_indexed_type<I>(indexed_types_for<Ts...>{});
	}

	template <size_t I, class TP>
	struct type_pack_element_getter;

	template <size_t I, class ... Ts>
	struct type_pack_element_getter<I, type_pack<Ts...> >
	{
		using type = typename decltype(get<I>(type_pack<Ts...>()))::type;
	};

	template <size_t POS, typename TP>
	using type_pack_element = typename type_pack_element_getter<POS, TP>::type;

}

#endif /* SEQUENCE_H_ */
