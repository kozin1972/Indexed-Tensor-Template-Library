/*
 * expression.h
 *
 *  Created on: 27 мая 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef TEXPRESS_H_
#define TEXPRESS_H_

namespace tpp
{

	template <typename T, typename ST>
	class tensor;

	template <typename T, typename ST>
	struct number_tensor_c
	{
		T number;
		const tensor<T, ST>& t;
		number_tensor_c(T number, const tensor<T, ST>& t): number(number), t(t) {}
	};
	template <typename T, typename ST>
	struct number_tensor
	{
		T number;
		tensor<T, ST>& t;
		number_tensor(T number, const tensor<T, ST>& t): number(number), t(t) {}
	};
	template <typename T, typename ST>
	number_tensor_c<T,ST> operator*(T number, const tensor<T,ST>& t) { return number_tensor_c<T, ST>(number, t); }
	template <typename T, typename ST>
	number_tensor_c<T,ST> operator*(const tensor<T,ST>& t, T number) { return number_tensor_c<T, ST>(number, t); }
	template <typename T, typename ST>
	number_tensor_c<T,ST> operator/(const tensor<T,ST>& t, T number) { return number_tensor_c<T, ST>(1.0/number, t); }
	template <typename T, typename ST>
	number_tensor<T,ST> operator*(T number, tensor<T,ST>& t) { return number_tensor<T, ST>(number, t); }
	template <typename T, typename ST>
	number_tensor<T,ST> operator*(tensor<T,ST>& t, T number) { return number_tensor<T, ST>(number, t); }
	template <typename T, typename ST>
	number_tensor<T,ST> operator/(tensor<T,ST>& t, T number) { return number_tensor<T, ST>(1.0/number, t); }
	template <typename T, typename ST>
	number_tensor_c<T,ST> operator*(int number, const tensor<T,ST>& t) { return number_tensor_c<T, ST>(number, t); }
	template <typename T, typename ST>
	number_tensor_c<T,ST> operator*(const tensor<T,ST>& t, int number) { return number_tensor_c<T, ST>(number, t); }
	template <typename T, typename ST>
	number_tensor_c<T,ST> operator/(const tensor<T,ST>& t, int number) { return number_tensor_c<T, ST>(1.0/(T)number, t); }
	template <typename T, typename ST>
	number_tensor<T,ST> operator*(int number, tensor<T,ST>& t) { return number_tensor<T, ST>(number, t); }
	template <typename T, typename ST>
	number_tensor<T,ST> operator*(tensor<T,ST>& t, int number) { return number_tensor<T, ST>(number, t); }
	template <typename T, typename ST>
	number_tensor<T,ST> operator/(tensor<T,ST>& t, int number) { return number_tensor<T, ST>(1.0/(T)number, t); }
	template <typename T, typename ST0, typename ST1>
	struct number_tt_c
	{
		T number;
		const tensor<T, ST0>& t0;
		const tensor<T, ST1>& t1;
		number_tt_c(T number, const tensor<T, ST0>& t0, const tensor<T, ST1>& t1): number(number), t0(t0), t1(t1) {}
	};
	template <typename T, typename ST0, typename ST1>
	number_tt_c<T,ST0,ST1> operator*(const tensor<T, ST0>& t0, const tensor<T, ST1>& t1) { return number_tt_c<T,ST0,ST1>(1.0,t0,t1); }
	template <typename T, typename ST0, typename ST1>
	number_tt_c<T,ST0,ST1> operator*(tensor<T, ST0>&& t0, tensor<T, ST1>&& t1) { return number_tt_c<T,ST0,ST1>(1.0,t0,t1); }
	template <typename T, typename ST0, typename ST1>
	number_tt_c<T,ST0,ST1> operator*(T number, const number_tt_c<T,ST0,ST1>& ncc) { return number_tt_c<T,ST0,ST1>(ncc.number*number,ncc.t0,ncc.t1); }
	template <typename T, typename ST0, typename ST1>
	number_tt_c<T,ST0,ST1> operator*(const number_tt_c<T,ST0,ST1>& ncc, T number) { return number_tt_c<T,ST0,ST1>(ncc.number*number,ncc.t0,ncc.t1); }
	template <typename T, typename ST0, typename ST1>
	number_tt_c<T,ST0,ST1> operator/(const number_tt_c<T,ST0,ST1>& ncc, T number) { return number_tt_c<T,ST0,ST1>(ncc.number/number,ncc.t0,ncc.t1); }
	template <typename T, typename ST0, typename ST1>
	number_tt_c<T,ST0,ST1> operator*(int number, const number_tt_c<T,ST0,ST1>& ncc) { return number_tt_c<T,ST0,ST1>(ncc.number*(T)number,ncc.t0,ncc.t1); }
	template <typename T, typename ST0, typename ST1>
	number_tt_c<T,ST0,ST1> operator*(const number_tt_c<T,ST0,ST1>& ncc, int number) { return number_tt_c<T,ST0,ST1>(ncc.number*(T)number,ncc.t0,ncc.t1); }
	template <typename T, typename ST0, typename ST1>
	number_tt_c<T,ST0,ST1> operator/(const number_tt_c<T,ST0,ST1>& ncc, int number) { return number_tt_c<T,ST0,ST1>(ncc.number/(T)number,ncc.t0,ncc.t1); }
	template <typename T, typename ST0, typename ST1>
	number_tt_c<T,ST0,ST1> operator*(const tensor<T, ST0>& t, const number_tensor_c<T,ST1>& nt) { return number_tt_c<T,ST0,ST1>(nt.number,t,nt.t); }
	template <typename T, typename ST0, typename ST1>
	number_tt_c<T,ST0,ST1> operator*(const number_tensor_c<T,ST1>& nt, const tensor<T, ST0>& t) { return number_tt_c<T,ST0,ST1>(nt.number,t,nt.t); }
	template <typename T, typename ST0, typename ST1>
	number_tt_c<T,ST0,ST1> operator*(const number_tensor_c<T,ST0>& nt0, const number_tensor_c<T,ST1>& nt1) { return number_tt_c<T,ST0,ST1>(nt0.number*nt1.number,nt0.t,nt1.t); }
};

#endif /* TEXPRESS_H_ */
