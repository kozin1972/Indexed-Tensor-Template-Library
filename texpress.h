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

namespace iTTL
{

	template <typename T, typename ST>
	class tensor;

	template <typename T, typename ST>
	struct number_tensor_c
	{
		T number;
		const tensor<T, ST>& t;
		constexpr number_tensor_c(T number, const tensor<T, ST>& t): number(number), t(t) {}
	};
	template <typename T, bool IS_INDEXED>
	constexpr T operator-(const tensor<T, stuple<0,0,IS_INDEXED,void> >& t) { return -*t.data_ptr(); }
	template <typename T, typename ST>
	constexpr number_tensor_c<T, ST> operator-(const tensor<T, ST>& t) { return {-1.0, t}; }
	template <typename T, typename ST>
	constexpr number_tensor_c<T, ST> operator-(const number_tensor_c<T, ST>& t) { return {-t.number, t.t}; }
	template <typename T, typename ST>
	constexpr number_tensor_c<T,ST> operator*(T number, const tensor<T,ST>& t) { return {number, t}; }
	template <typename T, typename ST>
	constexpr number_tensor_c<T,ST> operator*(const tensor<T,ST>& t, T number) { return {number, t}; }
	template <typename T, typename ST>
	constexpr number_tensor_c<T,ST> operator/(const tensor<T,ST>& t, T number) { return {1.0/number, t}; }
	template <typename T, typename ST>
	constexpr number_tensor_c<T,ST> operator*(int number, const tensor<T,ST>& t) { return {(T)number, t}; }
	template <typename T, typename ST>
	constexpr number_tensor_c<T,ST> operator*(const tensor<T,ST>& t, int number) { return {(T)number, t}; }
	template <typename T, typename ST>
	constexpr number_tensor_c<T,ST> operator/(const tensor<T,ST>& t, int number) { return {1.0/(T)number, t}; }
	template <typename T, typename ST0, typename ST1>
	struct number_tt_c
	{
		T number;
		const tensor<T, ST0>& t0;
		const tensor<T, ST1>& t1;
		constexpr number_tt_c(T number, const tensor<T, ST0>& t0, const tensor<T, ST1>& t1): number(number), t0(t0), t1(t1) {}
	};
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator-(const number_tt_c<T,ST0,ST1>& t) { return {-t.number, t.t0, t.t1}; }
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator*(const tensor<T, ST0>& t0, const tensor<T, ST1>& t1) { return {1.0, t0, t1}; }
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator*(tensor<T, ST0>&& t0, tensor<T, ST1>&& t1) { return {1.0, t0, t1}; }
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator*(T number, const number_tt_c<T,ST0,ST1>& ncc) { return {ncc.number*number, ncc.t0, ncc.t1}; }
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator*(const number_tt_c<T,ST0,ST1>& ncc, T number) { return {ncc.number*number, ncc.t0, ncc.t1}; }
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator/(const number_tt_c<T,ST0,ST1>& ncc, T number) { return {ncc.number/number, ncc.t0, ncc.t1}; }
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator*(int number, const number_tt_c<T,ST0,ST1>& ncc) { return {ncc.number*(T)number, ncc.t0, ncc.t1}; }
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator*(const number_tt_c<T,ST0,ST1>& ncc, int number) { return {ncc.number*(T)number, ncc.t0, ncc.t1}; }
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator/(const number_tt_c<T,ST0,ST1>& ncc, int number) { return {ncc.number/(T)number, ncc.t0, ncc.t1}; }
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator*(const tensor<T, ST0>& t, const number_tensor_c<T,ST1>& nt) { return {nt.number, t, nt.t}; }
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator*(const number_tensor_c<T,ST1>& nt, const tensor<T, ST0>& t) { return {nt.number, t, nt.t}; }
	template <typename T, typename ST0, typename ST1>
	constexpr number_tt_c<T,ST0,ST1> operator*(const number_tensor_c<T,ST0>& nt0, const number_tensor_c<T,ST1>& nt1) { return {nt0.number*nt1.number, nt0.t, nt1.t}; }
	template <typename T, typename T2>
	struct sum_num_struct
	{
		T summand;
		T multiplier;
		const T2& t2;
		constexpr sum_num_struct(T summand, T multiplier, const T2& t2): summand(summand), multiplier(multiplier), t2(t2) {}
	};
	template <typename T, typename T1, typename T2>
	struct sum_struct
	{
		T m1;
		const T1& t1;
		T m2;
		const T2& t2;
		constexpr sum_struct(T m1, const T1& t1, T m2, const T2& t2): m1(m1), t1(t1), m2(m2), t2(t2) {}
	};
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator-(const sum_num_struct<T, T2>& sn) { return {-sn.summand, -sn.multiplier, sn.t2}; }
	template <typename T, bool IS_INDEXED>
	constexpr T operator+(T number, const tensor<T, stuple<0,0,IS_INDEXED,void> >& t) { return number+*t.data_ptr(); }
	template <typename T, bool IS_INDEXED>
	constexpr T operator+(const tensor<T, stuple<0,0,IS_INDEXED,void> >& t, T number) { return number+*t.data_ptr(); }
	template <typename T, bool IS_INDEXED>
	constexpr T operator-(T number, const tensor<T, stuple<0,0,IS_INDEXED,void> >& t) { return number-*t.data_ptr(); }
	template <typename T, bool IS_INDEXED>
	constexpr T operator+(int number, const tensor<T, stuple<0,0,IS_INDEXED,void> >& t) { return (T)number+*t.data_ptr(); }
	template <typename T, bool IS_INDEXED>
	constexpr T operator+(const tensor<T, stuple<0,0,IS_INDEXED,void> >& t, int number) { return (T)number+*t.data_ptr(); }
	template <typename T, bool IS_INDEXED>
	constexpr T operator-(int number, const tensor<T, stuple<0,0,IS_INDEXED,void> >& t) { return (T)number-*t.data_ptr(); }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator+(T number, const tensor<T, ST>& t) { return {number, 1.0, t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator-(T number, const tensor<T, ST>& t) { return {number, -1.0, t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator+(const tensor<T, ST>& t, T number) { return {number, 1.0, t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator-(const tensor<T, ST>& t, T number) { return {-number, 1.0, t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator+(int number, const tensor<T, ST>& t) { return {(T)number, 1.0, t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator-(int number, const tensor<T, ST>& t) { return {(T)number, -1.0, t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator+(const tensor<T, ST>& t, int number) { return {(T)number, 1.0, t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator-(const tensor<T, ST>& t, int number) { return {-(T)number, 1.0, t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator+(T number, const number_tensor_c<T, ST>& t) { return {number, t.number, t.t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator-(T number, const number_tensor_c<T, ST>& t) { return {number, -t.number, t.t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator+(const number_tensor_c<T, ST>& t, T number) { return {number, t.number, t.t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator-(const number_tensor_c<T, ST>& t, T number) { return {-number, t.number, t.t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator+(int number, const number_tensor_c<T, ST>& t) { return {(T)number, t.number, t.t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator-(int number, const number_tensor_c<T, ST>& t) { return {(T)number, -t.number, t.t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator+(const number_tensor_c<T, ST>& t, int number) { return {(T)number, t.number, t.t}; }
	template <typename T, typename ST>
	constexpr sum_num_struct<T, tensor<T, ST> > operator-(const number_tensor_c<T, ST>& t, int number) { return {-(T)number, t.number, t.t}; }
	template <typename T, typename ST0, typename ST1>
	constexpr sum_num_struct<T, number_tt_c<T, ST0, ST1> > operator+(T number, const number_tt_c<T, ST0, ST1>& t) { return {number, 1.0, t}; }
	template <typename T, typename ST0, typename ST1>
	constexpr sum_num_struct<T, number_tt_c<T, ST0, ST1> > operator-(T number, const number_tt_c<T, ST0, ST1>& t) { return {number, -1.0, t}; }
	template <typename T, typename ST0, typename ST1>
	constexpr sum_num_struct<T, number_tt_c<T, ST0, ST1> > operator+(const number_tt_c<T, ST0, ST1>& t, T number) { return {number, 1.0, t}; }
	template <typename T, typename ST0, typename ST1>
	constexpr sum_num_struct<T, number_tt_c<T, ST0, ST1> > operator-(const number_tt_c<T, ST0, ST1>& t, T number) { return {-number, 1.0, t}; }
	template <typename T, typename ST0, typename ST1>
	constexpr sum_num_struct<T, number_tt_c<T, ST0, ST1> > operator+(int number, const number_tt_c<T, ST0, ST1>& t) { return {(T)number, 1.0, t}; }
	template <typename T, typename ST0, typename ST1>
	constexpr sum_num_struct<T, number_tt_c<T, ST0, ST1> > operator-(int number, const number_tt_c<T, ST0, ST1>& t) { return {(T)number, -1.0, t}; }
	template <typename T, typename ST0, typename ST1>
	constexpr sum_num_struct<T, number_tt_c<T, ST0, ST1> > operator+(const number_tt_c<T, ST0, ST1>& t, int number) { return {(T)number, 1.0, t}; }
	template <typename T, typename ST0, typename ST1>
	constexpr sum_num_struct<T, number_tt_c<T, ST0, ST1> > operator-(const number_tt_c<T, ST0, ST1>& t, int number) { return {-(T)number, 1.0, t}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator+(T number, const sum_num_struct<T, T2>& t) { return {number+t.summand, t.multiplier, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator-(T number, const sum_num_struct<T, T2>& t) { return {number-t.summand, -t.multiplier, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator+(const sum_num_struct<T, T2>& t, T number) { return {number+t.summand, t.multiplier, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator-(const sum_num_struct<T, T2>& t, T number) { return {-number+t.summand, t.multiplier, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator+(int number, const sum_num_struct<T, T2>& t) { return {(T)number+t.summand, t.multiplier, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator-(int number, const sum_num_struct<T, T2>& t) { return {(T)number-t.summand, -t.multiplier, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator+(const sum_num_struct<T, T2>& t, int number) { return {(T)number+t.summand, t.multiplier, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator-(const sum_num_struct<T, T2>& t, int number) { return {-(T)number+t.summand, t.multiplier, t.t2}; }

	template <typename T, typename T1, typename T2>
	constexpr sum_struct<T, T1, T2> operator-(const sum_struct<T, T1, T2>& ss) { return {-ss.m1, ss.t1, -ss.m2, ss.t2}; }
	template <typename T, bool IS_INDEXED0, bool IS_INDEXED1>
	constexpr T operator+(const tensor<T, stuple<0,0,IS_INDEXED0,void> >& t0, const tensor<T, stuple<0,0,IS_INDEXED1,void> >& t1) { return *t0.data_ptr()+*t1.data_ptr(); }
	template <typename T, bool IS_INDEXED0, bool IS_INDEXED1>
	constexpr T operator-(const tensor<T, stuple<0,0,IS_INDEXED0,void> >& t0, const tensor<T, stuple<0,0,IS_INDEXED1,void> >& t1) { return *t0.data_ptr()-*t1.data_ptr(); }

	template <typename T, typename T0, typename T1>
	constexpr sum_num_struct<T, sum_struct<T, T0, T1> > operator+(T number, const sum_struct<T, T0, T1>& t) { return {number, 1.0, t}; }
	template <typename T, typename T0, typename T1>
	constexpr sum_num_struct<T, sum_struct<T, T0, T1> > operator-(T number, const sum_struct<T, T0, T1>& t) { return {number, -1.0, t}; }
	template <typename T, typename T0, typename T1>
	constexpr sum_num_struct<T, sum_struct<T, T0, T1> > operator+(const sum_struct<T, T0, T1>& t, T number) { return {number, 1.0, t}; }
	template <typename T, typename T0, typename T1>
	constexpr sum_num_struct<T, sum_struct<T, T0, T1> > operator-(const sum_struct<T, T0, T1>& t, T number) { return {-number, 1.0, t}; }
	template <typename T, typename T0, typename T1>
	constexpr sum_num_struct<T, sum_struct<T, T0, T1> > operator+(int number, const sum_struct<T, T0, T1>& t) { return {(T)number, 1.0, t}; }
	template <typename T, typename T0, typename T1>
	constexpr sum_num_struct<T, sum_struct<T, T0, T1> > operator-(int number, const sum_struct<T, T0, T1>& t) { return {(T)number, -1.0, t}; }
	template <typename T, typename T0, typename T1>
	constexpr sum_num_struct<T, sum_struct<T, T0, T1> > operator+(const sum_struct<T, T0, T1>& t, int number) { return {(T)number, 1.0, t}; }
	template <typename T, typename T0, typename T1>
	constexpr sum_num_struct<T, sum_struct<T, T0, T1> > operator-(const sum_struct<T, T0, T1>& t, int number) { return {-(T)number, 1.0, t}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator*(const sum_num_struct<T, T2>& t, T number) { return {number*t.summand, number*t.multiplier, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator/(const sum_num_struct<T, T2>& t, T number) { return {t.summand/number, t.multiplier/number, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator*(T number, const sum_num_struct<T, T2>& t) { return {number*t.summand, number*t.multiplier, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator*(const sum_num_struct<T, T2>& t, int number) { return {(T)number*t.summand, (T)number*t.multiplier, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator/(const sum_num_struct<T, T2>& t, int number) { return {t.summand/(T)number, t.multiplier/(T)number, t.t2}; }
	template <typename T, typename T2>
	constexpr sum_num_struct<T, T2> operator*(int number, const sum_num_struct<T, T2>& t) { return {(T)number*t.summand, (T)number*t.multiplier, t.t2}; }

	template <typename T, typename ST, typename T2>
	constexpr sum_struct<T, tensor<T, ST>, T2 > operator+(const tensor<T, ST>& t0, const T2& t2) { return {1.0, t0, 1.0, t2}; }
	template <typename T, typename ST, typename T2>
	constexpr sum_struct<T, tensor<T, ST>, T2 > operator-(const tensor<T, ST>& t0, const T2& t2) { return {1.0, t0, -1.0, t2}; }
	template <typename T, typename ST, typename T2>
	constexpr sum_struct<T, tensor<T, ST>, T2 > operator+(const number_tensor_c<T, ST>& t0, const T2& t2) { return {t0.number, t0.t, 1.0, t2}; }
	template <typename T, typename ST, typename T2>
	constexpr sum_struct<T, tensor<T, ST>, T2 > operator-(const number_tensor_c<T, ST>& t0, const T2& t2) { return {t0.number, t0.t, -1.0, t2}; }
	template <typename T, typename ST0, typename ST1, typename T2>
	constexpr sum_struct<T, number_tt_c<T, ST0, ST1>, T2 > operator+(const number_tt_c<T, ST0, ST1>& t0, const T2& t2) { return {1.0, t0, 1.0, t2}; }
	template <typename T, typename ST0, typename ST1, typename T2>
	constexpr sum_struct<T, number_tt_c<T, ST0, ST1>, T2 > operator-(const number_tt_c<T, ST0, ST1>& t0, const T2& t2) { return {1.0, t0, -1.0, t2}; }
	template <typename T, typename T1, typename T2>
	constexpr sum_struct<T, sum_num_struct<T, T1>, T2 > operator+(const sum_num_struct<T, T1>& sn, const T2& t2) { return {1.0, sn, 1.0, t2}; }
	template <typename T, typename T1, typename T2>
	constexpr sum_struct<T, sum_num_struct<T, T1>, T2 > operator-(const sum_num_struct<T, T1>& sn, const T2& t2) { return {1.0, sn, -1.0, t2}; }
	template <typename T, typename T1, typename T2, typename T3>
	constexpr sum_struct<T, sum_struct<T, T1, T2>, T3 > operator+(const sum_struct<T, T1, T2>& ss, const T3& t3) { return {1.0, ss, 1.0, t3}; }
	template <typename T, typename T1, typename T2, typename T3>
	constexpr sum_struct<T, sum_struct<T, T1, T2>, T3 > operator-(const sum_struct<T, T1, T2>& ss, const T3& t3) { return {1.0, ss, -1.0, t3}; }
	template <typename T, typename T1, typename T2>
	constexpr sum_struct<T, T1, T2> operator*(const sum_struct<T, T1, T2>& ss, T number) { return {ss.m1*number, ss.t1, ss.m2*number, ss.t2}; }
	template <typename T, typename T1, typename T2>
	constexpr sum_struct<T, T1, T2> operator/(const sum_struct<T, T1, T2>& ss, T number) { return {ss.m1/number, ss.t1, ss.m2/number, ss.t2}; }
	template <typename T, typename T1, typename T2>
	constexpr sum_struct<T, T1, T2> operator*(T number, const sum_struct<T, T1, T2>& ss) { return {ss.m1*number, ss.t1, ss.m2*number, ss.t2}; }
	template <typename T, typename T1, typename T2>
	constexpr sum_struct<T, T1, T2> operator*(const sum_struct<T, T1, T2>& ss, int number) { return {ss.m1*(T)number, ss.t1, ss.m2*(T)number, ss.t2}; }
	template <typename T, typename T1, typename T2>
	constexpr sum_struct<T, T1, T2> operator/(const sum_struct<T, T1, T2>& ss, int number) { return {ss.m1/(T)number, ss.t1, ss.m2/(T)number, ss.t2}; }
	template <typename T, typename T1, typename T2>
	constexpr sum_struct<T, T1, T2> operator*(int number, const sum_struct<T, T1, T2>& ss) { return {ss.m1*(T)number, ss.t1, ss.m2*(T)number, ss.t2}; }

};

#endif /* TEXPRESS_H_ */
