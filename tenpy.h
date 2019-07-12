/*
 * tenpy.h
 *
 *  Created on: 11 июл. 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 *
 */

#ifndef TENPY_H_
#define TENPY_H_

#include <Python.h>

/* Disable old Numpy API */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/npy_math.h>
#include "numpy/arrayobject.h"

#include <tensor.h>

namespace iTTL
{

	class UnsupportedElementType: public exception
	{
		char *UnsupportedElementType_mes(int pn, const char type)
		{
			static char s[100];
			sprintf(s,"Error. Parameter #%d has unsupported element type '%c'", pn+1, type);
			return s;
		}
		char *UnsupportedElementType_mes(int pn, const char *type)
		{
			static char s[100];
			sprintf(s,"Error. Parameter #%d has unsupported element type '%s'", pn+1, type);
			return s;
		}
	public:
		UnsupportedElementType(int pn, const char type):exception(UnsupportedElementType_mes(pn, type)) {}
		UnsupportedElementType(int pn, const char *type):exception(UnsupportedElementType_mes(pn, type)) {}
	};
	class UnexpectedNumberOfDimensions: public exception
	{
		char *UnexpectedNumberOfDimensions_mes(int pn, int expected, int real)
		{
			static char s[200];
			sprintf(s,"Error. Parameter #%d has %d number of dimensions but %d is expected", pn+1, real, expected);
			return s;
		}
	public:
		UnexpectedNumberOfDimensions(int pn, int expected, int real):exception(UnexpectedNumberOfDimensions_mes(pn, expected, real)) {}
	};
	class UnsupportedStride: public exception
	{
		char *UnsupportedStride_mes(int pn, int dn, const char *type, size_t size, int stride)
		{
			static char s[200];
			sprintf(s,"Error. Parameter #%d has '%s' type of elements. Size of %s is %d. The PyArray_STRIDE returns %d for dimension %d. This stride is not supported", pn+1, type, type, (int)size, stride, dn);
			return s;
		}
	public:
		UnsupportedStride(int pn, int dn, const char *type, size_t size, int stride):exception(UnsupportedStride_mes(pn, dn, type, size, stride)) {}
	};

	static const int max_cont_type=5;
	template <typename T>
	int get_cont_type(PyArrayObject *A)
	{
		int ndim=PyArray_NDIM(A);
		npy_intp *d=PyArray_SHAPE(A);
		npy_intp *st=PyArray_STRIDES(A);
		size_t expectedStride=sizeof(T);
		int cont=0;
		for (int i=ndim-1;i>=0;i--)
		{
			if ((size_t)st[i]!=expectedStride)
				break;
			cont++;
			expectedStride*=d[i];
		}
//		if (cont>=ndim)
//			return max_cont_type;
		if (cont>max_cont_type)
			cont=max_cont_type;
		return cont;
	}
	template <>
	int get_cont_type<void>(PyArrayObject *A)
	{
		return 0;
	}

	template <int... V_TYPE>
	struct tensor_py
	{
		PyArrayObject *A;
	};

	template <int... V_TYPE>
	tensor_py<V_TYPE...> ind_tenpy(PyArrayObject *A, defaultIndex<V_TYPE>... ind)
	{
		return {A};
	}

	template <typename T, typename ST, int C_TYPE, int... V_TYPE>
	struct c_type_stuple;

	template <typename T, typename ... ST, int V_TYPE0, int... V_TYPE>
	struct c_type_stuple<T, type_sequence<ST...>, max_cont_type, V_TYPE0, V_TYPE...>:
		public c_type_stuple<T, type_sequence<ST..., segment<V_TYPE0, USAGE_FULL, 0, 0> >, max_cont_type, V_TYPE...>
	{
		typedef c_type_stuple<T, type_sequence<ST..., segment<V_TYPE0, USAGE_FULL, 0, 0> >, max_cont_type, V_TYPE...> base;
		static typename base::type wrap(int pn, PyArrayObject *A, npy_intp *d, npy_intp *st)
		{
			typename base::type ret=base::wrap(pn, A, d, st);
		    size_t step=st[sizeof...(ST)]/sizeof(T);
		    if (step*sizeof(T)!=(size_t)st[sizeof...(ST)])
			    	throw UnsupportedStride(pn, sizeof...(ST), typeid(T).name(), sizeof(T), st[sizeof...(ST)]);
//		    printf("Initializing %d array %d dimension with %d,%d\n",pn,(int)sizeof...(ST),(int)d[sizeof...(ST)],(int)step);
		    iTTL::get<sizeof...(ST)>(ret).init(d[sizeof...(ST)],step);
		    return ret;
		}
	};
	template <typename T, typename ... ST>
	struct c_type_stuple<T, type_sequence<ST...>, max_cont_type>
	{
		typedef tensor<T, stuple<sizeof...(ST), sizeof...(ST), true, void, ST...> > type;
		static type wrap(int pn, PyArrayObject *A, npy_intp *d, npy_intp *st)
		{
			int ndim=PyArray_NDIM(A);
			if ((size_t)ndim!=sizeof...(ST))
				throw UnexpectedNumberOfDimensions(pn, sizeof...(ST), ndim);
			return type(stuple<sizeof...(ST), sizeof...(ST), true, void, ST...>{}, (T *)PyArray_DATA(A));
		}
	};
	template <typename T, typename ... ST, int C_TYPE, int V_TYPE0, int... V_TYPE>
	struct c_type_stuple<T, type_sequence<ST...>, C_TYPE, V_TYPE0, V_TYPE...>:
		public std::conditional<(C_TYPE<=sizeof...(V_TYPE)),
					c_type_stuple<T, type_sequence<ST..., segment<V_TYPE0, USAGE_ENUM, 0, 0> >, C_TYPE, V_TYPE...>,
					typename std::conditional<(C_TYPE==sizeof...(V_TYPE)+1),
						c_type_stuple<T, type_sequence<ST..., segment<V_TYPE0, USAGE_CONT, 0, 0> >, C_TYPE, V_TYPE...>,
						c_type_stuple<T, type_sequence<ST..., segment<V_TYPE0, USAGE_FULL, 0, 0> >, C_TYPE, V_TYPE...>
					>::type
				>::type
	{
		typedef typename std::conditional<(C_TYPE<=sizeof...(V_TYPE)),
				c_type_stuple<T, type_sequence<ST..., segment<V_TYPE0, USAGE_ENUM, 0, 0> >, C_TYPE, V_TYPE...>,
				typename std::conditional<(C_TYPE==sizeof...(V_TYPE)+1),
					c_type_stuple<T, type_sequence<ST..., segment<V_TYPE0, USAGE_CONT, 0, 0> >, C_TYPE, V_TYPE...>,
					c_type_stuple<T, type_sequence<ST..., segment<V_TYPE0, USAGE_FULL, 0, 0> >, C_TYPE, V_TYPE...>
				>::type
			>::type base;
		static typename base::type wrap(int pn, PyArrayObject *A, npy_intp *d, npy_intp *st)
		{
			typename base::type ret=base::wrap(pn, A, d, st);
		    size_t step=st[sizeof...(ST)]/sizeof(T);
		    if (step*sizeof(T)!=(size_t)st[sizeof...(ST)])
			    	throw UnsupportedStride(pn, sizeof...(ST), typeid(T).name(), sizeof(T), st[sizeof...(ST)]);
//		    printf("Initializing %d array %d dimension with %d,%d\n",pn,(int)sizeof...(ST),(int)d[sizeof...(ST)],(int)step);
		    iTTL::get<sizeof...(ST)>(ret).init(d[sizeof...(ST)],step);
		    return ret;
		}
	};
	template <typename T, typename ... ST, int C_TYPE>
	struct c_type_stuple<T, type_sequence<ST...>, C_TYPE>
	{
		typedef tensor<T, stuple<sizeof...(ST), sizeof...(ST), true, void, ST...> > type;
		static type wrap(int pn, PyArrayObject *A, npy_intp *d, npy_intp *st)
		{
			int ndim=PyArray_NDIM(A);
			if ((size_t)ndim!=sizeof...(ST))
				throw UnexpectedNumberOfDimensions(pn, sizeof...(ST), ndim);
			return type(stuple<sizeof...(ST), sizeof...(ST), true, void, ST...>{}, (T *)PyArray_DATA(A));
		}
	};

	template <template <typename...> class clTmpl, typename CT, typename AT, typename ... args>
	struct objectCreatorPy;

	template <int MAXCONT, template <typename...> class clTmpl, typename T, typename CT, typename AT, typename ... args>
	struct objectCreatorPyMC;

	template <template <typename...> class clTmpl, typename T, typename ... CT, typename ... AT, int ... V_TYPE, typename ... Ts>
	struct objectCreatorPyMC<5, clTmpl, T, type_sequence<CT...>, type_sequence<AT...>, tensor_py<V_TYPE...>, Ts...>
	{
		inline static void create_object_py(AT... c_args, const tensor_py<V_TYPE...>& tpy, Ts... args)
		{
			int c_type=get_cont_type<T>(tpy.A);
//			printf("c_type for %d is %d\n", (int)sizeof...(AT), c_type);
			typedef c_type_stuple<T, type_sequence<>, 0, V_TYPE...> TT0;
			typedef c_type_stuple<T, type_sequence<>, 1, V_TYPE...> TT1;
			typedef c_type_stuple<T, type_sequence<>, 2, V_TYPE...> TT2;
			typedef c_type_stuple<T, type_sequence<>, 3, V_TYPE...> TT3;
			typedef c_type_stuple<T, type_sequence<>, 4, V_TYPE...> TT4;
			typedef c_type_stuple<T, type_sequence<>, 5, V_TYPE...> TT5;
			switch (c_type)
			{
			case 0:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT0::type>, type_sequence<AT..., typename TT0::type>, Ts...>::create_object_py(c_args...,
						TT0::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 1:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT1::type>, type_sequence<AT..., typename TT1::type>, Ts...>::create_object_py(c_args...,
						TT1::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 2:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT2::type>, type_sequence<AT..., typename TT2::type>, Ts...>::create_object_py(c_args...,
						TT2::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 3:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT3::type>, type_sequence<AT..., typename TT3::type>, Ts...>::create_object_py(c_args...,
						TT3::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 4:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT4::type>, type_sequence<AT..., typename TT4::type>, Ts...>::create_object_py(c_args...,
						TT4::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 5:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT5::type>, type_sequence<AT..., typename TT5::type>, Ts...>::create_object_py(c_args...,
						TT5::wrap(sizeof...(CT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			}
		}
	};

	template <template <typename...> class clTmpl, typename T, typename ... CT, typename ... AT, int ... V_TYPE, typename ... Ts>
	struct objectCreatorPyMC<4, clTmpl, T, type_sequence<CT...>, type_sequence<AT...>, tensor_py<V_TYPE...>, Ts...>
	{
		inline static void create_object_py(AT... c_args, const tensor_py<V_TYPE...>& tpy, Ts... args)
		{
			int c_type=get_cont_type<T>(tpy.A);
//			printf("c_type for %d is %d\n", (int)sizeof...(AT), c_type);
			typedef c_type_stuple<T, type_sequence<>, 0, V_TYPE...> TT0;
			typedef c_type_stuple<T, type_sequence<>, 1, V_TYPE...> TT1;
			typedef c_type_stuple<T, type_sequence<>, 2, V_TYPE...> TT2;
			typedef c_type_stuple<T, type_sequence<>, 3, V_TYPE...> TT3;
			typedef c_type_stuple<T, type_sequence<>, 4, V_TYPE...> TT4;
			switch (c_type)
			{
			case 0:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT0::type>, type_sequence<AT..., typename TT0::type>, Ts...>::create_object_py(c_args...,
						TT0::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 1:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT1::type>, type_sequence<AT..., typename TT1::type>, Ts...>::create_object_py(c_args...,
						TT1::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 2:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT2::type>, type_sequence<AT..., typename TT2::type>, Ts...>::create_object_py(c_args...,
						TT2::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 3:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT3::type>, type_sequence<AT..., typename TT3::type>, Ts...>::create_object_py(c_args...,
						TT3::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 4:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT4::type>, type_sequence<AT..., typename TT4::type>, Ts...>::create_object_py(c_args...,
						TT4::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			}
		}
	};

	template <template <typename...> class clTmpl, typename T, typename ... CT, typename ... AT, int ... V_TYPE, typename ... Ts>
	struct objectCreatorPyMC<3, clTmpl, T, type_sequence<CT...>, type_sequence<AT...>, tensor_py<V_TYPE...>, Ts...>
	{
		inline static void create_object_py(AT... c_args, const tensor_py<V_TYPE...>& tpy, Ts... args)
		{
			int c_type=get_cont_type<T>(tpy.A);
//			printf("c_type for %d is %d\n", (int)sizeof...(AT), c_type);
			typedef c_type_stuple<T, type_sequence<>, 0, V_TYPE...> TT0;
			typedef c_type_stuple<T, type_sequence<>, 1, V_TYPE...> TT1;
			typedef c_type_stuple<T, type_sequence<>, 2, V_TYPE...> TT2;
			typedef c_type_stuple<T, type_sequence<>, 3, V_TYPE...> TT3;
			switch (c_type)
			{
			case 0:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT0::type>, type_sequence<AT..., typename TT0::type>, Ts...>::create_object_py(c_args...,
						TT0::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 1:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT1::type>, type_sequence<AT..., typename TT1::type>, Ts...>::create_object_py(c_args...,
						TT1::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 2:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT2::type>, type_sequence<AT..., typename TT2::type>, Ts...>::create_object_py(c_args...,
						TT2::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 3:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT3::type>, type_sequence<AT..., typename TT3::type>, Ts...>::create_object_py(c_args...,
						TT3::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			}
		}
	};

	template <template <typename...> class clTmpl, typename T, typename ... CT, typename ... AT, int ... V_TYPE, typename ... Ts>
	struct objectCreatorPyMC<2, clTmpl, T, type_sequence<CT...>, type_sequence<AT...>, tensor_py<V_TYPE...>, Ts...>
	{
		inline static void create_object_py(AT... c_args, const tensor_py<V_TYPE...>& tpy, Ts... args)
		{
			int c_type=get_cont_type<T>(tpy.A);
//			printf("c_type for %d is %d\n", (int)sizeof...(AT), c_type);
			typedef c_type_stuple<T, type_sequence<>, 0, V_TYPE...> TT0;
			typedef c_type_stuple<T, type_sequence<>, 1, V_TYPE...> TT1;
			typedef c_type_stuple<T, type_sequence<>, 2, V_TYPE...> TT2;
			switch (c_type)
			{
			case 0:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT0::type>, type_sequence<AT..., typename TT0::type>, Ts...>::create_object_py(c_args...,
						TT0::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 1:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT1::type>, type_sequence<AT..., typename TT1::type>, Ts...>::create_object_py(c_args...,
						TT1::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 2:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT2::type>, type_sequence<AT..., typename TT2::type>, Ts...>::create_object_py(c_args...,
						TT2::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			}
		}
	};

	template <template <typename...> class clTmpl, typename T, typename ... CT, typename ... AT, int ... V_TYPE, typename ... Ts>
	struct objectCreatorPyMC<1, clTmpl, T, type_sequence<CT...>, type_sequence<AT...>, tensor_py<V_TYPE...>, Ts...>
	{
		inline static void create_object_py(AT... c_args, const tensor_py<V_TYPE...>& tpy, Ts... args)
		{
			int c_type=get_cont_type<T>(tpy.A);
//			printf("c_type for %d is %d\n", (int)sizeof...(AT), c_type);
			typedef c_type_stuple<T, type_sequence<>, 0, V_TYPE...> TT0;
			typedef c_type_stuple<T, type_sequence<>, 1, V_TYPE...> TT1;
			switch (c_type)
			{
			case 0:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT0::type>, type_sequence<AT..., typename TT0::type>, Ts...>::create_object_py(c_args...,
						TT0::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			case 1:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT1::type>, type_sequence<AT..., typename TT1::type>, Ts...>::create_object_py(c_args...,
						TT1::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			}
		}
	};

	template <template <typename...> class clTmpl, typename T, typename ... CT, typename ... AT, int ... V_TYPE, typename ... Ts>
	struct objectCreatorPyMC<0, clTmpl, T, type_sequence<CT...>, type_sequence<AT...>, tensor_py<V_TYPE...>, Ts...>
	{
		inline static void create_object_py(AT... c_args, const tensor_py<V_TYPE...>& tpy, Ts... args)
		{
			int c_type=get_cont_type<T>(tpy.A);
//			printf("c_type for %d is %d\n", (int)sizeof...(AT), c_type);
			typedef c_type_stuple<T, type_sequence<>, 0, V_TYPE...> TT0;
			switch (c_type)
			{
			case 0:
				objectCreatorPy<clTmpl, type_sequence<CT..., typename TT0::type>, type_sequence<AT..., typename TT0::type>, Ts...>::create_object_py(c_args...,
						TT0::wrap(sizeof...(AT), tpy.A, PyArray_SHAPE(tpy.A), PyArray_STRIDES(tpy.A)), args...);
				break;
			}
		}
	};

	template <template <typename...> class clTmpl, typename ... CT, typename ... AT, int ... V_TYPE, typename ... Ts>
	struct objectCreatorPy<clTmpl, type_sequence<CT...>, type_sequence<AT...>, tensor_py<V_TYPE...>, Ts...>
	{
		inline static void create_object_py(AT... c_args, const tensor_py<V_TYPE...>& tpy, Ts... args)
		{
			npy_intp datatype=PyArray_TYPE(tpy.A);
			switch (datatype)
			{
			case NPY_DOUBLE:
				objectCreatorPyMC<(sizeof...(V_TYPE)>max_cont_type?max_cont_type:(int)sizeof...(V_TYPE)),clTmpl, double, type_sequence<CT...>, type_sequence<AT...>, tensor_py<V_TYPE...>, Ts...>::create_object_py(c_args..., tpy, args...);
				return;
				break;
			case NPY_FLOAT:
				objectCreatorPyMC<(sizeof...(V_TYPE)>max_cont_type?max_cont_type:(int)sizeof...(V_TYPE)),clTmpl, float, type_sequence<CT...>, type_sequence<AT...>, tensor_py<V_TYPE...>, Ts...>::create_object_py(c_args..., tpy, args...);
				return;
				break;
			default:
				PyArray_Descr *descr=PyArray_DESCR(tpy.A);
				throw UnsupportedElementType(sizeof...(AT), descr->typeobj->tp_name);
			}
		}
	};


	template <template <typename...> class clTmpl, typename ... CT, typename ... AT, typename T0, typename ... Ts>
	struct objectCreatorPy<clTmpl, type_sequence<CT...>, type_sequence<AT...>, T0, Ts...>
	{
		inline static void create_object_py(AT... c_args, T0 arg0, Ts... args)
		{
			objectCreatorPy<clTmpl, type_sequence<CT...>, type_sequence<AT..., T0>, Ts...>::create_object_py(c_args..., arg0, args...);
		}
	};
	template <template <typename...> class clTmpl, typename ... CT, typename ... AT>
	struct objectCreatorPy<clTmpl, type_sequence<CT...>, type_sequence<AT...> >
	{
		inline static void create_object_py(AT... c_args)
		{
			clTmpl<CT...> obj(c_args...);
		}
	};
//	template <template <typename...> class clTmpl, typename ... CT, typename ... AT>
//	struct objectCreatorPy<clTmpl, void, type_sequence<CT...>, type_sequence<AT...> >
//	{
//		inline static void create_object_py(AT... c_args)
//		{
//		}
//	};

	template <template <typename...> class clTmpl, typename ... Ts>
	void createObjectPy(Ts... args)
	{
		objectCreatorPy<clTmpl, iTTL::type_sequence<>, iTTL::type_sequence<>, Ts...>::create_object_py(args...);
	}
};

#endif /* TENPY_H_ */
