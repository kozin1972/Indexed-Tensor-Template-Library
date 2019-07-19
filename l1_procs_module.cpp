/*
 * l1_procs_module.cpp
 *
 *  Created on: 10 июл. 2019 г.
 *      Author: Alexey Kozin
 *
 *  Python extension module l1_procs
 *  BR_solve function is exported to python
 *  BR_solve(y,A,x,ridge=0.0)
 *  y - input vector
 *  A - input matrix
 *  x - solution vector
 *  ridge - optional input coefficient
 *  Barrodale-Roberts method.
 *  Finds x to optimize
 *  sum[ sum(abs( y-Ax )) + abs(ridge*x) ] --> min
 *  or
 *  sum[ sum(abs( y-(A**T)x )) + abs(ridge*x) ] --> min
 *
 *  The length of y is usually greater than the length of x.
 *  The method transposes A if and only if Ax cannot be calculated but (A**T)x can be calculated.
 *
 *  Element data types must be either double for all numpy arrays or float for all numpy arrays
 *  The type of ridge is double
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#include <tenpy.h>

#include <l1_procs.h>

//template <typename TTY, typename TTA, typename TTX>
//class BR_solver_runtime;
//
////template <typename TY, typename STY, typename TA, typename STA, typename TX, typename STX>
////class BR_solver_runtime<iTTL::tensor<TY, STY>, iTTL::tensor<TA, STA>, iTTL::tensor<TX, STX> >
////{
////public:
////	BR_solver_runtime(const iTTL::tensor<TY, STY>& y, const iTTL::tensor<TA, STA>& A, iTTL::tensor<TX, STX>& x, double ridge=0.0)
////	{
////		if (typeid(TY)!=typeid(TA))
////			throw iTTL::DifferentTypesOfElements("BR_solver parameters #1 (Y) and #2 (A)", typeid(TY).name(), typeid(TA).name());
////		if (typeid(TX)!=typeid(TA))
////			throw iTTL::DifferentTypesOfElements("BR_solver parameters #2 (A) and #3 (X)", typeid(TA).name(), typeid(TX).name());
////	}
////};
//
//template <typename T, typename STY, typename STA, typename STX>
//class BR_solver_runtime<iTTL::tensor<T, STY>, iTTL::tensor<T, STA>, iTTL::tensor<T, STX> >:
//	public BR_solver<iTTL::tensor<T, STY>, iTTL::tensor<T, STA>, iTTL::tensor<T, STX> >
//{
//public:
//	inline BR_solver_runtime(const iTTL::tensor<T, STY>& y, const iTTL::tensor<T, STA>& A, iTTL::tensor<T, STX>& x, T ridge=0.0): BR_solver<iTTL::tensor<T, STY>, iTTL::tensor<T, STA>, iTTL::tensor<T, STX> >(y, A, x, ridge) {}
//};

/* Errors */
static PyObject *L1_Error;

/*
 * Enter point from Python
 */
static PyObject*
BR_solve(PyObject* self, PyObject* args)
{
    PyArrayObject *Y=NULL;
    PyArrayObject *A=NULL;
    PyArrayObject *X=NULL;
    int ndY = 0, ndA = 0, ndX = 0;
    double ridge=0.0;

    if (!PyArg_ParseTuple(args, "O!O!O!|d", &PyArray_Type, &Y, &PyArray_Type, &A, &PyArray_Type, &X, &ridge))
    {
//    	printf("Parse failed\n");
    	return NULL;
    }

// Get number of dimensions
    ndY = PyArray_NDIM(Y);
    ndA = PyArray_NDIM(A);
    ndX = PyArray_NDIM(X);
    if (ndY != 1)
    {
        PyErr_SetString(L1_Error, "Error. The first argument (input vector) must have 1 dimension");
        return NULL;
    }
    if (ndA != 2)
    {
        PyErr_SetString(L1_Error, "Error. The second argument (matrix) must have 2 dimension");
        return NULL;
    }
    if (ndX != 1)
    {
        PyErr_SetString(L1_Error, "Error. The third argument (output vector) must have 1 dimension");
        return NULL;
    }

//    PyArray_Descr *Yt=PyArray_DESCR(Y);
//    PyArray_Descr *At=PyArray_DESCR(A);
//    PyArray_Descr *Xt=PyArray_DESCR(X);
//	if (Yt->type_num!=At->type_num)
//	{
//		char errText[200];
//		sprintf(errText,"Error. Different element types '%s' and '%s' are passed to BR_solver parameters #1 (Y) and #2 (A)", Yt->typeobj->tp_name, At->typeobj->tp_name);
//        PyErr_SetString(L1_Error, errText);
//        return NULL;
//	}
//	if (At->type_num!=Xt->type_num)
//	{
//		char errText[200];
//		sprintf(errText,"Error. Different element types '%s' and '%s' are passed to BR_solver parameters #2 (A) and #3 (X)", At->typeobj->tp_name, Xt->typeobj->tp_name);
//        PyErr_SetString(L1_Error, errText);
//        return NULL;
//	}

//    npy_intp datatype=PyArray_TYPE(Y);
//    if (datatype!=PyArray_TYPE(A) || datatype!=PyArray_TYPE(X))
//    {
//        PyErr_SetString(L1_Error, "All arrays must have the same data type of elements");
//        return NULL;
//    }
//    int element_size=0;
//    switch (datatype)
//    {
//    case NPY_DOUBLE:
//    	element_size=8;
//    	break;
//    }
//    if (element_size==0)
//    {
//        PyErr_SetString(L1_Error, "This type of array elements is not supported currently");
//        return NULL;
//    }


// Get shapes of arrays
    npy_intp *dY = PyArray_SHAPE(Y);
    npy_intp *dA = PyArray_SHAPE(A);
    npy_intp *dX = PyArray_SHAPE(X);

    if (dY[0]<dX[0])
    {
    	PyErr_SetString(L1_Error, "Error. Source vector must be longer or equal to the destination vector");
    	return NULL;
    }

//    npy_intp *stY=PyArray_STRIDES(Y);
//    npy_intp *stA=PyArray_STRIDES(A);
//    npy_intp *stX=PyArray_STRIDES(X);
//
//    size_t stepY[1];
//    stepY[0]=stY[0]/element_size;
//    size_t stepA[2];
//    stepA[0]=stA[0]/element_size;
//    stepA[1]=stA[1]/element_size;
//    size_t stepX[1];
//    stepX[0]=stX[0]/element_size;
//    if (stepY[0]*element_size!=(size_t)stY[0])
//    {
//        PyErr_SetString(L1_Error, "The first argument (input vector) has bad stride");
//        return NULL;
//    }
//    if (stepA[0]*element_size!=(size_t)stA[0] || stepA[1]*element_size!=(size_t)stA[1])
//    {
//        PyErr_SetString(L1_Error, "The second argument (matrix) has bad stride");
//        return NULL;
//    }
//    if (stepX[0]*element_size!=(size_t)stX[0])
//    {
//        PyErr_SetString(L1_Error, "The third argument (output vector) has bad stride");
//        return NULL;
//    }

    if (dY[0]==dA[0] && dX[0]==dA[1])
    {
    	try
    	{
    		DECLARE_defaultIndex(I);
    		DECLARE_defaultIndex(J);
    		auto yi=iTTL::ind_tenpy(Y,I);
    		auto uij=iTTL::ind_tenpy(A,I,J);
    		auto vj=iTTL::ind_tenpy(X,J);
//    		printf("Normal A\n");

    		iTTL::createObjectPySimple<BR_solver>(yi,uij,vj,ridge);
//    		iTTL::classCreator<BR_solver, void, iTTL::type_sequence<>, iTTL::type_sequence<>, decltype(yi), decltype(uij), decltype(vj), double>::create_class(yi,uij,vj,ridge);
//			typedef iTTL::segment<1, iTTL::USAGE_CONT, 0, 0> sYt;
//			typedef iTTL::segment<1, iTTL::USAGE_CONT, 0, 0> sA0t;
//			typedef iTTL::segment<2, iTTL::USAGE_CONT, 0, 0> sA1t;
//			typedef iTTL::segment<2, iTTL::USAGE_CONT, 0, 0> sXt;
//			iTTL::stuple<1, 0, true, void, sYt> stplY;
//			iTTL::get<0>(stplY).init(dY[0],stepY[0]);
//			iTTL::stuple<2, 0, true, void, sA0t, sA1t> stplA;
//			iTTL::get<0>(stplA).init(dA[0],stepA[0]);
//			iTTL::get<1>(stplA).init(dA[1],stepA[1]);
//			iTTL::stuple<1, 0, true, void, sXt> stplX;
//			iTTL::get<0>(stplX).init(dX[0],stepX[0]);
//			iTTL::tensor<double, iTTL::stuple<1, 0, true, void, sYt> > tY(stplY, (double *)PyArray_DATA(Y));
//			iTTL::tensor<double, iTTL::stuple<2, 0, true, void, sA0t, sA1t> > tA(stplA, (double *)PyArray_DATA(A));
//			iTTL::tensor<double, iTTL::stuple<1, 0, true, void, sXt> > tX(stplX, (double *)PyArray_DATA(X));
//			BR_solve_one(tY, tA, tX, ridge);
    	}
    	catch (std::exception& ex)
    	{
        	PyErr_SetString(L1_Error, ex.what());
        	return NULL;
    	}
    }
    else if (dY[0]==dA[1] && dX[0]==dA[0])
    {
    	try
    	{
    		DECLARE_defaultIndex(I);
    		DECLARE_defaultIndex(J);
    		auto yi=iTTL::ind_tenpy(Y,I);
    		auto uji=iTTL::ind_tenpy(A,J,I);
    		auto vj=iTTL::ind_tenpy(X,J);
//    		printf("Transposed A\n");
    		iTTL::createObjectPySimple<BR_solver>(yi,uji,vj,ridge);
//
//    		typedef iTTL::segment<1, iTTL::USAGE_CONT, 0, 0> sYt;
//			typedef iTTL::segment<2, iTTL::USAGE_CONT, 0, 0> sA0t;
//			typedef iTTL::segment<1, iTTL::USAGE_CONT, 0, 0> sA1t;
//			typedef iTTL::segment<2, iTTL::USAGE_CONT, 0, 0> sXt;
//			iTTL::stuple<1, 0, true, void, sYt> stplY;
//			iTTL::get<0>(stplY).init(dY[0],stepY[0]);
//			iTTL::stuple<2, 0, true, void, sA0t, sA1t> stplA;
//			iTTL::get<0>(stplA).init(dA[0],stepA[0]);
//			iTTL::get<1>(stplA).init(dA[1],stepA[1]);
//			iTTL::stuple<1, 0, true, void, sXt> stplX;
//			iTTL::get<0>(stplX).init(dX[0],stepX[0]);
//			iTTL::tensor<double, iTTL::stuple<1, 0, true, void, sYt> > tY(stplY, (double *)PyArray_DATA(Y));
//			iTTL::tensor<double, iTTL::stuple<2, 0, true, void, sA0t, sA1t> > tA(stplA, (double *)PyArray_DATA(A));
//			iTTL::tensor<double, iTTL::stuple<1, 0, true, void, sXt> > tX(stplX, (double *)PyArray_DATA(X));
//			BR_solve_one(tY, tA, tX, ridge);
    	}
    	catch (std::exception& ex)
    	{
        	PyErr_SetString(L1_Error, ex.what());
        	return NULL;
    	}
    }
    else
    {
    	PyErr_SetString(L1_Error, "Error. Incompatible sizes of arrays");
    	return NULL;
    }
    Py_RETURN_NONE;
}

// Array with methods

static PyMethodDef l1_methods[] =
{
// name from python, name in C-file, ..., __doc__ string of method
     {"BR_solve", BR_solve, METH_VARARGS, "L1 regression."},
     {NULL, NULL, 0, NULL}
};


// Init our module in Python

PyMODINIT_FUNC
initl1_procs(void)
{
    PyObject *m=Py_InitModule("l1_procs", l1_methods);
    static char errname[]="l1_procs.error";
    import_array(); // Don't forget this line!
    L1_Error = PyErr_NewException(errname, NULL, NULL);
    Py_INCREF(L1_Error); // Increment reference count for object
    PyModule_AddObject(m, "error", L1_Error);
}

