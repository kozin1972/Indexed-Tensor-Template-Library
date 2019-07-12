/*
 * gem_sample.cpp
 *
 *  Created on: 20 мая 2019 г.
 *      Author: Alexey Kozin
 *
 *  The provided code illustrates the logic of gem procedure.
 *  gem is General Multiplication for iTTL::TENSOR objects
 *
 *  Possible output:
 *
 *  Starting normal gem...
 *  ...Done 7ms
 *  Starting nested loops...
 *  ...Done 5242ms
 *  Verification of normal gem...
 *  ...Done
 *  Starting unoptimized gem...
 *  ...Done 137ms
 *  Verification of unoptimized gem...
 *  ...Done
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#include <chrono>
#include <tensor.h>

void gem_logic()
{
//  Define dimensions.
//	Dimensions are constants because we are going to show the gem logic with a help of traditional C-style arrays
	const size_t n_size=2;//2
	const size_t m_size=150; //150
	const size_t p_size=5; //5
	const size_t q_size=6; //6
	const size_t i_size=12; //12
	const size_t j_size=13; //13
	const size_t k_size=200; //200
	const size_t x_size=3;//3
	const size_t y_size=4;//4
//	Define TENSOR indices. All of them are default. It is for clear understanding.
	DECLARE_defaultIndex(I);
	DECLARE_defaultIndex(J);
	DECLARE_defaultIndex(K);
	DECLARE_defaultIndex(M);
	DECLARE_defaultIndex(N);
	DECLARE_defaultIndex(P);
	DECLARE_defaultIndex(Q);
	DECLARE_defaultIndex(X);
	DECLARE_defaultIndex(Y);
//	Allocate memory for tensors
	iTTL::TENSOR<7> R({n_size,p_size,q_size,n_size,m_size,i_size,j_size});
	iTTL::TENSOR<4> A({n_size,x_size,k_size,m_size});
	iTTL::TENSOR<5> B({y_size,n_size,k_size,i_size,j_size});
//	Allocate memory for C-style copy of tensors
	double (*r)[p_size][q_size][n_size][m_size][i_size][j_size]=new double[n_size][p_size][q_size][n_size][m_size][i_size][j_size];
	double (*a)[x_size][k_size][m_size]=new double[n_size][x_size][k_size][m_size];
	double (*b)[n_size][k_size][i_size][j_size]=new double[y_size][n_size][k_size][i_size][j_size];
//	Initialize data with the same numbers
	double *pd=(double *)r;
	double *Pd=R.data_ptr();
	size_t r_size=R.size();
	size_t a_size=A.size();
	size_t b_size=B.size();
	for (size_t i=0;i<r_size;i++)
	{
		pd[i]=i;
		Pd[i]=i;
	}
	pd=(double *)a;
	Pd=A.data_ptr();
	for (size_t i=0;i<a_size;i++)
	{
		pd[i]=i;
		Pd[i]=i;
	}
	pd=(double *)b;
	Pd=B.data_ptr();
	for (size_t i=0;i<b_size;i++)
	{
		pd[i]=i;
		Pd[i]=i;
	}
//	Two numeric parameters for gem
	double alpha=13.0; 	// default value is 1.0
	double beta=17.0;	// default value is 0.0
//	The gem call being explained ( R = alpha*A*B + beta*R )
	printf("Starting normal gem...\n");
	std::chrono::steady_clock::time_point tp0=std::chrono::steady_clock::now();
	R(N,P,Q,N,M,I,J).gem(A(N,X,K,M),B(Y,N,K,I,J),alpha,beta);
	std::chrono::milliseconds dt0=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-tp0);
	printf("...Done %ldms\nStarting nested loops...\n",dt0.count());
	std::chrono::steady_clock::time_point tp1=std::chrono::steady_clock::now();
// 	The following nested loops do logically the same as the gem call above
	for (size_t n=0;n<n_size;n++)
		for (size_t p=0;p<p_size;p++)
			for (size_t q=0;q<q_size;q++)
				for (size_t m=0;m<m_size;m++)
					for (size_t i=0;i<i_size;i++)
						for (size_t j=0;j<j_size;j++)
						{
							r[n][p][q][n][m][i][j]*=beta;
//  Indices X,Y,K are absent in the result R(N,P,Q,N,M,I,J)
							for (size_t x=0;x<x_size;x++)
								for (size_t y=0;y<y_size;y++)
									for (size_t k=0;k<k_size;k++)
										r[n][p][q][n][m][i][j]+=a[n][x][k][m]*b[y][n][k][i][j]*alpha;
						}
// Verification
	std::chrono::milliseconds dt1=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-tp1);
	printf("...Done %ldms\nVerification of normal gem...\n",dt1.count());
	pd=(double *)r;
	Pd=R.data_ptr();
	int errors=0;
	for (size_t i=0;i<r_size;i++)
	{
		if (pd[i]!=Pd[i])
			errors++;
	}
	if (errors)
		printf("gem_logic leads to %d errors\n",errors);
	printf("...Done\n");
	Pd=R.data_ptr();
	for (size_t i=0;i<r_size;i++)
	{
		Pd[i]=i;
	}

	DECLARE_simpleIndex(FI);
	DECLARE_simpleIndex(FJ);
	DECLARE_simpleIndex(FK);
	DECLARE_simpleIndex(FM);
	DECLARE_simpleIndex(FN);
	DECLARE_simpleIndex(FP);
	DECLARE_simpleIndex(FQ);
	DECLARE_simpleIndex(FX);
	DECLARE_simpleIndex(FY);
	printf("Starting unoptimized gem...\n");
	std::chrono::steady_clock::time_point tp2=std::chrono::steady_clock::now();
	R(FN,FP,FQ,FN,FM,FI,FJ).gem(A(FN,FX,FK,FM),B(FY,FN,FK,FI,FJ),alpha,beta);
	std::chrono::milliseconds dt2=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-tp2);
	printf("...Done %ldms\nVerification of unoptimized gem...\n",dt2.count());
// Verification of simpleIndex
	pd=(double *)r;
	Pd=R.data_ptr();
	errors=0;
	for (size_t i=0;i<r_size;i++)
	{
		if (pd[i]!=Pd[i])
			errors++;
	}
	if (errors)
		printf("gem_logic on simpleIndex leads to %d errors\n",errors);
	printf("...Done\n");

//	Free C-style arrays
	delete[] b;
	delete[] a;
	delete[] r;
}


