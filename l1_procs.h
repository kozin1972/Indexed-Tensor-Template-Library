/*
 * l1_procs.h
 *
 *  Created on: 28 мая 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef L1_PROCS_H_
#define L1_PROCS_H_

#include <tensor.h>

template <typename T>
void BR_iterations(tpp::MATRIX<T>& A, size_t oq, size_t vq, std::vector<size_t>& idx, std::vector<T>& ratio, tpp::VECTOR<T>& AP, tpp::VECTOR<T>& AQ, tpp::VECTOR<T>& S, tpp::shared_container<std::vector<size_t> >& soZ, size_t *q)
{
//    Barrodale-Roberts iterations
//    with Bloomfield and Steiger modification

	DECLARE_defaultIndex(ni);
	DECLARE_defaultIndex(nj);
	DECLARE_segmentIndex(nv,vq,0);
	DECLARE_segmentIndex(no,oq,0);
//	T *Avq=&(T&)(A(vq,0));
	auto Avq=A(vq,ni);
	tpp::VECTOR<T> hm({vq});
	tpp::VECTOR<T> g({vq});
	tpp::VECTOR<T> BS_coef({vq});
	DECLARE_containerIndex(Z,soZ,0);

	while(1)
	{
//		printf("Current matrix:\n");
//		for (size_t j=0;j<oq;j++)
//		{
//			for (size_t i=0;i<vq;i++)
//				printf("%lf; ",(T&)A(i,j));
//			printf("%lf\n",(T&)A(vq,j));
//		}

		soZ.cont.clear();

		S(no).sign(A(vq,no));
		for (size_t j=0;j<oq;j++)
		{
			if (Avq(j)==0.0)
				soZ.cont.push_back(j);
		}

		T minL=1.0;
		size_t p=-1;
		g(nv).gem(A(nv,no),S(no));
		hm(nv).asum(g(nv));
		g(nv).asum(A(nv,Z));
		g(nv).axpy(hm(nv),-1);
		BS_coef(nv).asum(A(nv,no));
		g(nv).div(BS_coef(nv));
		for (size_t i=0;i<vq;i++)
		{
			T l=g(i);
			if (l<minL)
			{
				minL=l;
				p=i;
			}
		}
		if (minL>0.0)
		{
			break;
		}


// Solving subproblem for row p and destination in the last row. Finding q.

		auto X=A(p,ni);
		auto Y=A(vq,ni);
		idx.clear();
		for (size_t i=0;i<oq;i++)
		{
			if (tpp::abs((T&)X(i))>0.0)
			{
				ratio[i]=(T&)Y(i)/X(i);
				idx.push_back(i);
			}
		}
		size_t pqty=idx.size();
		if (pqty==0)
		{
			return;
		}
		sort(idx.begin(), idx.end(), [&ratio](size_t i1, size_t i2) {return ratio[i1] < ratio[i2];});
		size_t minq;
		if (pqty==1)
		{
			minq=idx[0];
		}
		else
		{
			size_t *pli=idx.data();
			size_t *pri=&(*idx.rbegin());
			T suml=tpp::abs((T&)X(*pli));
			T sumr=tpp::abs((T&)X(*pri));
			T *x=&(T&)X(0);
			while (pri>pli+1)
			{
				if (suml<sumr)
				{
					pli++;
					suml+=tpp::abs<T>(x[*pli]);
				}
				else
				{
					pri--;
					sumr+=tpp::abs<T>(x[*pri]);
				}
			}
			if (suml<sumr)
				pli++;
			else
				pri--;
			minq=*pli;
		}
		if (q[p]==minq)
			break;
		q[p]=minq;

//		printf("zqty=%llu, p=%llu, q[p]=%llu\n",(unsigned long long)soZ.cont.size(),(unsigned long long)p,(unsigned long long)q[p]);

		T sc=1.0/A(p,q[p]);
		A(p,ni).scal(sc);
		AP(ni)=A(p,ni);
		AQ(ni)=A(ni,q[p]);
		AQ(p)=0;
		A(ni,nj).gem(AP(nj),AQ(ni),-1,1);
		A(ni,q[p])=0;
		A(p,q[p])=1;
	}
}

template <typename T, typename STX, typename STU, typename STV>
void BR_solve_one(const tpp::tensor<T, STX>& x, const tpp::tensor<T, STU>& U, tpp::tensor<T, STV>&& v, T ridge=0.0)
{
//    Barrodale-Roberts method. Finds v to optimize sum(abs( x-Uv )) + sum(abs(ridge*v)) --> min
//
//	  ridge is an optional RIDGE-like non-negative coefficient.
//	  The physical units of ridge is the same as the dimension of U(i,j).
//	  So, if U(i,j) contains meters, ridge contains meters also.
//	  If coefficients of U have different physical units, only 0 value for ridge makes sense
//	  You should normalize U correctly if U coefficients are heterogeneous before using positive value for ridge

//  Test if lengths of correspondent indices are equal
	tpp::test_shape_length<STX, STU, STV>::test(std::tuple<STX, STU, STV>(x, U, v));
//  Test if U is a matrix
	static_assert(STU::snum==2,"Matrix of coefficients for BR_solve_one should have 2 dimensions (should be a matrix)");
// 	valence_parser<>::vd_by_mask is a sequence of 8 sequence of valence_data structures
//  Each of 8 sequences corresponds to a specific mask of valence
//  If valence is present in a tensor number N, the (1<<N) bit of the mask is set
//  Example: we have two tensors with shapes STA and STB. We need valences of only STB shape.
//  The result is in 2 sequence of valences if we write valence_parser<STA, STB>::vd_by_mask
//  The result is in 1 sequence of valences if we write valence_parser<STB, STA>::vd_by_mask
//  The tseq_element<>::type returns the element of a type_sequence
//  Each type sequence has static const element size
	typedef typename tpp::valence_parser<STX, STU, STV>::vd_by_mask vd_by_mask;
	static_assert(tpp::size(typename decltype(tpp::get<1>(vd_by_mask()))::type())==0,"x should not have free indices");
	static_assert(tpp::size(typename decltype(tpp::get<2>(vd_by_mask()))::type())==0,"U should not have free indices");
	static_assert(tpp::size(typename decltype(tpp::get<3>(vd_by_mask()))::type())==1,"x should be linked with U by one index");
	static_assert(tpp::size(typename decltype(tpp::get<4>(vd_by_mask()))::type())==0,"v should not have free indices");
//	static_assert(tpp::size(typename decltype(tpp::get<5>(vd_by_mask()))::type())<=1,"x and v can be linked by at most 1 index");
	static_assert(tpp::size(typename decltype(tpp::get<5>(vd_by_mask()))::type())==0,"multiple solution is not supported");
	static_assert(tpp::size(typename decltype(tpp::get<6>(vd_by_mask()))::type())==1,"v should be linked with U by one index");
	static_assert(tpp::size(typename decltype(tpp::get<7>(vd_by_mask()))::type())==0,"Through indices are not supported yet");
//  Non-empty type sequences has a head type defined
//  pos0, pos1, pos2 are position numbers of shape (dimension) of tensors included to the mask of the valence
	size_t xsize=tpp::get<decltype(tpp::head(tpp::type_pack_element<3,vd_by_mask>()))::type::pos0, STX>(x).length();
	size_t vsize=tpp::get<decltype(tpp::head(tpp::type_pack_element<6,vd_by_mask>()))::type::pos1, STV>(v).length();
	size_t Asize=xsize+vsize;
	size_t msize=xsize;
    if (ridge!=0.0)
    {
    	msize+=vsize;
    	Asize+=vsize;
    }

	tpp::MATRIX<> A({vsize+1, Asize});
//  Valences of our tensors are known. They are stored as static members in valence_data type
//  We create indices of known valences using usual C++ syntax
	tpp::segmentIndex<decltype(tpp::head(tpp::type_pack_element<6,vd_by_mask>()))::type::v_type> vi(vsize,0);
	tpp::segmentIndex<decltype(tpp::head(tpp::type_pack_element<3,vd_by_mask>()))::type::v_type> xi(xsize,0);
//  Now we can copy the matrix U to the beginning of A
	A(vi,xi)=U;
	A(vsize,xi)=x;
	DECLARE_segmentIndex(i,vsize,0);
	DECLARE_segmentIndex(j,vsize+1,0);
	if (ridge!=0.0)
	{
		A(j,i+xsize)=0.0;
		A(i,i+xsize)=ridge;
	}
	A(j,i+msize)=0.0;
	A(i,i+msize)=1.0;

	tpp::VECTOR<> AP({Asize});
	tpp::VECTOR<> AQ({vsize+1});
	tpp::VECTOR<> S({msize});
	tpp::shared_container<std::vector<size_t> > soZ(msize);
	std::vector<T> ratio(Asize);
	std::vector<size_t> idx(Asize);


	size_t *q=new size_t[vsize];
	for (size_t i=0;i<vsize;i++)
		q[i]=i+msize;
	BR_iterations(A,msize,vsize,idx,ratio,AP,AQ,S,soZ,q);

//  The problem is currently solved. However errors are accumulated due to iterations
//  Now we use the found solution completely defined by q as a starting point for the new short optimization cycle
//  The new optimization usually just checks the result, however if errors are large the result could be corrected

//	printf("Result matrix:\n");
//	for (size_t j=0;j<Asize;j++)
//	{
//		for (size_t i=0;i<vsize;i++)
//			printf("%lf; ",(T&)A(i,j));
//		printf("%lf\n",(T&)A(vsize,j));
//	}

	tpp::MATRIX<> C({vsize,vsize});
	tpp::VECTOR<> CX({vsize});
//  We create an alias of U where x-dimension precedes v-dimension. So, mU is an optionally transposed U
	auto mU=U.template order_indices<decltype(tpp::head(tpp::type_pack_element<3,vd_by_mask>()))::type::v_type, decltype(tpp::head(tpp::type_pack_element<6,vd_by_mask>()))::type::v_type>();
//	copy the identity matrix to the end of A
	A(j,i+msize)=0.0;
	A(i,i+msize)=1.0;
//	fills C
	for (size_t j=0;j<vsize;j++)
		if (q[j]<xsize)
			C(vi,j)=mU(q[j],vi);
		else
		{
			C(vi,j)=0.0;
			if (q[j]-msize==j)
				C(j,j)=1.0;
			else
				if (q[j]-xsize==j)
					C(j,j)=ridge;
				else
					C(j,j)=0.0;
		}
//	Calculate inverse of C at the end of A
	DECLARE_segmentIndex(vi2,vsize,0);
	DECLARE_segmentIndex(vi3,vsize,0);

	C(vi,vi2).gesv(A(vi3,vi2+msize));
//	Calculate the main part of A
	A(vi2,xi).gem(U,A(vi2,vi+msize));
//  Fill the rest of A
	if (ridge!=0.0)
	{
		A(vi2,i+xsize)=A(vi2,i+msize);
		A(vi2,i+xsize).scal(ridge);
		A(vsize,i+xsize)=0;
	}
	DECLARE_segmentIndex (vxi,Asize,0);
	for (size_t i=0;i<vsize;i++)
		if (q[i]<xsize)
			CX(i)=x(q[i]);
		else
			CX(i)=0.0;
	A(vsize,xi)=x;
	A(vsize,i+msize)=0;
	A(vsize,vxi).gem(A(vi,vxi),CX(vi),-1.0,1.0);
	for (size_t i=0;i<vsize;i++)
	{
		A(j,q[i])=0.0;
		A(i,q[i])=1.0;
	}

//	printf("\nRecalculated matrix:\n");
//	for (size_t j=0;j<Asize;j++)
//	{
//		for (size_t i=0;i<vsize;i++)
//			printf("%lf; ",(T&)A(i,j));
//		printf("%lf\n",(T&)A(vsize,j));
//	}

	BR_iterations(A,msize,vsize,idx,ratio,AP,AQ,S,soZ,q);

//	for (size_t j=0;j<vsize;j++)
//		if (q[j]<msize)
//			C(vi,j)=mU(q[j],vi);
//		else
//			C(vi,j)=0.0;

	v=A(vsize,vi+msize);
	v.scal(-1.0);
	delete[] q;
//	return C;
}


#endif /* L1_PROCS_H_ */
