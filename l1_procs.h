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

template <typename T, typename STX, typename STU, typename STV>
class BR_solver
{
//	tpp::test_shape_length<STX, STU, STV>::test(std::tuple<STX, STU, STV>(x, U, v));
//  Test if U is a matrix
	static_assert(STU::snum==2,"Matrix of coefficients for BR_solve_one should have 2 dimensions (should be a matrix)");
// 	valence_parser<>::type is a sequence of 8 sequence of valence_data structures
//  Each of 8 sequences corresponds to a specific mask of valence
//  If valence is present in a tensor number N, the (1<<N) bit of the mask is set
//  Example: we have two tensors with shapes STA and STB. We need valences of only STB shape.
//  The result is in 2 sequence of valences if we write valence_parser<STA, STB>::type
//  The result is in 1 sequence of valences if we write valence_parser<STB, STA>::type
//  The tseq_element<> returns the element of a type_sequence
//  Each type sequence has static const element size
//  Non-empty type sequences has a head type defined
//  pos0, pos1, pos2 are position numbers of shape (dimension) of tensors included to the mask of the valence
	typedef typename tpp::valence_parser<false, STX, STU, STV>::type vd_by_mask;
	static_assert(tpp::tseq_element<1,vd_by_mask>::size==0,"x should not have free indices");
	static_assert(tpp::tseq_element<2,vd_by_mask>::size==0,"U should not have free indices");
	static_assert(tpp::tseq_element<3,vd_by_mask>::size==1,"x should be linked with U by one index");
	static_assert(tpp::tseq_element<4,vd_by_mask>::size==0,"v should not have free indices");
//	static_assert(tpp::tseq_element<5,vd_by_mask>::size<=1,"x and v can be linked by at most 1 index");
	static_assert(tpp::tseq_element<5,vd_by_mask>::size==0,"multiple solution is not supported");
	static_assert(tpp::tseq_element<6,vd_by_mask>::size==1,"v should be linked with U by one index");
	static_assert(tpp::tseq_element<7,vd_by_mask>::size==0,"Through indices are not supported yet");
	size_t xsize;
	size_t vsize;
	size_t msize;
	size_t Asize;
	tpp::MATRIX<> A;
//  Valences of our tensors are known. They are stored as static members in valence_data type
//  We create indices of known valences using usual C++ syntax
	tpp::segmentIndex<tpp::tseq_element<6,vd_by_mask>::head::v_type> vi;
	tpp::segmentIndex<tpp::tseq_element<3,vd_by_mask>::head::v_type> xi;
	tpp::defaultIndex<tpp::tseq_element<6,vd_by_mask>::head::v_type> I;
	tpp::defaultIndex<tpp::tseq_element<3,vd_by_mask>::head::v_type> J;
	tpp::VECTOR<> AP;
	tpp::VECTOR<> AQ;
	tpp::VECTOR<> S;
	tpp::shared_container<std::vector<size_t> > soZ;
	std::vector<T> ratio;
	std::vector<size_t> idx;
	tpp::shared_container<std::vector<size_t> > QC;
	std::vector<size_t>& q;
	tpp::VECTOR<T> hm;
	tpp::VECTOR<T> g;
	tpp::VECTOR<T> BS_coef;
	tpp::segmentIndex<tpp::tseq_element<3,vd_by_mask>::head::v_type> mi;
	void BR_iterations()
	{
//  Barrodale-Roberts iterations
//  with Bloomfield and Steiger modification

		auto Avq=A(vsize,J);
		tpp::container<std::vector<size_t> >:: template index<tpp::tseq_element<3,vd_by_mask>::head::v_type> Z(soZ,0);

		while(1)
		{
//			printf("Current matrix:\n");
//			for (size_t j=0;j<msize;j++)
//			{
//				for (size_t i=0;i<vsize;i++)
//					printf("%lf; ",(T&)A(i,j));
//				printf("%lf\n",(T&)A(vsize,j));
//			}

			soZ.cont.clear();

			S(mi).sign(A(vsize,mi));
			for (size_t j=0;j<msize;j++)
			{
				if (Avq(j)==0.0)
					soZ.cont.push_back(j);
			}

			T minL=1.0;
			size_t p=-1;
			g(vi).gem(A(vi,mi),S(mi));
			hm(I).asum(g(I));
			g(vi).asum(A(vi,Z));
			g(I).axpy(hm(I),-1);
			BS_coef(vi).asum(A(vi,mi));
			g(I).div(BS_coef(I));
			for (size_t i=0;i<vsize;i++)
			{
				T l=g(i);
				if (l<minL)
				{
					minL=l;
					p=i;
				}
			}
			if (minL>=0.0)
			{
				break;
			}


// Solving subproblem for row p and destination in the last row. Finding q.

			auto X=A(p,J);
			auto Y=A(vsize,J);
			idx.clear();
			for (size_t i=0;i<msize;i++)
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
			sort(idx.begin(), idx.end(), [&](size_t i1, size_t i2) {return ratio[i1] < ratio[i2];});
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
			q[p]=minq;

//		printf("zqty=%llu, p=%llu, q[p]=%llu\n",(unsigned long long)soZ.cont.size(),(unsigned long long)p,(unsigned long long)q[p]);

			T sc=1.0/A(p,q[p]);
			A(p,J).scal(sc);
			AP(J)=A(p,J);
			AQ(I)=A(I,q[p]);
			AQ(p)=0;
			A(I,J).gem(AP(J),AQ(I),-1,1);
			A(J,q[p])=0;
			A(p,q[p])=1;
		}
	}
public:
	BR_solver(const tpp::tensor<T, STX>& x, const tpp::tensor<T, STU>& U, tpp::tensor<T, STV>& v, T ridge=0.0):
	xsize(tpp::get<tpp::tseq_element<3,vd_by_mask>::head::pos0>(x).length()),
	vsize(tpp::get<tpp::tseq_element<6,vd_by_mask>::head::pos1>(v).length()),
	msize(ridge==0.0?xsize:xsize+vsize),
	Asize(msize+vsize),
	A({vsize+1, Asize}),
	vi(vsize,0),
	xi(xsize,0),
//	vi2(vsize,0),
	AP({Asize}),
	AQ({vsize+1}),
	S({msize}),
	soZ(msize),
	ratio(Asize),
	idx(Asize),
	QC(vsize),
	q(QC.cont),
	hm({vsize}),
	g({vsize}),
	BS_coef({vsize}),
	mi(msize,0)
	{
//  Barrodale-Roberts method. Finds v to optimize sum(abs( x-Uv )) + sum(abs(ridge*v)) --> min
//
//	ridge is an optional RIDGE-like non-negative coefficient.
//	The physical units of ridge is the same as the dimension of U(i,j).
//	So, if U(i,j) contains meters, ridge contains meters also.
//	If coefficients of U have different physical units, only 0 value for ridge makes sense
//	You should normalize U correctly if U coefficients are heterogeneous before using positive value for ridge
//  Test if lengths of correspondent indices are equal
		tpp::check_shape_length<vd_by_mask>(x, U, v);

//      Now we copy the matrix U and x to the beginning of A
		A(vi,xi)=U;
		A(vsize,xi)=x;
		if (ridge!=0.0)
		{
			A(J,vi+xsize)=0.0;
			A(vi,vi+xsize)=ridge;
		}
		A(J,vi+msize)=0.0;
		A(vi,vi+msize)=1.0;

		for (size_t i=0;i<vsize;i++)
			q[i]=i+msize;
		BR_iterations();
		v=A(vsize,vi+msize);
		v.scal(-1.0);

//  The problem is currently solved. However errors are accumulated due to iterations
//  Now we use the found solution completely defined by q as a starting point for the new short optimization cycle
//  The new optimization usually just checks the result, however if errors are large the result could be corrected

//		printf("Result matrix:\n");
//		for (size_t j=0;j<Asize;j++)
//		{
//			for (size_t i=0;i<vsize;i++)
//				printf("%lf; ",(T&)A(i,j));
//			printf("%lf\n",(T&)A(vsize,j));
//		}

		tpp::VECTOR<> CX({vsize});

		A(vi,xi)=U;
		A(vsize,xi)=x;
		if (ridge!=0.0)
		{
			A(J,vi+xsize)=0.0;
			A(vi,vi+xsize)=ridge;
		}
		A(J,vi+msize)=0.0;
		A(vi,vi+msize)=1.0;

		DECLARE_segmentIndex(vi2, vsize, 0);
		auto Qi=tpp::index_creator<tpp::container<std::vector<size_t>>::index>::create(vi2,QC,0);
		auto lu=A(vi,Qi).lu();
		if (lu.solve(A(vi,J)))
			return;

		CX(vi2)=A(vsize,Qi);
		A(vsize,J).gem(A(vi,J),CX(vi),-1.0,1.0);
		A(I, Qi)=0.0;
		A(vi2, Qi)=1.0;

	//	printf("\nCorrected matrix:\n");
	//	for (size_t j=0;j<Asize;j++)
	//	{
	//		for (size_t i=0;i<vsize;i++)
	//			printf("%lf; ",(T&)A(i,j));
	//		printf("%lf\n",(T&)A(vsize,j));
	//	}

		BR_iterations();

		v=A(vsize,vi+msize);
		v.scal(-1.0);
	}
};

template <typename T, typename STX, typename STU, typename STV>
void BR_solve_one(const tpp::tensor<T, STX>& x, const tpp::tensor<T, STU>& U, tpp::tensor<T, STV>&& v, T ridge=0.0)
{
//    Barrodale-Roberts method. Finds v to optimize sum[ sum(abs( x-Uv )) + abs(ridge*v) ] --> min
	BR_solver<T, STX, STU, STV>(x, U, v, ridge);
}


#endif /* L1_PROCS_H_ */
