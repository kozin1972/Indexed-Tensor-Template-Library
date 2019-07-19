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


template <typename TTY, typename TTA, typename TTX>
class BR_solver;

//template <typename TY, typename STY, typename TA, typename STA, typename TX, typename STX>
//class BR_solver<iTTL::tensor<TY, STY>, iTTL::tensor<TA, STA>, iTTL::tensor<TX, STX> >
//{
//	static_assert(typeid(TY)==typeid(TA),"BR_solver Y and A parameters have different element types");
//	static_assert(typeid(TX)==typeid(TA),"BR_solver X and A parameters have different element types");
//};
//
//template <typename T, typename STY, typename STA, typename STX>
template <typename TY, typename STY, typename TA, typename STA, typename TX, typename STX>
class BR_solver<iTTL::tensor<TY, STY>, iTTL::tensor<TA, STA>, iTTL::tensor<TX, STX> >
{
	typedef typename iTTL::large_type<TY, TA, TX>::type T;
//	iTTL::test_shape_length<STY, STA, STX>::test(std::tuple<STY, STA, STX>(y, A, x));
//  Test if A is a matrix
	static_assert(STA::snum==2,"Matrix of coefficients for BR_solve_one should have 2 dimensions (should be a matrix)");
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
	typedef typename iTTL::valence_parser<false, STY, STA, STX>::type vd_by_mask;
	static_assert(iTTL::tseq_element<1,vd_by_mask>::size==0,"y should not have free indices");
	static_assert(iTTL::tseq_element<2,vd_by_mask>::size==0,"A should not have free indices");
	static_assert(iTTL::tseq_element<3,vd_by_mask>::size==1,"y should be linked with A by one index");
	static_assert(iTTL::tseq_element<4,vd_by_mask>::size==0,"x should not have free indices");
//	static_assert(iTTL::tseq_element<5,vd_by_mask>::size<=1,"y and x can be linked by at most 1 index");
	static_assert(iTTL::tseq_element<5,vd_by_mask>::size==0,"multiple solution is not supported");
	static_assert(iTTL::tseq_element<6,vd_by_mask>::size==1,"x should be linked with A by one index");
	static_assert(iTTL::tseq_element<7,vd_by_mask>::size==0,"Through indices are not supported yet");

	size_t ysize;
	size_t xsize;
	size_t msize;
	size_t Wsize;
	iTTL::MATRIX<T> W;
//  Valences of our tensors are known. They are stored as static members in valence_data type
//  We create indices of known valences using usual C++ syntax
	iTTL::segmentIndex<iTTL::tseq_element<6,vd_by_mask>::head::v_type> xi;
	iTTL::segmentIndex<iTTL::tseq_element<3,vd_by_mask>::head::v_type> yi;
	iTTL::defaultIndex<iTTL::tseq_element<6,vd_by_mask>::head::v_type> I;
	iTTL::defaultIndex<iTTL::tseq_element<3,vd_by_mask>::head::v_type> J;
	iTTL::VECTOR<T> WP;
	iTTL::VECTOR<T> WQ;
	iTTL::VECTOR<T> S;
	iTTL::shared_container<std::vector<size_t> > soZ;
	std::vector<T> ratio;
	std::vector<size_t> idx;
	iTTL::shared_container<std::vector<size_t> > QC;
	std::vector<size_t>& q;
	iTTL::VECTOR<T> hm;
	iTTL::VECTOR<T> g;
	iTTL::VECTOR<T> BS_coef;
	iTTL::segmentIndex<iTTL::tseq_element<3,vd_by_mask>::head::v_type> mi;
	void BR_iterations()
	{
//  Barrodale-Roberts iterations
//  with Bloomfield and Steiger modification

		auto Wxq=W(xsize,J);
		iTTL::container<std::vector<size_t> >:: template index<iTTL::tseq_element<3,vd_by_mask>::head::v_type> Z(soZ,0);

		while(1)
		{
//			printf("Current matrix:\n");
//			for (size_t j=0;j<msize;j++)
//			{
//				for (size_t i=0;i<xsize;i++)
//					printf("%lf; ",(T&)W(i,j));
//				printf("%lf\n",(T&)W(xsize,j));
//			}

			soZ.cont.clear();

			S(mi).sign(W(xsize,mi));
			for (size_t j=0;j<msize;j++)
			{
				if (Wxq(j)==0.0)
					soZ.cont.push_back(j);
			}

			T minL=1.0;
			size_t p=-1;
			g(xi).gem(W(xi,mi),S(mi));
			hm(I).asum(g(I));
			g(xi).asum(W(xi,Z));
			g(I).axpy(hm(I),-1);
			BS_coef(xi).asum(W(xi,mi));
			g(I).div(BS_coef(I));
			for (size_t i=0;i<xsize;i++)
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

			auto Y=W(p,J);
			auto D=W(xsize,J);
			idx.clear();
			for (size_t i=0;i<msize;i++)
			{
				if (iTTL::abs((T&)Y(i))>0.0)
				{
					ratio[i]=(T&)D(i)/Y(i);
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
				T suml=iTTL::abs((T&)Y(*pli));
				T sumr=iTTL::abs((T&)Y(*pri));
				T *y=&(T&)Y(0);
				while (pri>pli+1)
				{
					if (suml<sumr)
					{
						pli++;
						suml+=iTTL::abs<T>(y[*pli]);
					}
					else
					{
						pri--;
						sumr+=iTTL::abs<T>(y[*pri]);
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

			T sc=1.0/W(p,q[p]);
			W(p,J).scal(sc);
			WP(J)=W(p,J);
			WQ(I)=W(I,q[p]);
			WQ(p)=0;
			W(I,J).gem(WP(J),WQ(I),-1,1);
			W(J,q[p])=0;
			W(p,q[p])=1;
		}
	}
public:
	BR_solver(const iTTL::tensor<TY, STY>& y, const iTTL::tensor<TA, STA>& A, iTTL::tensor<TX, STX>& x, T ridge=0.0):
	ysize(iTTL::get<iTTL::tseq_element<3,vd_by_mask>::head::pos0>(y).length()),
	xsize(iTTL::get<iTTL::tseq_element<6,vd_by_mask>::head::pos1>(x).length()),
	msize(ridge==0.0?ysize:ysize+xsize),
	Wsize(msize+xsize),
	W({xsize+1, Wsize}),
	xi(xsize,0),
	yi(ysize,0),
//	vi2(xsize,0),
	WP({Wsize}),
	WQ({xsize+1}),
	S({msize}),
	soZ(msize),
	ratio(Wsize),
	idx(Wsize),
	QC(xsize),
	q(QC.cont),
	hm({xsize}),
	g({xsize}),
	BS_coef({xsize}),
	mi(msize,0)
	{
//  Barrodale-Roberts method. Finds x to optimize sum(abs( y-Ax )) + sum(abs(ridge*x)) --> min
//
//	ridge is an optional RIDGE-like non-negative coefficient.
//	The physical units of ridge is the same as the dimension of A(i,j).
//	So, if A(i,j) contains meters, ridge contains meters also.
//	If coefficients of A have different physical units, only 0 value for ridge makes sense
//	You should normalize A correctly if A coefficients are heterogeneous before using positive value for ridge
//  Test if lengths of correspondent indices are equal
//		static_assert(std::is_same<T,TA>::value,"Type of elements are different for Y and A");
//		static_assert(std::is_same<T,TX>::value,"Type of elements are different for Y and X");
//		int u0=(int)iTTL::check_shape<typename STA::template element_type<0>>::usage;
//		int u1=(int)iTTL::check_shape<typename STA::template element_type<1>>::usage;
//		printf("A-USAGE[0]=%d, A-USAGE[1]=%d\n", u0, u1);

		iTTL::check_shape_length<vd_by_mask>(y, A, x);

//      Now we copy the matrix A and y to the beginning of W
		W(xi,yi)=A;
		W(xsize,yi)=y;
		if (ridge!=0.0)
		{
			W(J,xi+ysize)=0.0;
			W(xi,xi+ysize)=ridge;
		}
		W(J,xi+msize)=0.0;
		W(xi,xi+msize)=1.0;

		for (size_t i=0;i<xsize;i++)
			q[i]=i+msize;
		BR_iterations();
		x=W(xsize,xi+msize);
		x.scal(-1.0);

//  The problem is currently solved. However errors are accumulated due to iterations
//  Now we use the found solution completely defined by q as a starting point for the new short optimization cycle
//  The new optimization usually just checks the result, however if errors are large the result could be corrected

//		printf("Result matrix:\n");
//		for (size_t j=0;j<Wsize;j++)
//		{
//			for (size_t i=0;i<xsize;i++)
//				printf("%lf; ",(T&)W(i,j));
//			printf("%lf\n",(T&)W(xsize,j));
//		}

		iTTL::VECTOR<T> CY({xsize});

		W(xi,yi)=A;
		W(xsize,yi)=y;
		if (ridge!=0.0)
		{
			W(J,xi+ysize)=0.0;
			W(xi,xi+ysize)=ridge;
		}
		W(J,xi+msize)=0.0;
		W(xi,xi+msize)=1.0;

		DECLARE_segmentIndex(xi2, xsize, 0);
		auto Qi=iTTL::index_creator<iTTL::container<std::vector<size_t>>::index>::create(xi2,QC,0);
		auto lu=W(xi,Qi).lu();
		if (lu.solve(W(xi,J)))
			return;

		CY(xi2)=W(xsize,Qi);
		W(xsize,J).gem(W(xi,J),CY(xi),-1.0,1.0);
		W(I, Qi)=0.0;
		W(xi2, Qi)=1.0;

	//	printf("\nCorrected matrix:\n");
	//	for (size_t j=0;j<Wsize;j++)
	//	{
	//		for (size_t i=0;i<xsize;i++)
	//			printf("%lf; ",(T&)W(i,j));
	//		printf("%lf\n",(T&)W(xsize,j));
	//	}

		BR_iterations();

		x=W(xsize,xi+msize);
		x.scal(-1.0);
	}
};

//template <typename T, typename STY, typename STA, typename STX>
template <typename TY, typename STY, typename TA, typename STA, typename TX, typename STX>
void BR_solve_one(const iTTL::tensor<TY, STY>& y, const iTTL::tensor<TA, STA>& A, iTTL::tensor<TX, STX>&& x, typename iTTL::large_type<TY, TA, TX>::type ridge=0.0)
{
//    Barrodale-Roberts method. Finds x to optimize sum[ sum(abs( y-Ax )) + abs(ridge*x) ] --> min
	BR_solver<iTTL::tensor<TY, STY>, iTTL::tensor<TA, STA>, iTTL::tensor<TX, STX> >(y, A, x, ridge);
}

//template <typename T, typename STY, typename STA, typename STX>
template <typename TY, typename STY, typename TA, typename STA, typename TX, typename STX>
void BR_solve_one(const iTTL::tensor<TY, STY>& y, const iTTL::tensor<TA, STA>& A, iTTL::tensor<TX, STX>& x, typename iTTL::large_type<TY, TA, TX>::type ridge=0.0)
{
//    Barrodale-Roberts method. Finds x to optimize sum[ sum(abs( y-Ax )) + abs(ridge*x) ] --> min
	BR_solver<iTTL::tensor<TY, STY>, iTTL::tensor<TA, STA>, iTTL::tensor<TX, STX> >(y, A, x, ridge);
}

#endif /* L1_PROCS_H_ */
