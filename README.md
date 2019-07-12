
                                     Indexed Tensor Template Library

Indexed Tensor Template Library is driven by indices and based on BLAS. It works at compile time. Procedures made with Tensor Indexed Template Library can be used for Python extensions. See iTTL.docx for details.

Examples:

M(I,I)=1;			// Set diagonal elements of matrix M to 1

A(I,J)=B(I,K)*C(K,J);		// Matrix multiplication. A=BC

a(I)=b(I)*c(I);			// Element-wise vector multiplication

a(I).asum(M(I,J));		// Store sum of absolute values of each row of M in a

tr=M(I,I).sum();		// Trace of matrix



	Indexed Tensor Template Library is designed for C++ developers who implements heavy optimizations 
including statistical models, data mining, big data analysis. 
Power of BLAS libraries gives a possibility to process high-dimension data. 
Power of C++ gives a possibility to create fast non-trivial algorithms. 
Power of Indexed Tensor Template Library makes code readable and decreases development time.

	Indexed Tensor Template Library works at compile time. 
It chooses a proper BLAS subroutine and joins loops if possible. 
Suppose you have a matrix M of 1000000 rows and 2 columns. 
Assignment M=0 is translated by Tensor Template Library into two steps. 
In the first step the total size of M is calculated. The second step is a call to dcopy subroutine of BLAS.
