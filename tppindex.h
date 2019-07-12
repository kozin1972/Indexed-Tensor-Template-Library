/*
 * tppindex.h
 *
 *  Created on: 7 апр. 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef TPPINDEX_H_
#define TPPINDEX_H_

#include <cstddef>
#include <memory>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <sstream>

namespace iTTL
{
	class exception: public std::exception
	{
	protected:
		std::string w;
		exception(const char *s)
		{
			w=s;
		}
		exception(const std::string& s)
		{
			w=s;
		}
	public:
		virtual const char *what() const noexcept
		{
			return w.c_str();
		}
	};
	class outOfBounds: public exception
	{
		char *indexValue(unsigned long long value, unsigned long long i, unsigned long long length) noexcept
		{
			static char s[200];
			if (length)
				sprintf(s,"Out of bounds error. Incorrect value %llu of index #%llu is applied. Correct index values for the underlying shape are [0..%llu]",value,i,length-1);
			else
				sprintf(s,"Out of bounds error. Incorrect value %llu of index #%llu is applied. The underlying shape is not initialized and cannot be used",value,i);
			return s;
		}
		char *newIndex(const char *min_max, unsigned long long i, unsigned long long sample, unsigned long long length) noexcept
		{
			static char s[200];
			if (length)
				sprintf(s,"Out of bounds error. The %s possible value of index #%llu is %llu. Correct index values for the underlying shape are [0..%llu]",min_max,i,sample,length-1);
			else
				sprintf(s,"Out of bounds error. The %s possible value of index #%llu is %llu. The underlying shape is not initialized and cannot be used",min_max,i,sample);
			return s;
		}
	public:
		outOfBounds(): exception("Attempt to access non-allocated tensor data") {}
		outOfBounds(unsigned long long value, unsigned long long i, unsigned long long length): exception(indexValue(value,i,length)) {}
		outOfBounds(const char *min_max, unsigned long long i, unsigned long long sample, unsigned long long length): exception(newIndex(min_max,i,sample,length)) {}
	};
	class incompatibleIndices: public exception
	{
		char *incompatible(unsigned long long i0, unsigned long long i1, int v_type, unsigned long long length0, unsigned long long length1) noexcept
		{
			static char s[200];
			sprintf(s,"Incompatible indices error. Indices #%llu and #%llu have the same valence type #%d and different lengths: %llu and %llu",i0,i1,v_type,length0,length1);
			return s;
		}
		char *incompatible(unsigned long long p0, int s0, unsigned long long p1, int s1, int v_type, unsigned long long length0, unsigned long long length1) noexcept
		{
			static char s[200];
			sprintf(s,"Incompatible dimensions error. Dimension #%llu of tensor #%d and #%llu of tensor #%d have the same valence type #%d and different lengths: %llu and %llu",p0,s0,p1,s1,v_type,length0,length1);
			return s;
		}
	public:
		incompatibleIndices(unsigned long long i0, unsigned long long i1, int v_type, unsigned long long length0, unsigned long long length1): exception(incompatible(i0,i1,v_type,length0,length1)) {}
		incompatibleIndices(unsigned long long p0, int s0, unsigned long long p1, int s1, int v_type, unsigned long long length0, unsigned long long length1): exception(incompatible(p0,s0,p1,s1,v_type,length0,length1)) {}
	};
	class reshapeException: public exception
	{
		template <size_t ODIM, size_t NDIM>
		std::string reshapeEx(const size_t (&oaw)[ODIM], const size_t (&naw)[NDIM]) noexcept
		{
			size_t tot;
			std::stringstream res;
			res<<"Reshape error. Total size of desired shape is greater then existing one. Initial shape: [";
			if (ODIM>0)
			{
				res<<oaw[0];
				tot=oaw[0];
				for (size_t i=1;i<ODIM;i++)
				{
					res<<"x";
					res<<oaw[i];
					tot*=oaw[i];
				}
				res<<"]=";
				res<<tot;
			}
			else
				res<<"]=1";
			res<<" < desired shape: [";

			if (NDIM>0)
			{
				res<<naw[0];
				tot=naw[0];
				for (size_t i=1;i<NDIM;i++)
				{
					res<<"x";
					res<<naw[i];
					tot*=naw[i];
				}
				res<<"]=";
				res<<tot;
			}
			else
				res<<"]=1";
			return res.str();
		}

	public:
		template <size_t ODIM, size_t NDIM>
		reshapeException(const size_t (&oaw)[ODIM], const size_t (&naw)[NDIM]): exception(reshapeEx(oaw,naw)) {}
	};
	class MatrixIsNotSquareException: public exception
	{
		char *message(size_t (&shape)[2])
		{
			static char s[200];
			sprintf(s,"Dimension size of matrix are different: [%llu x %llu]. Matrix should be square",(unsigned long long)shape[0],(unsigned long long)shape[1]);
			return s;
		}
	public:
		MatrixIsNotSquareException(size_t (&shape)[2]): exception(message(shape)) {}
	};
	class DifferentTypesOfElements: public exception
	{
		char *DifferentTypesOfElements_mes(const char *method, const char *t0, const char *t1)
		{
			static char s[100];
			sprintf(s,"Error. Different element types '%s' and '%s' are passed to '%s'\n", t0, t1, method);
			return s;
		}
	public:
		DifferentTypesOfElements(const char *method, const char *t0, const char *t1):exception(DifferentTypesOfElements_mes(method, t0, t1)) {}
	};
	template <typename T>
	class shared_array
	{
		struct lock_info
		{
			size_t shared_qty;
		};
		static const size_t extra_element_qty=(sizeof(lock_info)+sizeof(T)-1)/sizeof(T);
	public:
		inline static T *alloc(size_t size)
		{
			if (size==0)
				return NULL;
			T *p=new T[size+extra_element_qty];
			((lock_info *)p)->shared_qty=1;
			return p+extra_element_qty;
		}
		inline static void share(const T *p) noexcept
		{
			if (!p)
				return;
			p-=extra_element_qty;
			lock_info *pli=(lock_info *)p;
			pli->shared_qty++;
		}
		inline static void free(const T *p)
		{
			if (!p)
				return;
			p-=extra_element_qty;
			lock_info *pad=(lock_info *)p;
			pad->shared_qty--;
			if (!pad->shared_qty)
				delete[] p;
		}
	};

	template <typename T, typename _ = void>
	class shared_container
	{
		struct lock_info: public T
		{
			size_t shared_qty;
			template <typename ... Ts>
			lock_info(Ts...args): T(args...), shared_qty(1) {}
			static lock_info& unobj(T& object) noexcept
			{
				return *(static_cast<lock_info*>(&object));
			}
		};
	public:
		T& cont;
		size_t shared_qty() const noexcept { return lock_info::unobj(cont).shared_qty; }
		typedef T base_container_type;
		template <typename...Ts>
		shared_container(Ts...args): cont(*(new lock_info(args...))) {}
		shared_container(std::initializer_list<typename T::value_type> il, const typename T::allocator_type& __a = typename T::allocator_type()): cont(*(new lock_info(il))) {}
		shared_container(const shared_container& src) noexcept :cont(src.cont)
		{
			lock_info& li=lock_info::unobj(src.cont);
			li.shared_qty++;
		}
		~shared_container()
		{
			lock_info& li=lock_info::unobj(cont);
			if ((--li.shared_qty)==0)
				delete &li;
		}
		operator T&() { return cont; }
	};

	enum dSU
	{
		USAGE_FULL=0,
		USAGE_CONT=1,
		USAGE_ENUM=2,
		USAGE_SLAVE=3
	};

	template <int V_TYPE>
	class defaultIndex;

	template <int V_TYPE>
	class segmentIndex;

	template <int V_TYPE>
	class forwardIndex;

	template <int V_TYPE>
	class reverseIndex;

	template <int V_TYPE>
	class simpleIndex;


	template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	class segment;

	struct segment_data
	{
		size_t l;
		size_t s;
	};
	template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO>
	class segment<V_TYPE,USAGE,USE_ALSO,0>: protected segment_data
	{
		template <int, enum dSU, size_t, size_t> friend struct segment;
		template <typename, typename> friend struct segment_gluer;
	public:
		inline segment() = default;
		template <typename PARENT_SHAPE>
		inline segment(const defaultIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=p.length(); this->s=p.step(); }
		template <typename PARENT_SHAPE>
		inline segment(const segmentIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=i.length(); this->s=p.step(); }
		template <int O_V_TYPE, enum dSU O_USAGE>
		inline segment(const segment<O_V_TYPE, O_USAGE, USE_ALSO, 0>& src) noexcept { this->l=src.length(); this->s=src.step(); }
		inline segment(size_t l, size_t s) { this->l=l; this->s=s; }
		inline size_t length() const noexcept { return l; }
		inline size_t step() const noexcept { return s; }
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) noexcept
		{
			p+=o*s;
			return o;
		}
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) const noexcept
		{
			p+=o*s;
			return o;
		}
		template <typename T>
		class iterator
		{
			size_t s;
			size_t c;
			size_t end;
		public:
			inline iterator(const segment& i, const T * const p) noexcept: s(i.s), c(0), end(i.l) {}
			inline void move_one(T *&p) noexcept { p+=s; c++; }
			inline void move_one(const T *&p) noexcept { p+=s; c++; }
			inline void move(T *&p, ptrdiff_t shift) noexcept { p+=s*shift; c+=shift; }
			inline void move(const T *&p, ptrdiff_t shift) noexcept { p+=s*shift; c+=shift; }
			inline bool not_end(const T *p) const noexcept { return c<end; }
		};
	};

	template <int V_TYPE, enum dSU USAGE>
	class segment<V_TYPE,USAGE,0,0>: protected segment_data
	{
		template <int, enum dSU, size_t, size_t> friend struct segment;
		template <typename, typename> friend struct segment_gluer;
	public:
		inline void init(size_t l,size_t s) noexcept { this->l=l; this->s=s; }
		inline segment() noexcept { this->l=0; this->s=0; };
		template <typename PARENT_SHAPE>
		inline segment(const defaultIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=p.length(); this->s=p.step(); }
		template <typename PARENT_SHAPE>
		inline segment(const segmentIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=i.length(); this->s=p.step(); }
		template <int O_V_TYPE, enum dSU O_USAGE>
		inline segment(const segment<O_V_TYPE, O_USAGE, 0, 0>& src) noexcept { this->l=src.length(); this->s=src.step(); }
		inline segment(size_t l, size_t s) { this->l=l; this->s=s; }
		inline size_t length() const noexcept { return l; }
		inline size_t step() const noexcept { return s; }
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) noexcept
		{
			p+=o*s;
			return o;
		}
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) const noexcept
		{
			p+=o*s;
			return o;
		}
		template <typename T>
		class iterator
		{
			size_t s;
			const T * const end;
		public:
			inline iterator(const segment& i, const T * const p) noexcept: s(i.s), end(p+i.l*i.s) {}
			inline void move_one(T *&p) noexcept { p+=s; }
			inline void move_one(const T *&p) noexcept { p+=s; }
			inline void move(T *&p, ptrdiff_t shift) noexcept { p+=s*shift; }
			inline void move(const T *&p, ptrdiff_t shift) noexcept { p+=s*shift; }
			inline bool not_end(const T *p) const noexcept { return p<end; }
		};
		template <typename T>
		class cont_iterator
		{
			const T * const end;
		public:
			inline cont_iterator(const segment& i, const T * const p) noexcept: end(p+i.l) {}
			inline void move_one(T *&p) noexcept { p++; }
			inline void move_one(const T *&p) noexcept { p++; }
			inline void move(T *&p, ptrdiff_t shift) noexcept { p+=shift; }
			inline void move(const T *&p, ptrdiff_t shift) noexcept { p+=shift; }
			inline bool not_end(const T *p) const noexcept { return p<end; }
		};
	};

	struct general_segment_data
	{
		size_t l;
	};
	template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	class segment: protected general_segment_data
	{
	public:
		inline segment() = default;
		template <typename PARENT_SHAPE>
		inline segment(const defaultIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=p.length(); }
		template <typename PARENT_SHAPE>
		inline segment(const segmentIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=i.length(); }
		template <int OLD_V_TYPE>
		inline segment(const segment<OLD_V_TYPE, USAGE, USE_ALSO, PARENT>& src) noexcept { this->l=src.length(); }
		size_t length() const noexcept { return this->l; }
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) noexcept
		{
			return o;
		}
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) const noexcept
		{
			return o;
		}
		template <typename T>
		class iterator
		{
			size_t c;
			size_t end;
		public:
			inline iterator(const segment& i, const T * const p) noexcept: c(0), end(i.l)
			{
			}
			inline ptrdiff_t move_one(T *&p) noexcept { c++; return 1; }
			inline ptrdiff_t move_one(const T *&p) noexcept { c++; return 1; }
			inline ptrdiff_t move(T *&p, ptrdiff_t shift) noexcept { c+=shift; return shift; }
			inline ptrdiff_t move(const T *&p, ptrdiff_t shift) noexcept { c+=shift; return shift; }
			inline bool not_end(const T *p) const noexcept { return c<end; }
		};
	};

	template <int V_TYPE>
	class segmentIndex
	{
		size_t l;
		size_t o;
	public:
		inline void check_bounds(size_t inum, size_t len) const
		{
			if (l+o>len)
				throw outOfBounds("maximum", inum, l+o-1, len);
			if (o>len)
				throw outOfBounds("minimum", inum, o-1, len);
		}
		inline segmentIndex(size_t l, size_t o) noexcept:l(l),o(o){}
		inline segmentIndex operator+(ptrdiff_t offset) const noexcept { return segmentIndex(this->l,o+offset); }
		inline segmentIndex operator-(ptrdiff_t offset) const noexcept { return segmentIndex(l,o-offset); }
		inline reverseIndex<V_TYPE> operator-() const noexcept { return reverseIndex<V_TYPE>(l,o); }
		inline size_t min_element() const noexcept { return o; }
		inline size_t first() const noexcept { return o; }
		inline size_t max_element() const noexcept { return o+l-1; }
		inline size_t length() const noexcept { return l; }

		template <typename PARENT_SHAPE, size_t P_POS>
		struct applyer;

		template <template <int, enum dSU, size_t, size_t> class shapeClass, int P_V_TYPE, enum dSU P_USAGE, size_t P_USE_ALSO, size_t P_PARENT, size_t P_POS>
		struct applyer<shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>, P_POS>
		{
			using parent_shape=shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>;
			using new_shape = segment<V_TYPE, (P_USAGE<USAGE_CONT?USAGE_CONT:P_USAGE),0,P_POS+1>;
		};

		template <int P_V_TYPE, enum dSU P_USAGE, size_t P_POS>
		struct applyer<segment<P_V_TYPE, P_USAGE, 0, 0>, P_POS>
		{
			using parent_shape=segment<P_V_TYPE, P_USAGE, 0, 0>;
			using new_shape = segment<V_TYPE, (P_USAGE<USAGE_CONT?USAGE_CONT:P_USAGE) ,0, 0>;
		};
	};
#define DECLARE_segmentIndex(I,length,offset) iTTL::segmentIndex<__COUNTER__> I(length,offset);

	template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	class forward;

	template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO>
	class forward<V_TYPE,USAGE,USE_ALSO,0>: protected segment_data
	{
		template <int, enum dSU, size_t, size_t> friend struct forward;
	public:
		inline forward() = default;
		template <typename PARENT_SHAPE>
		inline forward(const defaultIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=p.length(); this->s=p.step(); }
		template <typename PARENT_SHAPE>
		inline forward(const simpleIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=p.length(); this->s=p.step(); }
		template <typename PARENT_SHAPE>
		inline forward(const forwardIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=i.length(); this->s=p.step(); }
		template <int O_V_TYPE, enum dSU O_USAGE>
		inline forward(const forward<O_V_TYPE, O_USAGE, USE_ALSO, 0>& src) noexcept { this->l=src.length(); this->s=src.step(); }
		inline forward(size_t l, size_t s) { this->l=l; this->s=s; }
		inline size_t length() const noexcept { return l; }
		inline size_t step() const noexcept { return s; }
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) noexcept
		{
			p+=o*s;
			return o;
		}
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) const noexcept
		{
			p+=o*s;
			return o;
		}
		template <typename T>
		class iterator
		{
			size_t s;
			size_t c;
			size_t end;
		public:
			inline iterator(const forward& i, const T * const p) noexcept: s(i.s), c(0), end(i.l) {}
			inline void move_one(T *&p) noexcept { p+=s; c++; }
			inline void move_one(const T *&p) noexcept { p+=s; c++; }
			inline void move(T *&p, ptrdiff_t shift) noexcept { p+=s*shift; c+=shift; }
			inline void move(const T *&p, ptrdiff_t shift) noexcept { p+=s*shift; c+=shift; }
			inline bool not_end(const T *p) const noexcept { return c<end; }
		};
	};

	template <int V_TYPE, enum dSU USAGE>
	class forward<V_TYPE,USAGE,0,0>: protected segment_data
	{
		template <int, enum dSU, size_t, size_t> friend struct forward;
	public:
		inline void init(size_t l,size_t s) noexcept { this->l=l; this->s=s; }
		inline forward() noexcept { this->l=0; this->s=0; };
		template <typename PARENT_SHAPE>
		inline forward(const defaultIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=p.length(); this->s=p.step(); }
		template <typename PARENT_SHAPE>
		inline forward(const simpleIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=p.length(); this->s=p.step(); }
		template <typename PARENT_SHAPE>
		inline forward(const forwardIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=i.length(); this->s=p.step(); }
		template <int O_V_TYPE, enum dSU O_USAGE>
		inline forward(const forward<O_V_TYPE, O_USAGE, 0, 0>& src) noexcept { this->l=src.length(); this->s=src.step(); }
		inline forward(size_t l, size_t s) { this->l=l; this->s=s; }
		inline size_t length() const noexcept { return l; }
		inline size_t step() const noexcept { return s; }
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) noexcept
		{
			p+=o*s;
			return o;
		}
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) const noexcept
		{
			p+=o*s;
			return o;
		}
		template <typename T>
		class iterator
		{
			size_t s;
			const T * const end;
		public:
			inline iterator(const forward& i, const T * const p) noexcept: s(i.s), end(p+i.l*i.s) {}
			inline void move_one(T *&p) noexcept { p+=s; }
			inline void move_one(const T *&p) noexcept { p+=s; }
			inline void move(T *&p, ptrdiff_t shift) noexcept { p+=s*shift; }
			inline void move(const T *&p, ptrdiff_t shift) noexcept { p+=s*shift; }
			inline bool not_end(const T *p) const noexcept { return p<end; }
		};
	};

	template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	class forward: protected general_segment_data
	{
	public:
		inline forward() = default;
		template <typename PARENT_SHAPE>
		inline forward(const defaultIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=p.length(); }
		template <typename PARENT_SHAPE>
		inline forward(const simpleIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=p.length(); }
		template <typename PARENT_SHAPE>
		inline forward(const forwardIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=i.length(); }
		template <int OLD_V_TYPE>
		inline forward(const forward<OLD_V_TYPE, USAGE, USE_ALSO, PARENT>& src) noexcept { this->l=src.length(); }
		size_t length() const noexcept { return this->l; }
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) noexcept
		{
			return o;
		}
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) const noexcept
		{
			return o;
		}
		template <typename T>
		class iterator
		{
			size_t c;
			size_t end;
		public:
			inline iterator(const forward& i, const T * const p) noexcept: c(0), end(i.l)
			{
			}
			inline ptrdiff_t move_one(T *&p) noexcept { c++; return 1; }
			inline ptrdiff_t move_one(const T *&p) noexcept { c++; return 1; }
			inline ptrdiff_t move(T *&p, ptrdiff_t shift) noexcept { c+=shift; return shift; }
			inline ptrdiff_t move(const T *&p, ptrdiff_t shift) noexcept { c+=shift; return shift; }
			inline bool not_end(const T *p) const noexcept { return c<end; }
		};
	};

	template <int V_TYPE>
	class forwardIndex
	{
		size_t l;
		size_t o;
	public:
		inline void check_bounds(size_t inum, size_t len) const
		{
			if (l+o>len)
				throw outOfBounds("maximum", inum, l+o-1, len);
			if (o>len)
				throw outOfBounds("minimum", inum, o-1, len);
		}
		inline forwardIndex(size_t l, size_t o) noexcept:l(l),o(o){}
		inline forwardIndex operator+(ptrdiff_t offset) const noexcept { return forwardIndex(this->l,o+offset); }
		inline forwardIndex operator-(ptrdiff_t offset) const noexcept { return forwardIndex(l,o-offset); }
		inline reverseIndex<V_TYPE> operator-() const noexcept { return reverseIndex<V_TYPE>(l,o); }
		inline size_t min_element() const noexcept { return o; }
		inline size_t first() const noexcept { return o; }
		inline size_t max_element() const noexcept { return o+l-1; }
		inline size_t length() const noexcept { return l; }

		template <typename PARENT_SHAPE, size_t P_POS>
		struct applyer;

		template <template <int, enum dSU, size_t, size_t> class shapeClass, int P_V_TYPE, enum dSU P_USAGE, size_t P_USE_ALSO, size_t P_PARENT, size_t P_POS>
		struct applyer<shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>, P_POS>
		{
			using parent_shape=shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>;
			using new_shape = forward<V_TYPE, (P_USAGE<USAGE_ENUM?USAGE_ENUM:P_USAGE),0,P_POS+1>;
		};

		template <int P_V_TYPE, enum dSU P_USAGE, size_t P_POS>
		struct applyer<forward<P_V_TYPE, P_USAGE, 0, 0>, P_POS>
		{
			using parent_shape=forward<P_V_TYPE, P_USAGE, 0, 0>;
			using new_shape = forward<V_TYPE, (P_USAGE<USAGE_ENUM?USAGE_ENUM:P_USAGE) ,0, 0>;
		};
	};
#define DECLARE_forwardIndex(I,length,offset) iTTL::forwardIndex<__COUNTER__> I(length,offset);



	template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	class reverse;

	template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO>
	class reverse<V_TYPE,USAGE,USE_ALSO,0>: protected segment_data
	{
		template <int, enum dSU, size_t, size_t> friend struct segment;
	public:
		inline reverse() = default;
		template <int P_V_TYPE, enum dSU P_USAGE, size_t P_USE_ALSO>
		inline reverse(const reverseIndex<V_TYPE>& i, const segment<P_V_TYPE,P_USAGE,P_USE_ALSO,0>& p) noexcept { this->l=i.length(); this->s=p.step(); }
		template <int OLD_V_TYPE>
		inline reverse(const reverse<OLD_V_TYPE, USAGE, USE_ALSO, 0>& src) noexcept { this->l=src.length(); this->s=src.step(); }
		inline size_t length() const  noexcept{ return l; }
		inline size_t step() const noexcept { return s; }
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) noexcept
		{
			p-=o*s;
			return -o;
		}
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) const noexcept
		{
			p-=o*s;
			return -o;
		}
		template <typename T>
		class iterator
		{
			size_t s;
			size_t c;
			size_t end;
		public:
			inline iterator(const reverse& i, const T * const p) noexcept: s(i.s), c(0), end(i.l) {}
			inline void move_one(T *&p) noexcept { p-=s; c++; }
			inline void move_one(const T *&p) noexcept { p-=s; c++; }
			inline void move(T *&p, ptrdiff_t shift) noexcept { p-=s*shift; c+=shift; }
			inline void move(const T *&p, ptrdiff_t shift) noexcept { p-=s*shift; c+=shift; }
			inline bool not_end(const T *p) const noexcept { return c<end; }
		};
	};

	template <int V_TYPE, enum dSU USAGE>
	class reverse<V_TYPE,USAGE,0,0>: protected segment_data
	{
		template <int, enum dSU, size_t, size_t> friend struct segment;
	public:
		inline reverse() = default;
		template <int P_V_TYPE, enum dSU P_USAGE, size_t P_USE_ALSO>
		inline reverse(const reverseIndex<V_TYPE>& i, const segment<P_V_TYPE,P_USAGE,P_USE_ALSO,0>& p) noexcept { this->l=i.length(); this->s=p.step(); }
		template <int OLD_V_TYPE>
		inline reverse(const reverse<OLD_V_TYPE, USAGE, 0, 0>& src) noexcept { this->l=src.length(); this->s=src.step(); }
		inline size_t length() const noexcept { return l; }
		inline size_t step() const noexcept { return s; }
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) noexcept
		{
			p-=o*s;
			return -o;
		}
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) const noexcept
		{
			p-=o*s;
			return -o;
		}
		template <typename T>
		class iterator
		{
			size_t s;
			const T * const end;
		public:
			inline iterator(const reverse& i, const T * const p) noexcept: s(i.s), end(p-i.l*i.s) {}
			inline void move_one(T *&p) noexcept { p-=s; }
			inline void move_one(const T *&p) noexcept { p-=s; }
			inline void move(T *&p, ptrdiff_t shift) noexcept { p-=s*shift; }
			inline void move(const T *&p, ptrdiff_t shift) noexcept { p-=s*shift; }
			inline bool not_end(const T *p) const noexcept { return p>end; }
		};
	};

	template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	class reverse: protected general_segment_data
	{
		template <int, enum dSU, size_t, typename> friend struct g_segment;
	public:
		inline reverse() = default;
		template <typename PARENT_SHAPE>
		inline reverse(const reverseIndex<V_TYPE>& i, const PARENT_SHAPE& p) noexcept { this->l=i.length(); }
		template <int OLD_V_TYPE>
		inline reverse(const reverse<OLD_V_TYPE, USAGE, USE_ALSO,PARENT>& src) noexcept { this->l=src.length(); }
		inline size_t length() const noexcept { return this->l; }
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) noexcept
		{
			return -o;
		}
		template <typename T>
		inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) const noexcept
		{
			return -o;
		}
		template <typename T>
		class iterator
		{
			size_t c;
			size_t end;
		public:
			inline iterator(const reverse& i, const T * const p) noexcept: c(0), end(i.l) {}
			inline ptrdiff_t move_one(T *&p) noexcept { c++; return -1; }
			inline ptrdiff_t move_one(const T *&p) noexcept { c++; return -1; }
			inline ptrdiff_t move(T *&p, ptrdiff_t shift) noexcept { c+=shift; return -shift; }
			inline ptrdiff_t move(const T *&p, ptrdiff_t shift) noexcept { c+=shift; return -shift; }
			inline bool not_end(const T *p) const noexcept { return c<end; }
		};
	};

	template <int V_TYPE>
	class reverseIndex
	{
		size_t l;
		size_t o;
	public:
		inline void check_bounds(size_t inum, size_t len) const
		{
			if (l+o>len)
				throw outOfBounds("maximum", inum, l+o-1, len);
			if (o>len)
				throw outOfBounds("minimum", inum, o-1, len);
		}
		inline reverseIndex(size_t l, size_t o) noexcept:l(l),o(o){}
		inline reverseIndex operator+(ptrdiff_t offset) const noexcept { return reverseIndex(this->l,o+offset); }
		inline reverseIndex operator-(ptrdiff_t offset) const noexcept { return reverseIndex(l,o-offset); }
		inline forwardIndex<V_TYPE> operator-() const noexcept { return forwardIndex<V_TYPE>(l,o); }
		inline size_t min_element() const noexcept { return o; }
		inline size_t first() const noexcept { return o+l-1; }
		inline size_t max_element() const noexcept { return o+l-1; }
		inline size_t length() const noexcept { return l; }

		template <typename PARENT_SHAPE, size_t P_POS>
		struct applyer;

		template <template <int, enum dSU, size_t, size_t> class shapeClass, int P_V_TYPE, enum dSU P_USAGE, size_t P_USE_ALSO, size_t P_PARENT, size_t P_POS>
		struct applyer<shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>, P_POS>
		{
			using parent_shape=shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>;
			using new_shape = reverse<V_TYPE, (P_USAGE<USAGE_ENUM?USAGE_ENUM:P_USAGE),0,P_POS+1>;
		};

		template <int P_V_TYPE, enum dSU P_USAGE, size_t P_POS>
		struct applyer<segment<P_V_TYPE, P_USAGE, 0, 0>, P_POS>
		{
			using parent_shape=segment<P_V_TYPE, P_USAGE, 0, 0>;
			using new_shape = reverse<V_TYPE, (P_USAGE<USAGE_ENUM?USAGE_ENUM:P_USAGE) ,0, 0>;
		};
	};
#define DECLARE_reverseIndex(I,length,offset) iTTL::reverseIndex<__COUNTER__> I(length,offset);


	template <int V_TYPE>
	class defaultIndex
	{
	public:
		inline void check_bounds(size_t inum, size_t len) const noexcept {}
		inline size_t first() const noexcept { return 0; }
		inline size_t min_element() const noexcept { return 0; }
		inline size_t max_element() const noexcept { return 0; }
		template <typename PARENT_SHAPE, size_t P_POS>
		struct applyer;

		template <template <int, enum dSU, size_t, size_t> class shapeClass, int P_V_TYPE, enum dSU P_USAGE, size_t P_USE_ALSO, size_t P_PARENT, size_t P_POS>
		struct applyer<shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>, P_POS>
		{
			using parent_shape=shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>;
			using new_shape = segment<V_TYPE, P_USAGE,0,P_POS+1>;
		};
		template <template <int, enum dSU, size_t, size_t> class shapeClass, int P_V_TYPE, enum dSU P_USAGE, size_t P_POS>
		struct applyer<shapeClass<P_V_TYPE, P_USAGE, 0, 0>, P_POS>
		{
			using parent_shape=shapeClass<P_V_TYPE, P_USAGE, 0, 0>;
			using new_shape=shapeClass<V_TYPE, P_USAGE, 0, 0>;
		};
	};
#define DECLARE_defaultIndex(I) iTTL::defaultIndex<__COUNTER__> I;


	template <int V_TYPE>
	class simpleIndex
	{
	public:
		inline void check_bounds(size_t inum, size_t len) const noexcept {}
		inline size_t first() const noexcept { return 0; }
		inline size_t min_element() const noexcept { return 0; }
		inline size_t max_element() const noexcept { return 0; }
		template <typename PARENT_SHAPE, size_t P_POS>
		struct applyer;

		template <template <int, enum dSU, size_t, size_t> class shapeClass, int P_V_TYPE, enum dSU P_USAGE, size_t P_USE_ALSO, size_t P_PARENT, size_t P_POS>
		struct applyer<shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>, P_POS>
		{
			using parent_shape=shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>;
			using new_shape = forward<V_TYPE,(P_USAGE>USAGE_ENUM?P_USAGE:USAGE_ENUM),0,P_POS+1>;
		};
		template <template <int, enum dSU, size_t, size_t> class shapeClass, int P_V_TYPE, enum dSU P_USAGE, size_t P_POS>
		struct applyer<shapeClass<P_V_TYPE, P_USAGE, 0, 0>, P_POS>
		{
			using parent_shape=shapeClass<P_V_TYPE, P_USAGE, 0, 0>;
			using new_shape=forward<V_TYPE, (P_USAGE>USAGE_ENUM?P_USAGE:USAGE_ENUM), 0, 0>;
		};
	};
#define DECLARE_simpleIndex(I) iTTL::simpleIndex<__COUNTER__> I;


	template <typename CONTAINER>
	class container
	{
	public:
		template <int V_TYPE>
		class index;
		template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
		class shape
		{
			shared_container<CONTAINER> so;
			typename CONTAINER::iterator start_it;
			template <int> friend class container<CONTAINER>::index;
			template <int, enum dSU, size_t, size_t> friend class container<CONTAINER>::shape;
		public:
			template <typename PARENT_SHAPE>
			inline shape(const index<V_TYPE>& i, const PARENT_SHAPE& p): so(i.sc()), start_it(so.cont.begin()) {}
			template <int OLD_V_TYPE>
			inline shape(const shape<OLD_V_TYPE, USAGE, USE_ALSO,PARENT>& src): so(src.so), start_it(src.start_it) {}
			size_t length() const { return this->so.cont.size(); }
			template <typename T>
			inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o)
			{
				typename CONTAINER::iterator bit=this->start_it;
				if (o<0)
					for (ptrdiff_t i=0;i>o;i--) this->start_it--;
				else
					for (ptrdiff_t i=0;i<o;i++) this->start_it++;
				return (*(start_it)-*bit);
			}
			template <typename T>
			inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) const
			{
				typename CONTAINER::iterator bit=this->start_it;
				if (o<0)
					for (ptrdiff_t i=0;i>o;i--) bit--;
				else
					for (ptrdiff_t i=0;i<o;i++) bit++;
				return (*bit-*(start_it));
			}
			template <typename T>
			class iterator
			{
				typename CONTAINER::iterator it;
				typename CONTAINER::iterator end;
				size_t o;
			public:
				inline iterator(const shape& i, const T * const p) : it(i.so.cont.begin()), end(i.so.cont.end()), o(*it) {}
				inline ptrdiff_t move_one(T *&p) { ptrdiff_t d=*(++it)-o; o+=d; return d; }
				inline ptrdiff_t move_one(const T *&p) { ptrdiff_t d=*(++it)-o; o+=d; return d; }
				inline ptrdiff_t move(T *&p, ptrdiff_t shift)
				{
					ptrdiff_t n;
					std::advance(it,shift);
					ptrdiff_t d=*it-o;
					o+=d;
					return d;
				}
				inline ptrdiff_t move(const T *&p, ptrdiff_t shift)
				{
					ptrdiff_t n;
					std::advance(it,shift);
					ptrdiff_t d=*it-o;
					o+=d;
					return d;
				}
				inline bool not_end(const T *p) const noexcept { return it!=end; } // noexcept???
			};
		};

		template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO>
		class shape<V_TYPE, USAGE, USE_ALSO, 0>
		{
			shared_container<CONTAINER> so;
			size_t s;
			typename CONTAINER::iterator start_it;
			template <int> friend class container<CONTAINER>::index;
			template <int, enum dSU, size_t, size_t> friend class container<CONTAINER>::shape;
		public:
			template <int P_V_TYPE, enum dSU P_USAGE>
			inline shape(const index<V_TYPE>& i, const segment<P_V_TYPE, P_USAGE, 0, 0>& p): so(i.sc()), s(p.step()), start_it(so.cont.begin()) {}
			template <int P_V_TYPE, enum dSU P_USAGE, size_t P_USE_ALSO, size_t P_PARENT>
			inline shape(const defaultIndex<V_TYPE>& i, const shape<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>& p): so(p.so), s(p.s), start_it(so.cont.begin()) {}
			template <int OLD_V_TYPE>
			inline shape(const shape<OLD_V_TYPE, USAGE, USE_ALSO,0>& src): so(src.so), s(src.s), start_it(src.start_it) {}
			size_t length() const noexcept(noexcept(so.cont.size())) { return this->so.cont.size(); }
			template <typename T>
			inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o)
			{
				typename CONTAINER::iterator bit=this->start_it;
				if (o<0)
					for (ptrdiff_t i=0;i>o;i--) this->start_it--;
				else
					for (ptrdiff_t i=0;i<o;i++) this->start_it++;
				p+=(*(start_it)-*bit)*s;
				return (*(start_it)-*bit);
			}
			template <typename T>
			inline ptrdiff_t apply_numeric_index(T *&p, ptrdiff_t o) const
			{
				typename CONTAINER::iterator bit=this->start_it;
				if (o<0)
					for (ptrdiff_t i=0;i>o;i--) bit--;
				else
					for (ptrdiff_t i=0;i<o;i++) bit++;
				p+=(*bit-*(start_it))*s;
				return (*bit-*(start_it));
			}
			template <typename T>
			class iterator
			{
				typename CONTAINER::iterator it;
				typename CONTAINER::iterator end;
				size_t s;
				size_t o;
			public:
				inline iterator(const shape& i, const T * const p): it(i.start_it), end(i.so.cont.end()), s(i.s), o(*it) {}
				inline void move_one(T *&p) { ++it; if (it==end) return; size_t n=*it; p+=(n-o)*s; o=n; }
				inline void move_one(const T *&p) { ++it; if (it==end) return; size_t n=*it; p+=(n-o)*s; o=n; }
				inline void move(T *&p, ptrdiff_t shift)
				{
					ptrdiff_t n;
					std::advance(it,shift);
					n=*it;
					p+=(n-o)*s;
					o=n;
				}
				inline void move(const T *&p, ptrdiff_t shift)
				{
					ptrdiff_t n;
					std::advance(it,shift);
					n=*it;
					p+=(n-o)*s;
					o=n;
				}
				inline bool not_end(const T *p) const noexcept { return it!=end; } // noexcept???
			};
		};

		template <int V_TYPE>
		class index: public shared_container<CONTAINER>
		{
			ptrdiff_t o;
			template <int, enum dSU, size_t, size_t> friend class shape;
		public:
			inline void check_bounds(size_t inum, size_t len) const
			{
				if (max_element()>=len)
					throw outOfBounds("maximum", inum, max_element(), len);
				if (min_element()>=len)
					throw outOfBounds("minimum", inum, min_element(), len);
			}
			inline index(const shared_container<CONTAINER>& so, ptrdiff_t offset) noexcept: shared_container<CONTAINER>(so), o(offset) {}
			template <int V_TYPE_SRC>
			inline index(const index<V_TYPE_SRC>& i, ptrdiff_t offset) noexcept: shared_container<CONTAINER>(i.sc()), o(i.o+offset) {}
			inline size_t first() const { return this->cont.size()==0?o:*(this->cont.begin())+o; }
			inline ptrdiff_t offset() const noexcept { return o; }
			inline size_t min_element() const { return o+(this->cont.begin()==this->cont.end()?0:*(std::min_element(this->cont.begin(),this->cont.end()))); }
			inline size_t max_element() const { return o+(this->cont.begin()==this->cont.end()?0:*(std::max_element(this->cont.begin(),this->cont.end()))); }
			inline const shared_container<CONTAINER>& sc() const noexcept { return *this; }
			inline index operator+(ptrdiff_t offset) const noexcept { return index(*this,offset); }
			inline index operator-(ptrdiff_t offset) const noexcept { return index(*this,-offset); }

			template <typename PARENT_SHAPE, size_t P_POS>
			struct applyer;

			template <template <int, enum dSU, size_t, size_t> class shapeClass, int P_V_TYPE, enum dSU P_USAGE, size_t P_USE_ALSO, size_t P_PARENT, size_t P_POS>
			struct applyer<shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>, P_POS>
			{
				using parent_shape=shapeClass<P_V_TYPE, P_USAGE, P_USE_ALSO, P_PARENT>;
				using new_shape = shape<V_TYPE, P_USAGE<USAGE_ENUM?USAGE_ENUM:P_USAGE,0,P_POS+1>;
			};
			template <int P_V_TYPE, enum dSU P_USAGE, size_t P_POS>
			struct applyer<segment<P_V_TYPE, P_USAGE, 0, 0>, P_POS>
			{
				using parent_shape=segment<P_V_TYPE, P_USAGE, 0, 0>;
				using new_shape=shape<V_TYPE, P_USAGE<USAGE_ENUM?USAGE_ENUM:P_USAGE, 0, 0>;
			};
		};
	};
#define DECLARE_containerIndex(I,cont,offset) typename iTTL::container<typename std::remove_reference<decltype(cont)>::type::base_container_type>::template index<__COUNTER__> I(cont,offset);

	template<template <int> class indexClass>
	struct index_creator
	{
		template <template <int> class srcIndex, int V_TYPE, typename ... Ts>
		static indexClass<V_TYPE> create(const srcIndex<V_TYPE>&, Ts...args) { return indexClass<V_TYPE>(args...); }
	};

	template <typename MASTER, typename SLAVE>
	struct segment_gluer;

	template <int M_V_TYPE, enum dSU M_USAGE, size_t M_USE_ALSO, int V_TYPE, enum dSU USAGE>
	struct segment_gluer<segment<M_V_TYPE, M_USAGE, M_USE_ALSO, 0>, segment<V_TYPE, USAGE, 0, 0> >
	{
		using MASTER=segment<M_V_TYPE, M_USAGE, M_USE_ALSO, 0>;
		static void update(MASTER& msh, const segment<V_TYPE, USAGE, 0, 0>& nsh) noexcept
		{
//			msh.l=nsh.l;
			msh.s+=nsh.s;
		}
	};

	template <typename SHAPE>
	struct check_shape;

	template <template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	struct check_shape<shapeClass<V_TYPE,USAGE,USE_ALSO,PARENT> >
	{
		static const bool need_parent=true;
		static const bool need_glue=false;
		static const bool full_segment=false;
		static const int v_type=V_TYPE;
		static const enum dSU usage=USAGE;
	};

	template <template <int, enum dSU, size_t, size_t> class shapeClass, int V_TYPE, enum dSU USAGE, size_t USE_ALSO>
	struct check_shape<shapeClass<V_TYPE,USAGE,USE_ALSO,0> >
	{
		static const bool need_parent=false;
		static const bool need_glue=false;
		static const bool full_segment=false;
		static const int v_type=V_TYPE;
		static const enum dSU usage=USAGE;
	};

	template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO, size_t PARENT>
	struct check_shape<segment<V_TYPE,USAGE,USE_ALSO,PARENT> >
	{
		static const bool need_parent=true;
		static const bool need_glue=false;
		static const bool full_segment=false;
		static const int v_type=V_TYPE;
		static const enum dSU usage=USAGE;
	};

	template <int V_TYPE, enum dSU USAGE, size_t USE_ALSO>
	struct check_shape<segment<V_TYPE,USAGE,USE_ALSO,0> >
	{
		static const bool need_parent=false;
		static const bool need_glue=true;
		static const bool full_segment=false;
		static const int v_type=V_TYPE;
		static const enum dSU usage=USAGE;
	};

	template <int V_TYPE, size_t USE_ALSO>
	struct check_shape<segment<V_TYPE,USAGE_FULL,USE_ALSO,0> >
	{
		static const bool need_parent=false;
		static const bool need_glue=true;
		static const bool full_segment=true;
		static const int v_type=V_TYPE;
		static const enum dSU usage=USAGE_FULL;
	};
};

#endif /* TPPINDEX_H_ */
