/*
 * dim_init.h
 *
 *  Created on: 4 июн. 2019 г.
 *      Author: Alexey Kozin
 *
 *  Please report bugs to tpptensor@mail.ru
 */

#ifndef DIM_INIT_H_
#define DIM_INIT_H_

namespace tpp
{

	struct aw_getter
	{

		template <typename T0, typename T1, typename ... Ts>
		static bool get_aw(size_t (&aw)[1], T0 w0, T1 w1, Ts ... wr) { aw[0]=sizeof...(Ts)+2; return true; }

		template <typename T, typename T0, typename ... Ts>
		static void init_data(size_t (&aw)[1], T *data, T0 v0, Ts ...  args)
		{
			data[0]=v0;
			init_data(aw, data+1, args...);
		}

		template <typename T>
		static void init_data(size_t (&aw)[1], T *data)
		{
		}

		template <typename T, size_t NUM0, size_t NUM1>
		static void init_data(size_t (&aw)[2], T *data, T arr[NUM0][NUM1])
		{
			for (size_t i=0;i<NUM0;i++)
				for (size_t j=0;j<NUM1;j++)
					data[i*NUM1+j]=arr[i][j];
		}

		template <typename T, size_t NUM0, size_t NUM1>
		static bool get_aw(size_t (&aw)[2], T (&v)[NUM0][NUM1]) { aw[0]=NUM0; aw[1]=NUM1; return true; }

		template <typename T>
		static void init_data(size_t (&aw)[2], T *data, std::initializer_list<std::initializer_list<T> > list)
		{
			size_t row=0;
			for (auto it=list.begin();it!=list.end();it++,row++)
			{
				std::copy(it->begin(),it->end(),data+row*aw[1]);
				if (it->size()<aw[1])
				{
					for (size_t i=it->size();i<aw[1];i++)
						data[row*aw[1]+i]=0;
				}
			}
		}

		template <size_t DIM, typename T, typename ...Ts>
		inline static void init_data(size_t (&aw)[DIM], T *data, Ts ... args)
		{
		}

		template <typename T>
		static bool get_aw(size_t (&aw)[2], std::initializer_list<std::initializer_list<T> > list)
		{
			aw[0]=list.size();
			aw[1]=0;
			for (auto it=list.begin();it!=list.end();it++)
				if (it->size()>aw[1])
					aw[1]=it->size();
			return true;
		}

//		template <typename T, size_t NUM0, size_t NUM1>


//		template <size_t DIM, typename ... Ts>
//		init_data(size_t (&aw)[DIM], T *data, Ts ...  args)

//		bool get_aw(size_t (&aw)[0]) { return false; }
//
////		bool get_aw(size_t (&aw)[1], int w0) { aw[0]=w0; return false; }
//
//		bool get_aw(size_t (&aw)[1], size_t w0) { aw[0]=w0; return false; }
//
//		bool get_aw(size_t (&aw)[1], int (&awp)[1]) { aw[0]=awp[0]; return false; }
//
//		bool get_aw(size_t (&aw)[1], size_t (&awp)[1]) { aw[0]=awp[0]; return false; }
//
//		template <typename T, size_t NUM>
//		bool get_aw(size_t (&aw)[1], T (&v)[NUM]) { aw[0]=NUM; return true; }

////	bool get_aw(size_t (&aw)[2], int w0, int w1) { aw[0]=w0; aw[1]=w1; return false; }
//
//	bool get_aw(size_t (&aw)[2], size_t w0, size_t w1) { aw[0]=w0; aw[1]=w1; return false; }
//
//	bool get_aw(size_t (&aw)[2], int (&awp)[2]) { aw[0]=awp[0]; aw[1]=awp[1]; return false; }
//
//	bool get_aw(size_t (&aw)[2], size_t (&awp)[2]) { aw[0]=awp[0]; aw[1]=awp[1]; return false; }
//
//
////	bool get_aw(size_t (&aw)[3], int w0, int w1, int w2) { aw[0]=w0; aw[1]=w1; aw[2]=w2; return false; }
//
//	bool get_aw(size_t (&aw)[3], size_t w0, size_t w1, size_t w2) { aw[0]=w0; aw[1]=w1; aw[2]=w2; return false; }
//
//	bool get_aw(size_t (&aw)[3], int (&awp)[3]) { aw[0]=awp[0]; aw[1]=awp[1]; aw[2]=awp[2]; return false; }
//
//	bool get_aw(size_t (&aw)[3], size_t (&awp)[3]) { aw[0]=awp[0]; aw[1]=awp[1]; aw[2]=awp[2]; return false; }
//
//	template <typename T, size_t NUM0, size_t NUM1, size_t NUM2>
//	bool get_aw(size_t (&aw)[3], T (&v)[NUM0][NUM1][NUM2]) { aw[0]=NUM0; aw[1]=NUM1; aw[2]=NUM2; return true; }

		template <size_t DIM, typename ... Ts>
		static bool get_aw_plain(size_t (&aw)[DIM], int w0, Ts ... wr)
		{
			aw[0]=w0;
			get_aw_plain(*((size_t (*)[DIM-1])&aw[1]), wr...);
			return false;
		}
		static bool get_aw_plain(int (&aw)[0])
		{
			return false;
		}
		static bool get_aw_plain(size_t (&aw)[0])
		{
			return false;
		}
		template <size_t DIM, typename ... Ts>
		static bool get_aw_plain(size_t (&aw)[DIM], size_t w0, Ts ... wr)
		{
			aw[0]=w0;
			get_aw_plain(*((size_t (*)[DIM-1])&aw[1]), wr...);
			return false;
		}

		template <size_t DIM, typename ... Ts>
		static bool get_aw(size_t (&aw)[DIM], int w0, Ts ... wr)
		{
			static_assert(sizeof...(Ts)+1==DIM,"Wrong number of initializers");
			aw[0]=w0;
			get_aw_plain(*((size_t (*)[DIM-1])&aw[1]), wr...);
			return false;
		}
		template <size_t DIM, typename ... Ts>
		static bool get_aw(size_t (&aw)[DIM], size_t w0, Ts ... wr)
		{
			static_assert(sizeof...(Ts)+1==DIM,"Wrong number of initializers");
			aw[0]=w0;
			get_aw_plain(*((size_t (*)[DIM-1])&aw[1]), wr...);
			return false;
		}

		template <size_t DIM>
		static bool get_aw(size_t (&aw)[DIM], int (&awp)[DIM]) { for (size_t i=0; i<DIM; i++) aw[i]=awp[i]; return false; }

		template <size_t DIM>
		static bool get_aw(size_t (&aw)[DIM], size_t (&awp)[DIM]) { for (size_t i=0; i<DIM; i++) aw[i]=awp[i]; return false; }

	};
};

#endif /* DIM_INIT_H_ */
