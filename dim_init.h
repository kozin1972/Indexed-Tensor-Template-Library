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
		static bool get_aw(size_t (&aw)[1], T0 w0, T1 w1, Ts ... wr) noexcept { aw[0]=sizeof...(Ts)+2; return true; }

		template <typename T, typename T0, typename ... Ts>
		static void init_data(size_t (&aw)[1], T *data, T0 v0, Ts ...  args) noexcept
		{
			data[0]=v0;
			init_data(aw, data+1, args...);
		}

		template <typename T>
		static void init_data(size_t (&aw)[1], T *data) noexcept
		{
		}

		template <typename T, size_t NUM0, size_t NUM1>
		static void init_data(size_t (&aw)[2], T *data, T arr[NUM0][NUM1]) noexcept
		{
			for (size_t i=0;i<NUM0;i++)
				for (size_t j=0;j<NUM1;j++)
					data[i*NUM1+j]=arr[i][j];
		}

		template <typename T, size_t NUM0, size_t NUM1>
		static bool get_aw(size_t (&aw)[2], T (&v)[NUM0][NUM1]) noexcept { aw[0]=NUM0; aw[1]=NUM1; return true; }

		template <typename T>
		static void init_data(size_t (&aw)[2], T *data, std::initializer_list<std::initializer_list<T> > list) noexcept
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
		inline static void init_data(size_t (&aw)[DIM], T *data, Ts ... args) noexcept
		{
		}

		template <typename T>
		static bool get_aw(size_t (&aw)[2], std::initializer_list<std::initializer_list<T> > list) noexcept
		{
			aw[0]=list.size();
			aw[1]=0;
			for (auto it=list.begin();it!=list.end();it++)
				if (it->size()>aw[1])
					aw[1]=it->size();
			return true;
		}

		template <size_t DIM, typename ... Ts>
		static bool get_aw_plain(size_t (&aw)[DIM], int w0, Ts ... wr) noexcept
		{
			aw[0]=w0;
			get_aw_plain(*((size_t (*)[DIM-1])&aw[1]), wr...);
			return false;
		}
		static bool get_aw_plain(int (&aw)[0]) noexcept
		{
			return false;
		}
		static bool get_aw_plain(size_t (&aw)[0]) noexcept
		{
			return false;
		}
		template <size_t DIM, typename ... Ts>
		static bool get_aw_plain(size_t (&aw)[DIM], size_t w0, Ts ... wr) noexcept
		{
			aw[0]=w0;
			get_aw_plain(*((size_t (*)[DIM-1])&aw[1]), wr...);
			return false;
		}

		template <size_t DIM, typename ... Ts>
		static bool get_aw(size_t (&aw)[DIM], int w0, Ts ... wr) noexcept
		{
			static_assert(sizeof...(Ts)+1==DIM,"Wrong number of initializers");
			aw[0]=w0;
			get_aw_plain(*((size_t (*)[DIM-1])&aw[1]), wr...);
			return false;
		}
		template <size_t DIM, typename ... Ts>
		static bool get_aw(size_t (&aw)[DIM], size_t w0, Ts ... wr) noexcept
		{
			static_assert(sizeof...(Ts)+1==DIM,"Wrong number of initializers");
			aw[0]=w0;
			get_aw_plain(*((size_t (*)[DIM-1])&aw[1]), wr...);
			return false;
		}

		template <size_t DIM>
		static bool get_aw(size_t (&aw)[DIM], int (&awp)[DIM]) noexcept { for (size_t i=0; i<DIM; i++) aw[i]=awp[i]; return false; }

		template <size_t DIM>
		static bool get_aw(size_t (&aw)[DIM], size_t (&awp)[DIM]) noexcept { for (size_t i=0; i<DIM; i++) aw[i]=awp[i]; return false; }

	};
};

#endif /* DIM_INIT_H_ */
