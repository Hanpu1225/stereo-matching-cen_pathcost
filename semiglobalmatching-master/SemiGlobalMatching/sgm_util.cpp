/* -*-c++-*- SemiGlobalMatching - Copyright (C) 2020.
* Author	: Yingsong Li(Ethan Li) <ethan.li.whu@gmail.com>
* https://github.com/ethan-li-coding/SemiGlobalMatching
* Describe	: implement of sgm_util
*/

#include "stdafx.h"
#include "sgm_util.h"
#include <algorithm>
#include <cassert>
#include <vector>
#include <queue>
#include<numeric>

//void sgm_util::census_transform_5x5(const uint8* source, uint32* census, const sint32& width,
//	const sint32& height)
//{
//	float32 sum = 0u;
//	float32 aver = 0u;
//	float32 beita = 0u;
//	const float32 size = width * height;
//	std::vector<float32>wd_data;
//	wd_data.reserve(size);
//
//	if (source == nullptr || census == nullptr || width <= 5 || height <= 5)
//	{
//		return;
//	}
//	// 逐像素计算census值
//	for (sint32 i = 2; i < height - 2; i++)
//	{
//		for (sint32 j = 2; j < width - 2; j++)
//		{
//			wd_data.clear();
//
//			for (sint32 r = -2; r <= 2; r++)
//			{
//				for (sint32 c = -2; c <= 2; c++)
//				{
//					const uint8 gray = source[(i + r) * width + j + c];
//					wd_data.push_back(gray);
//				}
//			}
//
//			/**********************************左上角9个数据******************************************/
//			std::vector<float32>left_vector(wd_data.begin(), wd_data.begin() + 13);//前面13个
//			std::vector<float32>left_vector1;//存储左上角的9个数据
//
//			left_vector1.clear();
//			for (int i = 0; i < left_vector.size(); i++)
//			{
//				if (i != 3 && i != 4 && i != 8 && i != 9)
//				{
//					left_vector1.push_back(left_vector[i]);
//				}
//			}
//			float32 leftUpsum = std::accumulate(left_vector1.begin(), left_vector1.end(), 0.0f);//9个数据的和
//			float32 leftUpAver = leftUpsum / left_vector1.size();//9个数据的平均值
//			std::sort(left_vector1.begin(), left_vector1.end());//对九个数据进行排序
//			float32 leftUpMid = left_vector1[left_vector1.size() / 2];//9个数据的中值
//
//			float32 leftUpSUM = 0;
//			for (int i = 0; i < left_vector1.size(); i++)
//			{
//				float32 v = abs((leftUpAver - left_vector1[i]) * (leftUpAver - left_vector1[i]));
//				leftUpSUM += v;
//			}
//
//			/**********************************右上角9个数据******************************************/
//
//			std::vector<float32>right_vector(wd_data.begin() + 2, wd_data.begin() + 15);
//			std::vector<float32>right_vector1;
//
//			right_vector1.clear();
//			for (int i = 0; i < right_vector.size(); i++)
//			{
//				if (i != 3 && i != 4 && i != 8 && i != 9)
//				{
//					right_vector1.push_back(right_vector[i]);
//				}
//
//			}
//			float32 rightUpsum = std::accumulate(right_vector1.begin(), right_vector1.end(), 0.0f);
//			float32 rightUpAver = rightUpsum / right_vector1.size();
//			std::sort(right_vector1.begin(), right_vector1.end());
//			float32 rightUpMid = right_vector1[right_vector1.size() / 2];
//
//			float32 rightUpSUM = 0;
//			for (int i = 0; i < right_vector1.size(); i++)
//			{
//				float32 v = abs((rightUpAver - right_vector1[i]) * (rightUpAver - right_vector1[i]));
//				rightUpSUM += v;
//			}
//			/**********************************左下角9个数据******************************************/
//			std::vector<float32>leftDown_vector(wd_data.begin() + 10, wd_data.end() - 2);
//			std::vector<float32>leftDown_vector1;
//
//			leftDown_vector1.clear();
//			for (int i = 0; i < leftDown_vector.size(); i++)
//			{
//				if (i != 3 && i != 4 && i != 8 && i != 9)
//				{
//					leftDown_vector1.push_back(leftDown_vector[i]);
//				}
//
//			}
//
//			float32 leftDownsum = std::accumulate(leftDown_vector1.begin(), leftDown_vector1.end(), 0.0f);
//			float32 leftDownAver = leftDownsum / leftDown_vector1.size();
//			std::sort(leftDown_vector1.begin(), leftDown_vector1.end());
//			float32 leftDownMid = leftDown_vector1[leftDown_vector1.size() / 2];
//
//			float32 leftDownSUM = 0;
//			for (int i = 0; i < leftDown_vector1.size(); i++)
//			{
//				float32 v = abs((leftDownAver - leftDown_vector1[i]) * (leftDownAver - leftDown_vector1[i]));
//				leftDownSUM += v;
//			}
//			/**********************************右下角9个数据******************************************/
//			std::vector<int> rightDown_vector(wd_data.begin() + 12, wd_data.end());
//			std::vector<int>rightDown_vector1;
//
//			rightDown_vector1.clear();
//			for (int i = 0; i < rightDown_vector.size(); i++)
//			{
//				if (i != 3 && i != 4 && i != 8 && i != 9)
//				{
//					rightDown_vector1.push_back(rightDown_vector[i]);
//				}
//
//			}
//			float32 rightDownsum = std::accumulate(rightDown_vector1.begin(), rightDown_vector1.end(), 0.0f);
//			float32 rightDownAver = rightDownsum / rightDown_vector1.size();
//			std::sort(rightDown_vector1.begin(), rightDown_vector1.end());
//			float32 rightDownMid = rightDown_vector1[rightDown_vector1.size() / 2];
//
//			float32 rightDownSUM = 0;
//			for (int i = 0; i < rightDown_vector1.size(); i++)
//			{
//				float32 v = abs((rightDownAver - rightDown_vector1[i]) * (rightDownAver - rightDown_vector1[i]));
//				rightDownSUM += v;
//			}
//
//			/**********************************中间窗口9个数据******************************************/
//			std::vector<int> Middle_vector(wd_data.begin() + 6, wd_data.end() - 6);
//			std::vector<int>Middle_vector1;
//
//			Middle_vector1.clear();
//			for (int i = 0; i < Middle_vector.size(); i++)
//			{
//				if (i != 3 && i != 4 && i != 8 && i != 9)
//				{
//					Middle_vector1.push_back(Middle_vector[i]);
//				}
//
//			}
//
//			float32 Midsum = std::accumulate(Middle_vector1.begin(), Middle_vector1.end(), 0.0f);
//			float32 Midaver = Midsum / Middle_vector1.size();
//			std::sort(Middle_vector1.begin(), Middle_vector1.end());
//			float32 Midmid = Middle_vector1[Middle_vector1.size() / 2];
//
//			float32 MidSUM = 0;
//			for (int i = 0; i < Middle_vector1.size(); i++)
//			{
//				float32 v = abs((Midaver - Middle_vector1[i]) * (Midaver - Middle_vector1[i]));
//				MidSUM += v;
//			}
//
//			/**********************************计算最小灰度均匀分布******************************************/
//
//			float32 Aver = 0;
//			float32 min = leftUpSUM;
//			auto maxPostion = *max_element(left_vector1.begin(), left_vector1.end());
//			auto minPostion = *min_element(left_vector1.begin(), left_vector1.end());
//			float32 leftUpbeita = maxPostion - minPostion;
//			if (leftUpbeita > 10)
//			{
//				Aver = leftUpMid;
//			}
//			else
//			{
//				Aver = leftUpAver;
//			}
//
//			if (min > rightUpSUM)
//			{
//				min = rightUpSUM;
//				auto maxPostion = *max_element(right_vector1.begin(), right_vector1.end());
//				auto minPostion = *min_element(right_vector1.begin(), right_vector1.end());
//				int rightUpbeita = maxPostion - minPostion;
//				if (rightUpbeita > 10)
//				{
//					Aver = rightUpMid;
//				}
//				else
//				{
//					Aver = rightUpAver;
//				}
//			}
//
//			if (min > leftDownSUM)
//			{
//				min = leftDownSUM;
//				auto maxPostion = *max_element(leftDown_vector1.begin(), leftDown_vector1.end());
//				auto minPostion = *min_element(leftDown_vector1.begin(), leftDown_vector1.end());
//				int leftDownbeita = maxPostion - minPostion;
//				if (leftDownbeita > 10)
//				{
//					Aver = leftDownMid;
//				}
//				else
//				{
//					Aver = leftDownAver;
//				}
//			}
//
//			if (min > rightDownSUM)
//			{
//				min = rightDownSUM;
//				auto maxPostion = *max_element(rightDown_vector1.begin(), rightDown_vector1.end());
//				auto minPostion = *min_element(rightDown_vector1.begin(), rightDown_vector1.end());
//				int rightDownbeita = maxPostion - minPostion;
//				if (rightDownbeita > 10)
//				{
//					Aver = rightDownMid;
//				}
//				else
//				{
//					Aver = rightDownAver;
//				}
//			}
//
//			if (min > MidSUM)
//			{
//				min = MidSUM;
//				auto maxPostion = *max_element(Middle_vector1.begin(), Middle_vector1.end());
//				auto minPostion = *min_element(Middle_vector1.begin(), Middle_vector1.end());
//				int Middlebeita = maxPostion - minPostion;
//				if (Middlebeita > 10)
//				{
//					Aver = Midmid;
//				}
//				else
//				{
//					Aver = Midaver;
//				}
//			}
//
//			uint32 census_val = 0u;
//			for (auto iter = wd_data.cbegin(); iter < wd_data.cend(); iter++)
//			{
//				census_val <<= 1;
//				if (*iter < Aver)
//				{
//					census_val += 1;
//				}
//			}
//			// 中心像素的census值
//			census[i * width + j] = census_val;
//		}
//	}
//}

//void sgm_util::census_transform_5x5(const uint8* source, uint32* census, const sint32& width,
//	const sint32& height)
//{
//	float32 sum = 0u;
//	float32 aver = 0u;
//	float32 beita = 0u;
//	const float32 size = width * height;
//	std::vector<float32>wd_data;
//	wd_data.reserve(size);
//
//	if (source == nullptr || census == nullptr || width <= 5 || height <= 5)
//	{
//		return;
//	}
//	// 逐像素计算census值
//	for (sint32 i = 2; i < height - 2; i++)
//	{
//		for (sint32 j = 2; j < width - 2; j++)
//		{
//			wd_data.clear();
//
//			for (sint32 r = -2; r <= 2; r++)
//			{
//				for (sint32 c = -2; c <= 2; c++)
//				{
//					const uint8 gray = source[(i + r) * width + j + c];
//					wd_data.push_back(gray);
//				}
//			}
//
//			std::vector<float32>::iterator i1 = wd_data.begin(), i2 = wd_data.begin() + 15;
//			std::vector<float32>::iterator i3 = wd_data.end() - 15, i4 = wd_data.end();
//			std::vector<float32>::iterator i5 = wd_data.begin() + 5, i6 = wd_data.end() - 5;
//			std::vector<float32>wd_data1(i1, i2);//前面15个
//			std::vector<float32>wd_data2(i3, i4);//后面15个
//			std::vector<float32>wd_data3(i5, i6);//中间15个
//
//			sum = std::accumulate(wd_data3.begin(), wd_data3.end(), 0.0f);//中间15个的灰度均值
//			aver = sum / wd_data3.size();
//			std::sort(wd_data3.begin(), wd_data3.end());
//			float32 mid = wd_data3[wd_data3.size() / 2];
//
//			float32 SUM = 0.0f;//计算窗口内全部灰度均匀分布度
//			for (sint32 i = 0; i < wd_data3.size(); i++)
//			{
//				//V = SUM(f(i,j)-fp(i,j));
//				float32 v = abs((aver - wd_data3[i]) * (aver - wd_data3[i]));
//				SUM += v;
//			}
//
//
//			float32 sum1 = std::accumulate(wd_data1.begin(), wd_data1.end(), 0.0f);
//			float32 aver1 = sum1 / wd_data1.size();//前13像素的灰度值的平均值
//			std::sort(wd_data1.begin(), wd_data1.end());
//			float32 mid1 = wd_data1[wd_data1.size() / 2];
//
//
//			float32 SUM1 = 0.0f;//计算前面13个灰度均匀分布度
//			for (sint32 i = 0; i < wd_data1.size(); i++)
//			{
//				//V = SUM(f(i,j)-fp(i,j));
//				float32 v = abs((aver1 - wd_data1[i]) * (aver1 - wd_data1[i]));
//				SUM1 += v;
//			}
//
//			float32 sum2 = std::accumulate(wd_data2.begin(), wd_data2.end(), 0.0f);//后面13个像素灰度平均值
//			float32 aver2 = sum2 / wd_data2.size();
//			std::sort(wd_data2.begin(), wd_data2.end());
//			float32 mid2 = wd_data2[wd_data2.size() / 2];
//
//			float32 SUM2 = 0.0f;//计算后面13个灰度均匀分布度
//			for (sint32 i = 0; i < wd_data2.size(); i++)
//			{
//				float32 v = abs((aver2 - wd_data2[i]) * (aver2 - wd_data2[i]));
//				SUM2 += v;
//			}
//
//			float32 Aver = 0.0f;
//			float32 min = SUM;
//			auto maxPostion = *max_element(wd_data3.begin(), wd_data3.end());
//			auto minPostion = *min_element(wd_data3.begin(), wd_data3.end());
//			beita = maxPostion - minPostion;
//			if (beita > 10)
//			{
//				Aver = mid;
//			}
//			else
//			{
//				Aver = aver;
//			}
//
//			if (min > SUM1)
//			{
//				min = SUM1;
//				auto maxPostion = *max_element(wd_data1.begin(), wd_data1.end());
//				auto minPostion = *min_element(wd_data1.begin(), wd_data1.end());
//				beita = maxPostion - minPostion;
//				if (beita > 10)
//				{
//					Aver = mid1;
//				}
//				else
//				{
//					Aver = aver1;
//				}
//				//Aver = aver1;
//
//			}
//			if (min > SUM2)
//			{
//				min = SUM2;
//				auto maxPostion = *max_element(wd_data2.begin(), wd_data2.end());
//				auto minPostion = *min_element(wd_data2.begin(), wd_data2.end());
//				beita = maxPostion - minPostion;
//				if (beita > 10)
//				{
//					Aver = mid2;
//				}
//				else
//				{
//					Aver = aver2;
//				}
//				//Aver = aver2;
//			}
//
//			uint32 census_val = 0u;
//			for (auto iter = wd_data.cbegin(); iter < wd_data.cend(); iter++)
//			{
//				census_val <<= 1;
//				if (*iter < Aver)
//				{
//					census_val += 1;
//				}
//			}
//			// 中心像素的census值
//			census[i * width + j] = census_val;
//		}
//	}
//}


//void sgm_util::census_transform_5x5(const uint8* source, uint32* census, const sint32& width,
//	const sint32& height)
//{
//	if (source == nullptr || census == nullptr || width <= 5 || height <= 5) {
//		return;
//	}
//
//	float32 sum = 0u;
//	float32 aver = 0u;
//	float32 beita = 0u;
//	const sint32 size = 25u;
//	std::vector<float32>wd_data;
//	wd_data.reserve(size);
//	// 逐像素计算census值
//	for (sint32 i = 2; i < height - 2; i++) 
//	{
//		for (sint32 j = 2; j < width - 2; j++) 
//		{
//			wd_data.clear();
//			// 中心像素值
//			const uint8 gray_center = source[i * width + j];
//			//uint32 census_val = 0u;
//			for (sint32 r = -2; r <= 2; r++) 
//			{
//				for (sint32 c = -2; c <= 2; c++) 
//				{
//					//census_val <<= 1;
//					const uint8 gray = source[(i + r) * width + j + c];
//					wd_data.push_back(gray);
//
//				}
//			}
//
//			auto maxPostion = *max_element(wd_data.begin(), wd_data.end());
//			auto minPostion = *min_element(wd_data.begin(), wd_data.end());
//			beita = maxPostion - minPostion;
//			if (beita > 10)
//			{
//				std::sort(wd_data.begin(), wd_data.end());
//				aver = wd_data[wd_data.size() / 2];
//			}
//			else
//			{
//				sum = std::accumulate(wd_data.begin(), wd_data.end(), 0.0f);
//				aver = sum / wd_data.size();
//			}
//			uint32 census_val = 0u;
//			for (auto iter = wd_data.cbegin(); iter < wd_data.cend(); iter++)
//			{
//				census_val <<= 1;
//				if (*iter < aver)
//				{
//					census_val += 1;
//				}
//			}
//			// 中心像素的census值
//			census[i * width + j] = census_val;
//		}
//	}
//}

//void sgm_util::census_transform_5x5(const uint8* source, uint32* census, const sint32& width,
//	const sint32& height)
//{
//	if (source == nullptr || census == nullptr || width <= 5 || height <= 5) {
//		return;
//	}
//
//	// 逐像素计算census值
//	for (sint32 i = 2; i < height - 2; i++) {
//		for (sint32 j = 2; j < width - 2; j++) {
//
//			// 中心像素值
//			const uint8 gray_center = source[i * width + j];
//
//			// 遍历大小为5x5的窗口内邻域像素，逐一比较像素值与中心像素值的的大小，计算census值
//			uint32 census_val = 0u;
//			for (sint32 r = -2; r <= 2; r++) {
//				for (sint32 c = -2; c <= 2; c++) {
//					census_val <<= 1;
//					const uint8 gray = source[(i + r) * width + j + c];
//					if (gray < gray_center) {
//						census_val += 1;
//					}
//				}
//			}
//
//			// 中心像素的census值
//			census[i * width + j] = census_val;
//		}
//	}
//}

void sgm_util::census_transform_9x7(const uint8* source, uint64* census, const sint32& width,
	const sint32& height, const sint32 wnd_size)
{
	/*if (source == nullptr || census == nullptr || width <= 9 || height <= 7) {
		return;
	}*/

	const sint32 radius = wnd_size / 2;
	const sint32 size = wnd_size * wnd_size;

	const double PI = 4.0 * atan(1.0);//圆周率pi赋值
		//const double PI = 3.14;
		// 存储局部窗口内的数据
	std::vector<double> wnd_data;
	wnd_data.reserve(size);
	std::vector<double> wnd_weight;

	// 逐像素计算census值
	for (sint32 i = radius; i < height - radius; i++) {
		for (sint32 j = radius; j < width - radius; j++) {

			wnd_data.clear();
			wnd_weight.clear();
			// 中心像素值
			const float32 gray_center = source[i * width + j];
			double weight = 0.0f;
			double iSum = 0.0f;

			// 遍历大小为5x5的窗口内邻域像素，逐一比较像素值与中心像素值的的大小，计算census
			for (sint32 r = -radius; r <= radius; r++)
			{
				for (sint32 c = -radius; c <= radius; c++)
				{
					double distance_Square = r * r + c * c;
					double weight = (1 / (PI * 2 * 1.5 * 1.5) * exp(-1 * distance_Square / (2 * 1.5 * 1.5)));
					wnd_weight.push_back(weight);
					const float32 gray = source[(i + r) * width + j + c];
					iSum += weight;//总权重
					wnd_data.push_back(gray);
				}
			}

			double Sum = 0.0f;
			for (int i = 0; i < wnd_weight.size(); i++)
			{
				wnd_weight[i] /= iSum;
				//每个权重/总的权重
				double Gass = wnd_weight[i] * wnd_data[i];//灰度值*权重值
				Sum += Gass;
			}
			double I = 0.0f;
			if (abs(gray_center - Sum) <= 18)
			{
				I = gray_center;
			}
			else
			{
				I = Sum;
			}

			uint64 census_val = 0u;
			for (int i = 0; i < wnd_data.size(); i++)
			{
				census_val <<= 1;
				if (wnd_data[i] < I)
				{
					census_val += 1;
				}
			}

			// 中心像素的census值
			census[i * width + j] = census_val;
		}
	}
}

uint8 sgm_util::Hamming32(const uint32& x, const uint32& y)
{
	uint32 dist = 0, val = x ^ y;

	// Count the number of set bits
	while (val) {
		++dist;
		val &= val - 1;
	}

	return static_cast<uint8>(dist);
}

uint8 sgm_util::Hamming64(const uint64& x, const uint64& y)
{
	uint64 dist = 0, val = x ^ y;

	// Count the number of set bits
	while (val) {
		++dist;
		val &= val - 1;
	}

	return static_cast<uint8>(dist);
}


PGradient sgm_util::GetGradient(const PGradient* grad_data, const sint32& width, const sint32& x, const sint32& y)
{
	return grad_data[y * width + x];
}

void sgm_util::MedianFilter(const float32* in, float32* out, const sint32& width, const sint32& height,
	const sint32 wnd_size)
{
	const sint32 radius = wnd_size / 2;
	const sint32 size = wnd_size * wnd_size;

	// 存储局部窗口内的数据
	std::vector<float32> wnd_data;
	wnd_data.reserve(size);

	for (sint32 i = 0; i < height; i++) {
		for (sint32 j = 0; j < width; j++) {
			wnd_data.clear();

			// 获取局部窗口数据
			for (sint32 r = -radius; r <= radius; r++) {
				for (sint32 c = -radius; c <= radius; c++) {
					const sint32 row = i + r;
					const sint32 col = j + c;
					if (row >= 0 && row < height && col >= 0 && col < width) {
						wnd_data.push_back(in[row * width + col]);
					}
				}
			}

			// 排序
			std::sort(wnd_data.begin(), wnd_data.end());
			// 取中值
			out[i * width + j] = wnd_data[wnd_data.size() / 2];
		}
	}
}

void sgm_util::AdaptiveMedianFilter(const float32* in, float32* out, const sint32& width, const sint32& height,
	sint32 wnd_size, int maxSize)
{
	const sint32 radius = wnd_size / 2;
	const sint32 size = wnd_size * wnd_size;

	// 存储局部窗口内的数据
	std::vector<float32> wnd_data;
	wnd_data.reserve(size);

	for (sint32 i = 0; i < height; i++)
	{
		for (sint32 j = 0; j < width; j++)
		{
			wnd_data.clear();
			const float gray_center = in[i * width + j];//窗口中心点 *

			// 获取局部窗口数据
			for (sint32 r = -radius; r <= radius; r++)
			{
				for (sint32 c = -radius; c <= radius; c++)
				{
					const sint32 row = i + r;
					const sint32 col = j + c;
					if (row >= 0 && row < height && col >= 0 && col < width)
					{
						wnd_data.push_back(in[row * width + col]);
					}
				}
			}

			// 排序
			std::sort(wnd_data.begin(), wnd_data.end());
			// 取中值
			//out[i * width + j] = wnd_data[wnd_data.size() / 2];

			auto min = wnd_data[0];//最小灰度值
			auto max = wnd_data[wnd_data.size() - 1];//最大的灰度值
			auto med = wnd_data[wnd_data.size() / 2];//中值
			auto zxy = gray_center;//未排序前的中心像素

			if (med > min && med < max)
			{
				if (zxy > min && zxy < max)

					out[i * width + j] = zxy;
				else
					out[i * width + j] = wnd_data[wnd_data.size() / 2];

			}
			else
			{
				wnd_size += 2;
				if (wnd_size <= maxSize)
				{
					AdaptiveMedianFilter(in, out, width, height, wnd_size, maxSize);
				}
				else
				{
					out[i * width + j] = wnd_data[wnd_data.size() / 2];
				}
			}
		}
	}
}

/*
在自适应中值滤波算法中，A步骤里面会先判断是否满足Zmin<Zmed<Zmax。
这一步骤实质是判断当前区域的中值点是否是噪声点，通常来说是满足Zmin<Zmed<Zmax这个条件的，
此时中值点不是噪声点，跳转到B；考虑一些特殊情况，如果Zmed=Zmin或者Zmed=Zmax，则认为是噪声点，应该扩大窗口尺寸，
在一个更大的范围内寻找一个合适的非噪声点，随后再跳转到B，否则输出的中值点是噪声点；
接下来考虑跳转到B之后的情况：判断中心点的像素值是否是噪声点，判断条件为Zmin<Zxy<Zmax，
原理同上，因为如果Zxy=Zmin或者Zxy=Zmax，则认为是噪声点。如果不是噪声点，我们可以保留当前像素点的灰度值；
如果是噪声点，则使用中值替代原始灰度值，滤去噪声。
*/