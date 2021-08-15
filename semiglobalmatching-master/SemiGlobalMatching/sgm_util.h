/* -*-c++-*- SemiGlobalMatching - Copyright (C) 2020.
* Author	: Yingsong Li(Ethan Li) <ethan.li.whu@gmail.com>
* https://github.com/ethan-li-coding/SemiGlobalMatching
* Describe	: implement of sgm_util
*/

#pragma once
#include "sgm_types.h"

#ifndef SAFE_DELETE
#define SAFE_DELETE(P) {if(P) delete[](P);(P)=nullptr;}
#endif

namespace sgm_util
{
	//······ census工具集
	// census变换

	/**
	 * \brief census变换
	 * \param source	输入，影像数据
	 * \param census	输出，census值数组
	 * \param width		输入，影像宽
	 * \param height	输入，影像高
	 */
	//void census_transform_5x5(const uint8* source, uint32* census, const sint32& width, const sint32& height);
	void census_transform_9x7(const uint8* source, uint64* census, const sint32& width, const sint32& height, const sint32 wnd_size);
	// Hamming距离
	uint8 Hamming32(const uint32& x, const uint32& y);
	uint8 Hamming64(const uint64& x, const uint64& y);

	PGradient GetGradient(const PGradient* grad_data, const sint32& width_, const sint32& x, const sint32& y);
	
	void MedianFilter(const float32* in, float32* out, const sint32& width, const sint32& height, const sint32 wnd_size);

	void AdaptiveMedianFilter(const float32* in, float32* out, const sint32& width, const sint32& height,
		sint32 wnd_size, int maxSize);

}