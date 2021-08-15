/* -*-c++-*- SemiGlobalMatching - Copyright (C) 2020.
* Author	: Yingsong Li(Ethan Li) <ethan.li.whu@gmail.com>
* https://github.com/ethan-li-coding/SemiGlobalMatching
* Describe	: header of sgm_types
*/

#pragma once

#include <cstdint>
#include <limits>

/** \brief float无效值 */
constexpr auto Invalid_Float = std::numeric_limits<float>::infinity();

constexpr auto Large_Float = 99999.0f;
constexpr auto Small_Float = -99999.0f;
/** \brief 基础类型别名 */
typedef int8_t			sint8;		// 有符号8位整数
typedef uint8_t			uint8;		// 无符号8位整数
typedef int16_t			sint16;		// 有符号16位整数
typedef uint16_t		uint16;		// 无符号16位整数
typedef int32_t			sint32;		// 有符号32位整数
typedef uint32_t		uint32;		// 无符号32位整数
typedef int64_t			sint64;		// 有符号64位整数
typedef uint64_t		uint64;		// 无符号64位整数
typedef float			float32;	// 单精度浮点
typedef double			float64;	// 双精度浮点


/**
* \brief 颜色结构体
*/
struct ADColor {
	uint8 r, g, b;
	ADColor() : r(0), g(0), b(0) {}
	ADColor(uint8 _b, uint8 _g, uint8 _r) {
		r = _r; g = _g; b = _b;
	}
};

/**
 * \brief 梯度结构体
 */
struct PGradient {
	sint16 x, y;
	PGradient() : x(0), y(0) {}
	PGradient(sint16 _x, sint16 _y) {
		x = _x; y = _y;
	}
};