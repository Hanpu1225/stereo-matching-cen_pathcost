/* -*-c++-*- SemiGlobalMatching - Copyright (C) 2020.
* Author	: Yingsong Li(Ethan Li) <ethan.li.whu@gmail.com>
* https://github.com/ethan-li-coding/SemiGlobalMatching
* Describe	: header of sgm_types
*/

#pragma once

#include <cstdint>
#include <limits>

/** \brief float��Чֵ */
constexpr auto Invalid_Float = std::numeric_limits<float>::infinity();

constexpr auto Large_Float = 99999.0f;
constexpr auto Small_Float = -99999.0f;
/** \brief �������ͱ��� */
typedef int8_t			sint8;		// �з���8λ����
typedef uint8_t			uint8;		// �޷���8λ����
typedef int16_t			sint16;		// �з���16λ����
typedef uint16_t		uint16;		// �޷���16λ����
typedef int32_t			sint32;		// �з���32λ����
typedef uint32_t		uint32;		// �޷���32λ����
typedef int64_t			sint64;		// �з���64λ����
typedef uint64_t		uint64;		// �޷���64λ����
typedef float			float32;	// �����ȸ���
typedef double			float64;	// ˫���ȸ���


/**
* \brief ��ɫ�ṹ��
*/
struct ADColor {
	uint8 r, g, b;
	ADColor() : r(0), g(0), b(0) {}
	ADColor(uint8 _b, uint8 _g, uint8 _r) {
		r = _r; g = _g; b = _b;
	}
};

/**
 * \brief �ݶȽṹ��
 */
struct PGradient {
	sint16 x, y;
	PGradient() : x(0), y(0) {}
	PGradient(sint16 _x, sint16 _y) {
		x = _x; y = _y;
	}
};