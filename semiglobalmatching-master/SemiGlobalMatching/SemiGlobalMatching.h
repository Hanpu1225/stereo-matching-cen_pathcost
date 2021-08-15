/* -*-c++-*- SemiGlobalMatching - Copyright (C) 2020.
* Author	: Yingsong Li(Ethan Li) <ethan.li.whu@gmail.com>
* https://github.com/ethan-li-coding/SemiGlobalMatching
* Describe	: header of semi-global matching class
*/

#pragma once

#include "sgm_types.h"
#include <vector>

/**
 * \brief SemiGlobalMatching�ࣨGeneral implementation of Semi-Global Matching��
 */
class SemiGlobalMatching
{
public:
	SemiGlobalMatching();
	~SemiGlobalMatching();


	/** \brief Census���ڳߴ����� */
	enum CensusSize {
		Census5x5 = 0,
		Census9x7
	};

	/** \brief SGM�����ṹ�� */
	struct SGMOption {
		uint8	num_paths;			// �ۺ�·���� 4 and 8
		sint32  min_disparity;		// ��С�Ӳ�
		sint32	max_disparity;		// ����Ӳ�

		sint32	lambda_ad;			// ����AD����ֵ�Ĳ���
		sint32	lambda_census;		// ����Census����ֵ�Ĳ���
		sint32  lambda_dc;

		CensusSize census_size;		// census���ڳߴ�

		bool	is_check_unique;	// �Ƿ���Ψһ��
		float32	uniqueness_ratio;	// Ψһ��Լ����ֵ ����С����-����С����)/��С���� > ��ֵ Ϊ��Ч����

		bool	is_check_lr;		// �Ƿ�������һ����
		float32	lrcheck_thres;		// ����һ����Լ����ֵ

		bool	is_remove_speckles;	// �Ƿ��Ƴ�С����ͨ��
		int		min_speckle_aera;	// ��С����ͨ���������������

		bool	is_fill_holes;		// �Ƿ�����Ӳ�ն�

		// P1,P2 
		// P2 = P2_init / (Ip-Iq)
		sint32  p1;				// �ͷ������P1
		sint32  p2_init;		// �ͷ������P2

		SGMOption(): num_paths(8), min_disparity(0), max_disparity(64), census_size(Census5x5), 
			/*
			*�Լ��Ĳ�����42 29 15
			* AD+Census��AD:10 Census:30
			*Census+Grad:�Ҷ�25 Census15
			*
			*/
			         lambda_ad(42), lambda_census(29), lambda_dc(15),
		             is_check_unique(true), uniqueness_ratio(0.95f),
		             is_check_lr(true), lrcheck_thres(1.0f),
		             is_remove_speckles(true), min_speckle_aera(20),
		             is_fill_holes(true),
		             p1(10), p2_init(150) { }
	};
public:
	/**
	 * \brief ��ĳ�ʼ�������һЩ�ڴ��Ԥ���䡢������Ԥ���õ�
	 * \param width		���룬�������Ӱ���
	 * \param height	���룬�������Ӱ���
	 * \param option	���룬SemiGlobalMatching����
	 */
	bool Initialize(const sint32& width, const sint32& height, const SGMOption& option);

	/**
	 * \brief ִ��ƥ��
	 * \param img_left	���룬��Ӱ������ָ�� 
	 * \param img_right	���룬��Ӱ������ָ��
	 * \param disp_left	�������Ӱ���Ӳ�ͼָ�룬Ԥ�ȷ����Ӱ��ȳߴ���ڴ�ռ�
	 */
	bool Match(const uint8* img_left, const uint8* img_right, float32* disp_left);

	/**
	 * \brief ����
	 * \param width		���룬�������Ӱ���
	 * \param height	���룬�������Ӱ���
	 * \param option	���룬SemiGlobalMatching����
	 */
	bool Reset(const uint32& width, const uint32& height, const SGMOption& option);

private:

	/** \brief ����Ҷ����� */
	void ComputeGray();

	/** \brief �����ݶ����� */
	void ComputeGradient();

	/** \brief Census�任 */
	void CensusTransform() const;

	/** \brief ���ۼ���	 */
	void ComputeCost() const;


	/** \brief �Ӳ����	 */
	void ComputeDisparity() const;

	/** \brief �Ӳ����	 */
	void ComputeDisparityRight() const;


	/** \brief �ڴ��ͷ�	 */
	void Release();

private:
	/** \brief SGM����	 */
	SGMOption option_;

	/** \brief Ӱ���	 */
	sint32 width_;

	/** \brief Ӱ���	 */
	sint32 height_;

	/** \brief ��Ӱ������	 */
	const uint8* img_left_;

	/** \brief ��Ӱ������	 */
	const uint8* img_right_;
	
	/** \brief ��Ӱ��Ҷ�����	 */
	uint8* gray_left_;
	/** \brief ��Ӱ��Ҷ�����	 */
	uint8* gray_right_;

	/** \brief ��Ӱ���ݶ�����	 */
	PGradient* grad_left_;
	/** \brief ��Ӱ���ݶ�����	 */
	PGradient* grad_right_;

	/** \brief ��Ӱ��censusֵ	*/
	void* census_left_;
	
	/** \brief ��Ӱ��censusֵ	*/
	void* census_right_;


	
	/** \brief ��ʼƥ�����	*/
	float32* cost_init_;
	
	/** \brief �ۺ�ƥ�����	*/
	float32* cost_aggr_;


	/** \brief ��Ӱ���Ӳ�ͼ	*/
	float32* disp_left_;
	/** \brief ��Ӱ���Ӳ�ͼ	*/
	float32* disp_right_;

	/** \brief �Ƿ��ʼ����־	*/
	bool is_initialized_;

	/** \brief �ڵ������ؼ�	*/
	std::vector<std::pair<int, int>> occlusions_;
	/** \brief ��ƥ�������ؼ�	*/
	std::vector<std::pair<int, int>> mismatches_;
};

