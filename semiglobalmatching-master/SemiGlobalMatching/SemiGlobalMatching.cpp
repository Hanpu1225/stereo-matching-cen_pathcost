/* -*-c++-*- SemiGlobalMatching - Copyright (C) 2020.
* Author	: Yingsong Li(Ethan Li) <ethan.li.whu@gmail.com>
* https://github.com/ethan-li-coding/SemiGlobalMatching
* Describe	: implement of semi-global matching class
*/

#include "stdafx.h"
#include "SemiGlobalMatching.h"
#include "sgm_util.h"
#include <algorithm>
#include <vector>
#include <cassert>
#include <chrono>
using namespace std::chrono;

SemiGlobalMatching::SemiGlobalMatching(): width_(0), height_(0), img_left_(nullptr), img_right_(nullptr),
                                          gray_left_(nullptr), gray_right_(nullptr),
                                          grad_left_(nullptr), grad_right_(nullptr),
                                          census_left_(nullptr), census_right_(nullptr),
                                          cost_init_(nullptr), cost_aggr_(nullptr),                                       
                                          disp_left_(nullptr), disp_right_(nullptr),
                                          is_initialized_(false)
{
}


SemiGlobalMatching::~SemiGlobalMatching()
{
    Release();
    is_initialized_ = false;
}

bool SemiGlobalMatching::Initialize(const sint32& width, const sint32& height, const SGMOption& option)
{
    // ··· 赋值
    
	// 影像尺寸
    width_ = width;
    height_ = height;
    // SGM参数
    option_ = option;

    if(width == 0 || height == 0) {
        return false;
    }

    //··· 开辟内存空间

    // census值（左右影像）
    const sint32 img_size = width * height;

    gray_left_ = new uint8[img_size]();
    gray_right_ = new uint8[img_size]();

    // 梯度数据
    grad_left_ = new PGradient[img_size]();
    grad_right_ = new PGradient[img_size]();

    
    census_left_ = new uint64[img_size]();
    census_right_ = new uint64[img_size]();
   

    // 视差范围
    const sint32 disp_range = option.max_disparity - option.min_disparity;
    if (disp_range <= 0) {
        return false;
    }

    // 匹配代价（初始/聚合）
    const sint32 size = width * height * disp_range;
    cost_init_  = new float32[size]();
    cost_aggr_  = new float32[size]();


    // 视差图
    disp_left_ = new float32[img_size]();
    disp_right_ = new float32[img_size]();

    is_initialized_ =gray_left_&&gray_right_&& grad_left_ && grad_right_&& census_left_ && census_right_ && cost_init_ && cost_aggr_ && disp_left_;

    return is_initialized_;
}


void SemiGlobalMatching::Release()
{
    // 释放内存
    SAFE_DELETE(gray_left_);
    SAFE_DELETE(gray_right_);
    SAFE_DELETE(grad_left_);
    SAFE_DELETE(grad_right_);

    SAFE_DELETE(census_left_);
    SAFE_DELETE(census_right_);
    SAFE_DELETE(cost_init_);
    SAFE_DELETE(cost_aggr_);
    SAFE_DELETE(disp_left_);
    SAFE_DELETE(disp_right_);
}

bool SemiGlobalMatching::Match(const uint8* img_left, const uint8* img_right, float32* disp_left)
{
    if(!is_initialized_) {
        return false;
    }
    if (img_left == nullptr || img_right == nullptr) {
        return false;
    }

    img_left_ = img_left;
    img_right_ = img_right;

    auto start = std::chrono::steady_clock::now();

    ComputeGray();

    ComputeGradient();
    // census变换
    CensusTransform();

    // 代价计算
    ComputeCost();

    auto end = steady_clock::now();
    auto tt = duration_cast<milliseconds>(end - start);
    printf("computing cost! timing :	%lf s\n", tt.count() / 1000.0);
    start = steady_clock::now();

    // 视差计算
    ComputeDisparity();

    end = steady_clock::now();
    tt = duration_cast<milliseconds>(end - start);
    printf("computing disparities! timing :	%lf s\n", tt.count() / 1000.0);
    start = steady_clock::now();

    // 中值滤波
     sgm_util::MedianFilter(disp_left_, disp_left_, width_, height_, 3);
     //自适应中值滤波
   //sgm_util::AdaptiveMedianFilter(disp_left_, disp_left_, width_, height_, 3, 7);

    // 输出视差图
    memcpy(disp_left, disp_left_, height_ * width_ * sizeof(float32));

	return true;
}

bool SemiGlobalMatching::Reset(const uint32& width, const uint32& height, const SGMOption& option)
{
    // 释放内存
    Release();

    // 重置初始化标记
    is_initialized_ = false;

    // 初始化
    return Initialize(width, height, option);
}

void SemiGlobalMatching::ComputeGray()
{
    // 彩色转灰度
    for (sint32 n = 0; n < 2; n++) {
        const auto color = (n == 0) ? img_left_ : img_right_;
        auto& gray = (n == 0) ? gray_left_ : gray_right_;
        for (sint32 y = 0; y < height_; y++) {
            for (sint32 x = 0; x < width_; x++) {
                const auto b = color[y * width_ * 3 + 3 * x];
                const auto g = color[y * width_ * 3 + 3 * x + 1];
                const auto r = color[y * width_ * 3 + 3 * x + 2];
                gray[y * width_ + x] = uint8(r * 0.299 + g * 0.587 + b * 0.114);
            }
        }
    }
}

void SemiGlobalMatching::ComputeGradient()
{
    const sint32 width = width_;
    const sint32 height = height_;
    if (width <= 0 || height <= 0 ||
        grad_left_ == nullptr || grad_right_ == nullptr ||
        gray_left_ == nullptr || gray_right_ == nullptr) {
        return;
    }

    // Sobel梯度算子
    for (sint32 n = 0; n < 2; n++) {
        auto* gray = (n == 0) ? gray_left_ : gray_right_;
        auto* grad = (n == 0) ? grad_left_ : grad_right_;
        for (int y = 1; y < height - 1; y++) {
            for (int x = 1; x < width - 1; x++) {
                const auto grad_x = (-gray[(y - 1) * width + x - 1] + gray[(y - 1) * width + x + 1]) +
                    (-2 * gray[y * width + x - 1] + 2 * gray[y * width + x + 1]) +
                    (-gray[(y + 1) * width + x - 1] + gray[(y + 1) * width + x + 1]);
                const auto grad_y = (-gray[(y - 1) * width + x - 1] - 2 * gray[(y - 1) * width + x] - gray[(y - 1) * width + x + 1]) +
                    (gray[(y + 1) * width + x - 1] + 2 * gray[(y + 1) * width + x] + gray[(y + 1) * width + x + 1]);
                grad[y * width + x].x = grad_x / 8;
                grad[y * width + x].y = grad_y / 8;
            }
        }
    }
}

void SemiGlobalMatching::CensusTransform() const
{
	// 左右影像census变换   
   /* sgm_util::census_transform_5x5(gray_left_, static_cast<uint32*>(census_left_), width_, height_);
    sgm_util::census_transform_5x5(gray_right_, static_cast<uint32*>(census_right_), width_, height_);*/

    sgm_util::census_transform_9x7(gray_left_, static_cast<uint64*>(census_left_), width_, height_, 9);
    sgm_util::census_transform_9x7(gray_right_, static_cast<uint64*>(census_right_), width_, height_, 9);
    
}

void SemiGlobalMatching::ComputeCost() const
{
    const sint32& min_disparity = option_.min_disparity;
    const sint32& max_disparity = option_.max_disparity;
    const sint32 disp_range = max_disparity - min_disparity;
    if (disp_range <= 0) {
        return;
    }

    const sint32& lambda_ad = option_.lambda_ad;
    const sint32& lambda_census = option_.lambda_census;
    const sint32& lambda_dc = option_.lambda_dc;

    for (sint32 y = 0; y < height_; y++)
    {
        for (sint32 x = 0; x < width_; x++)
        {

            const auto bl = img_left_[y * width_ * 3 + 3 * x];
            const auto gl = img_left_[y * width_ * 3 + 3 * x + 1];
            const auto rl = img_left_[y * width_ * 3 + 3 * x + 2];

            const auto& grad_p = sgm_util::GetGradient(grad_left_, width_, x, y);

            // 逐视差计算代价值
            for (sint32 d = min_disparity; d < max_disparity; d++) {
                auto& cost = cost_init_[y * width_ * disp_range + x * disp_range + (d - min_disparity)];
                const sint32 xr = x - d;
                if (x - d < 0 || x - d >= width_) {
                    cost = UINT8_MAX;//1.0f
                    continue;
                }
                // ad代价
                const auto br = img_right_[y * width_ * 3 + 3 * xr];
                const auto gr = img_right_[y * width_ * 3 + 3 * xr + 1];
                const auto rr = img_right_[y * width_ * 3 + 3 * xr + 2];
                const float32 cost_ad = (abs(bl - br) + abs(gl - gr) + abs(rl - rr)) / 3.0f;

                const auto& grad_q = sgm_util::GetGradient(grad_right_, width_, xr, y);
                const float32 cost_dc = (abs(grad_p.x - grad_q.x)) + (abs(grad_p.y - grad_q.y));

                const auto& census_val_l = static_cast<uint64*>(census_left_)[y * width_ + x];
                const auto& census_val_r = static_cast<uint64*>(census_right_)[y * width_ + xr];
                const float32 cost_census = sgm_util::Hamming64(census_val_l, census_val_r);

                //cost = sgm_util::Hamming32(census_val_l, census_val_r);

                
               cost = 1 - exp(-cost_ad / lambda_ad) + 1 - exp(-cost_census / lambda_census) + 1 - exp(-cost_dc / lambda_dc);

                //cost = 1 - exp(-cost_ad / lambda_ad) + 1 - exp(-cost_census / lambda_census);
                //cost = 1 - exp(-cost_dc / lambda_dc) + 1 - exp(-cost_census / lambda_census);
            }
        }
    } 


}

void SemiGlobalMatching::ComputeDisparity() const
{
    const sint32& min_disparity = option_.min_disparity;
    const sint32& max_disparity = option_.max_disparity;
    const sint32 disp_range = max_disparity - min_disparity;
    if (disp_range <= 0) {
        return;
    }

    // 左影像视差图
    const auto disparity = disp_left_;
    // 左影像聚合代价数组
    const auto cost_ptr = cost_init_;

    const sint32 width = width_;
    const sint32 height = height_;
    const bool is_check_unique = option_.is_check_unique;
    const float32 uniqueness_ratio = option_.uniqueness_ratio;

    // 为了加快读取效率，把单个像素的所有代价值存储到局部数组里
    std::vector<float32> cost_local(disp_range);

    // ---逐像素计算最优视差
    for (sint32 i = 0; i < height; i++) {
        for (sint32 j = 0; j < width; j++) {
            float32 min_cost = Large_Float;
            float32 sec_min_cost = Large_Float;
            sint32 best_disparity = 0;

            // ---遍历视差范围内的所有代价值，输出最小代价值及对应的视差值
            for (sint32 d = min_disparity; d < max_disparity; d++) {
                const sint32 d_idx = d - min_disparity;
                const auto& cost = cost_local[d_idx] = cost_ptr[i * width * disp_range + j * disp_range + d_idx];
                if (min_cost > cost) {
                    min_cost = cost;
                    best_disparity = d;
                }
            }

            if (is_check_unique) {
                // 再遍历一次，输出次最小代价值
                for (sint32 d = min_disparity; d < max_disparity; d++) {
                    if (d == best_disparity) {
                        // 跳过最小代价值
                        continue;
                    }
                    const auto& cost = cost_local[d - min_disparity];
                    sec_min_cost = std::min(sec_min_cost, cost);
                }

                // 判断唯一性约束
                // 若(min-sec)/min < min*(1-uniquness)，则为无效估计
                if (sec_min_cost - min_cost <= static_cast<float32>(min_cost * (1 - uniqueness_ratio))) {
                    disparity[i * width + j] = Invalid_Float;
                    continue;
                }
            }

            // ---子像素拟合
            if (best_disparity == min_disparity || best_disparity == max_disparity - 1) {
                disparity[i * width + j] = Invalid_Float;
                continue;
            }
            // 最优视差前一个视差的代价值cost_1，后一个视差的代价值cost_2
            const sint32 idx_1 = best_disparity - 1 - min_disparity;
            const sint32 idx_2 = best_disparity + 1 - min_disparity;
            const float32 cost_1 = cost_local[idx_1];
            const float32 cost_2 = cost_local[idx_2];
            // 解一元二次曲线极值
            const float32 denom = cost_1 + cost_2 - 2 * min_cost;
            disparity[i * width + j] = static_cast<float32>(best_disparity) + static_cast<float32>(cost_1 - cost_2) / (denom * 2.0f);
        }
    }
}

void SemiGlobalMatching::ComputeDisparityRight() const
{
    const sint32& min_disparity = option_.min_disparity;
    const sint32& max_disparity = option_.max_disparity;
    const sint32 disp_range = max_disparity - min_disparity;
    if (disp_range <= 0) {
        return;
    }

    // 右影像视差图
    const auto disparity = disp_right_;
    // 左影像聚合代价数组
    const auto cost_ptr = cost_init_;

    const sint32 width = width_;
    const sint32 height = height_;
    const bool is_check_unique = option_.is_check_unique;
    const float32 uniqueness_ratio = option_.uniqueness_ratio;

    // 为了加快读取效率，把单个像素的所有代价值存储到局部数组里
    std::vector<float32> cost_local(disp_range);

    // ---逐像素计算最优视差
    // 通过左影像的代价，获取右影像的代价
    // 右cost(xr,yr,d) = 左cost(xr+d,yl,d)
    for (sint32 i = 0; i < height; i++) {
        for (sint32 j = 0; j < width; j++) {
            float32 min_cost = Large_Float;
            float32 sec_min_cost = Large_Float;
            sint32 best_disparity = 0;

            // ---统计候选视差下的代价值
            for (sint32 d = min_disparity; d < max_disparity; d++) {
                const sint32 d_idx = d - min_disparity;
                const sint32 col_left = j + d;
                if (col_left >= 0 && col_left < width) {
                    const auto& cost = cost_local[d_idx] = cost_ptr[i * width * disp_range + col_left * disp_range + d_idx];
                    if (min_cost > cost) {
                        min_cost = cost;
                        best_disparity = d;
                    }
                }
                else {
                    cost_local[d_idx] = Large_Float;
                }
            }

            if (is_check_unique) {
                // 再遍历一次，输出次最小代价值
                for (sint32 d = min_disparity; d < max_disparity; d++) {
                    if (d == best_disparity) {
                        // 跳过最小代价值
                        continue;
                    }
                    const auto& cost = cost_local[d - min_disparity];
                    sec_min_cost = std::min(sec_min_cost, cost);
                }

                // 判断唯一性约束
                // 若(min-sec)/min < min*(1-uniquness)，则为无效估计
                if (sec_min_cost - min_cost <= static_cast<float32>(min_cost * (1 - uniqueness_ratio))) {
                    disparity[i * width + j] = Invalid_Float;
                    continue;
                }
            }

            // ---子像素拟合
            if (best_disparity == min_disparity || best_disparity == max_disparity - 1) {
                disparity[i * width + j] = Invalid_Float;
                continue;
            }

            // 最优视差前一个视差的代价值cost_1，后一个视差的代价值cost_2
            const sint32 idx_1 = best_disparity - 1 - min_disparity;
            const sint32 idx_2 = best_disparity + 1 - min_disparity;
            const float32 cost_1 = cost_local[idx_1];
            const float32 cost_2 = cost_local[idx_2];
            // 解一元二次曲线极值
            const float32 denom = cost_1 + cost_2 - 2 * min_cost;
            disparity[i * width + j] = static_cast<float32>(best_disparity) + static_cast<float32>(cost_1 - cost_2) / (denom * 2.0f);
        }
    }
}
