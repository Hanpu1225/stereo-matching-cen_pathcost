/* -*-c++-*- SemiGlobalMatching - Copyright (C) 2020.
* Author	: Yingsong Li(Ethan Li) <ethan.li.whu@gmail.com>
* Describe	: main
*/

#include "stdafx.h"
#include "SemiGlobalMatching.h"
#include <chrono>
using namespace std::chrono;

// opencv library
#include <opencv2/opencv.hpp>
#ifdef _DEBUG
#pragma comment(lib,"opencv_world310d.lib")
#else
#pragma comment(lib,"opencv_world310.lib")
#endif

/**
 * \brief
 * \param argv 3
 * \param argc argc[1]:左影像路径 argc[2]: 右影像路径 argc[3]: 最小视差[可选，默认0] argc[4]: 最大视差[可选，默认64]
 * \param eg. ..\Data\cone\im2.png ..\Data\cone\im6.png 0 64
 * \param eg. ..\Data\Reindeer\view1.png ..\Data\Reindeer\view5.png 0 128
 * \return
 */
int main(int argv, char** argc)
{
    if (argv < 3) {
        std::cout << "参数过少，请至少指定左右影像路径！" << std::endl;
        return -1;
    }

    //···············································································//
    // 读取影像
    std::string path_left = argc[1];
    std::string path_right = argc[2];

    cv::Mat img_left_c = cv::imread(path_left, cv::IMREAD_COLOR);
    cv::Mat img_left = cv::imread(path_left, cv::IMREAD_COLOR);
    cv::Mat img_right = cv::imread(path_right, cv::IMREAD_COLOR);

    if (img_left.data == nullptr || img_right.data == nullptr) {
        std::cout << "读取影像失败！" << std::endl;
        return -1;
    }
    if (img_left.rows != img_right.rows || img_left.cols != img_right.cols) {
        std::cout << "左右影像尺寸不一致！" << std::endl;
        return -1;
    }


    //···············································································//
    const sint32 width = static_cast<uint32>(img_left.cols);
    const sint32 height = static_cast<uint32>(img_right.rows);

    // 左右影像的彩色数据
    auto bytes_left = new uint8[width * height * 3];
    auto bytes_right = new uint8[width * height * 3];
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            bytes_left[i * 3 * width + 3 * j] = img_left.at<cv::Vec3b>(i, j)[0];
            bytes_left[i * 3 * width + 3 * j + 1] = img_left.at<cv::Vec3b>(i, j)[1];
            bytes_left[i * 3 * width + 3 * j + 2] = img_left.at<cv::Vec3b>(i, j)[2];
            bytes_right[i * 3 * width + 3 * j] = img_right.at<cv::Vec3b>(i, j)[0];
            bytes_right[i * 3 * width + 3 * j + 1] = img_right.at<cv::Vec3b>(i, j)[1];
            bytes_right[i * 3 * width + 3 * j + 2] = img_right.at<cv::Vec3b>(i, j)[2];
        }
    }

    printf("Loading Views...Done!\n");

    // SGM匹配参数设计
    SemiGlobalMatching::SGMOption sgm_option;
    // 聚合路径数
    sgm_option.num_paths = 8;
    // 候选视差范围
    sgm_option.min_disparity = argv < 4 ? 0 : atoi(argc[3]);
    sgm_option.max_disparity = argv < 5 ? 64 : atoi(argc[4]);
    // census窗口类型
    sgm_option.census_size = SemiGlobalMatching::Census9x7;

    printf("w = %d, h = %d, d = [%d,%d]\n\n", width, height, sgm_option.min_disparity, sgm_option.max_disparity);

    // 定义SGM匹配类实例
    SemiGlobalMatching sgm;

    //···············································································//
    // 初始化
	printf("SGM Initializing...\n");
    auto start = std::chrono::steady_clock::now();
    if (!sgm.Initialize(width, height, sgm_option)) {
        std::cout << "SGM初始化失败！" << std::endl;
        return -2;
    }
    auto end = std::chrono::steady_clock::now();
    auto tt = duration_cast<std::chrono::milliseconds>(end - start);
    printf("SGM Initializing Done! Timing : %lf s\n\n", tt.count() / 1000.0);

    //···············································································//
    // 匹配
	printf("SGM Matching...\n");
    start = std::chrono::steady_clock::now();
    // disparity数组保存子像素的视差结果
    auto disparity = new float32[uint32(width * height)]();
    if (!sgm.Match(bytes_left, bytes_right, disparity)) {
        std::cout << "SGM匹配失败！" << std::endl;
        return -2;
    }
    end = std::chrono::steady_clock::now();
    tt = duration_cast<std::chrono::milliseconds>(end - start);
    printf("\nSGM Matching...Done! Timing :   %lf s\n", tt.count() / 1000.0);

    //···············································································//
	// 显示视差图
    // 注意，计算点云不能用disp_mat的数据，它是用来显示和保存结果用的。计算点云要用上面的disparity数组里的数据，是子像素浮点数
    cv::Mat disp_mat = cv::Mat(height, width, CV_8UC1);
    float min_disp = width, max_disp = -width;
    for (sint32 i = 0; i < height; i++) {
        for (sint32 j = 0; j < width; j++) {
            const float32 disp = disparity[i * width + j];
            if (disp != Invalid_Float) {
                min_disp = std::min(min_disp, disp);
                max_disp = std::max(max_disp, disp);
            }
        }
    }
    for (sint32 i = 0; i < height; i++) {
        for (sint32 j = 0; j < width; j++) {
            const float32 disp = disparity[i * width + j];
            if (disp == Invalid_Float) {
                disp_mat.data[i * width + j] = 0;
            }
            else {
                disp_mat.data[i * width + j] = static_cast<uchar>((disp - min_disp) / (max_disp - min_disp) * 255);//static_cast<uchar>((disp - min_disp) / (max_disp - min_disp) * 255)
            }
        }
    }

    cv::imshow("视差图", disp_mat);
    cv::Mat disp_color;
    applyColorMap(disp_mat, disp_color, cv::COLORMAP_JET);
    cv::imshow("视差图-伪彩", disp_color);

    // 保存结果
    std::string disp_map_path = argc[1]; disp_map_path += ".d.png";
    std::string disp_color_map_path = argc[1]; disp_color_map_path += ".c.png";
    cv::imwrite(disp_map_path, disp_mat);
    cv::imwrite(disp_color_map_path, disp_color);


    cv::waitKey(0);

    //···············································································//
    // 释放内存
    delete[] disparity;
    disparity = nullptr;
    delete[] bytes_left;
    bytes_left = nullptr;
    delete[] bytes_right;
    bytes_right = nullptr;

    system("pause");
    return 0;
}

