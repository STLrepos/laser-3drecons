#ifndef __LASER3DRECONS_HPP__ 
#define __LASER3DRECONS_HPP__ 

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include <fstream>

#include <cv.h>
#include <highgui.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>
#include <opencv2/core/eigen.hpp>


class Laser3DReconst 
{
public:
  
      Laser3DReconst();
      
      ~Laser3DReconst();
      
      int init();
      
      int run();

private:
      
      Eigen::Matrix3d Kl, Kr;
      Eigen::Matrix3d F; 
      Eigen::MatrixXd Tl, Tr, Pl, Pr; 
      Eigen::MatrixXd c2Tc1, c2Rc1, c2tc1;
      
      // the Least-square equ vars
      Eigen::MatrixXd A; 
      Eigen::VectorXd b;
      Eigen::Vector3d r1,r2,r3,k1,k2,k3,t; 
      Eigen::Vector4d p1l,p2l,p3l,p1r,p2r,p3r;
      
      //int linethinning(const cv::Mat& image, cv::Mat& edge);
      int linethinning(cv::Mat& img, cv::Mat& edge);
      int findCorrespondence(const Eigen::Vector2d& ptL, const Eigen::MatrixXd& locs, Eigen::Vector2d& ptR, Eigen::VectorXi& flags);
      int trianglulate(const Eigen::Vector2d& ptL, const Eigen::Vector2d& ptR, Eigen::Vector3d& pt3);
      int reprojCheck(const Eigen::Vector2d& ptL, const Eigen::Vector2d& ptR, const Eigen::Vector3d& pt3);
      int drawEpiLines(cv::Mat& cR, const Eigen::Vector2d& pL, const Eigen::Matrix3d& F, const cv::Scalar& color);
      int drawEpiLinesLelf(cv::Mat& cL, const Eigen::Vector2d& pR, const Eigen::Matrix3d& F, const cv::Scalar& color);
};
#endif