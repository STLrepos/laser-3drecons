#include "Laser3DReconst.hpp"

using namespace cv;
using namespace std;
using namespace Eigen;


Laser3DReconst::Laser3DReconst()
{
  
}
      
      
      
Laser3DReconst::~Laser3DReconst()
{

}




int Laser3DReconst::linethinning(cv::Mat& img, cv::Mat& edge)
{
  int ret = 0;
  
  cvtColor(img, img, CV_BGR2GRAY);
  
  threshold(img, img, 127, 255, cv::THRESH_BINARY);
//   imshow("thresholded", img);
//   waitKey(0);
  
  Mat skel(img.size(), CV_8UC1, cv::Scalar(0));
  Mat temp(img.size(), CV_8UC1);
  
  Mat element = cv::getStructuringElement(cv::MORPH_CROSS, cv::Size(3, 3));
   
  bool done;
  do
  {
    cv::morphologyEx(img, temp, cv::MORPH_OPEN, element);
    cv::bitwise_not(temp, temp);
    cv::bitwise_and(img, temp, temp);
    cv::bitwise_or(skel, temp, skel);
    cv::erode(img, img, element);
    
    double max;
    cv::minMaxLoc(img, 0, &max);
    done = (max == 0);
  } while (!done);
      
//   imshow("thinned edge", skel);
//   waitKey(0);
  
  skel.copyTo(edge);
  
  return ret; 
}     
      
      
      
int Laser3DReconst::reprojCheck(const Vector2d& ptL, const Vector2d& ptR, const Vector3d& pt3)
{
  int ret = 0;
  
  double tol = 2.;   // make it 2 pixel error tolerance!
  
  double errL, errR;
  Vector3d ptR3, ptL3; 
  Vector4d pt34;
  
  pt34 << pt3(0),
          pt3(1), 
          pt3(2), 
          1.;
  
  //cout << "pt34 L =" << endl << pt34 << endl;
    
  Vector3d p3 = Kl*Tl*pt34;
  p3 = p3.array() / p3(2);
  Vector2d p2; 
  p2 << p3(0), 
        p3(1);
        
  //cout << "p2 L =" << endl << p2 << endl;
  //cout << "pt L =" << endl << ptL << endl;
  
  errL = (ptL - p2).norm();
  if (errL > tol) 
  {
    //cout << "reproj err: " << errL << endl;
    return 1;
  }
  
  p3 = Kr*Tr*pt34;
  p3 = p3.array() / p3(2);
  
  p2 << p3(0), 
        p3(1);
        
  //cout << "p2 R =" << endl << p2 << endl;
  //cout << "pt R =" << endl << ptR << endl;
  
  errR = (ptR - p2).norm();
  if (errR > tol) 
  {
    //cout << "reproj err: " << errR << endl;
    return 1;
  }
  
  return ret;
}
     
     
     
int Laser3DReconst::trianglulate(const Vector2d& ptL, const Vector2d& ptR, Vector3d& pt3)
{
  int ret = 0;
  
  Matrix<double,4,4> AA;
  
  AA.row(0) = ptL(0)*p3l.transpose() - p1l.transpose();
  AA.row(1) = ptL(1)*p3l.transpose() - p2l.transpose();
  AA.row(2) = ptR(0)*p3r.transpose() - p1r.transpose();
  AA.row(3) = ptR(1)*p3r.transpose() - p2r.transpose();

  // solve SVD of P for b, and then K 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(AA, Eigen::ComputeFullV); 
  
  Eigen::MatrixXd V = svd.matrixV(); 
  MatrixXd X = V.col(V.cols()-1); 
  
  MatrixXd zeroVec(4,1);
  zeroVec.fill(0.);
//   cout << "AA = " << endl << AA << endl;
//   cout << "zero vec = " << endl << zeroVec << endl;
//   cout << "X = " << endl << X << endl;
//   cout << "AA*x = " << endl << AA*X << endl; 
  
  if ((AA*X).isApprox(zeroVec))  
  {
//     cout << "Solution of A1x=b1 is GOOD!!" << endl;
  }// else
//    cout << "A1x=b1 does not have any solution!!" << endl;
           
  X = X.array() / X(3);
    
//   cout << "X after " << endl << X << endl; 
  pt3 = X.block(0,0,3,1);
 
  return ret;
}     
     
      
int Laser3DReconst::drawEpiLinesLelf(Mat& cL, const Vector2d& pR, const Matrix3d& F, const Scalar& color)
{
  int ret = 0;  
  
  Vector3d ptR3;
  ptR3(0) = pR(0);
  ptR3(1) = pR(1);
  ptR3(2) = 1.;
  Vector3d ptRF = ptR3.transpose()*F;
  
  // Calculate 2 points for the line 
  double a = ptRF(0), b = ptRF(1), c = ptRF(2);
  double x_ub = 0.;
  double x_lb = 719.;
  double y_ub = -(a*x_ub + c)/b;
  double y_lb = -(a*x_lb + c)/b;
  
  line(cL, Point(x_ub, y_ub), Point(x_lb, y_lb), color, 1, 8);
  
  return ret; 
}


int Laser3DReconst::drawEpiLines(Mat& cR, const Vector2d& pL, const Matrix3d& F, const Scalar& color)
{
  int ret = 0;  
  
  Vector3d ptL3;
  ptL3(0) = pL(0);
  ptL3(1) = pL(1);
  ptL3(2) = 1.;
  Vector3d FptL = F*ptL3;
  
  // Calculate 2 points for the line 
  double a = FptL(0), b = FptL(1), c = FptL(2);
  double x_ub = 0.;
  double x_lb = 719.;
  double y_ub = -(a*x_ub + c)/b;
  double y_lb = -(a*x_lb + c)/b;
  
  line(cR, Point(x_ub, y_ub), Point(x_lb, y_lb), color, 1, 8);
  
  return ret; 
}


int Laser3DReconst::run()
{
  int ret = 0;
  
  ofstream myfile("model.xyz");
  if (!myfile.is_open())
  {
    cout << "File cannot be opened!" << endl;
    return 1;
  }

  char nameL[100]=""; 
  char* fnameL = &nameL[0];  
  
  char nameR[100]=""; 
  char* fnameR = &nameR[0];    
  
  int miss;

  double t = (double)getTickCount();
  for(int i=0; i<314; i++)
  {
      Mat frameL, frameR;
      
      // ########## Read frames
      sprintf(fnameL,"./L/L%04d.jpg",i);
      sprintf(fnameR,"./R/R%04d.jpg",i);
      
      frameL = imread(fnameL);
      frameR = imread(fnameR);
  //     imshow("read frame L", frameL);
  //     imshow("read frame R", frameR);
  //     waitKey(1);
      
      // ########## Get thin lines 
      Mat edL, edR; 
      linethinning(frameL, edL);
      linethinning(frameR, edR);
//       imshow("edge frame L", edL);
//       imshow("edge frame R", edR);
//       waitKey(1);
      
      Mat locL, locR;
      findNonZero(edL, locL);
      findNonZero(edR, locR);
      locL.convertTo(locL,CV_64FC1);
      locR.convertTo(locR,CV_64FC1);
      
      VectorXi flags(locR.rows);
      flags.fill(0);
      if (locL.rows != locR.rows) 
      {
//         cout << "locL = " << locL.rows << endl;
//         cout << "locR = " << locR.rows << endl;
//         cout << "[Error] correspondence problem!!!" << " at frame " << i <<  endl;
        //return 1; 
      }
      
//       Mat pt = Mat(edL.size(), CV_8UC3, Scalar(0,0,0));
//       Mat ptr = Mat(edL.size(), CV_8UC3, Scalar(0,0,0));
//       for (int i=0; i<locL.rows; i++)
//       {
//         circle(pt, Point(locL.at<Vec2d>(i)[0], locL.at<Vec2d>(i)[1]), 1, Scalar(255,0,255), 1, 8, 0);
//       }
//       for (int i=0; i<locR.rows; i++)
//       {
//         circle(ptr, Point(locR.at<Vec2d>(i)[0], locR.at<Vec2d>(i)[1]), 1, Scalar(255,0,255), 1, 8, 0);
//       }
//       resize(pt, pt, Size(), .75, .75);
//       resize(ptr, ptr, Size(), .75, .75);
//       imshow("loc L", pt);
//       imshow("loc R", ptr);
//       waitKey(0);
        
      Vector2d pL, pR;
      MatrixXd locL_(locL.rows,2), locR_(locR.rows,2);
      for (int i=0; i<locL_.rows(); i++)
      {
        locL_(i,0) = locL.at<Vec2d>(i)[0];
        locL_(i,1) = locL.at<Vec2d>(i)[1];
      }
      for (int i=0; i<locR_.rows(); i++)
      {
        locR_(i,0) = locR.at<Vec2d>(i)[0];
        locR_(i,1) = locR.at<Vec2d>(i)[1];
      }
      
     RNG rng(0); 
     int ans;
     miss = 0;
     for (int i=0; i<locL.rows; i++)
     {
        Mat cL = Mat(edL.size(), CV_8UC3, Scalar(0,0,0));
        Mat cR = Mat(edR.size(), CV_8UC3, Scalar(0,0,0));
//        cout << locL_.row(i) << endl;
       pL = locL_.row(i);
       
       ans = findCorrespondence(pL, locR_, pR, flags); 
       
//        Scalar color(rng(256),rng(256),rng(256));
//        circle(cL, Point(pL(0),pL(1)), 2, Scalar(255,0,255), 2, 8, 0);
//        circle(cR, Point(pR(0),pR(1)), 2, Scalar(255,255,255), 2, 8, 0);
//        drawEpiLines(cR, pL, F, color);
//        drawEpiLinesLelf(cL, pR, F, color);
//        resize(cL, cL, Size(), .75, .75);
//        resize(cR, cR, Size(), .75, .75);
//        imshow("cor L", cL);
//        imshow("cor R", cR);
//        waitKey(0);
       
       if (!ans)  
       {
          Vector3d p3d;
          //trianglulate(pL, pR, p3d);
          trianglulate(pL, pR, p3d);
          if (!reprojCheck(pL, pR, p3d))
          {
            myfile <<  p3d(0) << " " << p3d(1) << " " << p3d(2) << endl; 
          } else
          {
//             cout << "reprojCheck fails!!" << endl; 
          }
       } else
       {
         miss++;
//          cout << "cannot find correspondence!" << endl;
       }
     }
     
//     cout << "total miss points: " << miss << endl;
  }
  
  myfile.close();
  t = ((double)getTickCount() - t)/getTickFrequency();

  cout << "Times passed in seconds: " << t << endl;
  
  return ret;
}
      
      
      
int Laser3DReconst::findCorrespondence(const Vector2d& ptL, const MatrixXd& locs, Vector2d& ptR, VectorXi& flags)
{
  int ret = 0;
  
  Vector2d pt;
  double min = 99999999.;
  
  Vector3d ptL3;
  ptL3(0) = ptL(0);
  ptL3(1) = ptL(1);
  ptL3(2) = 1.;
  Vector3d FptL = F*ptL3;
  
  // Calculate the conditions
  double a = FptL(0), b = FptL(1), c = FptL(2);
  double y_ub = -(a*0. + c)/b;
  double y_lb = -(a*719. + c)/b;
  double ub = std::max(y_ub, y_lb);
  double lb = std::min(y_ub, y_lb);
  
//   cout << "F*ptL" << endl << FptL << endl;
//   getchar();
  int id = 35000;
  for (int i=0; i<locs.rows(); i++)
    if (locs(i,1) <= ub &&  locs(i,1) >= lb && !flags(i))
  {
      Vector2d tmp = locs.row(i);
      Vector3d tmp3;  
      tmp3(0) = tmp(0);
      tmp3(1) = tmp(1);
      tmp3(2) = 1.;
      
      double ptRFptL = tmp3.transpose()*FptL; 
      
      if (abs(ptRFptL) < min)
      {
        min = abs(ptRFptL);
        pt = tmp;
        id = i;
      }
  }
  
  if (id != 35000)
  {
    flags(id) = 1;
  } else ret = 1;
  
  ptR = pt;
  
  return ret; 
}
      
      
      
int Laser3DReconst::init()
{
  int ret = 0;
  
  Kl << 1035.945388223665, 0.000000000000, 431.179652212212,
        0.000000000000, 1036.888690003534, 597.721238839792,
        0.000000000000, 0.000000000000, 1.000000000000;
        
  Kr << 1033.735500089685, 0.000000000000, 437.982707803673,
        0.000000000000, 1034.937552611470, 692.998336913086,
        0.000000000000, 0.000000000000, 1.000000000000;
        
  Tl.setZero(3,4);
  Tl << 1.000000000000, 0.000000000000, 0.000000000000, 0.000000000000,
        0.000000000000, 1.000000000000, 0.000000000000, 0.000000000000,
        0.000000000000, 0.000000000000, 1.000000000000, 0.000000000000;
  
  Tr.setZero(3,4);
  Tr << 0.932698609231, 0.007504048235, 0.360578692657, -73.407836521112,
        0.006026325702, 0.999319672845, -0.036385091230, -0.284060685179,
        -0.360606416660, 0.036109288630, 0.932018847201, 7.321666107008;
  
  F << 0.000000025550, 0.000000679875, -0.000506474386,
       0.000002138945, -0.000000278696, -0.008631002360,
      -0.001390792745, 0.007992747953, 1.000000000000;
      
  Pl = Kl*Tl;
  Pr = Kr*Tr;
  
  p1l = Pl.row(0);
  p2l = Pl.row(1);
  p3l = Pl.row(2);
  
  p1r = Pr.row(0);
  p2r = Pr.row(1);
  p3r = Pr.row(2);
      
  c2tc1 = Tr.col(3);
  r1 = Tr.col(0);
  r2 = Tr.col(1);
  r3 = Tr.col(2);
  
  t = c2tc1;
  r1 = Tr.col(0);
  r2 = Tr.col(1);
  r3 = Tr.col(2);

  k1 = Kr.row(0);
  k2 = Kr.row(1);
  k3 = Kr.row(2);

  return ret;
}

