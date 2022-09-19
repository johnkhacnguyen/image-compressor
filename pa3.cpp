/**
 *  @file pa3.cpp
 *  @description basic test cases for CPSC 221 PA3
 *  @author CPSC 221
 *
 *  You are highly recommended to add your own test cases to this file
 * 
 *  THIS FILE WILL NOT BE SUBMITTED TO PRAIRIELEARN
 */

#include <iostream>

#include "imgtree.h"

using namespace cs221util;
using namespace std;

// Test function declarations
void TestStats();
void TestImgTreeBuildRender();
void TestImgTreeFlipHorizontal();
void TestImgTreeCountLeavesPrune();


// Program entry point
int main(void) {

  // add or uncomment calls to test functions as needed
  TestStats();
  TestImgTreeBuildRender();
  TestImgTreeFlipHorizontal();
  TestImgTreeCountLeavesPrune();

  return 0;
}


// Test function implementations
void TestStats() {
  cout << "Entered TestStats..." << endl;

  cout << "Creating a 3x3 PNG... ";
  // make a 3x3 PNG with all alpha = 1.0, and RGB values as:
  // 12   3   6
  //  8  14   1
  //  9   7  11
  PNG testpng;
  testpng.resize(3, 3);
  
  RGBAPixel* px;
  px = testpng.getPixel(0, 0);
  *px = RGBAPixel(12, 12, 12, 1.0);
  px = testpng.getPixel(1, 0);
  *px = RGBAPixel(3, 3, 3, 1.0);
  px = testpng.getPixel(2, 0);
  *px = RGBAPixel(6, 6, 6, 1.0);
  px = testpng.getPixel(0, 1);
  *px = RGBAPixel(8, 8, 8, 1.0);
  px = testpng.getPixel(1, 1);
  *px = RGBAPixel(14, 14, 14, 1.0);
  px = testpng.getPixel(2, 1);
  *px = RGBAPixel(1, 1, 1, 1.0);
  px = testpng.getPixel(0, 2);
  *px = RGBAPixel(9, 9, 9, 1.0);
  px = testpng.getPixel(1, 2);
  *px = RGBAPixel(7, 7, 7, 1.0);
  px = testpng.getPixel(2, 2);
  *px = RGBAPixel(11, 11, 11, 1.0);

  cout << "done." << endl;
  cout << "Creating Stats object... ";

  Stats st(testpng);

  cout << "done." << endl;

  cout << "Getting average color of the 2x2 square in the lower right corner..." << endl;
  cout << "Expected value:\t" << RGBAPixel(8, 8, 8, 1.0) << endl;
  cout << "Your value:\t" << st.GetAvg(1, 1, 2, 2) << endl;

  cout << "Getting total sum square difference from average color of the 2x2 square in the lower right corner..." << endl;
  cout << "Expected value:\t" << 284.25 << endl;
  cout << "Your value:\t" << st.GetSumSqDev(1, 1, 2, 2) << endl;

  cout << "Leaving TestStats...\n" << endl;
}

void TestImgTreeBuildRender() {
  cout << "Entered TestImgTreeBuildRender..." << endl;

  PNG inputimg;
  inputimg.readFromFile("originals/kkkk-kmnn-960x540.png");

  cout << "Constructing ImgTree from input image... ";
  ImgTree mytree(inputimg);
  cout << "done." << endl;

  cout << "Rendering ImgTree to output... ";
  PNG outputimg = mytree.Render();
  cout << "done." << endl;

  cout << "Saving rendered PNG to file...";
  outputimg.writeToFile("images/kkkk-kmnn-build-render-960x540.png");
  cout << "done." << endl;

  cout << "Leaving TestImgTreeBuildRender...\n" << endl;
}

void TestImgTreeFlipHorizontal() {
  cout << "Entered TestImgTreeFlipHorizontal..." << endl;

  PNG inputimg;
  inputimg.readFromFile("originals/kkkk-kmnn-960x540.png");

  cout << "Constructing ImgTree from input image... ";
  ImgTree mytree(inputimg);
  cout << "done." << endl;

  cout << "Calling FlipHorizontal... ";
  mytree.FlipHorizontal();
  cout << "done." << endl;

  cout << "Rendering flipped ImgTree to output... ";
  PNG outputimg = mytree.Render();
  cout << "done." << endl;

  cout << "Saving rendered PNG to file...";
  outputimg.writeToFile("images/kkkk-kmnn-fliphorizontal-960x540.png");
  cout << "done." << endl;

  cout << "Leaving TestImgTreeFlipHorizontal...\n" << endl;
}

void TestImgTreeCountLeavesPrune() {
  cout << "Entered TestImgTreePrune..." << endl;

  PNG inputimg;
  inputimg.readFromFile("originals/kkkk-kmnn-960x540.png");

  cout << "Constructing ImgTree from input image... ";
  ImgTree mytree(inputimg);
  cout << "done." << endl;

  cout << "Counting leaves in the tree... ";
  unsigned int numleaves = mytree.CountLeaves();
  cout << "Tree contains " << numleaves << " leaves." << endl;

  cout << "Calling Prune with pct=85 and tol=0.02... ";
  mytree.Prune(85, 0.02);
  cout << "done." << endl;

  cout << "Counting leaves in the pruned tree... ";
  unsigned int numleavespruned = mytree.CountLeaves();
  cout << "Tree contains " << numleavespruned << " leaves." << endl;

  cout << "Rendering pruned ImgTree to output... ";
  PNG outputimg = mytree.Render();
  cout << "done." << endl;

  cout << "Saving rendered PNG to file...";
  outputimg.writeToFile("images/kkkk-kmnn-prune-pct85-tol002-960x540.png");
  cout << "done." << endl;

  cout << "Leaving TestImgTreePrune...\n" << endl;
}