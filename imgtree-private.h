/**
 *  @file imgtree-private.h
 *  @description student-defined functions for binary tree class storing image data for CPSC 221 PA3
 *  @author CPSC 221
 *
 *  SUBMIT THIS FILE TO PRAIRIELEARN, EVEN IF YOU DO NOT ADD ANYTHING TO IT
 * 
 *  Usage: As this file is included entirely into imgtree.h by the #include statement on line __
 *  you only need to write the function signature,
 *  e.g.
 *  void MyOwnSpecialFunction(int fave_num, ImgTreeNode* nd);
 * 
 *  and add/complete your implementation in imgtree.cpp
 */

int getArea(int upr, int lft, int lwr, int rgt);
unsigned int getXcoordToSplit(Stats& s, unsigned int upr, unsigned int lft, unsigned int lwr, unsigned int rt);
unsigned int getYcoordToSplit(Stats& s, unsigned int upr, unsigned int lft, unsigned int lwr, unsigned int rt);
void Render(ImgTreeNode* node, const PNG& img) const; 
void Clear(ImgTreeNode* subRoot);
ImgTreeNode* Copy(const ImgTreeNode* otherSubRoot);
unsigned int CountLeaves(ImgTreeNode* subRoot) const;
void FlipHorizontal(ImgTreeNode* node, bool isLeft, int parentLeft, int parentRight);
void Prune(double pct, double tol, ImgTreeNode* subRoot);
int PruneCurrentSubTree(const RGBAPixel& rootAvg, double tol, ImgTreeNode* subRoot);
