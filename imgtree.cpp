/**
 *  @file imgtree.cpp
 *  @description implementation of a binary tree class storing image data for CPSC 221 PA3
 *  @author CPSC 221
 *
 *  SUBMIT THIS FILE TO PRAIRIELEARN
 */

#include "imgtree.h"
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>

// not necessary to include imgtree-private.h since it is already included in imgtree.h
using std::fabs;
using std::abs;
using std::vector;
using std::swap;

    /**
     *  Default constructor creates an empty tree
     */
ImgTree::ImgTree() : root(nullptr), imgwidth(0), imgheight(0) {}

/**
 *  Parameterized constructor creates a tree from an input image.
 *  Every leaf in the tree corresponds to a single pixel in the PNG.
 *  Every non-leaf node corresponds to a rectangle of pixels defined
 *  by upper, left, lower, and right image coordinates. A node will
 *  also store the average color over its defined rectangle. Note that
 *  this average MUST be freshly computed from the individual pixels in
 *  the image and NOT by computing a weighted average from the colors
 *  of its direct children, as this would lead to an accumulation of
 *  rounding errors due to the integer nature of RGB color channels.
 *
 *  A non-leaf node's children correspond to a partition of the node's
 *  rectangle into two smaller rectangles, either by a horizontal line
 *  creating an upper rectangle and a lower rectangle, or a vertical line
 *  creating a left rectangle and a right rectangle.
 *
 *  The split is determined as follows:
 *  1. If the current node's rectangle width is the same or larger than its height,
 *     then a vertical line will divide the rectangle into left and right rectangles.
 *     If the current node's rectangle width is smaller than its height,
 *     then a horizontal line will divide the rectangle into upper and lower rectangles.
 *  2. The coordinate of the dividing line is chosen such that combined sum squared
 *     deviations from the mean color in the left/upper and right/lower rectangles is minimal
 *     among all other potential dividing lines of the same orientation.
 *     e.g. for a region (0,0) to (3,2), the candidates are the vertical lines dividing the region into:
 *     - (0,0) to (0,2) and (1,0) to (3,2)
 *     - (0,0) to (1,2) and (2,0) to (3,2)
 *     - (0,0) to (2,2) and (3,0) to (3,2)
 *     The line which produces the minimal combined sum squared scores on the left
 *     and right will be used for the split.
 *  3. In the unlikely event that multiple candidates produce the same score, the one which
 *     most evenly splits the rectangular area will be chosen.
 *  4. In the even more unlikely even that two candidates produce the same score and produce
 *     the same minimal area difference, the one with the smaller coordinate will be chosen.
 */

ImgTree::ImgTree(const PNG& img) : imgwidth(img.width()), imgheight(img.height()) {
  Stats statObj(img);
  root = BuildNode(statObj, 0, 0, imgheight - 1, imgwidth - 1);
}

/**
 *  Copy constructor creates a new tree that is structurally the same as the input tree and
 *  contains the same image data.
 *  @param other - the ImgTree to be copied
 */
ImgTree::ImgTree(const ImgTree& other) {
  // This implementation has already been completed for you
  Copy(other);
}

/**
 *  Overloaded assignment operator allows statements such as mytree_a = some_other_tree;
 *  Releases any existing memory and reconstructs this tree as a copy of other.
 *  @param rhs - the right hand side of the assignment statement
 */
ImgTree& ImgTree::operator=(const ImgTree& rhs) {
  // This implementation has already been completed for you
  if (this != &rhs) { // this and rhs are physically different trees in memory
    Clear(); // release any existing heap memory for this tree
    Copy(rhs);  // rebuild this tree as a copy of rhs
  }
  return *this;
}

/**
 *  Destructor releases any heap memory associated with this tree.
 */
ImgTree::~ImgTree() {
  // This implementation has already been completed for you
  Clear();
}

/**
 *  Releases all heap memory associated with this tree, restoring it to an "empty tree" state.
 *  Will be useful to define a recursive helper function for this.
 */
void ImgTree::Clear() {
  Clear(root);
  root = nullptr;
}

/*
* Private recursive helper for clear() function
*
*/
void ImgTree::Clear(ImgTreeNode* subRoot) {
  if (subRoot == nullptr) { return; } // Base case 1: Empty node passed in
  
  Clear(subRoot->A); // Recursive step
  Clear(subRoot->B);

  delete subRoot; // Work step
  subRoot = nullptr;
}

/**
 *  Copies the supplied parameter tree into this tree. Does not free any memory.
 *  Called by the copy constructor and operator=.
 *  Will be useful to define a recursive helper function for this.
 *  HINT: the helper should allocate a new node, and return a pointer to the allocated node.
 *        See the documention for BuildNode - this should work similarly.
 */
void ImgTree::Copy(const ImgTree& other) {
  // complete your implementation below
  imgwidth = other.imgwidth;
  imgheight = other.imgheight;
  root = Copy(other.root);
}

/* Private Recursive helper for Copy
 * @param otherSubRoot current node to copy from
 * @return constructed node with children that has copied data from otherSubRoot and its children
 * @post this->root has the contents of the copy
 */

ImgTreeNode* ImgTree::Copy(const ImgTreeNode* otherSubRoot) {
    if (otherSubRoot == nullptr) { return nullptr; } // Base case

    unsigned int upper = otherSubRoot->upper;        // Work step
    unsigned int left = otherSubRoot->left;
    unsigned int lower = otherSubRoot->lower;
    unsigned int right = otherSubRoot->right;
    RGBAPixel average(otherSubRoot->avg);

    ImgTreeNode* curSubRootNode = new ImgTreeNode(upper, left, lower, right, average);

    curSubRootNode->A = Copy(otherSubRoot->A);       // Recursive step
    curSubRootNode->B = Copy(otherSubRoot->B);

    return curSubRootNode;
}

/**
 *  Recursive helper function for initial construction of the tree. Constructs a single
 *  node according to supplied Stats and the requirements specified by the constructor
 *  documentation, and returns a pointer to the completed node.
 *  @param s - populated Stats object for computing this node's attributes
 *  @param upr - y-coordinate of the upper edge of the node's rectangular region
 *  @param lft - x-coordinate of the left side of the node's rectangular region
 *  @param lwr - y-coordinate of the lower edge of the node's rectangular region
 *  @param rt - x-coordinate of the right side of the node's rectangular region
 *  @return - pointer to a (completed) newly-allocated node for the specified parameters.
 */

ImgTreeNode* ImgTree::BuildNode(Stats& s, unsigned int upr, unsigned int lft, unsigned int lwr, unsigned int rt) {
  assert("ERROR: img area is negative" && !(((int(lwr) - int(upr) + 1) < 0) || ((int(rt) - int(lft) + 1) < 0)));

  ImgTreeNode* currentNode = new ImgTreeNode(upr, lft, lwr, rt, s.GetAvg(upr, lft, lwr, rt));
  int imgArea = getArea(upr, lft, lwr, rt);

  if (imgArea <= 1) { // Base Case: We reached a pixel, so just return node object
    return currentNode;
  }

  // Figure out if we need vertical or horizontal split
  int height = lwr - upr + 1;
  int width = rt - lft + 1;
  bool isVerticalLine = width >= height;

  // Work and Recursive steps
  if (isVerticalLine) { // Vertical Line divides rectangle into left and right rectangles
    int xCoordinate = getXcoordToSplit(s, upr, lft, lwr, rt);

    currentNode->A = BuildNode(s, upr, lft, lwr, xCoordinate);
    currentNode->B = BuildNode(s, upr, xCoordinate + 1, lwr, rt);
  }
  else { // Horizontal Line will divide rectangle into upper and lower rectangles
    int yCoordinate = getYcoordToSplit(s, upr, lft, lwr, rt);

    currentNode->A = BuildNode(s, upr, lft, yCoordinate, rt);
    currentNode->B = BuildNode(s, yCoordinate + 1, lft, lwr, rt);
  }

  return currentNode;
}

int ImgTree::getArea(int upr, int lft, int lwr, int rt) {
  return (lwr - upr + 1) * (rt - lft + 1);
}

unsigned int ImgTree::getXcoordToSplit(Stats& s, unsigned int upr, unsigned int lft, unsigned int lwr, unsigned int rt) {
  std::vector<SplitInfo> splitInfoArr;

  // Generate variance scores and coordinate data for each of those variances for vertical splits
  for (unsigned int x = lft; x < rt; ++x) {
    SplitInfo currentSplit;
    currentSplit.coordinate = x;
    currentSplit.sumsqscore = s.GetSumSqDev(upr, lft, lwr, x) + s.GetSumSqDev(upr, x + 1, lwr, rt);
    splitInfoArr.push_back(currentSplit);
  }

  // Index with the minimum variance
  size_t minVarianceIndex = 0;

  // Find the index with minimum variance
  for (unsigned int i = 1; i < splitInfoArr.size(); ++i) {
    if (splitInfoArr[i].sumsqscore < splitInfoArr[minVarianceIndex].sumsqscore) {
      minVarianceIndex = i;
    }
  }

  // For edge cases, find the coordinate which better splits the area
  for (size_t i = 0; i < splitInfoArr.size(); ++i) {

    if (i != minVarianceIndex && fabs(splitInfoArr[i].sumsqscore - splitInfoArr[minVarianceIndex].sumsqscore) <= 0.0001) {

      int mVI_leftArea = getArea(upr, lft, lwr, splitInfoArr[minVarianceIndex].coordinate);
      int mVI_rightArea = getArea(upr, splitInfoArr[minVarianceIndex].coordinate + 1, lwr, rt); 

      int i_leftArea = getArea(upr, lft, lwr, splitInfoArr[i].coordinate);
      int i_rightArea = getArea(upr, splitInfoArr[i].coordinate + 1, lwr, rt); 

      if (abs(mVI_rightArea - mVI_leftArea) > abs(i_leftArea - i_rightArea)) {
        minVarianceIndex = i;
      } 
    }

  }

  return splitInfoArr[minVarianceIndex].coordinate;
}

unsigned int ImgTree::getYcoordToSplit(Stats& s, unsigned int upr, unsigned int lft, unsigned int lwr, unsigned int rt) {
  std::vector<SplitInfo> splitInfoArr;

  // Generate variance scores and coordinate data for each of those variances for horizontal splits
  for (unsigned int y = upr; y < lwr; ++y) {
    SplitInfo currentSplit;
    currentSplit.coordinate = y;
    currentSplit.sumsqscore = s.GetSumSqDev(upr, lft, y, rt) + s.GetSumSqDev(y + 1, lft, lwr, rt);
    splitInfoArr.push_back(currentSplit);
  }

  // Index with the minimum variance
  size_t minVarianceIndex = 0;

  // Find the index with minimum variance
  for (size_t i = 1; i < splitInfoArr.size(); ++i) {
    if (splitInfoArr[i].sumsqscore < splitInfoArr[minVarianceIndex].sumsqscore) {
      minVarianceIndex = i;
    }
  }

  // For edge cases, find the coordinate which better splits the area
  for (size_t i = 0; i < splitInfoArr.size(); ++i) {

    if (i != minVarianceIndex && fabs(splitInfoArr[i].sumsqscore - splitInfoArr[minVarianceIndex].sumsqscore) <= 0.0001) {

      int mVI_topArea = getArea(upr, lft, splitInfoArr[minVarianceIndex].coordinate, rt);
      int mVI_botArea = getArea(splitInfoArr[minVarianceIndex].coordinate + 1, lft, lwr, rt); 

      int i_topArea = getArea(upr, lft, splitInfoArr[i].coordinate, rt);
      int i_botArea = getArea(splitInfoArr[minVarianceIndex].coordinate + 1, lft, lwr, rt);  

      if (abs(mVI_topArea - mVI_botArea) > abs(i_topArea - i_botArea)) {
        minVarianceIndex = i;
      } 
    }

  }

  return splitInfoArr[minVarianceIndex].coordinate;
}


/**
 *  Produces a PNG of appropriate dimensions and paints every leaf node's rectangle
 *  into the appropriate area of the PNG.
 *  May be called on pruned trees.
 *  @return fully-colored PNG, painted from the tree's leaf node data
 */

PNG ImgTree::Render() const {
  PNG img(imgwidth, imgheight);
  Render(root, img);
  return img;
}

void ImgTree::Render(ImgTreeNode* node, const PNG& img) const { 
  if (node != nullptr) {  

    if (node->A == nullptr || node->B == nullptr) {
      for (unsigned int y = node->upper; y <= node->lower; ++y) {
        for (unsigned int x = node->left; x <= node->right; ++x) {
          RGBAPixel* pixel = img.getPixel(x, y);
          pixel->r = node->avg.r;
          pixel->g = node->avg.g;
          pixel->b = node->avg.b;
          pixel->a = node->avg.a;
        }
      }
    } 

    else {
      Render(node->A, img);
      Render(node->B, img);
    }

  }
}

/**
 *  Rearranges a tree's internal pointers and node content so that its image data
 *  appears flipped horizontally when rendered.
 *  Beware that the tree may or may not have been pruned!
 *  Will be useful to define a recursive helper function for this.
 */
void ImgTree::FlipHorizontal() {
  if (root == nullptr) return;
  FlipHorizontal(root->A, true, root->left, root->right);
  FlipHorizontal(root->B, false, root->left, root->right);
}

// If prune was recently called, then we need to call img on root copy
void ImgTree::FlipHorizontal(ImgTreeNode* node, bool isLeft, int parentLeft, int parentRight) {
  if (node == nullptr) return;

  if (isLeft) {
    unsigned int diff = node->right - node->left;
    node->left = parentRight - diff;
    node->right = parentRight;
    FlipHorizontal(node->A, true, node->left, node->right);
    FlipHorizontal(node->B, false, node->left, node->right);

  } else { // isRight
    unsigned int diff = node->right - node->left;
    node->right = parentLeft + diff;
    node->left = parentLeft; 
    FlipHorizontal(node->A, true, node->left, node->right);
    FlipHorizontal(node->B, false, node->left, node->right);

  }


}


// if (count >= 50) return;
//cerr << "\n Before flip A: (" << subRoot->A->left << ", " << subRoot->A->upper << ") to (" << subRoot->A->right << ", " << subRoot->A->lower << ") and RGBA: " << subRoot->A->avg << "\n";
//cerr << "Before flip B: (" << subRoot->B->left << ", " << subRoot->B->upper << ") to (" << subRoot->B->right << ", " << subRoot->B->lower << ") and RGBA: " << subRoot->B->avg << "\n";
// cerr << "\n After flip A: (" << subRoot->A->left << ", " << subRoot->A->upper << ") to (" << subRoot->A->right << ", " << subRoot->A->lower << ") and RGBA: " << subRoot->A->avg << "\n";
// cerr << "After flip B: (" << subRoot->B->left << ", " << subRoot->B->upper << ") to (" << subRoot->B->right << ", " << subRoot->B->lower << ") and RGBA: " << subRoot->B->avg << "\n";

/**
 *  Rearranges a tree's internal pointers and node content so that its image data
 *  appears flipped vertically when rendered.
 *  Beware that the tree may or may not have been pruned!
 *  Will be useful to define a recursive helper function for this.
 */
void ImgTree::FlipVertical() {
  // complete your implementation below
  PNG img(imgwidth, imgheight);
  Render(root, img);

  vector<vector<RGBAPixel>> pixels;
  pixels.resize(img.width());

  for (unsigned int x = 0; x < img.width(); x++) {
    pixels[x].resize(img.height());

    for (unsigned int y = 0; y < img.height(); y++) {
      pixels[x][y] = *(img.getPixel(x, y));
    }
  }

  for (unsigned int x = 0; x < pixels.size(); x++) {
    vector<RGBAPixel> flip(pixels[x].rbegin(), pixels[x].rend());
    pixels[x].swap(flip);

    for (unsigned int y = 0; y < img.height(); y++) {
      RGBAPixel* pix = img.getPixel(x,y);
      pix->r = pixels[x][y].r;
      pix->g = pixels[x][y].g;
      pix->b = pixels[x][y].b;
      pix->a = pixels[x][y].a;
    }
  }

  Stats statObj(img);
  Clear(root);
  root = nullptr;
  root = BuildNode(statObj, 0, 0, imgheight - 1, imgwidth - 1);


}

/**
 *  Trims subtrees as high as possible in the tree.
 *  A subtree is pruned (all decendants deallocated and child links set to null)
 *  if at least pct (out of 100) of its leaves are within tol of the average
 *  color in the subtree's root.
 *  Assume that this will only be called on trees which have not previously been pruned.
 *  Will be useful to define AT LEAST one recursive helper function for this.
 *  @pre pct is a valid value between 0 and 100
 *  @param pct percentage (out of 100) of leaf node descendants must be within the tolerance threshold
 *             of color difference in order to be pruned
 *  @param tol threshold color difference to qualify for pruning
 */
void ImgTree::Prune(double pct, double tol) {
  // complete your implementation below
  Prune(pct, tol, root);
}

void ImgTree::Prune(double pct, double tol, ImgTreeNode* subRoot) {
  if (subRoot == nullptr) { return; } // Base case: Reached nullptr

  // Work step
  int howManyInTolerance = PruneCurrentSubTree(subRoot->avg, tol, subRoot);
  int numberOfLeaves = CountLeaves(subRoot);
  long double pctInTolerance = ((long double)howManyInTolerance / (long double)numberOfLeaves) * 100.0;

  if (pctInTolerance >= pct) {
    Clear(subRoot->A);
    Clear(subRoot->B);
    subRoot->A = nullptr;
    subRoot->B = nullptr;
    return;
  } else {

  // Recursive step
  Prune(pct, tol, subRoot->A);
  Prune(pct, tol, subRoot->B);
  }

}

int ImgTree::PruneCurrentSubTree(const RGBAPixel& rootAvg, double tol, ImgTreeNode* subRoot) {
  if (subRoot == nullptr) { return 0; } // Base case 1: (NULL root passed in)

  if (subRoot->A == nullptr || subRoot->B == nullptr) { // Base case 2: (Leaf node reached)
    double dist = rootAvg.dist(subRoot->avg);        
    int withinTol = dist <= tol;
    return withinTol;
  }
  
  // Work + Recursive step
  return PruneCurrentSubTree(rootAvg, tol, subRoot->A) + PruneCurrentSubTree(rootAvg, tol, subRoot->B);
}

/**
 *  Counts the number of leaf nodes in the tree.
 *  Will be useful to define a recursive helper function for this.
 */
unsigned int ImgTree::CountLeaves() const {
  return CountLeaves(root);
}

unsigned int ImgTree::CountLeaves(ImgTreeNode* subRoot) const {
  if (subRoot == nullptr) return 0;
  if (subRoot->A == nullptr || subRoot->B == nullptr) return 1;
  return CountLeaves(subRoot->A) + CountLeaves(subRoot->B);
}