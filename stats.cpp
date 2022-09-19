/**
 *  @file stats.cpp
 *  @description implementation of a stats class for rapid calculation of color averages
 *   and total color differences in CPSC 221 PA3
 *  @author CPSC 221
 *
 *  SUBMIT THIS FILE TO PRAIRIELEARN
 */

#include "stats.h"
#include <cmath>

 /**
  *  Computes/retrieves the sum of a single color channel in a defined rectangular region
  *  @pre channel is a valid channel identifier
  *  @pre upper, left, lower, and right are valid image coordinates
  *  @param channel - one of 'r', 'g', or 'b'
  *  @param upper - y-coordinate of the upper edge of the rectangular region
  *  @param left - x-coordinate of the left side of the rectangular region
  *  @param lower - y-coordinate of the lower edge of the rectangular region
  *  @param right - x-coordinate of the right side of the rectangular region
  *  @return the sum of the appropriate color channel values in the defined rectangular area
  */
 //upper=0, left=0, lower=539, right=959
unsigned long Stats::GetColorSum(char channel, unsigned int upper, unsigned int left, unsigned int lower, unsigned int right) {
    unsigned long greenAreaSum = 0;
    unsigned long orangeAreaSum = 0;
    unsigned long blueAreaSum = 0;
    unsigned long purpleAreaSum = 0;
    unsigned long greyAreaSum = 0;
    unsigned long blackAreaSum = 0;
    unsigned long redAreaSum = 0;

    // Edge cases: upper == 0, left == 0
    /*
    if (upper == 0 && left == 0) {
      if (channel == 'r') { return sumR[lower][right]; }
      else if (channel == 'g') { return sumG[lower][right]; }
      else if (channel == 'b') { return sumB[lower][right]; }
    }
    */
   
    if (channel == 'r') {
      greenAreaSum = sumR[lower][right];          
      orangeAreaSum = (int(upper) - 1) < 0 ? 0 : sumR[upper - 1][right];
      purpleAreaSum = (int(left) - 1) < 0 ? 0 : sumR[lower][left - 1];      
      greyAreaSum = (int(upper) - 1) < 0 || (int(left) - 1) < 0 ? 0 : sumR[upper - 1][left - 1]; 

      if (((long long)greenAreaSum - (long long)orangeAreaSum) < 0) { cerr << "\nERROR: greenArea is smaller than orange area\n"; }

      if (orangeAreaSum == 0) {
        blueAreaSum = greenAreaSum - purpleAreaSum;
        return blueAreaSum;
      } 
      else {
        blueAreaSum = greenAreaSum - orangeAreaSum;
      }

      if (((long long)purpleAreaSum - (long long)greyAreaSum) < 0) { cerr << "\nERROR: purpleArea is smaller than grey area\n"; }

      blackAreaSum = purpleAreaSum - greyAreaSum;

      if (((long long)blueAreaSum - (long long)blackAreaSum) < 0) { cerr << "\nERROR: blueArea is smaller than black area\n"; }

      redAreaSum = blueAreaSum - blackAreaSum;
      return redAreaSum;
    } 

    else if (channel == 'g') {
      greenAreaSum = sumG[lower][right];          
      orangeAreaSum = (int(upper) - 1) < 0 ? 0 : sumG[upper - 1][right];
      purpleAreaSum = (int(left) - 1) < 0 ? 0 : sumG[lower][left - 1];      
      greyAreaSum = (int(upper) - 1) < 0 || (int(left) - 1) < 0 ? 0 : sumG[upper - 1][left - 1]; 

      if (((long long)greenAreaSum - (long long)orangeAreaSum) < 0) { cerr << "\nERROR: greenArea is smaller than orange area\n"; }

      if (orangeAreaSum == 0) {
        blueAreaSum = greenAreaSum - purpleAreaSum;
        return blueAreaSum;
      } else {
        blueAreaSum = greenAreaSum - orangeAreaSum;
      }

      if (((long long)purpleAreaSum - (long long)greyAreaSum) < 0) { cerr << "\nERROR: purpleArea is smaller than grey area\n"; }

      blackAreaSum = purpleAreaSum - greyAreaSum;

      if (((long long)blueAreaSum - (long long)blackAreaSum) < 0) { cerr << "\nERROR: blueArea is smaller than black area\n"; }

      redAreaSum = blueAreaSum - blackAreaSum;
      return redAreaSum;

    } 

    else if (channel == 'b') {
      greenAreaSum = sumB[lower][right];          
      orangeAreaSum = (int(upper) - 1) < 0 ? 0 : sumB[upper - 1][right];
      purpleAreaSum = (int(left) - 1) < 0 ? 0 : sumB[lower][left - 1];      
      greyAreaSum = (int(upper) - 1) < 0 || (int(left) - 1) < 0 ? 0 : sumB[upper - 1][left - 1]; 

      if (((long long)greenAreaSum - (long long)orangeAreaSum) < 0) { cerr << "\nERROR: greenArea is smaller than orange area\n"; }

      if (orangeAreaSum == 0) {
        blueAreaSum = greenAreaSum - purpleAreaSum;
        return blueAreaSum;
      } else {
        blueAreaSum = greenAreaSum - orangeAreaSum;
      }

      if (((long long)purpleAreaSum - (long long)greyAreaSum) < 0) { cerr << "\nERROR: purpleArea is smaller than grey area\n"; }

      blackAreaSum = purpleAreaSum - greyAreaSum;

      if (((long long)blueAreaSum - (long long)blackAreaSum) < 0) { cerr << "\nERROR: blueArea is smaller than black area\n"; }

      redAreaSum = blueAreaSum - blackAreaSum;
      return redAreaSum;
    }

    return redAreaSum;
}

/**
 *  Computes/retrieves the sum of alpha values in a defined rectangular region
 *  @pre upper, left, lower, and right are valid image coordinates
 *  @param upper - y-coordinate of the upper edge of the rectangular region
 *  @param left - x-coordinate of the left side of the rectangular region
 *  @param lower - y-coordinate of the lower edge of the rectangular region
 *  @param right - x-coordinate of the right side of the rectangular region
 *  @return the sum of the alpha values in the defined rectangular area
 */
double Stats::GetAlphaSum(unsigned int upper, unsigned int left, unsigned int lower, unsigned int right) {
  double greenAreaSum = 0.0;
  double orangeAreaSum = 0.0;
  double blueAreaSum = 0.0;
  double purpleAreaSum = 0.0;
  double greyAreaSum = 0.0;
  double blackAreaSum = 0.0;
  double redAreaSum = 0.0;

  greenAreaSum = sumA[lower][right];          
  orangeAreaSum = (int(upper) - 1) < 0 ? 0.0 : sumA[upper - 1][right];
  purpleAreaSum = (int(left) - 1) < 0 ? 0.0 : sumA[lower][left - 1];      
  greyAreaSum = (int(upper) - 1) < 0 || (int(left) - 1) < 0 ? 0.0 : sumA[upper - 1][left - 1]; 

  if ((greenAreaSum - orangeAreaSum) < 0) { cerr << "\nERROR: greenArea is smaller than orange area\n"; }

  if (orangeAreaSum <= 1e-4) { // check if orangeAreaSum is virtually 0
    blueAreaSum = greenAreaSum - purpleAreaSum;
    return blueAreaSum;
  } 
  else {
    blueAreaSum = greenAreaSum - orangeAreaSum;
  }

  if ((purpleAreaSum - greyAreaSum) < 0) { cerr << "\nERROR: purpleArea is smaller than grey area\n"; }

  blackAreaSum = purpleAreaSum - greyAreaSum;

  if ((blueAreaSum - blackAreaSum) < 0) { cerr << "\nERROR: blueArea is smaller than black area\n"; }

  redAreaSum = blueAreaSum - blackAreaSum; 

  return redAreaSum;

}

/**
 *  Computes/retrieves the squared sum of a single color channel in a defined rectangular region
 *  @pre channel is a valid channel identifier
 *  @pre upper, left, lower, and right are valid image coordinates
 *  @param channel - one of 'r', 'g', or 'b'
 *  @param upper - y-coordinate of the upper edge of the rectangular region
 *  @param left - x-coordinate of the left side of the rectangular region
 *  @param lower - y-coordinate of the lower edge of the rectangular region
 *  @param right - x-coordinate of the right side of the rectangular region
 *  @return the squared sum of the appropriate color channel values in the defined rectangular area
 */
unsigned long Stats::GetColorSumSq(char channel, unsigned int upper, unsigned int left, unsigned int lower, unsigned int right) {
    unsigned long greenAreaSum = 0;
    unsigned long orangeAreaSum = 0;
    unsigned long blueAreaSum = 0;
    unsigned long purpleAreaSum = 0;
    unsigned long greyAreaSum = 0;
    unsigned long blackAreaSum = 0;
    unsigned long redAreaSum = 0;

    if (channel == 'r') {
      greenAreaSum = sumSqR[lower][right];          
      orangeAreaSum = (int(upper) - 1) < 0 ? 0 : sumSqR[upper - 1][right];
      purpleAreaSum = (int(left) - 1) < 0 ? 0 : sumSqR[lower][left - 1];      
      greyAreaSum = (int(upper) - 1) < 0 || (int(left) - 1) < 0 ? 0 : sumSqR[upper - 1][left - 1]; 

      if (((long long)greenAreaSum - (long long)orangeAreaSum) < 0) { cerr << "\nERROR: greenArea is smaller than orange area\n"; }

      if (orangeAreaSum == 0) {
        blueAreaSum = greenAreaSum - purpleAreaSum;
        return blueAreaSum;
      } 
      else {
        blueAreaSum = greenAreaSum - orangeAreaSum;
      }

      if (((long long)purpleAreaSum - (long long)greyAreaSum) < 0) { cerr << "\nERROR: purpleArea is smaller than grey area\n"; }

      blackAreaSum = purpleAreaSum - greyAreaSum;

      if (((long long)blueAreaSum - (long long)blackAreaSum) < 0) { cerr << "\nERROR: blueArea is smaller than black area\n"; }

      redAreaSum = blueAreaSum - blackAreaSum;
      return redAreaSum;
    } 

    else if (channel == 'g') {
      greenAreaSum = sumSqG[lower][right];          
      orangeAreaSum = (int(upper) - 1) < 0 ? 0 : sumSqG[upper - 1][right];
      purpleAreaSum = (int(left) - 1) < 0 ? 0 : sumSqG[lower][left - 1];      
      greyAreaSum = (int(upper) - 1) < 0 || (int(left) - 1) < 0 ? 0 : sumSqG[upper - 1][left - 1]; 

      if (((long long)greenAreaSum - (long long)orangeAreaSum) < 0) { cerr << "\nERROR: greenArea is smaller than orange area\n"; }

      if (orangeAreaSum == 0) {
        blueAreaSum = greenAreaSum - purpleAreaSum;
        return blueAreaSum;
      } else {
        blueAreaSum = greenAreaSum - orangeAreaSum;
      }

      if (((long long)purpleAreaSum - (long long)greyAreaSum) < 0) { cerr << "\nERROR: purpleArea is smaller than grey area\n"; }

      blackAreaSum = purpleAreaSum - greyAreaSum;

      if (((long long)blueAreaSum - (long long)blackAreaSum) < 0) { cerr << "\nERROR: blueArea is smaller than black area\n"; }

      redAreaSum = blueAreaSum - blackAreaSum;
      return redAreaSum;

    } 

    else if (channel == 'b') {
      greenAreaSum = sumSqB[lower][right];          
      orangeAreaSum = (int(upper) - 1) < 0 ? 0 : sumSqB[upper - 1][right];
      purpleAreaSum = (int(left) - 1) < 0 ? 0 : sumSqB[lower][left - 1];      
      greyAreaSum = (int(upper) - 1) < 0 || (int(left) - 1) < 0 ? 0 : sumSqB[upper - 1][left - 1]; 

      if (((long long)greenAreaSum - (long long)orangeAreaSum) < 0) { cerr << "\nERROR: greenArea is smaller than orange area\n"; }

      if (orangeAreaSum == 0) {
        blueAreaSum = greenAreaSum - purpleAreaSum;
        return blueAreaSum;
      } else {
        blueAreaSum = greenAreaSum - orangeAreaSum;
      }

      if (((long long)purpleAreaSum - (long long)greyAreaSum) < 0) { cerr << "\nERROR: purpleArea is smaller than grey area\n"; }

      blackAreaSum = purpleAreaSum - greyAreaSum;

      if (((long long)blueAreaSum - (long long)blackAreaSum) < 0) { cerr << "\nERROR: blueArea is smaller than black area\n"; }

      redAreaSum = blueAreaSum - blackAreaSum;
      return redAreaSum;
    }

    return redAreaSum;
}

/**
 *  Computes/retrieves the squared sum of alpha values in a defined rectangular region
 *  @pre upper, left, lower, and right are valid image coordinates
 *  @param upper - y-coordinate of the upper edge of the rectangular region
 *  @param left - x-coordinate of the left side of the rectangular region
 *  @param lower - y-coordinate of the lower edge of the rectangular region
 *  @param right - x-coordinate of the right side of the rectangular region
 *  @return the squared sum of the alpha values in the defined rectangular area
 */
double Stats::GetAlphaSumSq(unsigned int upper, unsigned int left, unsigned int lower, unsigned int right) {
  double greenAreaSum = 0.0;
  double orangeAreaSum = 0.0;
  double blueAreaSum = 0.0;
  double purpleAreaSum = 0.0;
  double greyAreaSum = 0.0;
  double blackAreaSum = 0.0;
  double redAreaSum = 0.0;

  greenAreaSum = sumSqA[lower][right];          
  orangeAreaSum = (int(upper) - 1) < 0 ? 0.0 : sumSqA[upper - 1][right];
  purpleAreaSum = (int(left) - 1) < 0 ? 0.0 : sumSqA[lower][left - 1];      
  greyAreaSum = (int(upper) - 1) < 0 || (int(left) - 1) < 0 ? 0.0 : sumSqA[upper - 1][left - 1]; 

  if ((greenAreaSum - orangeAreaSum) < 0) { cerr << "\nERROR: greenArea is smaller than orange area\n"; }

  if (orangeAreaSum <= 1e-4) { // check if orangeAreaSum is virtually 0
    blueAreaSum = greenAreaSum - purpleAreaSum;
    return blueAreaSum;
  } 
  else {
    blueAreaSum = greenAreaSum - orangeAreaSum;
  }

  if ((purpleAreaSum - greyAreaSum) < 0) { cerr << "\nERROR: purpleArea is smaller than grey area\n"; }

  blackAreaSum = purpleAreaSum - greyAreaSum;

  if ((blueAreaSum - blackAreaSum) < 0) { cerr << "\nERROR: blueArea is smaller than black area\n"; }

  redAreaSum = blueAreaSum - blackAreaSum; 

  return redAreaSum;
}

/**
 *  Simple function to compute the number of pixels in a defined rectangular region
 *  @pre upper, left, lower, and right are valid image coordinates
 *  @param upper - y-coordinate of the upper edge of the rectangular region
 *  @param left - x-coordinate of the left side of the rectangular region
 *  @param lower - y-coordinate of the lower edge of the rectangular region
 *  @param right - x-coordinate of the right side of the rectangular region
 *  @return the area of the defined rectangular area, in pixels
 */
unsigned int Stats::GetRectangleArea(unsigned int upper, unsigned int left, unsigned int lower, unsigned int right) {
  return (lower - upper + 1) * (right - left + 1);
}

/**
 *  Parameterized constructor
 *  Builds the channel sum and squared sum vectors from the supplied input image.
 *  Each channel sum vector's entry (x,y) will contain the sum of their respective
 *  color channels of all pixels in the rectangular image region bounded by (0,0) and (x,y).
 *  Likewise, each channel squared sum vector's entry (x,y) will contain the squared sum of their
 *  respective color channels of all pixels in the rectangular image region bounded by (0,0) and (x,y).
 *
 *  ***DON'T FORGET TO PRE-MULTIPLY THE ALPHA CHANNEL***
 * 
 *  @param img - input image from which the channel sum vectors will be populated
 */
Stats::Stats(const PNG& img) {
  // complete your implementation below
  unsigned long redColorSum = 0;
  unsigned long blueColorSum = 0;
  unsigned long greenColorSum = 0;
  double alphaSum = 0.0;

  unsigned long redColorSumSq = 0;
  unsigned long blueColorSumSq = 0;
  unsigned long greenColorSumSq = 0;
  long double alphaSumSq = 0.0;

  sumR.resize(img.height());
  sumG.resize(img.height());
  sumB.resize(img.height());
  sumA.resize(img.height());

  sumSqR.resize(img.height());
  sumSqG.resize(img.height());
  sumSqB.resize(img.height());
  sumSqA.resize(img.height());

  for (unsigned int y = 0; y < img.height(); ++y) {
    redColorSum = 0;
    blueColorSum = 0;
    greenColorSum = 0;
    alphaSum = 0.0;

    redColorSumSq = 0;
    blueColorSumSq = 0;
    greenColorSumSq = 0;
    alphaSumSq = 0.0;

    sumR[y].resize(img.width());
    sumG[y].resize(img.width());
    sumB[y].resize(img.width());
    sumA[y].resize(img.width());

    sumSqR[y].resize(img.width());
    sumSqG[y].resize(img.width());
    sumSqB[y].resize(img.width());
    sumSqA[y].resize(img.width());

    for (unsigned int x = 0; x < img.width(); ++x) {
      RGBAPixel* curPixel = img.getPixel(x, y);

      redColorSum += curPixel->r;
      blueColorSum += curPixel->b;
      greenColorSum += curPixel->g;
      alphaSum += curPixel->a * 255.0;

      redColorSumSq += curPixel->r * curPixel->r;
      blueColorSumSq += curPixel->b * curPixel->b;
      greenColorSumSq += curPixel->g * curPixel->g;
      alphaSumSq += (curPixel->a * 255.0) * (curPixel->a * 255.0);

      unsigned long prevRed = 0;
      unsigned long prevBlue = 0;
      unsigned long prevGreen = 0;
      double prevAlpha = 0.0;

      unsigned long prevRedSq = 0;
      unsigned long prevBlueSq = 0;
      unsigned long prevGreenSq = 0;
      double prevAlphaSq = 0.0;

      if (y != 0) {
        prevRed = sumR[y - 1][x]; // used to be += but not necessary
        prevBlue = sumB[y - 1][x];
        prevGreen = sumG[y - 1][x];
        prevAlpha = sumA[y - 1][x];

        prevRedSq = sumSqR[y - 1][x];
        prevBlueSq = sumSqB[y - 1][x];
        prevGreenSq = sumSqG[y - 1][x];
        prevAlphaSq = sumSqA[y - 1][x];
      }

      sumR[y][x] = redColorSum + prevRed;
      sumG[y][x] = greenColorSum + prevGreen;
      sumB[y][x] = blueColorSum + prevBlue;
      sumA[y][x] = alphaSum + prevAlpha;

      sumSqR[y][x] = redColorSumSq + prevRedSq;
      sumSqG[y][x] = greenColorSumSq + prevGreenSq;
      sumSqB[y][x] = blueColorSumSq + prevBlueSq;
      sumSqA[y][x] = alphaSumSq + prevAlphaSq;
    }
  }
  
}

/**
 *  Computes/retrieves the average color of all pixels contained in the rectangle
 *  bounded by upper, left, lower, and right. Fractional values should be
 *  truncated/rounded down for assignment into integer variables.
 *  @pre upper, left, lower, and right are valid image coordinates
 *  @param upper - y-coordinate of the upper edge of the rectangular region
 *  @param left - x-coordinate of the left side of the rectangular region
 *  @param lower - y-coordinate of the lower edge of the rectangular region
 *  @param right - x-coordinate of the right side of the rectangular region
 *  @return pixel containing the average color of the pixels in the defined rectangle
 */
RGBAPixel Stats::GetAvg(unsigned int upper, unsigned int left, unsigned int lower, unsigned int right) {
  unsigned long red = GetColorSum('r', upper, left, lower, right);
  unsigned long blue = GetColorSum('b', upper, left, lower, right); 
  unsigned long green = GetColorSum('g', upper, left, lower, right);
  double alpha = GetAlphaSum(upper, left, lower, right);
  int area = GetRectangleArea(upper, left, lower, right);

  // std::floor(red / double(area))
  int redAvg = red / area;
  int blueAvg = blue / area;
  int greenAvg = green / area;
  double alphaAvg = (alpha / area) / 255.0;

  return RGBAPixel(redAvg, greenAvg, blueAvg, alphaAvg);
}

/**
 *  Computes the total sum squared difference from the mean, for the specified rectangle.
 *  Each channel's sum squared difference is computed separately, and then added to form the total.
 *
 *  Note that using the GetAvg function in computing the sum squared difference will result in
 *  accumulation of error especially with larger rectangles.
 *  You should use more precise computation of the average color for this function's intermediate steps.
 * 
 *  @pre upper, left, lower, and right are valid image coordinates
 *  @param upper - y-coordinate of the upper edge of the rectangular region
 *  @param left - x-coordinate of the left side of the rectangular region
 *  @param lower - y-coordinate of the lower edge of the rectangular region
 *  @param right - x-coordinate of the right side of the rectangular region
 *  @return total sum of squared deviations from the mean, over all color channels.
 */

double Stats::GetSumSqDev(unsigned int upper, unsigned int left, unsigned int lower, unsigned int right) {
  unsigned long redSq = GetColorSumSq('r', upper, left, lower, right);
  unsigned long blueSq = GetColorSumSq('b', upper, left, lower, right); 
  unsigned long greenSq = GetColorSumSq('g', upper, left, lower, right);
  double alphaSq = GetAlphaSumSq(upper, left, lower, right);
 
  unsigned long red = GetColorSum('r', upper, left, lower, right);
  unsigned long blue = GetColorSum('b', upper, left, lower, right); 
  unsigned long green = GetColorSum('g', upper, left, lower, right);
  double alpha = GetAlphaSum(upper, left, lower, right);

  int area = GetRectangleArea(upper, left, lower, right);

  long double redAvg = (double)redSq - ((red * red) / (double)area);
  long double blueAvg = (double)blueSq - ((blue * blue) / (double)area);
  long double greenAvg = (double)greenSq - ((green * green) / (double)area);
  double alphaAvg = alphaSq - ((alpha * alpha) / double(area)); // Potientially add in a divide by 255 here???
  
  return redAvg + blueAvg + greenAvg + alphaAvg;
}