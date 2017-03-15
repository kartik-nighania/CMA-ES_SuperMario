#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstddef>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <string>

#include <boost/random.hpp>

typedef boost::mt19937 RNGType;  // Choose random number generator.

namespace mlpack {
namespace neuro_cmaes 
{

// Return the sign of a number.
template <typename T> 
  int sgn(T val) 
{
    return (T(0) < val) - (val < T(0));
}

template<typename T>
T square(T d)
{
  return d*d;
}

template<typename T>
T maxElement(const T* rgd, int len)
{
  return *std::max_element(rgd, rgd + len);
}

template<typename T>
T minElement(const T* rgd, int len)
{
  return *std::min_element(rgd, rgd + len);
}

template<typename T>
int maxIndex(const T* rgd, int len)
{
  return std::max_element(rgd, rgd + len) - rgd;
}

template<typename T>
int minIndex(const T* rgd, int len)
{
  return std::min_element(rgd, rgd + len) - rgd;
}

/** sqrt(a^2 + b^2) numerically stable. */
template<typename T>
T myhypot(T a, T b)
{
  const register T fabsa = std::fabs(a), fabsb = std::fabs(b);
  if(fabsa > fabsb)
  {
    const register T r = b / a;
    return fabsa*std::sqrt(T(1)+r*r);
  }
  else if(b != T(0))
  {
    const register T r = a / b;
    return fabsb*std::sqrt(T(1)+r*r);
  }
  else
    return T(0);
}

inline double sigmoid(double x)
 { 
  return 1.0 / (1.0 + exp(-x));
 }

inline double relu(double x)
 {
  return (x > 0)? x:0;
 }

// Random number generator.
RNGType rng;

// Set seed for random number generator
void Seed(int seedVal)
 {
  rng.seed(seedVal);
 }

// Set seed by time for random number generator
void TimeSeed() 
{
  rng.seed(time(0));
}

// Returns randomly either 1 or -1
int RandPosNeg()
 {
  boost::random::uniform_int_distribution<> dist(0, 1);
  return dist(rng);
 }

// Returns a random integer between [x, y]
// in case of ( 0 .. 1 ) returns 0
int RandInt(int x, int y)
 {
	boost::random::uniform_int_distribution<> dist(x, y);
    return dist(rng);
 }

// Return a random float between [0, 1]
double RandFloat() 
{
  boost::random::uniform_01<> dist;
  return dist(rng); 
}

// Return a random float between [x, y]
double RandFloat(double x, double y) {
  boost::random::uniform_real_distribution<> dist(x, y);
  return dist(rng);
}



}  // namespace mlpack
}  // namespace neuro_cmaes

#endif  // UTILS_HPP
