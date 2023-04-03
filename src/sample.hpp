#ifndef SAMPLE_H
#define SAMPLE_H

#include <vector>
#include <random>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/xdouble.h>
#include <NTL/BasicThreadPool.h>

const int HELIB_GAUSS_TRUNC = 8;
// First 70 odd decimal places of Pi == 4.0 * atan(1).
const long double PI =
    3.1415926535897932384626433832795028841971693993751058209749445923078164L;

template <typename T>
inline long lsize(const std::vector<T> &v)
{
  return (long)v.size();
}

template <typename T>
inline long lsize(const NTL::Vec<T> &v)
{
  return v.length();
}

inline double RandomReal()
{
  NTL::ZZ num;
  NTL::RandomBits(num, NTL_DOUBLE_PRECISION);

  double denom = std::ldexp(1.0, NTL_DOUBLE_PRECISION);
  // 2^NTL_DOUBLE_PRECISION

  return NTL::conv<double>(num) / denom;
}

// Choose a vector of continuous Gaussians
inline void sampleGaussian(std::vector<double> &dvec, long n, double stdev)
{
  if (n <= 0)
    n = lsize(dvec);
  if (n <= 0)
    return;

  dvec.resize(n); // allocate space for n variables

  // Uses the Box-Muller method to get two Normal(0,stdev^2) variables
  for (long i = 0; i < n; i += 2)
  {
    // r1, r2 are "uniform in (0,1)"
    double r1 = RandomReal();
    double r2 = RandomReal();
    while (r2 == 0)
      r2 = RandomReal();
    double theta = 2.0 * PI * r1;
    double rr = std::sqrt(-2.0 * log(r2));
    if (rr > HELIB_GAUSS_TRUNC)
    {
      // sanity-check, truncate at HELIB_GAUSS_TRUNC standard deviations
      rr = HELIB_GAUSS_TRUNC;
    }

    // Generate two Gaussians RV's
    dvec[i] = stdev * rr * std::cos(theta);
    if (i + 1 < n)
      dvec[i + 1] = stdev * rr * std::sin(theta);
  }
}

inline void sampleGaussian(std::vector<NTL::xdouble> &dvec, long n, NTL::xdouble stdev)
{
  if (n <= 0)
    n = lsize(dvec);
  if (n <= 0)
    return;

  dvec.resize(n); // allocate space for n variables

  // Uses the Box-Muller method to get two Normal(0,stdev^2) variables
  for (long i = 0; i < n; i += 2)
  {
    // r1, r2 are "uniform in (0,1)"
    double r1 = RandomReal();
    double r2 = RandomReal();
    while (r2 == 0)
      r2 = RandomReal();
    double theta = 2.0 * PI * r1;
    double rr = std::sqrt(-2.0 * log(r2));
    if (rr > HELIB_GAUSS_TRUNC)
    {
      // sanity-check, truncate at HELIB_GAUSS_TRUNC standard deviations
      rr = HELIB_GAUSS_TRUNC;
    }

    // Generate two Gaussians RV's
    dvec[i] = stdev * rr * std::cos(theta);
    if (i + 1 < n)
      dvec[i + 1] = stdev * rr * std::sin(theta);
  }
}

inline void sampleGaussian(NTL::ZZX &poly, long n, NTL::xdouble stdev)
{
  if (n <= 0)
    return;
  std::vector<NTL::xdouble> dvec;
  sampleGaussian(dvec, n, stdev); // sample continuous Gaussians

  // round and copy to coefficients of poly
  poly.SetLength(n); // allocate space for degree-(n-1) polynomial
  for (long i = 0; i < n; i++)
    NTL::conv(poly[i], dvec[i] + 0.5); // round to nearest integer
  poly.normalize();
}

inline void sampleGaussian(NTL::ZZ_pX &poly, long n, NTL::xdouble stdev)
{
  if (n <= 0)
    return;
  std::vector<NTL::xdouble> dvec;
  sampleGaussian(dvec, n, stdev); // sample continuous Gaussians

  NTL::ZZ tempInt;
  // round and copy to coefficients of poly
  poly.SetLength(n); // allocate space for degree-(n-1) polynomial
  for (long i = 0; i < n; i++)
  {
    NTL::conv(tempInt, dvec[i] + 0.5); // round to nearest integer
    NTL::conv(poly[i], tempInt);
  }
  // std::cout << "B: " << NTL::deg(poly) << ", ";
  poly.normalize();
  // std::cout << NTL::deg(poly) << ", ";
}

inline void sampleUniform(NTL::ZZ_pX &poly, long n, int eta)
{
  if (n <= 0)
    return;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(-eta, eta);

  poly.SetLength(n); // allocate space for degree-(n-1) polynomial
  for (long i = 0; i < n; i++)
  {
    poly[i] = dist(gen);
  }
  // std::cout << "B: " << NTL::deg(poly) << ", ";
  poly.normalize();
  // std::cout << NTL::deg(poly) << ", ";
}

inline void sampleUniform(NTL::ZZX &poly, long n, int eta)
{
  if (n <= 0)
    return;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(-eta, eta);

  poly.SetLength(n); // allocate space for degree-(n-1) polynomial
  for (long i = 0; i < n; i++)
  {
    poly[i] = dist(gen);
  }
  // std::cout << "B: " << NTL::deg(poly) << ", ";
  poly.normalize();
  // std::cout << NTL::deg(poly) << ", ";
}

#endif
