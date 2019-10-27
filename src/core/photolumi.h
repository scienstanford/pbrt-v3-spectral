/*
 This is a class for all photoluminescent phenomenon.
 Currently we only have the fluorescence effect for tissues.
*/
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_PHOTOLUMI_H
#define PBRT_CORE_PHOTOLUMI_H

// core/photolumi.h*
#include "pbrt.h"
#include "stringprint.h"
#include <Eigen/Dense>
#include "spectrum.h"

namespace pbrt {

extern bool PhotoLumiSamplesSorted(const Float *lambda, Float **vals,
                                   int n);
extern void SortPhotoLumiSamples(Float *lambda, Float **vals, int n);
extern Float AveragePhotoLumiSamples(const Float *lambda, Float **vals,
                                     int n, Float lambdaStart, Float lambdaEnd,
                                     int j);
enum class PhotoLumiType { Reflectance, Illuminant, Display };
extern void Blackbody(const Float *lambda, int n, Float T, Float *Le);
extern void BlackBodyNormalized(const Float *lambda, int n, Float T, Float *vals);

class PhotoLumi : public Eigen::Matrix<
    Float, nSpectralSamples, nSpectralSamples> {
 public:
  static constexpr int nSamples = nSpectralSamples;
  PhotoLumi() : Eigen::Matrix<Float, nSpectralSamples, nSpectralSamples>(
      PhotoLumi::Zero()) { }

  template<typename OtherDerived>
  PhotoLumi(const Eigen::MatrixBase<OtherDerived>& other) :
      Eigen::Matrix<Float, nSpectralSamples, nSpectralSamples>(other) {
  }

  template<typename OtherDerived>
  PhotoLumi& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
    this->Eigen::Matrix<
        Float, nSpectralSamples, nSpectralSamples>::operator=(other);
    return *this;
  }

  static PhotoLumi FromSampled(const Float *lambda, Float **v,
      int n) {
    // Sort samples if unordered, use sorted for returned photolumi.
    // This part might need to be changed as we assume wavelength is sorted.
    if (!PhotoLumiSamplesSorted(lambda, v, n)) {
      std::vector<Float> slambda(&lambda[0], &lambda[n]);
      SortPhotoLumiSamples(&slambda[0], v, n);
      return FromSampled(&slambda[0], v, n);
    }
    PhotoLumi p;
    for (int i = 0; i < nSpectralSamples; ++i) {
      Float lambda0 = Lerp(Float(i) / Float(nSpectralSamples),
                           sampledLambdaStart, sampledLambdaEnd);
      Float lambda1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                           sampledLambdaStart, sampledLambdaEnd);
      for (int j = 0; j < nSpectralSamples; ++j)
        p(i, j) = AveragePhotoLumiSamples(lambda, v, n, lambda0, lambda1, i);
    }
    return p;
  }

  static void Init() {

  }

  bool IsBlack() const {
    return this->cwiseAbs().maxCoeff() == 0.0;
  }
};

}  // namespace pbrt

#endif // PBRT_CORE_PHOTOLUMI_H
