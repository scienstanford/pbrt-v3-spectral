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
#include <unsupported/Eigen/Splines>
#include "spectrum.h"

namespace pbrt {

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

  static PhotoLumi FromSampled(const std::vector<Float>& lambda,
      const std::vector<std::vector<Float>> &v) {
    DCHECK(lambda.size() == v.size());

    // Perform average interpolation along each direction independently.
    Eigen::Array<Float, Eigen::Dynamic, nSpectralSamples> intermediate_result(
        lambda.size(), nSpectralSamples);
    for (int i = 0; i < lambda.size(); ++i) {
      for (int j = 0; j < nSpectralSamples; ++j) {
        Float lambda0 = Lerp(Float(j) / Float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
        Float lambda1 = Lerp(Float(j + 1) / Float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
        intermediate_result(i, j) = AverageSpectrumSamples(
            lambda.data(), v[i].data(), lambda.size(), lambda0, lambda1);
      }
    }

    PhotoLumi result;
    for (int j = 0; j < nSpectralSamples; j++) {
      auto v_column =
          intermediate_result.col(j);

      for (int i = 0; i < nSpectralSamples; i++) {
        Float lambda0 = Lerp(Float(i) / Float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
        Float lambda1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
        result(i, j) = AverageSpectrumSamples(
            lambda.data(), v_column.data(), lambda.size(), lambda0, lambda1);
      }
    }
    return result;
  }

  static void Init() {

  }

  bool IsBlack() const {
    return this->cwiseAbs().maxCoeff() == 0.0;
  }
};

}  // namespace pbrt

#endif // PBRT_CORE_PHOTOLUMI_H
