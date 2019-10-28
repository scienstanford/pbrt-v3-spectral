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

#include "pbrt.h"
#include "stringprint.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/Splines>
#include "spectrum.h"

namespace pbrt {

enum class PhotoLumiType { Reflectance, Illuminant, Display };

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
    DCHECK(lambda.size() == v.size() && lambda.size() > 1);

    // The values in the input matrix represents a ratio between in and out
    // light. Thus, the wavelength interpolation should be similar to the a
    // normal interpolation and is different from the one in Spectrum class.
    // However, special care is needed and we cannot use vanilla 2d
    // interpolation:
    //   We would like non-fluorescent (diagonal) matrix to remain diagonal
    //   after interpolation. That is, we don't want to introduce non-existence
    //   fluorescent in the interpolation.
    // Thus, we perform fluorescent and normal reflection (i.e. the diagonal)
    // interpolated independently.

    // Store non-fluorescent part to be interpolated later.
    PhotoLumi result;
    Eigen::Array<Float, 1, Eigen::Dynamic> nonFluorescent(lambda.size());
    for (int i = 0; i < lambda.size(); i++) {
      nonFluorescent(i) = v[i][i];
    }

    // Compute scaled / normalized input and output wavelength.
    Eigen::Array<Float, 1, Eigen::Dynamic> scaled_lambda = (Eigen::Map<
        const Eigen::Array<Float, 1, Eigen::Dynamic>>(
            lambda.data(), lambda.size()) - lambda[0]) / (
                lambda.back() - lambda[0]);
    Eigen::Array<Float, 1, nSpectralSamples> scaled_output_lambda = (
        Eigen::Array<Float, 1, nSpectralSamples>::LinSpaced(
            sampledLambdaStart, sampledLambdaEnd, nSpectralSamples) -
            lambda[0]) / (lambda.back() - lambda[0]);
    scaled_output_lambda = scaled_output_lambda.max(0.0f).min(1.0f);

    // Interpolate for fluorescent part.
    Eigen::Array<Float, Eigen::Dynamic, nSpectralSamples> intermediate_result(
        lambda.size(), nSpectralSamples);
    for (int i = 0; i < lambda.size(); ++i) {
      Eigen::Array<Float, 1, Eigen::Dynamic> current_row = Eigen::Map<
          const Eigen::Array<Float, 1, Eigen::Dynamic>>(
              v[i].data(), lambda.size());
      // remove non-fluorescent part.
      // Maybe it's better to actual remove this data point in interpolation.
      // But assuming that fluorescent happens at far enough wavelength from
      // input light, it should be safe to just set it to zero.
      current_row[i] = 0.0f;
      Eigen::Spline<Float, 1> spline(Eigen::SplineFitting<
          Eigen::Spline<Float, 1>>::Interpolate(current_row, 2, scaled_lambda));
      for (int j = 0; j < nSpectralSamples; j++) {
        intermediate_result(i, j) = spline(scaled_output_lambda(j))(0);
      }
    }

    for (int j = 0; j < nSpectralSamples; j++) {
      Eigen::Spline<Float, 1> spline(Eigen::SplineFitting<
          Eigen::Spline<Float, 1>>::Interpolate(
              intermediate_result.col(j).transpose(), 2, scaled_lambda));
      for (int i = 0; i < nSpectralSamples; i++) {
        result(i, j) = spline(scaled_output_lambda(i))(0);
      }
    }

    // Interpolate for non-fluorescent part.
    Eigen::Spline<Float, 1> spline(Eigen::SplineFitting<
        Eigen::Spline<Float, 1>>::Interpolate(
            nonFluorescent, 2, scaled_lambda));
    for (int i = 0; i < nSpectralSamples; i++) {
      result(i, i) = spline(scaled_output_lambda(i))(0);
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
