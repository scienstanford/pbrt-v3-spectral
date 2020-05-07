
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_SPECTRUM_H
#define PBRT_CORE_SPECTRUM_H

// core/spectrum.h*
#include "pbrt.h"
#include "stringprint.h"

#include <Eigen/Dense>

namespace pbrt {

// Spectrum Utility Declarations
extern bool SpectrumSamplesSorted(
    const Float* lambda, const Float* vals,
    int n);

extern void SortSpectrumSamples(Float* lambda, Float* vals, int n);

extern Float AverageSpectrumSamples(
    const Float* lambda, const Float* vals,
    int n, Float lambdaStart, Float lambdaEnd);

inline void XYZToRGB(const Float xyz[3], Float rgb[3]) {
  rgb[0] = 3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
  rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
  rgb[2] = 0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];
}

inline void RGBToXYZ(const Float rgb[3], Float xyz[3]) {
  xyz[0] = 0.412453f * rgb[0] + 0.357580f * rgb[1] + 0.180423f * rgb[2];
  xyz[1] = 0.212671f * rgb[0] + 0.715160f * rgb[1] + 0.072169f * rgb[2];
  xyz[2] = 0.019334f * rgb[0] + 0.119193f * rgb[1] + 0.950227f * rgb[2];
}

enum class SpectrumType {
  Reflectance, Illuminant, Display, Mouth, Loaded
};

extern Float InterpolateSpectrumSamples(
    const Float* lambda, const Float* vals,
    int n, Float l);

extern void Blackbody(const Float* lambda, int n, Float T, Float* Le);

extern void BlackbodyNormalized(
    const Float* lambda, int n, Float T,
    Float* vals);

// Spectral Data Declarations
static const int nLCDSamples = 101; // Added by TL for LCD-Apple display
extern const Float lcdApple_r[nLCDSamples];
extern const Float lcdApple_g[nLCDSamples];
extern const Float lcdApple_b[nLCDSamples];
extern const Float lcdApple_lambda[nLCDSamples];

static const int nCIESamples = 471;
extern const Float CIE_X[nCIESamples];
extern const Float CIE_Y[nCIESamples];
extern const Float CIE_Z[nCIESamples];
extern const Float CIE_lambda[nCIESamples];
static const Float CIE_Y_integral = 106.856895;
static const int nRGB2SpectSamples = 32;
extern const Float RGB2SpectLambda[nRGB2SpectSamples];
extern const Float RGBRefl2SpectWhite[nRGB2SpectSamples];
extern const Float RGBRefl2SpectCyan[nRGB2SpectSamples];
extern const Float RGBRefl2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBRefl2SpectYellow[nRGB2SpectSamples];
extern const Float RGBRefl2SpectRed[nRGB2SpectSamples];
extern const Float RGBRefl2SpectGreen[nRGB2SpectSamples];
extern const Float RGBRefl2SpectBlue[nRGB2SpectSamples];
extern const Float RGBIllum2SpectWhite[nRGB2SpectSamples];
extern const Float RGBIllum2SpectCyan[nRGB2SpectSamples];
extern const Float RGBIllum2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBIllum2SpectYellow[nRGB2SpectSamples];
extern const Float RGBIllum2SpectRed[nRGB2SpectSamples];
extern const Float RGBIllum2SpectGreen[nRGB2SpectSamples];
extern const Float RGBIllum2SpectBlue[nRGB2SpectSamples];


class RGBSpectrum: public Eigen::Array<Float, 3, 1> {
 public:
  static constexpr int nSamples = 3;
  RGBSpectrum() : Eigen::Array<Float, 3, 1>(RGBSpectrum::Zero()) { }

  RGBSpectrum(Float v) : Eigen::Array<Float, 3, 1>(
      RGBSpectrum::Constant(v)) { }

  template<typename OtherDerived>
  RGBSpectrum(const Eigen::ArrayBase<OtherDerived>& other) :
      Eigen::Array<Float, 3, 1>(other) {
  }

  template<typename OtherDerived>
  RGBSpectrum& operator=(const Eigen::ArrayBase<OtherDerived>& other) {
    this->Eigen::Array<Float, 3, 1>::operator=(other);
    return *this;
  }

  RGBSpectrum(
      const RGBSpectrum& s,
      SpectrumType type = SpectrumType::Reflectance) {
    *this = s;
  }

  RGBSpectrum Clamp(Float low = 0, Float high = Infinity) const {
    return this->min(high).max(low);
  }

  bool IsBlack() const {
    return this->abs().maxCoeff() == 0.0f;
  }

  bool HasNaNs() const {
    for (int i = 0; i < 3; ++i)
      if (std::isnan((*this)[i])) return true;
    return false;
  }

  static RGBSpectrum FromRGB(
      const Float rgb[3],
      SpectrumType type = SpectrumType::Reflectance) {
    RGBSpectrum s;
    s[0] = rgb[0];
    s[1] = rgb[1];
    s[2] = rgb[2];
    DCHECK(!s.HasNaNs());
    return s;
  }

  void ToRGB(Float* rgb) const {
    rgb[0] = (*this)[0];
    rgb[1] = (*this)[1];
    rgb[2] = (*this)[2];
  }

  const RGBSpectrum& ToRGBSpectrum() const { return *this; }

  void ToXYZ(Float xyz[3]) const {
    Float rgb[3];
    ToRGB(rgb);
    RGBToXYZ(rgb, xyz);
  }

  static RGBSpectrum FromXYZ(
      const Float xyz[3],
      SpectrumType type = SpectrumType::Reflectance) {
    Float rgb[3];
    XYZToRGB(xyz, rgb);
    return FromRGB(rgb);
  }

  Float y() const {
    Eigen::Vector3f YWeight;
    YWeight << 0.212671f, 0.715160f, 0.072169f;
    return YWeight.dot(this->matrix());
  }

  static RGBSpectrum FromSampled(const Float* lambda, const Float* v, int n) {
    // Sort samples if unordered, use sorted for returned spectrum
    if (!SpectrumSamplesSorted(lambda, v, n)) {
      std::vector<Float> slambda(&lambda[0], &lambda[n]);
      std::vector<Float> sv(&v[0], &v[n]);
      SortSpectrumSamples(&slambda[0], &sv[0], n);
      return FromSampled(&slambda[0], &sv[0], n);
    }
    Float xyz[3] = {0, 0, 0};
    for (int i = 0; i < nCIESamples; ++i) {
      Float val = InterpolateSpectrumSamples(lambda, v, n, CIE_lambda[i]);
      xyz[0] += val * CIE_X[i];
      xyz[1] += val * CIE_Y[i];
      xyz[2] += val * CIE_Z[i];
    }
    Float scale = Float(CIE_lambda[nCIESamples - 1] - CIE_lambda[0]) /
                  Float(CIE_Y_integral * nCIESamples);
    xyz[0] *= scale;
    xyz[1] *= scale;
    xyz[2] *= scale;
    return FromXYZ(xyz);
  }
};


class Spectrum: public Eigen::Array<Float, nSpectralSamples, 1> {
 public:
  static constexpr int nSamples = nSpectralSamples;
  Spectrum() : Eigen::Array<Float, nSpectralSamples, 1>(Spectrum::Zero()) { }

  explicit Spectrum(Float v) : Eigen::Array<Float, nSpectralSamples, 1>(
      Spectrum::Constant(v)) { }

  template<typename OtherDerived>
  Spectrum(const Eigen::ArrayBase<OtherDerived>& other) :
      Eigen::Array<Float, nSpectralSamples, 1>(other) {
  }

  template<typename OtherDerived>
  Spectrum& operator=(const Eigen::ArrayBase<OtherDerived>& other) {
    this->Eigen::Array<Float, nSpectralSamples, 1>::operator=(other);
    return *this;
  }

  static Spectrum FromSampled(const Float* lambda, const Float* v, int n) {
    // Sort samples if unordered, use sorted for returned spectrum
    if (!SpectrumSamplesSorted(lambda, v, n)) {
      std::vector<Float> slambda(&lambda[0], &lambda[n]);
      std::vector<Float> sv(&v[0], &v[n]);
      SortSpectrumSamples(&slambda[0], &sv[0], n);
      return FromSampled(&slambda[0], &sv[0], n);
    }
    Spectrum r = Spectrum::Zero();
    
    // ZLY: If sampled wavelengths are exactly equal to the hardcoded wavelengths in refWave,
    // directly copy lambda to the returned value r
    for (int i = 0; i < nSpectralSamples; ++i) {
      // If current wavelength exactly equals to reference wavelength, directly copy the value
      if (Float(refWave[i]) == lambda[i]) r[i] = v[i];
      else{
          // Compute average value of given SPD over $i$th sample's range
          Float lambda0 = Lerp(Float(i) / Float(nSpectralSamples),
                               sampledLambdaStart, sampledLambdaEnd);
          Float lambda1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                               sampledLambdaStart, sampledLambdaEnd);
          r[i] = AverageSpectrumSamples(lambda, v, n, lambda0, lambda1);
      }
    }
    return r;
  }

  static void Init() {
    // Compute XYZ matching functions for _SampledSpectrum_
    for (int i = 0; i < nSpectralSamples; ++i) {
      Float wl0 = Lerp(Float(i) / Float(nSpectralSamples),
                       sampledLambdaStart, sampledLambdaEnd);
      Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                       sampledLambdaStart, sampledLambdaEnd);
      X[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples, wl0, wl1);
      Y[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples, wl0, wl1);
      Z[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples, wl0, wl1);
    }

    // Compute RGB to spectrum functions for _SampledSpectrum_
    for (int i = 0; i < nSpectralSamples; ++i) {
      Float wl0 = Lerp(Float(i) / Float(nSpectralSamples),
                       sampledLambdaStart, sampledLambdaEnd);
      Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                       sampledLambdaStart, sampledLambdaEnd);
      rgbRefl2SpectWhite[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectWhite,
                                 nRGB2SpectSamples, wl0, wl1);
      rgbRefl2SpectCyan[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectCyan,
                                 nRGB2SpectSamples, wl0, wl1);
      rgbRefl2SpectMagenta[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectMagenta,
                                 nRGB2SpectSamples, wl0, wl1);
      rgbRefl2SpectYellow[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectYellow,
                                 nRGB2SpectSamples, wl0, wl1);
      rgbRefl2SpectRed[i] = AverageSpectrumSamples(
          RGB2SpectLambda, RGBRefl2SpectRed, nRGB2SpectSamples, wl0, wl1);
      rgbRefl2SpectGreen[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectGreen,
                                 nRGB2SpectSamples, wl0, wl1);
      rgbRefl2SpectBlue[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectBlue,
                                 nRGB2SpectSamples, wl0, wl1);

      rgbIllum2SpectWhite[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectWhite,
                                 nRGB2SpectSamples, wl0, wl1);
      rgbIllum2SpectCyan[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectCyan,
                                 nRGB2SpectSamples, wl0, wl1);
      rgbIllum2SpectMagenta[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectMagenta,
                                 nRGB2SpectSamples, wl0, wl1);
      rgbIllum2SpectYellow[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectYellow,
                                 nRGB2SpectSamples, wl0, wl1);
      rgbIllum2SpectRed[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectRed,
                                 nRGB2SpectSamples, wl0, wl1);
      rgbIllum2SpectGreen[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectGreen,
                                 nRGB2SpectSamples, wl0, wl1);
      rgbIllum2SpectBlue[i] =
          AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectBlue,
                                 nRGB2SpectSamples, wl0, wl1);
    }

    // Added by TL
    // Compute LCD-Apple SPD for _SampledSpectrum_
    for (int i = 0; i < nSpectralSamples; ++i) {
      Float wl0 = Lerp(Float(i) / Float(nSpectralSamples),
                       sampledLambdaStart, sampledLambdaEnd);
      Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                       sampledLambdaStart, sampledLambdaEnd);
      R[i] = AverageSpectrumSamples(
          lcdApple_lambda, lcdApple_r, nLCDSamples, wl0, wl1);
      G[i] = AverageSpectrumSamples(
          lcdApple_lambda, lcdApple_g, nLCDSamples, wl0, wl1);
      B[i] = AverageSpectrumSamples(
          lcdApple_lambda, lcdApple_b, nLCDSamples, wl0, wl1);
    }
      
    // Added by ZLY
    // Prepared for using Mouth Basis function for _SampledSpectrum_
    for (int i = 0; i < nSpectralSamples; ++i) {
      Float wl0 = Lerp(Float(i) / Float(nSpectralSamples),
                       sampledLambdaStart, sampledLambdaEnd);
      Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                       sampledLambdaStart, sampledLambdaEnd);
      MouthR[i] = AverageSpectrumSamples(
          lcdApple_lambda, lcdApple_r, nLCDSamples, wl0, wl1);
      MouthG[i] = AverageSpectrumSamples(
          lcdApple_lambda, lcdApple_g, nLCDSamples, wl0, wl1);
      MouthB[i] = AverageSpectrumSamples(
          lcdApple_lambda, lcdApple_b, nLCDSamples, wl0, wl1);
    }
  }

  Spectrum Clamp(Float low = 0, Float high = Infinity) const {
    return this->min(high).max(low);
  }

  bool IsBlack() const {
    return this->abs().maxCoeff() == 0.0f;
  }

  bool HasNaNs() const {
    for (int i = 0; i < nSpectralSamples; ++i)
      if (std::isnan((*this)[i])) return true;
    return false;
  }

  RGBSpectrum ToRGBSpectrum() const;

  void ToXYZ(Float xyz[3]) const {
    xyz[0] = ((*this) * X).sum();
    xyz[1] = ((*this) * Y).sum();
    xyz[2] = ((*this) * Z).sum();
    Float scale = Float(sampledLambdaEnd - sampledLambdaStart) /
                  Float(CIE_Y_integral * nSpectralSamples);
    xyz[0] *= scale;
    xyz[1] *= scale;
    xyz[2] *= scale;
  }

  Float y() const {
    Float yy = ((*this) * Y).sum();
    yy = (yy < 0) ? 0 : yy;
    return yy * Float(sampledLambdaEnd - sampledLambdaStart) /
           Float(CIE_Y_integral * nSpectralSamples);
  }

  void ToRGB(Float rgb[3]) const {
    Float xyz[3];
    ToXYZ(xyz);
    XYZToRGB(xyz, rgb);
  }

  static Spectrum FromRGB(
      const Float rgb[3], SpectrumType type = SpectrumType::Illuminant,
        Spectrum basisOne = Spectrum(0.f),
        Spectrum basisTwo = Spectrum(0.f),
        Spectrum basisThree = Spectrum(0.f));

  explicit Spectrum(const RGBSpectrum &r,
      SpectrumType type = SpectrumType::Reflectance);

  static Spectrum FromXYZ(
      const Float xyz[3], SpectrumType type = SpectrumType::Reflectance) {
    Float rgb[3];
    XYZToRGB(xyz, rgb);
    return FromRGB(rgb, type);
  }

  std::string ToString() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
  }

  //Trisha Added (6-2016)
  // Get the value in the spectrum at a specific wavelength.
  void GetValueAtWavelength(Float wavelength, Float* output) const {
    Float w0;
    Float w1;
    Float t;
    *output = 0;

    // A rare case but let's catch it
    if (wavelength == sampledLambdaEnd) {
      *output = (*this)[nSpectralSamples - 1];
    }

    for (int i = 0; i < nSpectralSamples; i++) {

      w0 = Lerp(Float(i) / Float(nSpectralSamples),
                sampledLambdaStart, sampledLambdaEnd);
      w1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                sampledLambdaStart, sampledLambdaEnd);

      if ((wavelength >= w0) && (wavelength < w1)) {
        t = (wavelength - w0) / (w1 - w0);
        *output = Lerp(t, (*this)[i], (*this)[i + 1]);
        return;
      }
    }
    // Make sure we don't return 0. If we do, then we get NaN's!
    // I commented this out, because what if the spectrum is zero?
    //assert(*output != 0);
  }

 private:
  // SampledSpectrum Private Data
  static Spectrum X, Y, Z;
  static Spectrum R, G, B; //Added by TL
  static Spectrum MouthR, MouthG, MouthB; // Added by ZLY
  static Spectrum TempR, TempG, TempB;   // Added by ZLY
  static Spectrum rgbRefl2SpectWhite, rgbRefl2SpectCyan;
  static Spectrum rgbRefl2SpectMagenta, rgbRefl2SpectYellow;
  static Spectrum rgbRefl2SpectRed, rgbRefl2SpectGreen;
  static Spectrum rgbRefl2SpectBlue;
  static Spectrum rgbIllum2SpectWhite, rgbIllum2SpectCyan;
  static Spectrum rgbIllum2SpectMagenta, rgbIllum2SpectYellow;
  static Spectrum rgbIllum2SpectRed, rgbIllum2SpectGreen;
  static Spectrum rgbIllum2SpectBlue;
};


inline RGBSpectrum Lerp(Float t, const RGBSpectrum& s1, const RGBSpectrum& s2) {
  return (1 - t) * s1 + t * s2;
}

inline Spectrum Lerp(Float t, const Spectrum& s1, const Spectrum& s2) {
  return (1 - t) * s1 + t * s2;
}

void ResampleLinearSpectrum(
    const Float* lambdaIn, const Float* vIn, int nIn,
    Float lambdaMin, Float lambdaMax, int nOut,
    Float* vOut);

}  // namespace pbrt

#endif  // PBRT_CORE_SPECTRUM_H
