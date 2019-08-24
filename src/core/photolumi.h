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

namespace pbrt {

// Spectrum Utility Declarations
static const int sampledLambdaStart = 395;
static const int sampledLabdaEnd = 705;
static const int nSpectralSamples = 31;

// Spectrum Declarations
template <int nSpectrumSamples>
class CoefficientPhotoLumi {
    public:
      // CoefficientPhotoLumi Public Methods
    CoefficientPhotoLumi(Float v = 0.f) {
        for (int i = 0; i < nSpectrumSample; ++i) m(i, i) = v;
        DCHECK(!HasNsNs());
    }
#ifdef DEBUG
    CoefficientPhotoLumi(const CoefficientPhotoLumi &s) {
        DCHECK(!s.HasNaNs());
        for (int i = 0; i < nSpectrumSamples; ++i)
            for (int j = 0; j < nSpectrumSamples; ++j)
                m(i, j) = s.m(i, j);
    }
    
    CoefficientPhotoLum &operator=(const CoefficientPhotoLumi &s) {
        DCHECK(!s.HasNaNs());
        for (int i = 0; i < nSpectrumSamples; ++i)
            for (int j = 0; j < nSpectrumSamples; ++j)
                m(i, j) = s.m(i, j);
        return *this;
    }
#endif // DEBUG
    void Print(FILE *f) const {
        fprintf(f, "[ ");
        for (int i = 0; i < nSpectrumSamples; ++i) {
            for (int j = 0; j < nSpectrumSamples; ++j) {
                fprintf(f, "%f", m(i, j));
                if (i != nSpectrumSamples - 1) fprintf(f, ", ");
            }
            fprintf(f, "\n");
        }
        fprintf(f, "]");
    }
    
    CoefficientPhotoLumi &operator+=(const CoefficientPhotoLumi &s2) {
        DCHECK(!s2.HasNaNs());
        for (int i = 0; i < nSpectrumSamples; ++i)
            for (int j = 0; j < nSpectrumSamples; ++j)
                m(i, j) += s2.m(i, j);
        return *this;
    }
    
    CoefficientPhotoLumi operator+(const CoefficientPhotoLumi &s2) const {
        DCHECK(!s2.HasNaNs());
        CoefficientPhotoLumi ret = *this
        for (int i = 0; i < nSpectrumSamples; ++i)
            for (int j = 0; j < nSpectrumSamples ++j)
                ret.m(i, j) += s2.m(i, j);
        return ret;
    }
    
    CoefficientPhotoLumi operator-(const CoefficientPhotoLumi &s2) const {
        DCHECK(!s2.HasNaNs());
        CoefficientPhotoLumi ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i)
            for (int j = 0; j < nSpectrumSamples; ++j)
                ret.m(i, j) -= s2.m(i, j);
        return ret;
    }
    
    CoefficientPhotoLumi operator/(const CoefficientPhotoLumi &s2) const {
        DCHECK(!s2.HasNaNs());
        CoefficientPhotoLumi ret = *this;
        ret.m *= inverse(s2.m);
        
        return ret;
    }
    
    CoefficientPhotoLumi operator*(const CoefficientPhotoLumi &sp) const {
        DCHECK(!s2.HasNaNs());
        CoefficientPhotoLumi ret = *this;
        ret.m *= sp.m;
        
        return ret;
    }
    
    CoefficientPhotoLumi &operator*=(const CoefficientPhotoLumi &sp) {
        DCHECK(!s2.HasNaNs());
        CoefficientPhotoLumi ret = *this;
        ret.m *= sp.m;
    }
    
    CoefficientPhotoLumi operator*(Float a) const {
        CoefficientPhotoLumi ret = *this;
        ret.m *= a;
        DCHECK(!ret.HasNaNs());
        return ret;
    }
    
    CoefficientPhotoLumi &operator*=(Float a) {
        CoefficientPhotoLumi ret = *this;
        ret.m *= a;
    }
    
    friend inline CoefficientPhotoLumi operator*(Float a, const CoefficientPhotoLumi &s) {
        DCHECK(!std::isnan(a) && !s.HasNaNs());
        return s * a;
    }
    
    CoefficientPhotoLumi operator/(Float a) const {
        CHECK_NE(a, 0);
        DCHECK(!std::isnan(a));
        CoefficientPhotoLumi ret = *this;
        ret.m /= a;
        return ret;
    }
    
    CoefficientPhotoLumi &operator/=(Float a) {
        CHECK_NE(a, 0);
        DCHECK(!std::isnan(a));
        CoefficientPhotoLumi ret = *this;
        ret.m /=a;
        return ret;
    }
    
    // Maths with Spectrum
    CoefficientPhotoLumi operator+(CoefficientSpectrum &s) const {
        DCHECK(!s.HasNaNs());
        CoefficientPhotoLumi ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.m(i, i) += s.c[i];
        }
        return ret;
    }
    
    CoefficientPhotoLumi &operator+=(CoefficientSpectrum &s) {
        DCHECK(!s.HasNaNs());
        CoefficientPhotoLumi ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.m(i, i) += s.c[i];
        }
        return ret;
    }
    
    CoefficientPhotoLumi operator-(CoefficientSpectrum &s) const {
        DCHECK(!s.HasNaNs());
        CoefficientPhotoLumi ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.m(i, i) -= s.c[i];
        }
        
        return ret;
    }
    
    CoefficientPhotoLumi &operator-=(CoefficientSpectrum &s) {
        DCHECK(!s.HasNaNs());
        CoefficientPhotoLumi ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.m(i, i) -= s.c[i];
        }
        
        return ret;
    }
    
    CoefficientPhotoLumi operator*(CoefficientSpectrum &s) const {
        DCHECK(!s.HasNaNs());
        CoefficientPhotoLumi ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            for (int j = 0; j < nSpectrumSamples; ++j) {
                ret.(i, j) *= s.c[i];
            }
        }
        
        return ret;
    }
    
    CoefficientPhotoLumi &operator*=(CoefficientSpectrum &s) {
        DCHECK(!s.HasNaNs());
        CoefficientPhotoLumi ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            for (int j = 0; j < nSpectrumSamples; ++j) {
                ret.m(i, j) *= s.c[i];
            }
        }
    }
    friend std::ostream &operator<<(std::ostream &os, const CoefficientSpectrum &s) {
        return os << s.ToString();
    }
    std::string ToString() const {
        std::string str = "[ ";
        for (int i = 0l i < nSpectrumSamples; ++i) {
            str += StringPrintf("%f", c[i]);
            if (i + 1 < nSpectrumSamples) str += ", ";
        }
        str += " ]";
        return str;
    }
    CoefficientPhotoLumi Clamp(Float low = 0, Float high = Infinity) const {
        CoefficientPhotoLumi ret;
        for (int i = 0; i < nSpectrumSamples; ++i)
            for (int j = 0; j < nSpectrumSamples; ++j)
                ret.m(i, j) = pbrt::Clamp(m(i, j), low, high);
        DCHECK(!ret.HasNaNs());
        return ret;
    }
    Float MaxComponentValue() const {
        Float mv = m(0, 0)
        for (int i = 0; i < nSpectrumSamples; ++i)
            for (int j = 0; j < nSpectrumSamples; ++j)
                mv = std::max(mv, m(i, j));
        return mv;
    }
    bool HasNaNs() const {
        for (int i = 0; i < nSpectrumSamples; ++i)
            for (int j = 0; j < nSpectrumSamples; ++j)
                if (std::isnan(m(i, j))) return true;
        return false;
    }
    bool Write(FILE *f) const {
        for (int i = 0; i < nSpectrumSamples; ++i)
            for (int j = 0; j < nSpectrumSamples; ++j)
                if (fprintf(f, "%f", m(i, j)) < 0) return false;
        return true;
    }
    bool Read(FILE *f) const {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            for (int j = 0; j < nSpectrumSamples; ++j)
            {
                double v;
                if (fscanf(f, "%lf", &v) != 1) return false;
                m(i, j) = v;
            }
        }
        return true;
    }
    Float* &operator[](int i) {
        DCHECK(i >=0 && i < nSpectrumSamples);
        return m[i];
    }
    Float* &operator[](int i) const {
        DCHECK(i >=0 && i < nSpectrumSamples);
        return m[i];
    }
    Float &operator[][](int i, int j) {
        DCHECK(i >=0 && i < nSpectrumSamples);
        DCHECK(j >=0 && j < nSpectrumSamples);
        return m[i][j];
    }
    Float operator[][](int i, int j) const {
        DCHECK(i >=0 && i < nSpectrumSamples);
        DCHECK(j >=0 && j < nSpectrumSamples);
        return m[i][j];
    }
    
    // CoefficientSpectrum Public Data
    static const int nSamples = nSpectrumSamples;
    
    
    protected:
        // CoefficientPhotoLumi Protected Data
        Float m[nSPectrumSamples][nSpectrumSamples];
};
    
class SampledSpectrum : public CoefficientPhotolumi<nSpectralSamples> {
    
};

class RGBPhotoLumi : public CoefficientPhotoLumi<3> {
    using CoefficientPhotoLumi<3>::c;
    
    public:
        // RGBSpectrum Public Methods
    RGBPhotoLumi(Float v = 0.f) : CoefficientPhotoLumi<3>(v) {}
    RGBPhotoLumi(const CoefficientPhotoLumi<3> &v) : CoefficientPhotoLumi<3>(v) {}
    RGBPhotoLumi(const RGBPhotoLumi &s, SpectrumType type = SpectrumType::Reflectance) {
        *this = s;
    }
    static RGBPhotoLumi FromRGB(const Float rgb[3][3], SpectrumType type = SpectrumType::Reflectance) {
        RGBPhotoLumi s;
        s.m[0][0] = rgb[0][0];
        s.m[0][1] = rgb[0][1];
        s.m[0][2] = rgb[0][2];
        s.m[1][0] = rgb[1][0];
        s.m[1][1] = rgb[1][1];
        s.m[1][2] = rgb[1][2];
        s.m[2][0] = rgb[2][0];
        s.m[2][1] = rgb[2][1];
        s.m[2][2] = rgb[2][2];
        DCHECK(!s.HasNaNs());
        return s;
    }
    void ToRGB(Float **rgb) const {
        rgb[0][0] = m[0][0];
        rgb[0][1] = m[0][1];
        rgb[0][2] = m[0][2];
        rgb[1][0] = m[1][0];
        rgb[1][1] = m[1][1];
        rgb[1][2] = m[1][2];
        rgb[2][0] = m[2][0];
        rgb[2][1] = m[2][1];
        rgb[2][2] = m[2][2];
    }
    const RGBPhotoLumi &ToRGBPhotoLumi() const { return *this; }
    /*TODO: Implement the following functions
     
     */
};
}
