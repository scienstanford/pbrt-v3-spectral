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

// Spectrum Utility Declarations
//static const int sampledLambdaStart = 395;
//static const int sampledLabdaEnd = 705;
//static const int nSpectralSamples = 31;
//extern bool SpectrumSamplesSorted(const Float *lambda, const Float *vals,
//                                  int n);
//extern void SortSpectrumSamples(Float *lambda, Float *vals, int n);
//extern Float AverageSpectrumSamples(const Float *lambda, const Float *vals,
//                                    int n, Float lambdaStart, Float lambdaEnd);
extern bool PhotoLumiSamplesSorted(const Float *lambda, Float **vals,
                                   int n);
extern void SortPhotoLumiSamples(Float *lambda, Float **vals, int n);
extern Float AveragePhotoLumiSamples(const Float *lambda, Float **vals,
                                     int n, Float lambdaStart, Float lambdaEnd,
                                     int j);
enum class PhotoLumiType { Reflectance, Illuminant, Display };
//extern Float InterpolateSpectrumSamples(const Float *lambda, const Float *vals,
//                                        int n, Float l);
extern void Blackbody(const Float *lambda, int n, Float T, Float *Le);
extern void BlackBodyNormalized(const Float *lambda, int n, Float T, Float *vals);
// Spectrum Declarations
template <int nSpectrumSamples>
class CoefficientPhotoLumi {

    public:
      // CoefficientPhotoLumi Public Methods
        CoefficientPhotoLumi(Float v = 0.f) {
            for (int i = 0; i < nSpectrumSamples; ++i)
                for (int j = 0; j < nSpectrumSamples; ++j) {
                    if (i == j) m(i, j) = v;
                    else m(i, j) = 0.f;
                }
            DCHECK(!HasNaNs());
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
            CoefficientPhotoLumi ret = *this;
            for (int i = 0; i < nSpectrumSamples; ++i)
                for (int j = 0; j < nSpectrumSamples; ++j)
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
    
        CoefficientPhotoLumi operator*(const CoefficientPhotoLumi &s2) const {
            DCHECK(!s2.HasNaNs());
            CoefficientPhotoLumi ret = *this;
            ret.m *= s2.m;
            
            return ret;
        }
    
        CoefficientPhotoLumi &operator*=(const CoefficientPhotoLumi &s2) {
            DCHECK(!s2.HasNaNs());
//            CoefficientPhotoLumi ret = *this;
//            ret.m *= s2.m;
//            return ret;
            m *= s2.m;
            return *this;
        }
    
        CoefficientPhotoLumi operator*(Float a) const {
            CoefficientPhotoLumi ret = *this;
            ret.m *= a;
            DCHECK(!ret.HasNaNs());
            return ret;
        }
    
        CoefficientPhotoLumi &operator*=(Float a) {
//            CoefficientPhotoLumi ret = *this;
//            ret.m /= a;
//            return ret;
            m *= a;
            return *this;
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
//            CoefficientPhotoLumi ret = *this;
//            ret.m /= a;
//            return ret;
            m /= a;
            return *this;
        }
    
        friend std::ostream &operator<<(std::ostream &os, const CoefficientPhotoLumi &s) {
            return os << s.ToString();
        }
        std::string ToString() const {
            std::string str = "[ ";
            for (int i = 0; i < nSpectrumSamples; ++i) {
                for (int j = 0; j < nSpectrumSamples; ++j) {
                    str += StringPrintf("%f", m(i, j));
                    if (j + 1 < nSpectrumSamples) str += ", ";
                }
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
            Float mv = m(0, 0);
            for (int i = 0; i < nSpectrumSamples; ++i)
                for (int j = 0; j < nSpectrumSamples; ++j)
                    mv = std::max(mv, m(i, j));
            return mv;
        }
        bool HasNaNs() const {
            for (int i = 0; i < nSpectrumSamples; ++i)
                for (int j = 0; j < nSpectrumSamples; ++j) {
                    if (std::isnan(m(i, j))) return true;
                }

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

//        Float &operator()(int i, int j) {
//            DCHECK(i >=0 && i < nSpectrumSamples);
//            DCHECK(j >=0 && j < nSpectrumSamples);
//            return m(i, j);
//        }
        Float operator()(int i, int j) {
            DCHECK(i >=0 && i < nSpectrumSamples);
            DCHECK(j >=0 && j < nSpectrumSamples);
            return m(i, j);
        }
    
        // Maths with Spectrum
        CoefficientPhotoLumi operator+(Spectrum &s) const {
            DCHECK(!s.HasNaNs());
            CoefficientPhotoLumi ret = *this;
            for (int i = 0; i < nSpectralSamples; ++i) {
                ret.m(i, i) += s.GetValueAtIndex(i);
            }
            return ret;
        }
    
        CoefficientPhotoLumi &operator+=(Spectrum &s) {
            DCHECK(!s.HasNaNs());
            for (int i = 0; i < nSpectralSamples; ++i) {
                m(i, i) = m(i, i) + s.GetValueAtIndex(i);
            }
            return *this;
        }
    
        CoefficientPhotoLumi operator-(Spectrum &s) const {
            DCHECK(!s.HasNaNs());
            CoefficientPhotoLumi ret = *this;
            for (int i = 0; i < nSpectralSamples; ++i) {
                ret.m(i, i) -= s.GetValueAtIndex(i);
            }
            
            return ret;
        }
    
        CoefficientPhotoLumi &operator-=(Spectrum &s) {
            DCHECK(!s.HasNaNs());
            CoefficientPhotoLumi ret = *this;
            for (int i = 0; i < nSpectralSamples; ++i) {
                ret.m(i, i) -= s.GetValueAtIndex(i);
            }
            
            return ret;
        }
    
        CoefficientPhotoLumi operator*(Spectrum &s) const {
            DCHECK(!s.HasNaNs());
            CoefficientPhotoLumi ret = *this;
            for (int i = 0; i < nSpectralSamples; ++i) {
                for (int j = 0; j < nSpectralSamples; ++j) {
                    ret.m(i, j) *= s.GetValueAtIndex(i);
                }
            }
            
            return ret;
        }
    
        CoefficientPhotoLumi &operator*=(Spectrum &s) {
            DCHECK(!s.HasNaNs());
            CoefficientPhotoLumi &ret = *this;
            for (int i = 0; i < nSpectralSamples; ++i) {
                for (int j = 0; j < nSpectralSamples; ++j) {
                    ret.m(i, j) *= s.GetValueAtIndex(i);
                }
            }
            return ret;
        }
        Eigen::Matrix<Float, nSpectrumSamples, nSpectrumSamples> getMatrix() {
            return m;
        }
        
        // CoefficientSpectrum Public Data
        static const int nSamples = nSpectrumSamples;
    protected:
        // CoefficientPhotoLumi Protected Data
        Eigen::Matrix<Float, nSpectrumSamples, nSpectrumSamples> m;

};
    
class SampledPhotoLumi : public CoefficientPhotoLumi<nSpectralSamples> {
    public:
      // SampledPhotoLumi Public Methods
//        using CoefficientPhotoLumi::operator*=;
//        using CoefficientPhotoLumi::operator*;
//        using CoefficientPhotoLumi::operator+;
//        using CoefficientPhotoLumi::operator+=;
//        using CoefficientPhotoLumi::operator-;
//        using CoefficientPhotoLumi::operator-=;
        SampledPhotoLumi(Float v = 0.f) : CoefficientPhotoLumi(v) {}
        SampledPhotoLumi(const CoefficientPhotoLumi<nSpectralSamples> &v)
            : CoefficientPhotoLumi<nSpectralSamples>(v) {}
        static SampledPhotoLumi FromSampled(const Float *lambda, Float **v,
                                           int n) {
            // Sort samples if unordered, use sorted for returned photolumi
            // This part might need to be changed as we are assuming the wavelength is in order
            if (!PhotoLumiSamplesSorted(lambda, v, n)) {
                std::vector<Float> slambda(&lambda[0], &lambda[n]);
                SortPhotoLumiSamples(&slambda[0], v, n);
                return FromSampled(&slambda[0], v, n);
            }
            SampledPhotoLumi p;
            for (int i = 0; i < nSpectralSamples; ++i) {
                Float lambda0 = Lerp(Float(i) / Float(nSpectralSamples),
                                     sampledLambdaStart, sampledLambdaEnd);
                Float lambda1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                                     sampledLambdaStart, sampledLambdaEnd);
                for (int j = 0; j < nSpectralSamples; ++j)
                    p.m(i, j) = AveragePhotoLumiSamples(lambda, v, n, lambda0, lambda1, i);
            }
            return p;
//            SampledPhotoLumi p;
//            for (int i = 0; i < nSpectralSamples; ++i) {
//                for (int j = 0; j < nSpectralSamples; ++j){
//                    p.m(i, j) = v[i][j];
//                }
//            }
//            return p;
        }
        const Spectrum getSpectrum() const { return s; }
        void ToSpectrum() {
            for (int i = 0; i < nSpectralSamples; ++i) {
                Float sum = 0.f;
                for (int j = 0; j < nSpectralSamples; ++j) {
                    sum += m(j, i);
                }
                s.AssignValueAtIndex(i, sum);
            }
        }
    
        void Transpose() const {
            m.transpose();
        }
    
        int Size() const {
            return m.rows();
        }
    
        void SetValue(int i, int j, Float value) {
            m(i, j) = value;
        }

        static void Init() {
            
        }
    
    private:
        Spectrum s;
    
};

class RGBPhotoLumi : public CoefficientPhotoLumi<3> {
    using CoefficientPhotoLumi<3>::m;
    
    public:
        // RGBSpectrum Public Methods
    RGBPhotoLumi(Float v = 0.f) : CoefficientPhotoLumi<3>(v) {}
    RGBPhotoLumi(const CoefficientPhotoLumi<3> &v) : CoefficientPhotoLumi<3>(v) {}
    RGBPhotoLumi(const RGBPhotoLumi &s, PhotoLumiType type = PhotoLumiType::Reflectance) {
        *this = s;
    }
    static RGBPhotoLumi FromRGB(const Float rgb[3][3], PhotoLumiType type = PhotoLumiType::Reflectance) {
        RGBPhotoLumi s;
        s.m(0, 0) = rgb[0][0];
        s.m(0, 1) = rgb[0][1];
        s.m(0, 2) = rgb[0][2];
        s.m(1, 0) = rgb[1][0];
        s.m(1, 1) = rgb[1][1];
        s.m(1, 2) = rgb[1][2];
        s.m(2, 0) = rgb[2][0];
        s.m(2, 1) = rgb[2][1];
        s.m(2, 2) = rgb[2][2];
        DCHECK(!s.HasNaNs());
        return s;
    }
    void ToRGB(Float **rgb) const {
        rgb[0][0] = m(0, 0);
        rgb[0][1] = m(0, 1);
        rgb[0][2] = m(0, 2);
        rgb[1][0] = m(1, 0);
        rgb[1][1] = m(1, 1);
        rgb[1][2] = m(1, 2);
        rgb[2][0] = m(2, 0);
        rgb[2][1] = m(2, 1);
        rgb[2][2] = m(2, 2);
    }
    const RGBPhotoLumi &ToRGBPhotoLumi() const { return *this; }
    
    const Spectrum getSpectrum() const { return s; };
    void ToSpectrum() {
        for (int i = 0; i < 2; ++i) {
            Float sum = 0.f;
            for (int j = 0; j < 2; ++j) {
                sum += m(i, j);
            }
            s.AssignValueAtIndex(i, sum);
        }
    }
    
    void Transpose() const {
        m.transpose();
    }
    
    int Size() const {
        return m.rows();
    }
    
    void SetValue(int i, int j, Float value) {
        m(i, j) = value;
    }
    private:
    // RGBSpectrum Private Data
        Spectrum s;
    /*
     void ToXYZ(Float xyz[3][3])
     
     static RGBSpectrum FromXYZ(const Float xyz[3][3],
        PhotoLumiType type = PhotoLumiType::Reflectance)
     
     Float y() const {}
     
     static RGBPhotoLumi FromSampled(const Float *lambda, const Float *v, int n){}
     
     */
};
    
// PhotoLumi Inline Functions
template <int nSpectrumSamples>
inline CoefficientPhotoLumi<nSpectrumSamples> Pow(
    const CoefficientPhotoLumi<nSpectrumSamples> &s, Float e) {
    CoefficientPhotoLumi<nSpectrumSamples> ret;
    for (int i = 0; i < nSpectrumSamples; ++i)
        for (int j = 0; j < nSpectrumSamples; ++j)
            ret.m[i][j] = std::pow(s.m[i][j], e);
    DCHECK(!ret.HasNaNs());
    return ret;
}

inline RGBPhotoLumi Lerp(Float t, const RGBPhotoLumi &s1, const RGBPhotoLumi &s2) {
    return (1 - t) * s1 + t * s2;
}

inline SampledPhotoLumi Lerp(Float t, const SampledPhotoLumi &s1, const SampledPhotoLumi &s2) {
    return (1 - t) * s1 + t * s2;
}
//void ResampleLinearPhotoLumi(const Float *lamdaIn, const Float *vIn, int nIn,
//                            Float lambdaMin, Float lambdaMax, int nOut,
//                            Float *vOut);
}  // namespace pbrt

#endif // PBRT_CORE_PHOTOLUMI_H
