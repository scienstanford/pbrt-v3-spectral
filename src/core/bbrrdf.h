//
//  bbrrdf.h
//  PBRT-V3
//  BBRRDF stands for bispectral bidirectional reflection
//  and reradiation function.
//
//  Created by ZhengLyu on 6/29/19.
//

#if defined(_MSC_VER)
#define NOMINMAX
#progama once
#endif

#ifndef PBRT_CORE_BBRRDF_H
#define PBRT_CORE_BBRRDF_H

// core/bbrrdf.h*
#include "interaction.h"
#include "reflection.h"
#include "stats.h"

namespace pbrt {

// BBRRDF Declarations
class BBRRDF {
    public:
        // BBRRDF Public Methods
        explicit BBRRDF(PhotoLumi reRadMatrix) : reRadMatrix(std::move(reRadMatrix)) {}
        virtual ~BBRRDF() = default;
        virtual PhotoLumi f(const Vector3f &wo, const Vector3f &wi) const = 0;
        virtual PhotoLumi Sample_f(const Vector3f &wo, Vector3f *wi,
                           const Point2f &sample, Float *pdf, BxDFType type = BSDF_ALL,
                           BxDFType *sampledType = nullptr) const = 0;
        virtual PhotoLumi rho(const Vector3f &wo, int nSamples,
                                            const Point2f*samples) const = 0;
        virtual PhotoLumi rho(const int nSamples, const Point2f *samples1,
                              const Point2f *samples2) const = 0;
        virtual Float Pdf(const Vector3f &wo, const Vector3f &wi) const = 0;
        int NumComponents() const;

    
    protected:
        // BBRRDF Protected Data
        int nBBRRDFs = 0;
        static PBRT_CONSTEXPR int MaxBBRRDFs = 8;
        PhotoLumi reRadMatrix;

};

inline int BBRRDF::NumComponents() const {
    return nBBRRDFs;
}

//class SubSurfaceBBRRDF : public BBRRDF {
//    public:
//        // SubSurfaceBBRRDF Public Methods
//        SubSurfaceBBRRDF(Eigen::MatrixXd reRadMatrix, const SurfaceInteraction &po, Float eta)
//            : BBRRDF(reRadMatrix),
//              po(po),
//              eta(eta) {}
//
//        // The functions below are temporarily used for demo, which are not true
//        PhotoLumi S(const SurfaceInteraction &pi, const Vector3f &wi) const;
//        PhotoLumi Sample_S(const Scene &scene, Float u1, const Point2f &u2,
//                                   MemoryArena &arena, SurfaceInteraction *si,
//                                   Float *pdf) const;
//
//    protected:
//        // SubSurfaceBSSRDF Protected Data
//        const SurfaceInteraction &po;
//        Float eta;
//};

class SurfaceBBRRDF : public BBRRDF {
    public:
        // SurfaceBBRRDF Public Methods
        explicit SurfaceBBRRDF(PhotoLumi reRadMatrix)
            : BBRRDF(std::move(reRadMatrix)) {}
    
        void Add(SurfaceBBRRDF *b) {
            CHECK_LT(nBBRRDFs, MaxBBRRDFs);
            bbrrdfs[nBBRRDFs++] = b;
        }
        // The functions below are temporarily used for demo, which are not true
        PhotoLumi f(const Vector3f &wo, const Vector3f &wi) const override;
        virtual PhotoLumi rho(const Vector3f &wo, int nSamples,
                              const Point2f*samples) const override { return reRadMatrix; }; // Reflectance, seems meaningless here?
        virtual PhotoLumi rho(const int nSamples, const Point2f *samples1,
                              const Point2f *samples2) const override { return reRadMatrix; };
        PhotoLumi Sample_f(const Vector3f &wo, Vector3f *wi,
                           const Point2f &sample, Float *pdf, BxDFType type = BSDF_ALL,
                           BxDFType *sampledType = nullptr) const override;
        Float Pdf(const Vector3f &wo, const Vector3f &wi) const override;
    private:
        SurfaceBBRRDF *bbrrdfs[MaxBBRRDFs];
        friend class MixMaterial;


};

class ScaledSurfaceBBRRDF : public SurfaceBBRRDF {
  public:
    ScaledSurfaceBBRRDF(SurfaceBBRRDF *bbrrdf, const Spectrum &scale)
        : SurfaceBBRRDF(reRadMatrix), bbrrdf(bbrrdf), scale(scale) {}
    PhotoLumi f(const Vector3f &wo, const Vector3f &wi) const {
        return ((bbrrdf->f(wo, wi)).array().rowwise() * scale.transpose()).matrix();
    }
    PhotoLumi Sample_f(const Vector3f &wo, Vector3f *wi,
                       const Point2f &sample, Float *pdf, BxDFType type = BSDF_ALL,
                       BxDFType *sampledType = nullptr) const {
        PhotoLumi f = bbrrdf->Sample_f(wo, wi, sample, pdf, type, sampledType);
        return (f.array().rowwise() * scale.transpose()).matrix();
    }
    
  
  private:
    BBRRDF *bbrrdf;
    Spectrum scale;
};

} // pbrt namespace

#endif // PBRT_CORE_BBRRDF_H
