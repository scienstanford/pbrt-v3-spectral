
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


// materials/retroreflective.cpp*
#include "materials/retroreflective.h"
#include "textures/constant.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"
#include "spectrum.h"
#include "scene.h"
#include "sampling.h"
#include "sampler.h"
#include "stats.h"
#include <stdarg.h>
#include "interpolation.h"
#include "microfacet.h"
#include "shape.h"
#include "bbrrdf.h"


namespace pbrt {

// Create a look up table for Effective Retroreflective Area
    inline Float ERA(const Float incidentAngle)
    {const int nViewangleSamples = 25;
        const Float Viewangle[nViewangleSamples] = {
            1.0000, 0.9970, 0.9899, 0.9797, 0.9645, 0.9535, 0.9406, 0.9243, 0.9102, 0.8986,
            0.8867, 0.8757, 0.8627, 0.8491, 0.8365, 0.8234, 0.8065, 0.6166, 0.5472, 0.4357,
            0.3948, 0.2711, 0.1546, 0.0095, 0.0000};
        const Float Era[nViewangleSamples] = {
            0.6700, 0.6620, 0.6180, 0.6180, 0.5730, 0.5370, 0.4950, 0.4450, 0.4050, 0.3640,
            0.3260, 0.2890, 0.2500, 0.2100, 0.1730, 0.1340, 0.0930, 0.0930, 0.0930, 0.0930,
            0.0930, 0.0930, 0.0830, 0.0030, 0};
        Float ERA_area;
        if(incidentAngle >= Viewangle[0]) ERA_area = 0.67;
        for (int i = 0; i < nViewangleSamples; ++i)
        {if (incidentAngle >= Viewangle[i+1] && incidentAngle <= Viewangle[i])
        ERA_area = Era[i];}
        if(incidentAngle <= Viewangle[24]) ERA_area = 0; return ERA_area;}

    //////////////////////////////////////////////////////////////////////////////////////////////
// ClassDiffuse term
//    class RetroDiffuse : public BSSRDF{
//    public:
//        //
//    };
// Class Mirror term
  
class RetroSpecularReflection : public BxDF {
    public:
        // RetroSpecularReflection Public Methods
        RetroSpecularReflection(const Spectrum &R, Fresnel *fresnel)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)),
        R(R),
        fresnel(fresnel) {}
        Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
            return Spectrum(0.f);
        }
        Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample,
                          Float *pdf, BxDFType *sampledType) const;
        Float Pdf(const Vector3f &wo, const Vector3f &wi) const { return 0; }
        std::string ToString() const;
        
    private:
        // ReroSpecularReflection Private Data
        const Spectrum R;
        const Fresnel *fresnel;
};

// Class Retroreflective term
class RetroReflection : public BxDF {
    public:
        // MicrofacetReflection Public Methods
        RetroReflection(const Spectrum &R, Fresnel *fresnel)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
        R(R),
        fresnel(fresnel) {}
        Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
        Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                          Float *pdf, BxDFType *sampledType) const;
        Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
        std::string ToString() const;
        
    private:
        // MicrofacetReflection Private Data
        const Spectrum R;
        const MicrofacetDistribution *distribution;
        const Fresnel *fresnel;
};

//////////////////////////////////////////////////////////////////////////////////////////////
// Diffuse term
    
// Mirror term

Spectrum RetroSpecularReflection::Sample_f(const Vector3f &wo, Vector3f *wi,
                                          const Point2f &sample, Float *pdf,
                                         BxDFType *sampledType) const {
       // Compute perfect specular reflection direction
       *wi = Vector3f(wo.x, wo.y, wo.z);
        //    *wi = Vector3f(-wo.x, -wo.y, wo.z);
        *pdf = 1;
        return fresnel->Evaluate(CosTheta(wo)) * R / AbsCosTheta(*wi);// change fresnel *wi to wo
    }

    std::string RetroSpecularReflection::ToString() const {
       return std::string("[ RetroSpecularReflection R: ") + R.ToString() +
        std::string(" fresnel: ") + fresnel->ToString() + std::string(" ]");
}

// Retro term
Spectrum RetroReflection::f(const Vector3f &wo, const Vector3f &wi) const {
        Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
//        Vector3f wh = wo + wi;
        // Handle degenerate cases for microfacet reflection
        if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0.);
//        if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
//        wh = Normalize(wh);
        Spectrum Fi = fresnel->Evaluate(CosTheta(wi));
        Spectrum Fo = fresnel->Evaluate(CosTheta(wo));
        Spectrum unit = Spectrum::Ones();
        //  (1 / (1 + Lambda(wo) + Lambda(wi)));
    Vector3f n = Vector3f(0, 0, 1);
    Float sin2ThetaI = std::max(Float(0), Float(1 - cosThetaI * cosThetaI));
    Float sin2ThetaT = 0.5 * 0.5 * sin2ThetaI;
    
    // Handle total internal reflection for transmission
    if (sin2ThetaT >= 1) return Spectrum::Zero();
    Float cosThetaT = std::sqrt(1 - sin2ThetaT);
    Vector3f wt = 0.5 * -wo + (0.5 * cosThetaI - cosThetaT) * (n);
    Float incidentAngle = AbsCosTheta(wt);
//    Vector3f wt = Refract(wi, n, 1.8f);
        return R * (unit-Fi)*(unit-Fo)*ERA(incidentAngle) / cosThetaI;// function E(wt) needs to be added.
    }
    
    std::string RetroReflection::ToString() const {
        return std::string("[ RetroReflection R: ") + R.ToString() +
        std::string(" distribution: ") + distribution->ToString() +
        std::string(" fresnel: ") + fresnel->ToString() + std::string(" ]");
}
    
    
Spectrum RetroReflection::Sample_f(const Vector3f &wo, Vector3f *wi,
                                                 const Point2f &u, Float *pdf,
                                                 BxDFType *sampledType) const {
        // Sample microfacet orientation $\wh$ and reflected direction $\wi$
        if (wo.z == 0) return Spectrum::Zero();
       // Vector3f wh = distribution->Sample_wh(wo, u);
        *wi = wo;
        //if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);
        // Compute PDF of _wi_ for microfacet reflection
        //*pdf = 1;
        Normal3f n;
        *pdf = Pdf(wo, *wi);
        return f(wo, *wi);
    }
    
    Float RetroReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {
        //if (!SameHemisphere(wo, wi)) return 0;
        //Vector3f wh = Normalize(wo + wi);
        return 1; // needs to add more complex calculation
}
////////////////////////////////////////////////////////////////////////////////////////////////
//
// RetroreflectiveMaterial Method Definitions
//
void retroreflectiveMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                               MemoryArena &arena,
                                               TransportMode mode,
                                               bool allowMultipleLobes) const {
    // Perform bump mapping with _bumpMap_, if present
    if (bumpMap) Bump(bumpMap, si);
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, 1.8f);
    Fresnel *frMf =ARENA_ALLOC(arena, FresnelDielectric)(1., 1.8);
    Spectrum R = Kr->Evaluate(*si).Clamp(); // Add specular reflection term from mirror
//    Spectrum ks = Ks->Evaluate(*si).Clamp();
    Spectrum kd = Kd->Evaluate(*si).Clamp();
    if (!R.IsBlack())
        si->bsdf->Add(ARENA_ALLOC(arena, RetroSpecularReflection)(
            kd, ARENA_ALLOC(arena, FresnelNoOp))); // FresnelNoOp // Change R to kd
    si->bsdf->Add(ARENA_ALLOC(arena, RetroReflection)(kd, frMf)); // Change ks to kd
//    // subsurface diffuse
    Spectrum mfree = Spectrum::Ones();
    
//
    Spectrum sig_a, sig_s;
    Float sig_a_rgb[3] = {.0011f, .0024f, .014f},
    sig_s_rgb[3] = {2.55f, 3.21f, 3.77f};
    sig_a = Spectrum::FromRGB(sig_a_rgb),
    sig_s = Spectrum::FromRGB(sig_s_rgb);
    SubsurfaceFromDiffuse(table, kd, mfree, &sig_a, &sig_s);
    si->bssrdf = ARENA_ALLOC(arena, TabulatedBSSRDF)(*si, this, mode, 1.8f,
                                                     sig_a, sig_s, table);
//    if (!kd.IsBlack()) {
//        BxDF *diff = ARENA_ALLOC(arena, LambertianReflection)(kd);
//        si->bsdf->Add(diff);
//}
    if (fluorescence != nullptr) {
      si->bbrrdf = ARENA_ALLOC(arena, SurfaceBBRRDF)(
          fluorescence->Evaluate(*si) * concentration->Evaluate(*si));
      // This need to be replaced by other scattering function at some time.
      si->bbrrdf->Add(ARENA_ALLOC(arena, SurfaceBBRRDF)(
          fluorescence->Evaluate(*si) * concentration->Evaluate(*si)));
    }

}

    
retroreflectiveMaterial *CreateRetroreflectiveMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> Kd =
    mp.GetSpectrumTexture("Kd", Spectrum(0.25f));
    std::shared_ptr<Texture<Spectrum>> Kr =
    mp.GetSpectrumTexture("Kr", Spectrum(0.9f));
    std::shared_ptr<Texture<Spectrum>> Ks =
    mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    std::shared_ptr<Texture<Float>> sigma = mp.GetFloatTexture("sigma", 0.f);
    std::shared_ptr<Texture<Float>> roughness =
        mp.GetFloatTexture("roughness", .01f);
    std::shared_ptr<Texture<Float>> uRoughness =
        mp.GetFloatTextureOrNull("uroughness");
    std::shared_ptr<Texture<Float>> vRoughness =
        mp.GetFloatTextureOrNull("vroughness");
    std::shared_ptr<Texture<PhotoLumi>> fluorescence =
        mp.GetPhotoLumiTextureOrNull("fluorescence");
    std::shared_ptr<Texture<Float>> concentration =
        mp.GetFloatTexture("concentration", 1.f);
    std::shared_ptr<Texture<Float>> bumpMap =
        mp.GetFloatTextureOrNull("bumpmap");
    bool remapRoughness = mp.FindBool("remaproughness", true);
    return new retroreflectiveMaterial(Kd, Kr, Ks, sigma, roughness, uRoughness, vRoughness, fluorescence, concentration, bumpMap,remapRoughness);
}

}  // namespace pbrt
