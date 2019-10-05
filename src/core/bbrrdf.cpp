//
//  bbrrdf.cpp
//  bsdftest
//
//  Created by ZhengLyu on 7/5/19.
//

#include "reflection.h"
#include "spectrum.h"
#include "sampler.h"
#include "sampling.h"
#include "interpolation.h"
#include "scene.h"
#include "interaction.h"
#include "stats.h"
#include <stdarg.h>
#include "photolumi.h"
#include "bbrrdf.h"

namespace pbrt {

//PhotoLumi SubSurfaceBBRRDF::S(const SurfaceInteraction &pi, const Vector3f &wi) const {
////    // Adapted from Sw
////    Float c = 1 - 2 * FresnelMoment1(1 / eta);
////    Spectrum f = (1 - FrDielectric(CosTheta(w), 1, eta)) / (c * Pi);
////    return reRadMatrix * f;
//    Spectrum tmpSp = Spectrum(1.f);
//    return reRadMatrix * tmpSp;
//}

PhotoLumi SurfaceBBRRDF::f(const Vector3f &wo, const Vector3f &wi) const {
    Spectrum tmpSp = Spectrum(1.f);
    return reRadMatrix * tmpSp;
}

PhotoLumi SurfaceBBRRDF::Sample_f(const Vector3f &wo, Vector3f *wi,
                                  const Point2f &u, Float *pdf, BxDFType type,
                                  BxDFType *sampledType) const {
    // Cosine-sample the hemisphere, flipping the direction if necessary
    *wi = CosineSampleHemisphere(u);
    if (wo.z < 0) wi->z *= -1;
    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}

Float SurfaceBBRRDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
}
}
