//
//  bbrrdf.cpp
//  bsdftest
//
//  We would for now assume Bidirectional Bispectral reflection and reradiation
//  distribution function follows Lambertian distribution, as mentioned in multiple
//  articles.
//
//  ZhengLyu & Haomiao, 2020
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

PhotoLumi SurfaceBBRRDF::f(const Vector3f &wo, const Vector3f &wi) const {
    return reRadMatrix * InvPi;
}

PhotoLumi SurfaceBBRRDF::Sample_f(const Vector3f &wo, Vector3f *wi,
                                  const Point2f &u, Float *pdf, BxDFType type,
                                  BxDFType *sampledType) const {
    // Cosine-sample the hemisphere, flipping the direction if necessary
    *wi = CosineSampleHemisphere(u);
    if (wo.z < 0) wi->z *= -1;
    *pdf = Pdf(wo, *wi);
    
    PhotoLumi f = PhotoLumi::Zero();
    for (int i = 0; i < nBBRRDFs; ++i)
        f = f + bbrrdfs[i]->f(wo, *wi);
    
    return f;
}

Float SurfaceBBRRDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    // Bought from BxDF
    return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
}


}
