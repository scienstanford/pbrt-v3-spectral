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
#include "pbrt.h"
#include "geometry.h"
#include "microfacet.h"
#include "shape.h"
#include "spectrum.h"
#include "reflection.h"
#include "Eigen/Dense"

namespace pbrt {

// BBRRDF Declarations
class BBRRDF : public BxDF {
    public:
      // BBRRDF Public Methods
    BBRRDF(BxDF *bxdf, const Eigen::MatrixXd &reradMatrix)
            : BxDF(BxDFType(bxdf->type)), bxdf(bxdf),
                        reradMatrix(reradMatrix) {}
        Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
        std::string ToString() const;
    
    private:
        BxDF *bxdf;
        Eigen::MatrixXd reradMatrix;
};

} // pbrt namespace

#endif // PBRT_CORE_BBRRDF_H
