//
//  bbrrdf.cpp
//  bsdftest
//
//  Created by ZhengLyu on 7/5/19.
//

#include "bbrrdf.h"
#include "interpolation.h"
#include "parallel.h"
#include "scene.h"

namespace pbrt {
    Spectrum BBRRDF::f(const Vector3f &wo, const Vector3f &wi) const {
        Spectrum tempF(.5f);
        Spectrum result(0.f);
        for (int i = 0; i < tempF.nSamples; ++i) { // col
            const Float thisF = tempF[i];
            for (int j = 0; j < tempF.nSamples; ++j) { // row
                result[i] += thisF * reradMatrix(j, i);
            }
        }
        
        return result;
    }
}
