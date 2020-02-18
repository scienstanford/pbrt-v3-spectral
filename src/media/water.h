
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

#ifndef PBRT_MEDIA_WATER_H
#define PBRT_MEDIA_WATER_H

// media/homogeneous.h*
#include "medium.h"

namespace pbrt {


class KopelevichPhaseFunction : public PhaseFunction {
    public:
    // PhaseFunction Interface
    KopelevichPhaseFunction(float cSmall, float cLarge);
    virtual ~KopelevichPhaseFunction();
    virtual Spectrum p(const Vector3f &wo, const Vector3f &wi) const;
    virtual Float Sample_p(const Vector3f &wo, Vector3f *wi,
                           const Point2f &u, Float wavelength) const;
    virtual std::string ToString() const {
        return StringPrintf("[ KopelevichPhaseFunction cSmall=%.3f cLarge=%.3f]",cSmall,cLarge);
    };
    
    Spectrum getScatter();

    static const int nAngularSamples = 64;
    
    private:
        const Float cSmall;
        const Float cLarge;
        Float phaseAngle[nAngularSamples];        
        std::vector<Spectrum> CDF;
    
        Spectrum sigma_s;
    
};

// WaterMedium Declarations
class WaterMedium : public Medium {
  public:
    // WaterMedium Public Methods
    WaterMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, const PhaseFunction *phase)
        : sigma_a(sigma_a),
          sigma_s(sigma_s),
          sigma_t(sigma_s + sigma_a), ph(phase){
    };
    ~WaterMedium(){
        delete ph;
    }
    Spectrum Tr(const Ray &ray, Sampler &sampler) const;
    Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                    MediumInteraction *mi) const;

  private:
    // WaterMedium Private Data
    const Spectrum sigma_a, sigma_s, sigma_t;
    const PhaseFunction *ph;
};

WaterMedium* createWaterMedium(Float cPlankton, Float aCDOM440, Float aNAP400, Float cSmall, Float cLarge);


}  // namespace pbrt

#endif  // PBRT_MEDIA_HOMOGENEOUS_H
