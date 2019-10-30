
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

#ifndef PBRT_MATERIALS_Retroreflective_H
#define PBRT_MATERIALS_Retroreflective_H

// materials/Retroreflective.h*
#include "pbrt.h"
#include "material.h"
#include "spectrum.h"
#include "bssrdf.h"

namespace pbrt {

// RetroreflectiveMaterial Declarations
class retroreflectiveMaterial : public Material {
  public:
    retroreflectiveMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,
                            const std::shared_ptr<Texture<Spectrum>> &Kr,
                            const std::shared_ptr<Texture<Spectrum>> &Ks,
                            const std::shared_ptr<Texture<Float>> &sigma,
                            const std::shared_ptr<Texture<Float>> &roughness,
                            const std::shared_ptr<Texture<Float>> &uRoughness,
                            const std::shared_ptr<Texture<Float>> &vRoughness,
                            const std::shared_ptr<Texture<PhotoLumi>> &fluorescence,
                            const std::shared_ptr<Texture<Float>> &bumpMap,
                            bool remapRoughness)
   
    : Kd(Kd),
        Kr(Kr),
        Ks(Ks),
        sigma(sigma),
        roughness(roughness),
        uRoughness(uRoughness),
        vRoughness(vRoughness),
        fluorescence(fluorescence),
        bumpMap(bumpMap),
        remapRoughness(remapRoughness),
        table(100, 64) {
        ComputeBeamDiffusionBSSRDF(0.0f, 1.8f, &table);
    }
void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // RetroreflectiveMaterial Private Data
    std::shared_ptr<Texture<Spectrum>> Kd, Kr, Ks;
    std::shared_ptr<Texture<Float>> sigma,roughness, uRoughness, vRoughness;
    std::shared_ptr<Texture<PhotoLumi>> fluorescence;
    std::shared_ptr<Texture<Float>> bumpMap;
    bool remapRoughness;
    BSSRDFTable table;
};

retroreflectiveMaterial *CreateRetroreflectiveMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_Retroreflective_H
