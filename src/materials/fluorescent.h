//
//  fluorescent.h
//  PBRT-V3
//
//  Created by ZhengLyu on 7/9/19.
//

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_FLUORESCENT_H
#define PBRT_MATERIALS_FLUORESCENT_H

// materials/fluorescent.h
#include "pbrt.h"
#include "reflection.h"
#include "material.h"
#include "bbrrdf.h"

namespace pbrt {

// FluorescentMaterial Declarations
class FluorescentMaterial : public Material {
    public:
        // FluorescentMaterial Public Methods
        FluorescentMaterial(const std::shared_ptr<Texture<PhotoLumi>> &fluorescence,
                            const std::shared_ptr<Texture<Float>> &bumpMap,
                            const std::shared_ptr<Texture<Spectrum>> &Kd,
                            const std::shared_ptr<Texture<Float>> &sigma)
            : fluorescence(fluorescence),
              bumpMap(bumpMap),
              Kd(Kd),
              sigma(sigma){}
        void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                        TransportMode mode,
                                        bool allowMultipleLobes) const;
    
    private:
        // FluorescentMaterial Private Data
        std::shared_ptr<Texture<PhotoLumi>> fluorescence;
        std::shared_ptr<Texture<Float>> bumpMap;
        std::shared_ptr<Texture<Spectrum>> Kd;
    std::shared_ptr<Texture<Float>> sigma;
};
    

FluorescentMaterial *CreateFluorescentMaterial(const TextureParams &mp);
    
}


#endif /* fluorescent_h */
