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
        FluorescentMaterial(const std::shared_ptr<Texture<PhotoLumi>> &reRadMatrix,
                            const std::shared_ptr<Texture<Float>> &bumpMap)
            : reRadMatrix(reRadMatrix),
              bumpMap(bumpMap){}
        void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                        TransportMode mode,
                                        bool allowMultipleLobes) const;
    
    private:
        // FluorescentMaterial Private Data
        std::shared_ptr<Texture<PhotoLumi>> reRadMatrix;
        std::shared_ptr<Texture<Float>> bumpMap;
 
};
    

FluorescentMaterial *CreateFluorescentMaterial(const TextureParams &mp);
    
}


#endif /* fluorescent_h */
