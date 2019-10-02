//
//  fluorescent.cpp
//  bsdftest
//
//  Created by ZhengLyu on 7/9/19.
//

// materials/fluorescent.cpp*
#include "materials/fluorescent.h"
#include "textures/constant.h"
#include "spectrum.h"
#include "texture.h"
#include "paramset.h"
#include "interaction.h"

namespace pbrt {

// FluprescentMaterial Method Definitions
void FluorescentMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
    bool allowMultipleLobes) const {
    // Perform bump mapping with _bumpMap_, if present
    if (bumpMap) Bump(bumpMap, si);
    PhotoLumi reRadMtx = reRadMatrix->Evaluate(*si);
    si->bbrrdf = ARENA_ALLOC(arena, SurfaceBBRRDF)(reRadMtx);
}

FluorescentMaterial *CreateFluorescentmaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<PhotoLumi>> reRadMatrix = mp.GetPhotoLumiTexture("photolumi", PhotoLumi(1.f));
    std::shared_ptr<Texture<Float>> bumpMap =
        mp.GetFloatTextureOrNull("bumpmap");
    return new FluorescentMaterial(reRadMatrix, bumpMap);
}
}
