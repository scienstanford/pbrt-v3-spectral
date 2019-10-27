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
#include "photolumi.h"

namespace pbrt {

// FluprescentMaterial Method Definitions
void FluorescentMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
    bool allowMultipleLobes) const {
    
    // Perform bump mapping with _bumpMap_, if present
    if (bumpMap) Bump(bumpMap, si);
    
    // From matte material (To be removed)
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
    Spectrum r = Kd->Evaluate(*si).Clamp();
    Float sig = Clamp(sigma->Evaluate(*si), 0, 90);
    if (!r.IsBlack()) {
        if (sig == 0)
            si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(r));
        else
            si->bsdf->Add(ARENA_ALLOC(arena, OrenNayar)(r, sig));
    }
    //
    PhotoLumi reRadMtx = reRadMatrix->Evaluate(*si);
    int size = PhotoLumi::nSamples;
    for (int i = 0; i < size; ++i)
        reRadMtx(i, (size - 1) * 9 / 10) = 0.5;
    si->bbrrdf = ARENA_ALLOC(arena, SurfaceBBRRDF)(reRadMtx);
}

FluorescentMaterial *CreateFluorescentMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<PhotoLumi>> reRadMatrix = mp.GetPhotoLumiTexture(
        "photolumi", PhotoLumi::Identity());
    std::shared_ptr<Texture<Float>> bumpMap =
        mp.GetFloatTextureOrNull("bumpmap");
    std::shared_ptr<Texture<Spectrum>> Kd =
        mp.GetSpectrumTexture("Kd", Spectrum(0.5f));
    std::shared_ptr<Texture<Float>> sigma = mp.GetFloatTexture("sigma", 0.f);
    return new FluorescentMaterial(reRadMatrix, bumpMap, Kd, sigma);
}
}
