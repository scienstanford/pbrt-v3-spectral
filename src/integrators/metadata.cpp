
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

// Metadata integrator added by Trisha on 11/2017. This is a modification of the direct lighting integrator. Instead of returning a radiance value however, we return a value corresponding to the metadata that we want (e.g. depth, material index, etc.)


// integrators/metadata.cpp*
#include "integrators/metadata.h"
#include "interaction.h"
#include "paramset.h"
#include "camera.h"
#include "film.h"
#include "stats.h"

namespace pbrt {

// MetadataIntegrator Method Definitions
void MetadataIntegrator::Preprocess(const Scene &scene,
                                          Sampler &sampler) {
    // Do we need anything in here?
}

Spectrum MetadataIntegrator::Li(const RayDifferential &ray,
                                      const Scene &scene, Sampler &sampler,
                                      MemoryArena &arena, int depth) const {
    
    ProfilePhase p(Prof::SamplerIntegratorLi);
    Spectrum L(0.f);
    
    // Find closest ray intersection or return background radiance
    SurfaceInteraction isect;
    if (!scene.Intersect(ray, &isect)) {
        // Did not hit anything. Return 0.
        return L;
    }
    
    // Depending on the strategy, return a different value
    if(strategy == MetadataStrategy::depth){
        Vector3f toIntersect = isect.p - ray.o;
        L = Spectrum(toIntersect.Length());
        
    }else if(strategy == MetadataStrategy::material){
        L = Spectrum(isect.materialId);
        
    }else if(strategy == MetadataStrategy::mesh){
        L = Spectrum(isect.primitiveId);
        
    }else if(strategy == MetadataStrategy::coordinates){
        // Return world coordinates of intersection.
        L[0] = isect.p.x;
        L[1] = isect.p.y;
        L[2] = isect.p.z;
        
    }else{
        // We should never get here, but just in case:
        Error("Could not recognize metadata strategy.");
    }
    
    return L;
}

MetadataIntegrator *CreateMetadataIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera) {

    MetadataStrategy strategy;
    std::string st = params.FindOneString("strategy", "depth");
    if (st == "depth")
        strategy = MetadataStrategy::depth;
    else if (st == "material")
        strategy = MetadataStrategy::material;
    else if (st == "mesh")
        strategy = MetadataStrategy::mesh;
    else if (st == "coordinates")
        strategy = MetadataStrategy::coordinates;
    else {
        Warning(
            "Strategy \"%s\" for metadata unknown. "
            "Using \"depth\".",
            st.c_str());
        strategy = MetadataStrategy::depth;
    }
    
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }
    
    return new MetadataIntegrator(strategy, camera, sampler,
                                        pixelBounds);
}

}  // namespace pbrt
