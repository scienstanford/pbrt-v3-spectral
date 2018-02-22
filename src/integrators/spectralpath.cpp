
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

// integrators/path.cpp*
#include "integrators/spectralpath.h"
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"
#include "progressreporter.h"

namespace pbrt {
    
    
    STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);
    
    STAT_PERCENT("Integrator/Zero-radiance paths", zeroRadiancePaths, totalPaths);
    STAT_INT_DISTRIBUTION("Integrator/Path length", pathLength);
    
    // SpectralPathIntegrator Method Definitions
    SpectralPathIntegrator::SpectralPathIntegrator(int maxDepth,
                                                   std::shared_ptr<const Camera> camera,
                                                   std::shared_ptr<Sampler> sampler,
                                                   const Bounds2i &pixelBounds, Float rrThreshold,
                                                   const std::string &lightSampleStrategy,
                                                   int numCABands)
    : SamplerIntegrator(camera, sampler, pixelBounds),
    maxDepth(maxDepth),
    rrThreshold(rrThreshold),
    lightSampleStrategy(lightSampleStrategy),
    numCABands(numCABands){
    }
    
    void SpectralPathIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
        lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);
    }
    
    Spectrum SpectralPathIntegrator::Li(const RayDifferential &r, const Scene &scene,
                                        Sampler &sampler, MemoryArena &arena,
                                        int depth) const {
        ProfilePhase p(Prof::SamplerIntegratorLi);
        Spectrum L(0.f), beta(1.f);
        RayDifferential ray(r);
        bool specularBounce = false;
        int bounces;
        // Added after book publication: etaScale tracks the accumulated effect
        // of radiance scaling due to rays passing through refractive
        // boundaries (see the derivation on p. 527 of the third edition). We
        // track this value in order to remove it from beta when we apply
        // Russian roulette; this is worthwhile, since it lets us sometimes
        // avoid terminating refracted rays that are about to be refracted back
        // out of a medium and thus have their beta value increased.
        Float etaScale = 1;
        
        for (bounces = 0;; ++bounces) {
            // Find next path vertex and accumulate contribution
            VLOG(2) << "Path tracer bounce " << bounces << ", current L = " << L
            << ", beta = " << beta;
            
            // Intersect _ray_ with scene and store intersection in _isect_
            SurfaceInteraction isect;
            bool foundIntersection = scene.Intersect(ray, &isect);
            
            // Possibly add emitted light at intersection
            if (bounces == 0 || specularBounce) {
                // Add emitted light at path vertex or from the environment
                if (foundIntersection) {
                    L += beta * isect.Le(-ray.d);
                    VLOG(2) << "Added Le -> L = " << L;
                } else {
                    for (const auto &light : scene.infiniteLights)
                        L += beta * light->Le(ray);
                    VLOG(2) << "Added infinite area lights -> L = " << L;
                }
            }
            
            // Terminate path if ray escaped or _maxDepth_ was reached
            if (!foundIntersection || bounces >= maxDepth) break;
            
            // Compute scattering functions and skip over medium boundaries
            isect.ComputeScatteringFunctions(ray, arena, true);
            if (!isect.bsdf) {
                VLOG(2) << "Skipping intersection due to null bsdf";
                ray = isect.SpawnRay(ray.d);
                bounces--;
                continue;
            }
            
            const Distribution1D *distrib = lightDistribution->Lookup(isect.p);
            
            // Sample illumination from lights to find path contribution.
            // (But skip this for perfectly specular BSDFs.)
            if (isect.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) >
                0) {
                ++totalPaths;
                Spectrum Ld = beta * UniformSampleOneLight(isect, scene, arena,
                                                           sampler, false, distrib);
                VLOG(2) << "Sampled direct lighting Ld = " << Ld;
                if (Ld.IsBlack()) ++zeroRadiancePaths;
                CHECK_GE(Ld.y(), 0.f);
                L += Ld;
            }
            
            // Sample BSDF to get new path direction
            Vector3f wo = -ray.d, wi;
            Float pdf;
            BxDFType flags;
            Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
                                              BSDF_ALL, &flags);
            VLOG(2) << "Sampled BSDF, f = " << f << ", pdf = " << pdf;
            if (f.IsBlack() || pdf == 0.f) break;
            beta *= f * AbsDot(wi, isect.shading.n) / pdf;
            VLOG(2) << "Updated beta = " << beta;
            CHECK_GE(beta.y(), 0.f);
            DCHECK(!std::isinf(beta.y()));
            specularBounce = (flags & BSDF_SPECULAR) != 0;
            if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
                Float eta = isect.bsdf->eta;
                // Update the term that tracks radiance scaling for refraction
                // depending on whether the ray is entering or leaving the
                // medium.
                etaScale *= (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
            }
            ray = isect.SpawnRay(wi);
            
            // Account for subsurface scattering, if applicable
            if (isect.bssrdf && (flags & BSDF_TRANSMISSION)) {
                // Importance sample the BSSRDF
                SurfaceInteraction pi;
                Spectrum S = isect.bssrdf->Sample_S(
                                                    scene, sampler.Get1D(), sampler.Get2D(), arena, &pi, &pdf);
                DCHECK(!std::isinf(beta.y()));
                if (S.IsBlack() || pdf == 0) break;
                beta *= S / pdf;
                
                // Account for the direct subsurface scattering component
                L += beta * UniformSampleOneLight(pi, scene, arena, sampler, false,
                                                  lightDistribution->Lookup(pi.p));
                
                // Account for the indirect subsurface scattering component
                Spectrum f = pi.bsdf->Sample_f(pi.wo, &wi, sampler.Get2D(), &pdf,
                                               BSDF_ALL, &flags);
                if (f.IsBlack() || pdf == 0) break;
                beta *= f * AbsDot(wi, pi.shading.n) / pdf;
                DCHECK(!std::isinf(beta.y()));
                specularBounce = (flags & BSDF_SPECULAR) != 0;
                ray = pi.SpawnRay(wi);
            }
            
            // Possibly terminate the path with Russian roulette.
            // Factor out radiance scaling due to refraction in rrBeta.
            Spectrum rrBeta = beta * etaScale;
            if (rrBeta.MaxComponentValue() < rrThreshold && bounces > 3) {
                Float q = std::max((Float).05, 1 - rrBeta.MaxComponentValue());
                if (sampler.Get1D() < q) break;
                beta /= 1 - q;
                DCHECK(!std::isinf(beta.y()));
            }
        }
        ReportValue(pathLength, bounces);
        return L;
    }
    
    // Overwrite render method
    void SpectralPathIntegrator::Render(const Scene &scene) {
        
        Preprocess(scene, *sampler);
        // Render image tiles in parallel
        
        // Compute number of tiles, _nTiles_, to use for parallel rendering
        Bounds2i sampleBounds = camera->film->GetSampleBounds();
        Vector2i sampleExtent = sampleBounds.Diagonal();
        const int tileSize = 16;
        Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
                       (sampleExtent.y + tileSize - 1) / tileSize);
        ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
        {
            ParallelFor2D([&](Point2i tile) {
                // Render section of image corresponding to _tile_
                
                // Allocate _MemoryArena_ for tile
                MemoryArena arena;
                
                // Get sampler instance for tile
                int seed = tile.y * nTiles.x + tile.x;
                std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);
                
                // Compute sample bounds for tile
                int x0 = sampleBounds.pMin.x + tile.x * tileSize;
                int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
                int y0 = sampleBounds.pMin.y + tile.y * tileSize;
                int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
                Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
                LOG(INFO) << "Starting image tile " << tileBounds;
                
                // Get _FilmTile_ for tile
                std::unique_ptr<FilmTile> filmTile =
                camera->film->GetFilmTile(tileBounds);
                
                // Calculate corresponding index positions on sampled spectrum (e.g. if nSpectralSamples = 32 and nCABands = 3, we want to divide the indices into (1 to 11), (12 to 22), and (23 to 32.) This delta index defines the spacing.)
                int deltaIndex = round(nSpectralSamples/numCABands);
                float deltaWave = (sampledLambdaEnd-sampledLambdaStart)/numCABands;
                
                // Loop over pixels in tile to render them
                for (Point2i pixel : tileBounds) {
                    {
                        ProfilePhase pp(Prof::StartPixel);
                        tileSampler->StartPixel(pixel);
                    }
                    
                    // Do this check after the StartPixel() call; this keeps
                    // the usage of RNG values from (most) Samplers that use
                    // RNGs consistent, which improves reproducability /
                    // debugging.
                    if (!InsideExclusive(pixel, pixelBounds))
                        continue;
                    
                    do {
                        // Initialize _CameraSample_ for current sample
                        CameraSample cameraSample =
                        tileSampler->GetCameraSample(pixel);
                        
                        Spectrum L(0.f); // This will be the final radiance for this bundle of rays of different wavelength.
                        Float rayWeight;
                        
                        // For each sample, we loop through  all the CA bands and trace a new ray per wavelength. We then put all the returned values in a spectrum for the original sample.
                        for(int s = 0; s < numCABands; s++){
                            
                            // Generate camera ray for current sample
                            RayDifferential ray;
                            
                            // Attach a wavelength value
                            ray.wavelength = sampledLambdaStart + deltaWave * s + (deltaWave/2); // Use middle wavelength of the spectrum band
                            
                            Spectrum Ls(0.f);
                            
                            rayWeight = camera->GenerateRayDifferential(cameraSample, &ray);
                            ray.ScaleDifferentials(1 / std::sqrt((Float)tileSampler->samplesPerPixel));
                            ++nCameraRays;
                            
                            // Evaluate radiance along camera ray
                    
                            // This specific ray (with an assigned wavelength band) will go through the rest of the rendering pipeline in "Li". This includes going out through the lens (where it will be refracted according to its wavelength), reflecting off objects, and finally hitting a light source. The radiance is returned here. The radiance is returned as a full spectrum, but we only care about the value associated with the ray's assigned wavelength. This is because the direction the ray exited the lens is dependent on the wavelength.
                            if (rayWeight > 0) Ls = Li(ray, scene, *tileSampler, arena, 0);
                            
                            // Issue warning if unexpected radiance value returned
                            if (Ls.HasNaNs()) {
                                LOG(ERROR) << StringPrintf(
                                                           "Not-a-number radiance value returned "
                                                           "for pixel (%d, %d), sample %d. Setting to black.",
                                                           pixel.x, pixel.y,
                                                           (int)tileSampler->CurrentSampleNumber());
                                Ls = Spectrum(0.f);
                            } else if (Ls.y() < -1e-5) {
                                LOG(ERROR) << StringPrintf(
                                                           "Negative luminance value, %f, returned "
                                                           "for pixel (%d, %d), sample %d. Setting to black.",
                                                           Ls.y(), pixel.x, pixel.y,
                                                           (int)tileSampler->CurrentSampleNumber());
                                Ls = Spectrum(0.f);
                            } else if (std::isinf(Ls.y())) {
                                LOG(ERROR) << StringPrintf(
                                                           "Infinite luminance value returned "
                                                           "for pixel (%d, %d), sample %d. Setting to black.",
                                                           pixel.x, pixel.y,
                                                           (int)tileSampler->CurrentSampleNumber());
                                Ls = Spectrum(0.f);
                            }
                            VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
                            ray << " -> L = " << Ls;
                            
                            // Isolate the returned radiance for this wavelength
                            
                            float Ls_lambda;
                            Ls.GetValueAtWavelength(ray.wavelength, &Ls_lambda);
                            
                            // Assign the result to all wavelengths sampled around the target wavelength. For example, if nWaveBands = 3, then we would split the spectrum into three equal parts and assign the result from the first wavelength to the first third, the second wavelength to the second third, etc.
                            int bottomIndex = deltaIndex*s;
                            int topIndex = std::min(deltaIndex*(s+1),nSpectralSamples-1);
                            for(int waveIndex = bottomIndex; waveIndex < topIndex; waveIndex++){
                                L.AssignValueAtIndex(waveIndex, Ls_lambda);
                            }
                            

                        }
                        
                        // Add camera ray's contribution to image
                        filmTile->AddSample(cameraSample.pFilm, L, rayWeight);
                        
                        // Free _MemoryArena_ memory from computing image sample
                        // value
                        arena.Reset();
                    } while (tileSampler->StartNextSample());
                }
                LOG(INFO) << "Finished image tile " << tileBounds;
                
                // Merge image tile into _Film_
                camera->film->MergeFilmTile(std::move(filmTile));
                reporter.Update();
            }, nTiles);
            reporter.Done();
        }
        LOG(INFO) << "Rendering finished";
        
        // Save final image after rendering
        camera->film->WriteImage();
        
        
    }
    
    SpectralPathIntegrator *CreateSpectralPathIntegrator(const ParamSet &params,
                                                         std::shared_ptr<Sampler> sampler,
                                                         std::shared_ptr<const Camera> camera) {
        int maxDepth = params.FindOneInt("maxdepth", 5);
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
        Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
        std::string lightStrategy =
        params.FindOneString("lightsamplestrategy", "spatial");
        
        // Get number of wavelength dependent waves to generate per single ray
        int numCABands = params.FindOneInt("numCABands", 4);
        if(numCABands != 1){
            Warning("Using spectral rendering. For every pixel sample we will trace %dx more rays. Rendering will be %d times slower.",numCABands,numCABands);
        }
        return new SpectralPathIntegrator(maxDepth, camera, sampler, pixelBounds,
                                          rrThreshold, lightStrategy,numCABands);
    }
    
}  // namespace pbrt
