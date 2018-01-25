
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


// cameras/environment.cpp*
#include "cameras/environment.h"
#include "paramset.h"
#include "sampler.h"
#include "stats.h"

namespace pbrt {
    
    // EnvironmentCamera Method Definitions
    Float EnvironmentCamera::GenerateRay(const CameraSample &sample,
                                         Ray *ray) const {
        ProfilePhase prof(Prof::GenerateCameraRay);
        // Compute environment camera ray direction
        Float theta = Pi * sample.pFilm.y / film->fullResolution.y; // Elevation
        Float phi = 2 * Pi * sample.pFilm.x / film->fullResolution.x; // Azimuth
        Vector3f dir(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
                     std::cos(theta)); // convert to xyz, assuming x/y plane is spherical reference plane.
        // Flip Y and Z-axis in order to match PBRT conventions (camera is looking down the z-axis with y-axis being "up" and x-axis being "right".)
        dir = Vector3f(dir.x,dir.z,dir.y);
        *ray = Ray(Point3f(0, 0, 0), dir, Infinity,
                   Lerp(sample.time, shutterOpen, shutterClose));
        
        // ----------------------------------------
        // Added by Trisha to create omnistereo panoramas. Much of this is taken from Blender Cycles.
        // For more info on the ODS format, take a look at:
        //Peleg, Shmuel, Moshe Ben-Ezra, and Yael Pritch. "Omnistereo: Panoramic stereo imaging." IEEE Transactions on Pattern Analysis and Machine Intelligence 23.3 (2001): 279-290.
        
        Float interocular_offset = ipd;
        
        /* Interocular offset of zero means either non stereo, or stereo without
         * spherical stereo. */
        if(interocular_offset != 0.0f){
            
            if(poleMergeTo > 0.0f) {
                Float altitude = fabsf(asinf(ray->d.z));
                if(altitude > poleMergeTo) {
                    interocular_offset = 0.0f;
                }
                else if(altitude > poleMergeFrom) {
                    Float fac = (altitude - poleMergeFrom) / (poleMergeTo - poleMergeFrom);
                    Float fade = cosf(fac * PiOver2);
                    interocular_offset *= fade;
                }
            }   
            
            Vector3f up = Vector3f(0.0f, 1.0f, 0.0f);
            Vector3f side;
            if(ray->d == up){
                side = Vector3f(0.f,0.f,0.f);
            }
            else{
                side = Normalize(Cross(ray->d, up));
            }
            Vector3f stereo_offset = side * interocular_offset;
            
            ray->o += stereo_offset;
            
            /* Convergence distance is Infinity in the case of parallel convergence mode,
             * no need to modify direction in this case either. */
            if(convergenceDistance != Infinity)
            {
                Vector3f screen_offset = convergenceDistance * ray->d;
                ray->d = Normalize(screen_offset - stereo_offset);
            }
            
        }
        // ----------------------------------------
        
        ray->medium = medium;
        *ray = CameraToWorld(*ray);
        return 1;
    }
    
    EnvironmentCamera *CreateEnvironmentCamera(const ParamSet &params,
                                               const AnimatedTransform &cam2world,
                                               Film *film, const Medium *medium) {
        // Extract common camera parameters from _ParamSet_
        Float shutteropen = params.FindOneFloat("shutteropen", 0.f);
        Float shutterclose = params.FindOneFloat("shutterclose", 1.f);
        if (shutterclose < shutteropen) {
            Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                    shutterclose, shutteropen);
            std::swap(shutterclose, shutteropen);
        }
        Float lensradius = params.FindOneFloat("lensradius", 0.f);
        Float focaldistance = params.FindOneFloat("focaldistance", 1e30f);
        Float frame = params.FindOneFloat(
                                          "frameaspectratio",
                                          Float(film->fullResolution.x) / Float(film->fullResolution.y));
        Bounds2f screen;
        if (frame > 1.f) {
            screen.pMin.x = -frame;
            screen.pMax.x = frame;
            screen.pMin.y = -1.f;
            screen.pMax.y = 1.f;
        } else {
            screen.pMin.x = -1.f;
            screen.pMax.x = 1.f;
            screen.pMin.y = -1.f / frame;
            screen.pMax.y = 1.f / frame;
        }
        int swi;
        const Float *sw = params.FindFloat("screenwindow", &swi);
        if (sw) {
            if (swi == 4) {
                screen.pMin.x = sw[0];
                screen.pMax.x = sw[1];
                screen.pMin.y = sw[2];
                screen.pMax.y = sw[3];
            } else
                Error("\"screenwindow\" should have four values");
        }
        
        // Parameters added by Trisha to support omnistereo panoramas. Much of this is taken directly from Blender Cycles.
        Float ipd = params.FindOneFloat("ipd", 0.f);
        Float poleMergeTo = params.FindOneFloat("poleMergeAngleTo", 90.f);
        Float poleMergeFrom = params.FindOneFloat("poleMergeAngleFrom", 90.f);
        Float convergenceDistance = params.FindOneFloat("convergencedistance", Infinity);
        
        (void)lensradius;     // don't need this
        (void)focaldistance;  // don't need this
        
        return new EnvironmentCamera(cam2world, shutteropen, shutterclose, film,
                                     medium,ipd,poleMergeTo,poleMergeFrom,convergenceDistance);
    }
    
}  // namespace pbrt
