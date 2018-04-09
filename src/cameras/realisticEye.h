
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

#ifndef PBRT_CAMERAS_REALISTICEYE_H
#define PBRT_CAMERAS_REALISTICEYE_H

// cameras/realistic.h*
#include "pbrt.h"
#include "camera.h"
#include "film.h"
#include <gsl/gsl_randist.h>
#include "spectrum.h" // This is necessary to declare Spectrum class in this header file.

namespace pbrt {
    
    // This is slightly different than our normal lens class. We have more variables to be able to describe the elements of the eye including biconic surfaces and unique ocular mediums.
    struct LensElementEye{
        float radiusX;
        float radiusY;
        float thickness;
        float mediumIndex; // Corresponds to the mediumElement. Describes media directly "behind" (-z direction, toward the retina) the surface.
        float semiDiameter;
        float conicConstantX;
        float conicConstantY;
    };
    
    // RealisticCamera Declarations
    class RealisticEye : public Camera {
    public:
        // RealisticCamera Public Methods
        RealisticEye(const AnimatedTransform &CameraToWorld,
                     Float shutterOpen,
                     Float shutterClose,
                     bool simpleWeighting,
                     bool noWeighting,
                     Film *film,
                     const Medium *medium,
                     std::string specfile,
                     Float pupilDiameter,
                     Float retinaDistance,
                     Float retinaRadius,
                     Float retinaSemiDiam,
                     std::vector<Spectrum> iorSpectra,
                     bool flipRad,
                     bool mmUnits,
                     bool diffractionEnabled);
        
        Float GenerateRay(const CameraSample &sample, Ray *) const;
        
    private:
        
        const bool simpleWeighting;
        const bool noWeighting;
        
        // Lens information
        std::vector<LensElementEye> lensEls;
        Float effectiveFocalLength;
        
        // Specific parameters for the human eye
        Float pupilDiameter;
        Float retinaDistance;
        Float retinaRadius;
        Float retinaSemiDiam;
        Float retinaDiag; // This will take the place of "film->diag"
        std::vector<Spectrum> iorSpectra;
        
        // Flags for conventions
        bool diffractionEnabled;
        float lensScaling;
        
        // Private methods for tracing through lens
        bool IntersectLensElAspheric(const Ray &r, Float *tHit, LensElementEye currElement, Float zShift, Vector3f *n) const;
        void applySnellsLaw(Float n1, Float n2, Float lensRadius, Vector3f &normalVec, Ray * ray ) const;
        Float lookUpIOR(int mediumIndex, const Ray &ray) const;
        void diffractHURB(Point3f intersect, Float apertureRadius, const Float wavelength, const Vector3f oldDirection, Vector3f *newDirection) const;
        
        // Handy method to explicity solve for the z(x,y) at a given point (x,y), for the biconic SAG
        Float BiconicZ(Float x, Float y, LensElementEye currElement) const;
        
        // GSL seed(?) for random number generation
        gsl_rng * r;
        
    };
    
    RealisticEye *CreateRealisticEye(const ParamSet &params,
                                     const AnimatedTransform &cam2world,
                                     Film *film, const Medium *medium);
    
}  // namespace pbrt

#endif  // PBRT_CAMERAS_REALISTICEYE_H
