
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

// cameras/omni.cpp*
#include "cameras/omni.h"
#include "paramset.h"
#include "sampler.h"
#include "sampling.h"
#include "floatfile.h"
#include "imageio.h"
#include "reflection.h"
#include "stats.h"
#include "lowdiscrepancy.h"
#include "light.h"
#include "samplers/random.h"
#include "ext/json.hpp"
#include <array>
#include <fstream>

using json = nlohmann::json;

namespace pbrt {

STAT_PERCENT("Camera/Rays vignetted by lens system", vignettedRays, totalRays);

// OmniCamera Method Definitions
OmniCamera::OmniCamera(const AnimatedTransform &CameraToWorld, Float shutterOpen,
    Float shutterClose, Float apertureDiameter, Float filmdistance,
    Float focusDistance, bool simpleWeighting, bool noWeighting,
    bool caFlag, std::vector<Float> &lensData, std::vector<Float> &microlensData,
    Vector2i microlensDims, std::vector<Float> & microlensOffsets, float microlensSensorOffset,
    Film *film,
    const Medium *medium)
    : Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
      simpleWeighting(simpleWeighting), noWeighting(noWeighting), caFlag(caFlag) {
    static RNG rng;

    auto floatsToElementInterfaces = [](std::vector<Float>& lensData, std::vector<LensElementInterface>& elementInterfaces, float apertureDiameter) {
        for (int i = 0; i < (int)lensData.size(); i += 4) {
            if (lensData[i] == 0) {
                if (apertureDiameter > lensData[i + 3]) {
                    Warning(
                        "Specified aperture diameter %f is greater than maximum "
                        "possible %f.  Clamping it.",
                        apertureDiameter, lensData[i + 3]);
                }
                else {
                    lensData[i + 3] = apertureDiameter;
                }
            }
            Float cRadius = lensData[i] * (Float).001;
            Float aRadius = lensData[i + 3] * (Float).001 / Float(2.);
            Float thickness = lensData[i + 1] * (Float).001;
            Float ior = lensData[i + 2];
            LensElementInterface interface(
                cRadius, aRadius, thickness, ior);
            interface.transform = Translate(Vector3f(0.001, 0.001, 0.0)*rng.UniformFloat()-Vector3f(0.0005,0.0005,0.0));
            elementInterfaces.push_back(interface);
        }
    };

    floatsToElementInterfaces(lensData, elementInterfaces, apertureDiameter);
    
    if (microlensData.size() > 0) {
        floatsToElementInterfaces(microlensData, microlens.elementInterfaces, apertureDiameter);
        microlens.offsets.resize(microlensOffsets.size() / 2);
        for (int i = 0; i < microlens.offsets.size(); ++i) {
            microlens.offsets[i].x = microlensOffsets[2 * i + 0];
            microlens.offsets[i].y = microlensOffsets[2 * i + 1];
        }
        microlens.dimensions = microlensDims;
        microlens.offsetFromSensor = microlensSensorOffset;
    }


    // Compute lens--film distance for given focus distance
    // TL: If a film distance is given, hardset the focus distance. If not, use the focus distance given.
    if(filmdistance == 0){
        Float fb = FocusBinarySearch(focusDistance);
        LOG(INFO) << StringPrintf("Binary search focus: %f -> %f\n", fb,
                                  FocusDistance(fb));
        elementInterfaces.back().thickness = FocusThickLens(focusDistance);
        LOG(INFO) << StringPrintf("Thick lens focus: %f -> %f\n",
                                  elementInterfaces.back().thickness,
                                    FocusDistance(elementInterfaces.back().thickness));
    } else {
        // Use given film distance
        LOG(INFO) << StringPrintf("Focus distance hard set: %f -> %f\n",filmdistance,
                                  FocusDistance(filmdistance));
        elementInterfaces.back().thickness = filmdistance;
    }
    
    // Print out film distance into terminal
    std::cout << "Distance from film to back of lens: " << elementInterfaces.back().thickness << " m" << std::endl;
    std::cout << "Focus distance in scene: " << FocusDistance(elementInterfaces.back().thickness) << " m" << std::endl;
          
    // Compute exit pupil bounds at sampled points on the film
    int nSamples = 64;
    exitPupilBounds.resize(nSamples);
    ParallelFor([&](int64_t i) {
        Float r0 = (Float)i / nSamples * film->diagonal / 2;
        Float r1 = (Float)(i + 1) / nSamples * film->diagonal / 2;
        exitPupilBounds[i] = BoundExitPupil(r0, r1);
    }, nSamples);

    if (simpleWeighting)
        Warning("\"simpleweighting\" option with OmniCamera no longer "
                "necessarily matches regular camera images. Further, pixel "
                "values will vary a bit depending on the aperture size. See "
                "this discussion for details: "
                "https://github.com/mmp/pbrt-v3/issues/162#issuecomment-348625837");
}

bool OmniCamera::TraceLensesFromFilm(const Ray &rCamera, 
    const std::vector<LensElementInterface>& interfaces, Ray *rOut,
    const Transform CameraToLens = Scale(1, 1, -1)) const {
    Float elementZ = 0;
    // Transform _rCamera_ from camera to lens system space
    Ray rLens = CameraToLens(rCamera);
    
    // Added by Trisha to keep wavelength information
    if(rOut){
        rLens.wavelength = rOut->wavelength;
    } else {
        rLens.wavelength = 550;
    }
    
    for (int i = interfaces.size() - 1; i >= 0; --i) {
        const LensElementInterface &element = interfaces[i];
        // Update ray from film accounting for interaction with _element_
        elementZ -= element.thickness;

        // Compute intersection of ray with lens element
        Float t;
        Normal3f n;
        bool isStop = (element.curvatureRadius.x == 0);
        if (isStop) {
            // The refracted ray computed in the previous lens element
            // interface may be pointed towards film plane(+z) in some
            // extreme situations; in such cases, 't' becomes negative.
            if (rLens.d.z >= 0.0) return false;
            t = (elementZ - rLens.o.z) / rLens.d.z;
        } else {
            Float radius = element.curvatureRadius.x;
            Float zCenter = elementZ + element.curvatureRadius.x;
            rLens = element.transform(rLens);
            if (!IntersectSphericalElement(radius, zCenter, rLens, &t, &n))
                return false;
        }
        CHECK_GE(t, 0);

        auto invTransform = Inverse(element.transform);
        rLens = invTransform(rLens);
        n = invTransform(n);
        // Test intersection point against element aperture
        Point3f pHit = rLens(t);
        Float r2 = pHit.x * pHit.x + pHit.y * pHit.y;
        if (r2 > element.apertureRadius.x * element.apertureRadius.x) return false;
        rLens.o = pHit;

        // Update ray path for element interface interaction
        if (!isStop) {
            Vector3f w;
            Float etaI = element.eta;
            Float etaT = (i > 0 && elementInterfaces[i - 1].eta != 0)
                             ? elementInterfaces[i - 1].eta
                             : 1;
            // Added by Trisha and Zhenyi (5/18)
            if(caFlag && (rLens.wavelength >= 400) && (rLens.wavelength <= 700))
            {                    
                if (etaI != 1)
                    etaI = (rLens.wavelength - 550) * -.04/(300)  +  etaI;
                if (etaT != 1)
                    etaT = (rLens.wavelength - 550) * -.04/(300)  +  etaT;
            }
            
            if (!Refract(Normalize(-rLens.d), n, etaI / etaT, &w)) return false;
            rLens.d = w;
        }
    }
    // Transform _rLens_ from lens system space back to camera space
    if (rOut != nullptr) {
        static const Transform LensToCamera = Scale(1, 1, -1);
        *rOut = LensToCamera(rLens);
    }
    return true;
}

bool OmniCamera::IntersectSphericalElement(Float radius, Float zCenter,
                                                const Ray &ray, Float *t,
                                                Normal3f *n) {
    // Compute _t0_ and _t1_ for ray--element intersection
    Point3f o = ray.o - Vector3f(0, 0, zCenter);
    Float A = ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z;
    Float B = 2 * (ray.d.x * o.x + ray.d.y * o.y + ray.d.z * o.z);
    Float C = o.x * o.x + o.y * o.y + o.z * o.z - radius * radius;
    Float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1)) return false;

    // Select intersection $t$ based on ray direction and element curvature
    bool useCloserT = (ray.d.z > 0) ^ (radius < 0);
    *t = useCloserT ? std::min(t0, t1) : std::max(t0, t1);
    if (*t < 0) return false;

    // Compute surface normal of element at ray intersection point
    *n = Normal3f(Vector3f(o + *t * ray.d));
    *n = Faceforward(Normalize(*n), -ray.d);
    return true;
}

bool OmniCamera::TraceLensesFromScene(const Ray &rCamera,
                                           Ray *rOut) const {
    Float elementZ = -LensFrontZ();
    // Transform _rCamera_ from camera to lens system space
    static const Transform CameraToLens = Scale(1, 1, -1);
    Ray rLens = CameraToLens(rCamera);
    for (size_t i = 0; i < elementInterfaces.size(); ++i) {
        const LensElementInterface &element = elementInterfaces[i];
        // Compute intersection of ray with lens element
        Float t;
        Normal3f n;
        bool isStop = (element.curvatureRadius.x == 0);
        if (isStop)
            t = (elementZ - rLens.o.z) / rLens.d.z;
        else {
            Float radius = element.curvatureRadius.x;
            Float zCenter = elementZ + element.curvatureRadius.x;
            if (!IntersectSphericalElement(radius, zCenter, rLens, &t, &n))
                return false;
        }
        CHECK_GE(t, 0);

        // Test intersection point against element aperture
        Point3f pHit = rLens(t);
        Float r2 = pHit.x * pHit.x + pHit.y * pHit.y;
        if (r2 > element.apertureRadius.x * element.apertureRadius.x) return false;
        rLens.o = pHit;

        // Update ray path for from-scene element interface interaction
        if (!isStop) {
            Vector3f wt;
            Float etaI = (i == 0 || elementInterfaces[i - 1].eta == 0)
                             ? 1
                             : elementInterfaces[i - 1].eta;
            Float etaT =
                (elementInterfaces[i].eta != 0) ? elementInterfaces[i].eta : 1;
            if (!Refract(Normalize(-rLens.d), n, etaI / etaT, &wt))
                return false;
            rLens.d = wt;
        }
        elementZ += element.thickness;
    }
    // Transform _rLens_ from lens system space back to camera space
    if (rOut != nullptr) {
        static const Transform LensToCamera = Scale(1, 1, -1);
        *rOut = LensToCamera(rLens);
    }
    return true;
}

void OmniCamera::DrawLensSystem() const {
    Float sumz = -LensFrontZ();
    Float z = sumz;
    for (size_t i = 0; i < elementInterfaces.size(); ++i) {
        const LensElementInterface &element = elementInterfaces[i];
        Float r = element.curvatureRadius.x;
        if (r == 0) {
            // stop
            printf("{Thick, Line[{{%f, %f}, {%f, %f}}], ", z,
                   element.apertureRadius.x, z, 2 * element.apertureRadius.x);
            printf("Line[{{%f, %f}, {%f, %f}}]}, ", z, -element.apertureRadius.x,
                   z, -2 * element.apertureRadius.x);
        } else {
            Float theta = std::abs(std::asin(element.apertureRadius.x / r));
            if (r > 0) {
                // convex as seen from front of lens
                Float t0 = Pi - theta;
                Float t1 = Pi + theta;
                printf("Circle[{%f, 0}, %f, {%f, %f}], ", z + r, r, t0, t1);
            } else {
                // concave as seen from front of lens
                Float t0 = -theta;
                Float t1 = theta;
                printf("Circle[{%f, 0}, %f, {%f, %f}], ", z + r, -r, t0, t1);
            }
            if (element.eta != 0 && element.eta != 1) {
                // connect top/bottom to next element
                CHECK_LT(i + 1, elementInterfaces.size());
                Float nextApertureRadius =
                    elementInterfaces[i + 1].apertureRadius.x;
                Float h = std::max(element.apertureRadius.x, nextApertureRadius);
                Float hlow =
                    std::min(element.apertureRadius.x, nextApertureRadius);

                Float zp0, zp1;
                if (r > 0) {
                    zp0 = z + element.curvatureRadius.x -
                          element.apertureRadius.x / std::tan(theta);
                } else {
                    zp0 = z + element.curvatureRadius.x +
                          element.apertureRadius.x / std::tan(theta);
                }

                Float nextCurvatureRadius =
                    elementInterfaces[i + 1].curvatureRadius.x;
                Float nextTheta = std::abs(
                    std::asin(nextApertureRadius / nextCurvatureRadius));
                if (nextCurvatureRadius > 0) {
                    zp1 = z + element.thickness + nextCurvatureRadius -
                          nextApertureRadius / std::tan(nextTheta);
                } else {
                    zp1 = z + element.thickness + nextCurvatureRadius +
                          nextApertureRadius / std::tan(nextTheta);
                }

                // Connect tops
                printf("Line[{{%f, %f}, {%f, %f}}], ", zp0, h, zp1, h);
                printf("Line[{{%f, %f}, {%f, %f}}], ", zp0, -h, zp1, -h);

                // vertical lines when needed to close up the element profile
                if (element.apertureRadius.x < nextApertureRadius) {
                    printf("Line[{{%f, %f}, {%f, %f}}], ", zp0, h, zp0, hlow);
                    printf("Line[{{%f, %f}, {%f, %f}}], ", zp0, -h, zp0, -hlow);
                } else if (element.apertureRadius.x > nextApertureRadius) {
                    printf("Line[{{%f, %f}, {%f, %f}}], ", zp1, h, zp1, hlow);
                    printf("Line[{{%f, %f}, {%f, %f}}], ", zp1, -h, zp1, -hlow);
                }
            }
        }
        z += element.thickness;
    }

    // 24mm height for 35mm film
    printf("Line[{{0, -.012}, {0, .012}}], ");
    // optical axis
    printf("Line[{{0, 0}, {%f, 0}}] ", 1.2f * sumz);
}

void OmniCamera::DrawRayPathFromFilm(const Ray &r, bool arrow,
                                          bool toOpticalIntercept) const {
    Float elementZ = 0;
    // Transform _ray_ from camera to lens system space
    static const Transform CameraToLens = Scale(1, 1, -1);
    Ray ray = CameraToLens(r);
    printf("{ ");
    if (!TraceLensesFromFilm(r, elementInterfaces, nullptr)) printf("Dashed, ");
    for (int i = elementInterfaces.size() - 1; i >= 0; --i) {
        const LensElementInterface &element = elementInterfaces[i];
        elementZ -= element.thickness;
        bool isStop = (element.curvatureRadius.x == 0);
        // Compute intersection of ray with lens element
        Float t;
        Normal3f n;
        if (isStop)
            t = -(ray.o.z - elementZ) / ray.d.z;
        else {
            Float radius = element.curvatureRadius.x;
            Float zCenter = elementZ + element.curvatureRadius.x;
            if (!IntersectSphericalElement(radius, zCenter, ray, &t, &n))
                goto done;
        }
        CHECK_GE(t, 0);

        printf("Line[{{%f, %f}, {%f, %f}}],", ray.o.z, ray.o.x, ray(t).z,
               ray(t).x);

        // Test intersection point against element aperture
        Point3f pHit = ray(t);
        Float r2 = pHit.x * pHit.x + pHit.y * pHit.y;
        Float apertureRadius2 = element.apertureRadius.x * element.apertureRadius.x;
        if (r2 > apertureRadius2) goto done;
        ray.o = pHit;

        // Update ray path for element interface interaction
        if (!isStop) {
            Vector3f wt;
            Float etaI = element.eta;
            Float etaT = (i > 0 && elementInterfaces[i - 1].eta != 0)
                             ? elementInterfaces[i - 1].eta
                             : 1;
            if (!Refract(Normalize(-ray.d), n, etaI / etaT, &wt)) goto done;
            ray.d = wt;
        }
    }

    ray.d = Normalize(ray.d);
    {
        Float ta = std::abs(elementZ / 4);
        if (toOpticalIntercept) {
            ta = -ray.o.x / ray.d.x;
            printf("Point[{%f, %f}], ", ray(ta).z, ray(ta).x);
        }
        printf("%s[{{%f, %f}, {%f, %f}}]", arrow ? "Arrow" : "Line", ray.o.z,
               ray.o.x, ray(ta).z, ray(ta).x);

        // overdraw the optical axis if needed...
        if (toOpticalIntercept)
            printf(", Line[{{%f, 0}, {%f, 0}}]", ray.o.z, ray(ta).z * 1.05f);
    }

done:
    printf("}");
}

void OmniCamera::DrawRayPathFromScene(const Ray &r, bool arrow,
                                           bool toOpticalIntercept) const {
    Float elementZ = LensFrontZ() * -1;

    // Transform _ray_ from camera to lens system space
    static const Transform CameraToLens = Scale(1, 1, -1);
    Ray ray = CameraToLens(r);
    for (size_t i = 0; i < elementInterfaces.size(); ++i) {
        const LensElementInterface &element = elementInterfaces[i];
        bool isStop = (element.curvatureRadius.x == 0);
        // Compute intersection of ray with lens element
        Float t;
        Normal3f n;
        if (isStop)
            t = -(ray.o.z - elementZ) / ray.d.z;
        else {
            Float radius = element.curvatureRadius.x;
            Float zCenter = elementZ + element.curvatureRadius.x;
            if (!IntersectSphericalElement(radius, zCenter, ray, &t, &n))
                return;
        }
        CHECK_GE(t, 0.f);

        printf("Line[{{%f, %f}, {%f, %f}}],", ray.o.z, ray.o.x, ray(t).z,
               ray(t).x);

        // Test intersection point against element aperture
        Point3f pHit = ray(t);
        Float r2 = pHit.x * pHit.x + pHit.y * pHit.y;
        Float apertureRadius2 = element.apertureRadius.x * element.apertureRadius.x;
        if (r2 > apertureRadius2) return;
        ray.o = pHit;

        // Update ray path for from-scene element interface interaction
        if (!isStop) {
            Vector3f wt;
            Float etaI = (i == 0 || elementInterfaces[i - 1].eta == 0.f)
                             ? 1.f
                             : elementInterfaces[i - 1].eta;
            Float etaT = (elementInterfaces[i].eta != 0.f)
                             ? elementInterfaces[i].eta
                             : 1.f;
            if (!Refract(Normalize(-ray.d), n, etaI / etaT, &wt)) return;
            ray.d = wt;
        }
        elementZ += element.thickness;
    }

    // go to the film plane by default
    {
        Float ta = -ray.o.z / ray.d.z;
        if (toOpticalIntercept) {
            ta = -ray.o.x / ray.d.x;
            printf("Point[{%f, %f}], ", ray(ta).z, ray(ta).x);
        }
        printf("%s[{{%f, %f}, {%f, %f}}]", arrow ? "Arrow" : "Line", ray.o.z,
               ray.o.x, ray(ta).z, ray(ta).x);
    }
}

void OmniCamera::ComputeCardinalPoints(const Ray &rIn, const Ray &rOut,
                                            Float *pz, Float *fz) {
    Float tf = -rOut.o.x / rOut.d.x;
    *fz = -rOut(tf).z;
    Float tp = (rIn.o.x - rOut.o.x) / rOut.d.x;
    *pz = -rOut(tp).z;
}

void OmniCamera::ComputeThickLensApproximation(Float pz[2],
                                                    Float fz[2]) const {
    // Find height $x$ from optical axis for parallel rays
    Float x = .001 * film->diagonal;

    // Compute cardinal points for film side of lens system
    Ray rScene(Point3f(x, 0, LensFrontZ() + 1), Vector3f(0, 0, -1));
    Ray rFilm;
    CHECK(TraceLensesFromScene(rScene, &rFilm))
        << "Unable to trace ray from scene to film for thick lens "
           "approximation. Is aperture stop extremely small?";
    ComputeCardinalPoints(rScene, rFilm, &pz[0], &fz[0]);

    // Compute cardinal points for scene side of lens system
    rFilm = Ray(Point3f(x, 0, LensRearZ() - 1), Vector3f(0, 0, 1));
    CHECK(TraceLensesFromFilm(rFilm, elementInterfaces, &rScene))
        << "Unable to trace ray from film to scene for thick lens "
           "approximation. Is aperture stop extremely small?";
    ComputeCardinalPoints(rFilm, rScene, &pz[1], &fz[1]);
}

Float OmniCamera::FocusThickLens(Float focusDistance) {
    Float pz[2], fz[2];
    ComputeThickLensApproximation(pz, fz);
    LOG(INFO) << StringPrintf("Cardinal points: p' = %f f' = %f, p = %f f = %f.\n",
                              pz[0], fz[0], pz[1], fz[1]);
    LOG(INFO) << StringPrintf("Effective focal length %f\n", fz[0] - pz[0]);
    // Compute translation of lens, _delta_, to focus at _focusDistance_
    Float f = fz[0] - pz[0];
    Float z = -focusDistance;
    Float c = (pz[1] - z - pz[0]) * (pz[1] - z - 4 * f - pz[0]);
    CHECK_GT(c, 0) << "Coefficient must be positive. It looks focusDistance: " << focusDistance << " is too short for a given lenses configuration";
    Float delta =
        0.5f * (pz[1] - z + pz[0] - std::sqrt(c));
    return elementInterfaces.back().thickness + delta;
}

Float OmniCamera::FocusBinarySearch(Float focusDistance) {
    Float filmDistanceLower, filmDistanceUpper;
    // Find _filmDistanceLower_, _filmDistanceUpper_ that bound focus distance
    filmDistanceLower = filmDistanceUpper = FocusThickLens(focusDistance);
    while (FocusDistance(filmDistanceLower) > focusDistance)
        filmDistanceLower *= 1.005f;
    while (FocusDistance(filmDistanceUpper) < focusDistance)
        filmDistanceUpper /= 1.005f;

    // Do binary search on film distances to focus
    for (int i = 0; i < 20; ++i) {
        Float fmid = 0.5f * (filmDistanceLower + filmDistanceUpper);
        Float midFocus = FocusDistance(fmid);
        if (midFocus < focusDistance)
            filmDistanceLower = fmid;
        else
            filmDistanceUpper = fmid;
    }
    return 0.5f * (filmDistanceLower + filmDistanceUpper);
}

Float OmniCamera::FocusDistance(Float filmDistance) {
    // Find offset ray from film center through lens
    Bounds2f bounds = BoundExitPupil(0, .001 * film->diagonal);

    const std::array<Float, 3> scaleFactors = {0.1f, 0.01f, 0.001f};
    Float lu = 0.0f;

    Ray ray;

    // Try some different and decreasing scaling factor to find focus ray
    // more quickly when `aperturediameter` is too small.
    // (e.g. 2 [mm] for `aperturediameter` with wide.22mm.dat),
    bool foundFocusRay = false;
    for (Float scale : scaleFactors) {
        lu = scale * bounds.pMax[0];
        if (TraceLensesFromFilm(Ray(Point3f(0, 0, LensRearZ() - filmDistance),
                                    Vector3f(lu, 0, filmDistance)), elementInterfaces,
                                &ray)) {
            foundFocusRay = true;
            break;
        }
    }

    if (!foundFocusRay) {
        Error(
            "Focus ray at lens pos(%f,0) didn't make it through the lenses "
            "with film distance %f?!??\n",
            lu, filmDistance);
        return Infinity;
    }

    // Compute distance _zFocus_ where ray intersects the principal axis
    Float tFocus = -ray.o.x / ray.d.x;
    Float zFocus = ray(tFocus).z;
    if (zFocus < 0) zFocus = Infinity;
    return zFocus;
}


Bounds2f OmniCamera::BoundExitPupil(Float pFilmX0, Float pFilmX1) const {
    Bounds2f pupilBounds;
    // Sample a collection of points on the rear lens to find exit pupil
    const int nSamples = 1024 * 1024;
    int nExitingRays = 0;

    // Compute bounding box of projection of rear element on sampling plane
    Float rearRadius = RearElementRadius();
    Bounds2f projRearBounds(Point2f(-1.5f * rearRadius, -1.5f * rearRadius),
                            Point2f(1.5f * rearRadius, 1.5f * rearRadius));
    for (int i = 0; i < nSamples; ++i) {
        // Find location of sample points on $x$ segment and rear lens element
        Point3f pFilm(Lerp((i + 0.5f) / nSamples, pFilmX0, pFilmX1), 0, 0);
        Float u[2] = {RadicalInverse(0, i), RadicalInverse(1, i)};
        Point3f pRear(Lerp(u[0], projRearBounds.pMin.x, projRearBounds.pMax.x),
                      Lerp(u[1], projRearBounds.pMin.y, projRearBounds.pMax.y),
                      LensRearZ());

        // Expand pupil bounds if ray makes it through the lens system
        if (Inside(Point2f(pRear.x, pRear.y), pupilBounds) ||
            TraceLensesFromFilm(Ray(pFilm, pRear - pFilm), elementInterfaces, nullptr)) {
            pupilBounds = Union(pupilBounds, Point2f(pRear.x, pRear.y));
            ++nExitingRays;
        }
    }

    // Return entire element bounds if no rays made it through the lens system
    if (nExitingRays == 0) {
        LOG(INFO) << StringPrintf("Unable to find exit pupil in x = [%f,%f] on film.",
                                  pFilmX0, pFilmX1);
        return projRearBounds;
    }

    // Expand bounds to account for sample spacing
    pupilBounds = Expand(pupilBounds, 2 * projRearBounds.Diagonal().Length() /
                                          std::sqrt(nSamples));
    return pupilBounds;
}

void OmniCamera::RenderExitPupil(Float sx, Float sy,
                                      const char *filename) const {
    Point3f pFilm(sx, sy, 0);

    const int nSamples = 2048;
    Float *image = new Float[3 * nSamples * nSamples];
    Float *imagep = image;

    for (int y = 0; y < nSamples; ++y) {
        Float fy = (Float)y / (Float)(nSamples - 1);
        Float ly = Lerp(fy, -RearElementRadius(), RearElementRadius());
        for (int x = 0; x < nSamples; ++x) {
            Float fx = (Float)x / (Float)(nSamples - 1);
            Float lx = Lerp(fx, -RearElementRadius(), RearElementRadius());

            Point3f pRear(lx, ly, LensRearZ());

            if (lx * lx + ly * ly > RearElementRadius() * RearElementRadius()) {
                *imagep++ = 1;
                *imagep++ = 1;
                *imagep++ = 1;
            } else if (TraceLensesFromFilm(Ray(pFilm, pRear - pFilm),
                elementInterfaces, nullptr)) {
                *imagep++ = 0.5f;
                *imagep++ = 0.5f;
                *imagep++ = 0.5f;
            } else {
                *imagep++ = 0.f;
                *imagep++ = 0.f;
                *imagep++ = 0.f;
            }
        }
    }

    WriteImage(filename, image,
               Bounds2i(Point2i(0, 0), Point2i(nSamples, nSamples)),
               Point2i(nSamples, nSamples));
    delete[] image;
}


Point3f OmniCamera::SampleExitPupil(const Point2f &pFilm,
                                         const Point2f &lensSample,
                                         Float *sampleBoundsArea) const {
    
    // Find exit pupil bound for sample distance from film center
    Float rFilm = std::sqrt(pFilm.x * pFilm.x + pFilm.y * pFilm.y);
    int rIndex = rFilm / (film->diagonal / 2) * exitPupilBounds.size();
    rIndex = std::min((int)exitPupilBounds.size() - 1, rIndex);
    Bounds2f pupilBounds = exitPupilBounds[rIndex];
    if (sampleBoundsArea) *sampleBoundsArea = pupilBounds.Area();

    // Generate sample point inside exit pupil bound
    Point2f pLens = pupilBounds.Lerp(lensSample);

    // Return sample point rotated by angle of _pFilm_ with $+x$ axis
    Float sinTheta = (rFilm != 0) ? pFilm.y / rFilm : 0;
    Float cosTheta = (rFilm != 0) ? pFilm.x / rFilm : 1;
    return Point3f(cosTheta * pLens.x - sinTheta * pLens.y,
                   sinTheta * pLens.x + cosTheta * pLens.y, LensRearZ());
}

static Vector2f mapMul(Vector2f v0, Vector2f v1) {
    return Vector2f(v0.x*v1.x, v0.y*v1.y);
}
static Vector2f mapDiv(Vector2f v0, Vector2f v1) {
    return Vector2f(v0.x/v1.x, v0.y/v1.y);
}
static Point2f mapDiv(Point2f v0, Vector2f v1) {
    return Point2f(v0.x / v1.x, v0.y / v1.y);
}


Point3f OmniCamera::SampleMicrolensPupil(const Point2f &pFilm, const Point2f &lensSample,
    Float *sampleBoundsArea) const {
    Bounds2f extent = film->GetPhysicalExtent();
    Vector2f normFilm = mapDiv(pFilm - extent.pMin, extent.pMax - extent.pMin);
    Vector2i d = microlens.dimensions;
    Vector2f df = Vector2f(d.x, d.y);
    Vector2f lensSpace = mapMul(normFilm, df);
    Point2f lensIndex(floor(lensSpace.x), floor(lensSpace.y));
    Point2f sampledLensSpacePt = mapDiv(lensIndex+lensSample, df);
    // sample on microlens
    Point2f result2 = extent.Lerp(sampledLensSpacePt);

    return Point3f(result2.x, result2.y, microlens.offsetFromSensor);
}

void OmniCamera::TestExitPupilBounds() const {
    Float filmDiagonal = film->diagonal;

    static RNG rng;

    Float u = rng.UniformFloat();
    Point3f pFilm(u * filmDiagonal / 2, 0, 0);

    Float r = pFilm.x / (filmDiagonal / 2);
    int pupilIndex =
        std::min((int)exitPupilBounds.size() - 1,
                 (int)std::floor(r * (exitPupilBounds.size() - 1)));
    Bounds2f pupilBounds = exitPupilBounds[pupilIndex];
    if (pupilIndex + 1 < (int)exitPupilBounds.size())
        pupilBounds = Union(pupilBounds, exitPupilBounds[pupilIndex + 1]);

    // Now, randomly pick points on the aperture and see if any are outside
    // of pupil bounds...
    for (int i = 0; i < 1000; ++i) {
        Point2f u2{rng.UniformFloat(), rng.UniformFloat()};
        Point2f pd = ConcentricSampleDisk(u2);
        pd *= RearElementRadius();

        Ray testRay(pFilm, Point3f(pd.x, pd.y, 0.f) - pFilm);
        Ray testOut;
        if (!TraceLensesFromFilm(testRay, elementInterfaces, &testOut)) continue;

        if (!Inside(pd, pupilBounds)) {
            fprintf(stderr,
                    "Aha! (%f,%f) went through, but outside bounds (%f,%f) - "
                    "(%f,%f)\n",
                    pd.x, pd.y, pupilBounds.pMin[0], pupilBounds.pMin[1],
                    pupilBounds.pMax[0], pupilBounds.pMax[1]);
            RenderExitPupil(
                (Float)pupilIndex / exitPupilBounds.size() * filmDiagonal / 2.f,
                0.f, "low.exr");
            RenderExitPupil((Float)(pupilIndex + 1) / exitPupilBounds.size() *
                                filmDiagonal / 2.f,
                            0.f, "high.exr");
            RenderExitPupil(pFilm.x, 0.f, "mid.exr");
            exit(0);
        }
    }
    fprintf(stderr, ".");
}


Float OmniCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    ++totalRays;
    // Find point on film, _pFilm_, corresponding to _sample.pFilm_
    Point2f s(sample.pFilm.x / film->fullResolution.x,
        sample.pFilm.y / film->fullResolution.y);
    Point2f pFilm2 = film->GetPhysicalExtent().Lerp(s);
    Point3f pFilm(-pFilm2.x, pFilm2.y, 0);

    Float exitPupilBoundsArea;
    Point3f pRear;
    // Trace ray from _pFilm_ through lens system
    if (HasMicrolens()) {
        pRear = SampleMicrolensPupil(Point2f(pFilm.x, pFilm.y), sample.pLens, &exitPupilBoundsArea);
    } else {
        pRear = SampleExitPupil(Point2f(pFilm.x, pFilm.y), sample.pLens, &exitPupilBoundsArea);
    }
    Ray rFilm = Ray(pFilm, pRear - pFilm, Infinity,
        Lerp(sample.time, shutterOpen, shutterClose));

    if (HasMicrolens()) {
        if (!TraceLensesFromFilm(rFilm, microlens.elementInterfaces, &rFilm)) {
            ++vignettedRays;
            return 0;
        }
    }

    if (!TraceLensesFromFilm(rFilm, elementInterfaces, ray)) {
        ++vignettedRays;
        return 0;
    }

    // Finish initialization of _OmniCamera_ ray
    *ray = CameraToWorld(*ray);
    ray->d = Normalize(ray->d);
    ray->medium = medium;


    // Return weighting for _OmniCamera_ ray
    if (HasMicrolens()) {
        // TODO: Proper weighting
        Float cosTheta = Normalize(rFilm.d).z;
        Float cos4Theta = (cosTheta * cosTheta) * (cosTheta * cosTheta);
        if (simpleWeighting)
            return cos4Theta * exitPupilBoundsArea / exitPupilBounds[0].Area();
        else
            return (shutterClose - shutterOpen) *
            (cos4Theta * exitPupilBoundsArea) / (LensRearZ() * LensRearZ());
    } else {
        Float cosTheta = Normalize(rFilm.d).z;
        Float cos4Theta = (cosTheta * cosTheta) * (cosTheta * cosTheta);
        if (simpleWeighting)
            return cos4Theta * exitPupilBoundsArea / exitPupilBounds[0].Area();
        else
            return (shutterClose - shutterOpen) *
            (cos4Theta * exitPupilBoundsArea) / (LensRearZ() * LensRearZ());
    }
}


bool OmniCamera::HasMicrolens() const {
    return microlens.elementInterfaces.size() > 0;
}

OmniCamera *CreateOmniCamera(const ParamSet &params,
    const AnimatedTransform &cam2world,
    Film *film, const Medium *medium) {
    Float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    Float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
            shutterclose, shutteropen);
        std::swap(shutterclose, shutteropen);
    }

    // Omni camera-specific parameters
    std::string lensFile = params.FindOneFilename("lensfile", "");
    Float apertureDiameter = params.FindOneFloat("aperturediameter", 1.0);
    Float focusDistance = params.FindOneFloat("focusdistance", 10.0);
    bool simpleWeighting = params.FindOneBool("simpleweighting", true);
    bool noWeighting = params.FindOneBool("noweighting", false);

    Float microlensSensorOffset = params.FindOneFloat("microlenssensoroffset", 0.0);
    int microlensSimulationRadius = params.FindOneInt("microlenssimulationradius", 0);

    if (microlensSimulationRadius != 0) {
        Warning("We only currently support simulating one microlens per pixel, switching microlensSimulationRadius to 0");
    }

    if (lensFile == "") {
        Error("No lens description file supplied!");
        return nullptr;
    }
    // Load element data from lens description file
    std::vector<Float> lensData;

    // Load element data from lens description file
    std::vector<Float> microlensData;

    std::vector<Float> microlensOffsets;

    Vector2i microlensDims;

    auto endsWith = [](const std::string& str, const std::string& suffix) {
        return str.size() >= suffix.size() && 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
    };

    if (endsWith(lensFile, ".dat")) {
        if (!ReadFloatFile(lensFile.c_str(), &lensData)) {
            Error("Error reading lens specification file \"%s\".",
                lensFile.c_str());
            return nullptr;
        }
        if (lensData.size() % 4 != 0) {
            // Trisha: If the size has an extra value, it's possible this lens type was meant for pbrt-v2-spectral and has an extra focal length value at the top. In this case, let's automatically convert it by removing this extra value.
            if (lensData.size() % 4 == 1) {
                Warning("Extra value in lens specification file, this lens file may be for pbrt-v2-spectral. Removing extra value to make it compatible with pbrt-v3-spectral...");
                lensData.erase(lensData.begin());
            }
            else {
                Error(
                    "Excess values in lens specification file \"%s\"; "
                    "must be multiple-of-four values, read %d.",
                    lensFile.c_str(), (int)lensData.size());
                return nullptr;
            }
        }
    } else {
        if (!endsWith(lensFile, ".json")) {
            Error("Invalid format for lens specification file \"%s\".",
                lensFile.c_str());
            return nullptr;
        }
        // read a JSON file
        std::ifstream i(lensFile);
        json j;
        if (i && (i>>j)) {
            // assert(j.is_object())
            // j["name"]
            // j["description"]
            if (j["name"].is_string()) {
                LOG(INFO) << StringPrintf("Loading lens %s.\n",
                    j["name"].get<std::string>().c_str());
            }
            if (j["description"].is_string()) {
                LOG(INFO) << StringPrintf("%s\n",
                    j["description"].get<std::string>().c_str());
            }
            auto jsurfaces = j["surfaces"];
            std::vector<OmniCamera::LensElementInterface> surfaces;

            auto toVec2 = [](json val) {
                if (val.is_number()) {
                    return Vector2f{ (Float)val, (Float)val };
                } else if (val.is_array() && val.size() == 2) {
                    return Vector2f{ val[0], val[1] };
                } 
                return Vector2f(); // Default value
            };
            auto toVec2i = [](json val) {
                if (val.is_number()) {
                    return Vector2i{ (int)val, (int)val };
                }
                else if (val.is_array() && val.size() == 2) {
                    return Vector2i{ val[0], val[1] };
                }
                return Vector2i(); // Default value
            };

            auto toTransform = [lensFile](json t) {
                // Stored in columns in json, but pbrt stores matrices 
                // in row-major order.
                if (t.is_null()) { // Perfectly fine to have no transform
                    return Transform();
                }
                if (!(t.size() == 4 && t[0].size() == 3 && 
                    t[1].size() == 3 && t[2].size() == 3 && 
                    t[3].size() == 3)) {
                    Error("Invalid transform in lens specification file \"%s\", must be an array of 4 arrays of 3 floats each (column-major transform).", lensFile.c_str());
                }
                Float m[4][4];
                for (int c = 0; c < 4; ++c) {
                    for (int r = 0; r < 4; ++r) {
                        m[r][c] = (Float)t[c][r];
                    }
                }
                m[0][3] = m[1][3] = m[2][3] = (Float)0.0;
                m[3][3] = (Float)1.0;
                return Transform(m);
            };

            auto toLensElementInterface = [toVec2, toTransform](json surf) {
                OmniCamera::LensElementInterface result;
                result.apertureRadius   = toVec2(surf["semi_aperture"]);
                result.conicConstant    = toVec2(surf["conic_constant"]);
                result.curvatureRadius  = toVec2(surf["radius"]);
                result.eta              = Float(surf["ior"]);
                result.thickness        = Float(surf["thickness"]);
                result.transform        = toTransform(surf["transform"]);
                return result;
            };


            if (jsurfaces.is_array() && jsurfaces.size() > 0) {
                for (auto jsurf : jsurfaces) {
                    surfaces.push_back(toLensElementInterface(jsurf));
                }
            } else {
                Error("Error, lens specification file without a valid surface array \"%s\".",
                    lensFile.c_str());
                return nullptr;
            }
            auto surfacesToFloats = [lensFile](const std::vector<OmniCamera::LensElementInterface>& in, std::vector<Float>& out) {
                for (auto& s : in) {
                    if (s.curvatureRadius.x != s.curvatureRadius.y) {
                        Error("Error, lens specification file  \"%s\" contains aspheric element, NYI.",
                            lensFile.c_str());
                        return nullptr;
                    }
                    if (s.apertureRadius.x != s.apertureRadius.y) {
                        Error("Error, lens specification file  \"%s\" contains asymmetric aperture radius element, NYI.",
                            lensFile.c_str());
                        return nullptr;
                    }
                    if (s.conicConstant.x != s.conicConstant.y || s.conicConstant.x != (Float)0.0) {
                        Error("Error, lens specification file  \"%s\" contains conic element, NYI.",
                            lensFile.c_str());
                        return nullptr;
                    }
                    if (!s.transform.IsIdentity()) {
                        Error("Error, lens specification file  \"%s\" contains non-identity transform, NYI.",
                            lensFile.c_str());
                        return nullptr;
                    }
                    out.push_back(s.curvatureRadius.x);
                    out.push_back(s.thickness);
                    out.push_back(s.eta);
                    out.push_back(s.apertureRadius.x);
                }
            };
            surfacesToFloats(surfaces, lensData);
            auto microlens = j["microlens"];
            if (!microlens.is_null()) {
                microlensDims = toVec2i(microlens["dimensions"]);

                if (microlensDims.x <= 0 || microlensDims.y <= 0) {
                    Error("Error, microlens specification without valid dimensions in \"%s\".",
                        lensFile.c_str());
                    return nullptr;
                }

                auto mljOffsets = microlens["offsets"];
                if (mljOffsets.is_array() && mljOffsets.size() > 0) {
                    if (mljOffsets.size() != microlensDims.x*microlensDims.y) {
                        Error("Error, microlens dimensions (%d x %d ) mismatch number of offsets (%d) in \"%s\".",
                            microlensDims.x, microlensDims.y, (int)mljOffsets.size(), lensFile.c_str());
                        return nullptr;
                    }
                    for (auto offset : mljOffsets) {
                        microlensOffsets.push_back(Float(offset));
                    }
                }

                auto mljSurfaces = microlens["surfaces"];
                std::vector<OmniCamera::LensElementInterface> mlSurfaces;

                if (mljSurfaces.is_array() && mljSurfaces.size() > 0) {
                    for (auto jsurf : mljSurfaces) {
                        mlSurfaces.push_back(toLensElementInterface(jsurf));
                    }
                }
                else {
                    Error("Error, microlens specification without a valid surface array in \"%s\".",
                        lensFile.c_str());
                    return nullptr;
                }

                surfacesToFloats(mlSurfaces, microlensData);


            } else {
            }

        } else {
            Error("Error reading lens specification file \"%s\".",
                lensFile.c_str());
            return nullptr;
        }
    }
    
    // Add functionality to hard set the film distance
    Float filmDistance = params.FindOneFloat("filmdistance", 0);
    
    // Chromatic aberration flag
    bool caFlag = params.FindOneBool("chromaticAberrationEnabled", false);
    
    return new OmniCamera(cam2world, shutteropen, shutterclose,
                               apertureDiameter, filmDistance, focusDistance, simpleWeighting, noWeighting, caFlag,
                               lensData, microlensData, microlensDims, microlensOffsets, microlensSensorOffset, film, medium);
}

}  // namespace pbrt
