
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
    bool caFlag, const std::vector<OmniCamera::LensElementInterface> &lensInterfaceData,
    const std::vector<OmniCamera::LensElementInterface> &microlensData,
    Vector2i microlensDims, const std::vector<Vector2f> & microlensOffsets,
    float microlensSensorOffset, Film *film, const Medium *medium)
    : Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
      simpleWeighting(simpleWeighting), noWeighting(noWeighting), caFlag(caFlag) {
    
    elementInterfaces = lensInterfaceData;
    if (microlensData.size() > 0) {
        microlens.elementInterfaces = microlensData;
        microlens.offsets = microlensOffsets;
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

    // Print out a mathematica command that generates useful figures; could be moved into lenstool and heavily parameterized
    const bool generateMathematicaDrawing = false;
    if (generateMathematicaDrawing) {
        printf("Graphics[{");
        DrawLensSystem();
        printf(",");
        for (int i = 0; i < 15; ++i) {
            Point3f pFilm(0, 0, 0);
            Point2f pLens(i / 14.0f, 0.5f);
            Float exitPupilBoundsArea;
            Point3f pRear;
            // Trace ray from _pFilm_ through lens system
            if (HasMicrolens()) {
                pRear = SampleMicrolensPupil(Point2f(pFilm.x, pFilm.y), pLens, &exitPupilBoundsArea);
            } else {
                pRear = SampleExitPupil(Point2f(pFilm.x, pFilm.y), pLens, &exitPupilBoundsArea);
            }
            Ray rFilm = Ray(pFilm, pRear - pFilm, Infinity,0);
            DrawRayPathFromFilm(rFilm, false, false);
            printf(",");
        }
        printf("}]");
    }


    if (simpleWeighting)
        Warning("\"simpleweighting\" option with OmniCamera no longer "
                "necessarily matches regular camera images. Further, pixel "
                "values will vary a bit depending on the aperture size. See "
                "this discussion for details: "
                "https://github.com/mmp/pbrt-v3/issues/162#issuecomment-348625837");
}

OmniCamera::IntersectResult OmniCamera::TraceElement(const LensElementInterface &element, const Ray& rLens, 
    const Float& elementZ, Float& t, Normal3f& n, bool& isStop, 
    const ConvexQuadf& bounds = ConvexQuadf()) const {
    isStop = (element.curvatureRadius.x == 0);
    auto invTransform = Inverse(element.transform);
    Ray rElement = invTransform(rLens);
    if (isStop) {
        // The refracted ray computed in the previous lens element
        // interface may be pointed "backwards" in some
        // extreme situations; in such cases, 't' becomes negative.
        t = (elementZ - rElement.o.z) / rElement.d.z;
        if (rElement.d.z == 0.0 || t < 0) return MISS;
    } else {
        Float radius = element.curvatureRadius.x;
        Float zCenter = elementZ + element.curvatureRadius.x;
        if (!IntersectSphericalElement(radius, zCenter, rElement, &t, &n))
            return MISS;
    }
    CHECK_GE(t, 0);
    // Transform the normal back into the original space.
    n = element.transform(n);

    // Test intersection point against element aperture
    Point3f pHit = rElement(t);
    Float r2 = pHit.x * pHit.x + pHit.y * pHit.y;
    Float apertureRadius2 = element.apertureRadius.x * element.apertureRadius.x;
    if (r2 > apertureRadius2) return CULLED_BY_APERTURE;
    if (!bounds.Contains(Point2f(pHit))) {
        return CULLED_BY_APERTURE;
    }
    return HIT;
};

bool OmniCamera::TraceLensesFromFilm(const Ray &rCamera, 
    const std::vector<LensElementInterface>& interfaces, Ray *rOut,
    const Transform CameraToLens = Scale(1, 1, -1), const ConvexQuadf& bounds = ConvexQuadf()) const {
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
        // Compute intersection of ray with lens element
        Float t;
        Normal3f n;
        bool isStop;
        elementZ -= element.thickness;
        IntersectResult result = TraceElement(element, rLens, elementZ, t, n, isStop, bounds);
        if (result != HIT)
            return false;

        rLens.o = rLens(t);

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
        *rOut = Inverse(CameraToLens)(rLens);
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
        bool isStop;
        IntersectResult result = TraceElement(element, rLens, elementZ, t, n, isStop);
        if (result != HIT)
            return false;
        rLens.o = rLens(t);

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

    auto toZXCoord = [](const Point3f& p, const Transform& xForm) {
        const Point3f p3 = xForm(p);
        return Point2f(p3.z, p3.x);
    };

    for (size_t i = 0; i < elementInterfaces.size(); ++i) {
        const LensElementInterface &element = elementInterfaces[i];
        const Transform& xForm = element.transform;
        Float r = element.curvatureRadius.x;
        if (r == 0) {
            Point2f startU = toZXCoord(Point3f(element.apertureRadius.x,0,z), xForm);
            Point2f startL = toZXCoord(Point3f(-element.apertureRadius.x, 0, z), xForm);
            Point2f endU = toZXCoord(Point3f(2*element.apertureRadius.x, 0, z), xForm);
            Point2f endL = toZXCoord(Point3f(-2* element.apertureRadius.x, 0, z), xForm);
            // stop
            printf("{Thick, Line[{{%f, %f}, {%f, %f}}], ", startU.x, startU.y, endU.x, endU.y);
            printf("Line[{{%f, %f}, {%f, %f}}]}, ", startL.x, startL.y, endL.x, endL.y);
        } else {
            // TODO compute proper tilt
            Point2f C = toZXCoord(Point3f(0, 0, z + r), xForm);
            Float theta = std::abs(std::asin(element.apertureRadius.x / r));
            if (r > 0) {
                // convex as seen from front of lens
                Float t0 = Pi - theta;
                Float t1 = Pi + theta;
                printf("Circle[{%f, %f}, %f, {%f, %f}], ", C.x, C.y, r, t0, t1);
            } else {
                // concave as seen from front of lens
                Float t0 = -theta;
                Float t1 = theta;
                printf("Circle[{%f, %f}, %f, {%f, %f}], ", C.x, C.y, -r, t0, t1);
            }
            if (element.eta != 0 && element.eta != 1) {
                // TODO: re-enable
                /*
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
                */
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
    printf("{ ");
    if (!TraceLensesFromFilm(r, elementInterfaces, nullptr)) printf("Dashed, ");
    // Transform _ray_ from camera to lens system space
    static const Transform CameraToLens = Scale(1, 1, -1);
    Ray ray = CameraToLens(r);
    for (int i = elementInterfaces.size() - 1; i >= 0; --i) {
        const LensElementInterface &element = elementInterfaces[i];
        Float t;
        Normal3f n;
        bool isStop;
        elementZ -= element.thickness;
        IntersectResult result = TraceElement(element, ray, elementZ, t, n, isStop);
        if (result == MISS)
            goto done;

        printf("Line[{{%f, %f}, {%f, %f}}],", ray.o.z, ray.o.x, ray(t).z,
               ray(t).x);

        if (result == CULLED_BY_APERTURE)
            goto done;

        ray.o = ray(t);

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
        Float t;
        Normal3f n;
        bool isStop;
        IntersectResult result = TraceElement(element, ray, elementZ, t, n, isStop);
        if (result == MISS) return;

        printf("Line[{{%f, %f}, {%f, %f}}],", ray.o.z, ray.o.x, ray(t).z,
            ray(t).x);

        if (result == CULLED_BY_APERTURE) return;

        ray.o = ray(t);

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

void OmniCamera::ComputeThickLensApproximation(Float pz[2], Float fz[2]) const {
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

    Point3f finalElementTranslation = elementInterfaces[elementInterfaces.size() - 1].transform(Point3f(0,0,0));
    Point2f xy(finalElementTranslation.x, finalElementTranslation.y);
    Bounds2f projRearBounds(Point2f(-1.5f * rearRadius, -1.5f * rearRadius) + xy,
                            Point2f(1.5f * rearRadius, 1.5f * rearRadius) + xy);
    // TODO: More sophisticated handling of the exit pupil for microlens arrays
    if (HasMicrolens()) {
        return projRearBounds;
    }
    for (int i = 0; i < nSamples; ++i) {
        // Find location of sample points on $x$ segment and rear lens element
        Point3f pFilm(Lerp((i + 0.5f) / nSamples, pFilmX0, pFilmX1), 0, 0);
        Float u[2] = {RadicalInverse(0, i), RadicalInverse(1, i)};
        Point3f pRear(Lerp(u[0], projRearBounds.pMin.x, projRearBounds.pMax.x),
                      Lerp(u[1], projRearBounds.pMin.y, projRearBounds.pMax.y),
                      LensRearZ() + finalElementTranslation.z);

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
static Vector2f mapDiv(Vector2f v0, Vector2i v1) {
    return Vector2f(v0.x / (pbrt::Float)v1.x, v0.y / (pbrt::Float)v1.y);
}
static Point2f mapDiv(Point2f v0, Vector2f v1) {
    return Point2f(v0.x / v1.x, v0.y / v1.y);
}

Point2i OmniCamera::MicrolensIndex(const Point2f& p) const {
    Bounds2f extent = film->GetPhysicalExtent();
    Vector2f normFilm = mapDiv(p - extent.pMin, extent.pMax - extent.pMin);
    Vector2i d = microlens.dimensions;
    Vector2f df = Vector2f(d.x, d.y);
    Vector2f lensSpace = mapMul(normFilm, df);
    return Point2i(floor(lensSpace.x), floor(lensSpace.y));
}

Point3f OmniCamera::SampleMicrolensPupil(const Point2f &pFilm, const Point2f &lensSample,
    Float *sampleBoundsArea) const {
    Point2f lensIndex = Point2f(MicrolensIndex(pFilm));
    Vector2i d = microlens.dimensions;
    Vector2f df = Vector2f(d.x, d.y);
    Point2f sampledLensSpacePt = mapDiv(lensIndex+lensSample, df);
    // sample on microlens
    Point2f result2 = film->GetPhysicalExtent().Lerp(sampledLensSpacePt);
    if (sampleBoundsArea) *sampleBoundsArea = film->GetPhysicalExtent().Area() / (df.x*df.y);
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

Transform OmniCamera::MicrolensElement::ComputeCameraToMicrolens() const {
    return Transform();
}

Point2f OmniCamera::MicrolensCenterFromIndex(const Point2i& idx) const {
    Vector2f indexf(idx.x, idx.y);
    const Vector2i& d = microlens.dimensions;
    Vector2f normalizedLensCenter = mapDiv(indexf + Vector2f(0.5f, 0.5f), d);
    Point2f lensCenter = film->GetPhysicalExtent().Lerp(Point2f(normalizedLensCenter));
    if (idx.x >= 0 && idx.y >= 0 && idx.x < d.x && idx.y < d.y) {
        lensCenter += microlens.offsets[idx.y*d.x + idx.x];
    }
    return lensCenter;
}

OmniCamera::MicrolensElement OmniCamera::MicrolensElementFromIndex(const Point2i& idx) const {
    MicrolensElement element;
    element.index = idx;
    element.center = MicrolensCenterFromIndex(idx);
    const Point2i offsets[4] = { {0,0}, {0,1}, {1,0}, {1,1} };
    Point2f corners[4];
    for (int i = 0; i < 4; ++i) {
        corners[i] = Point2f(0, 0);
        Point2i cornerIdx = idx-Vector2i(offsets[i]);
        for (const auto& off : offsets) {
            Point2f neighborCenter = MicrolensCenterFromIndex(cornerIdx + off);
            corners[i] += neighborCenter;
        }
        corners[i] *= (Float)0.25;
        corners[i] += -element.center; // Center the corners
    }
    element.centeredBounds = ConvexQuadf(corners[0], corners[1], corners[2], corners[3]);
    return element;
}

OmniCamera::MicrolensElement OmniCamera::ComputeMicrolensElement(const Ray& filmRay) const {
    Point3f pointOnMicrolens(filmRay(microlens.offsetFromSensor / filmRay.d.z));
    Point2f pointOnMicrolens2(pointOnMicrolens);
    return MicrolensElementFromIndex(MicrolensIndex(pointOnMicrolens2));
}


bool OmniCamera::TraceFullLensSystemFromFilm(const Ray& rIn, Ray* rOut) const {
    if (HasMicrolens()) {
        MicrolensElement element = ComputeMicrolensElement(rIn);
        Transform CameraToMicrolens = Scale(1, 1, -1)*Translate({ -element.center.x, -element.center.y, 0.0f });
        Ray rAfterMicrolens;
        if (rOut) rAfterMicrolens = *rOut;
        if (!TraceLensesFromFilm(rIn, microlens.elementInterfaces, &rAfterMicrolens, CameraToMicrolens, element.centeredBounds)) {
            return false;
        }
        return TraceLensesFromFilm(rAfterMicrolens, elementInterfaces, rOut);
    } else {
        return TraceLensesFromFilm(rIn, elementInterfaces, rOut);
    }
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

    if (!TraceFullLensSystemFromFilm(rFilm, ray)) {
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
            return cos4Theta;
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

    Float microlensSensorOffset = params.FindOneFloat("microlenssensoroffset", 0.001);
    int microlensSimulationRadius = params.FindOneInt("microlenssimulationradius", 0);

    if (microlensSimulationRadius != 0) {
        Warning("We only currently support simulating one microlens per pixel, switching microlensSimulationRadius to 0");
    }

    if (lensFile == "") {
        Error("No lens description file supplied!");
        return nullptr;
    }
    // Load element data from lens description file
    std::vector<OmniCamera::LensElementInterface> lensInterfaceData;

    // Load element data from lens description file
    std::vector<OmniCamera::LensElementInterface> microlensData;

    std::vector<Vector2f> microlensOffsets;

    Vector2i microlensDims;

    auto endsWith = [](const std::string& str, const std::string& suffix) {
        return str.size() >= suffix.size() && 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
    };

    if (endsWith(lensFile, ".dat")) {
        Error(
            "Invalid lens file format for file \"%s\"," 
            "must use new json. To convert \"realistic\""
            "camera style .dat files, use the `lenstool` command:"
            "  lenstool convert [input .dat filename] [output .json filename]",
            lensFile.c_str());
        return nullptr;
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
                for (int r = 0; r < 3; ++r) {
                    for (int c = 0; c < 3; ++c) {
                        m[r][c] = (Float)t[c][r];
                    }
                    // Translation specified in mm, needs to be converted to m
                    m[r][3] = (Float)t[3][r] * (Float).001;
                }
                m[3][0] = m[3][1] = m[3][2] = (Float)0.0;
                m[3][3] = (Float)1.0;
                return Transform(m);
            };
            auto toIORSpectrum = [lensFile](json jiors) {
                // Stored in columns in json, but pbrt stores matrices 
                // in row-major order.
                if (jiors.is_number()) { // Perfectly fine to have no transform
                    return (Float)jiors;
                }
                if (!(jiors.is_array()) || (jiors.size() != 2) || (!jiors[0].is_array())
                    || (!jiors[1].is_array()) || (jiors[0].size() != jiors[1].size())
                    || (!(jiors[0][0].is_number())) || (!(jiors[1][0].is_number()))) {
                    Error("Invalid ior in lens specification file \"%s\","
                        " must be either a single float, or a pair of parallel arrays."
                        " The first array should be a list of wavelengths (in nm),"
                        " The second array should be the ior at those wavelengths.", lensFile.c_str());
                }
                size_t numSamples = jiors[0].size();
                std::vector<Float> wavelengths(numSamples);
                std::vector<Float> iors(numSamples);
                for (int i = 0; i < numSamples; ++i) {
                    wavelengths[i] = (Float)jiors[0][i];
                    iors[i] = (Float)jiors[1][i];
                }
                SampledSpectrum s = SampledSpectrum::FromSampled(wavelengths.data(), iors.data(), (int)numSamples);
                for (int i = 0; i < numSamples-1; ++i) {
                    if (std::abs(s[i] - s[i + 1]) > 0.001) {
                        Error("Invalid ior in lens specification file \"%s\","
                            " spectrum must be constant (wavelength-varying ior NYI), back-to-back (%f-%f = %f)", 
                            lensFile.c_str(), s[i], s[i+1], std::abs(s[i] - s[i + 1]));
                    }
                }
                return s[0];
            };

            auto toLensElementInterface = [toVec2, toTransform, toIORSpectrum, apertureDiameter](json surf) {
                OmniCamera::LensElementInterface result;
                // Convert mm to m
                result.apertureRadius   = toVec2(surf["semi_aperture"]) * (Float).001;
                result.conicConstant    = toVec2(surf["conic_constant"]) * (Float).001;
                result.curvatureRadius  = toVec2(surf["radius"]) * (Float).001;
                result.eta              = toIORSpectrum(surf["ior"]);
                result.thickness        = Float(surf["thickness"]) * (Float).001;
                result.transform        = toTransform(surf["transform"]);
                if (result.curvatureRadius.x == 0.0f) {
                    Float apertureRadius = apertureDiameter * (Float).001 / Float(2.);
                    if (apertureRadius > result.apertureRadius.x) {
                        Warning(
                            "Specified aperture radius %f is greater than maximum "
                            "possible %f.  Clamping it.",
                            apertureRadius, result.apertureRadius.x);
                    } else {
                        result.apertureRadius.x = apertureRadius;
                        result.apertureRadius.y = apertureRadius;
                    }
                }
                return result;
            };

            if (jsurfaces.is_array() && jsurfaces.size() > 0) {
                for (auto jsurf : jsurfaces) {
                    lensInterfaceData.push_back(toLensElementInterface(jsurf));
                }
            } else {
                Error("Error, lens specification file without a valid surface array \"%s\".",
                    lensFile.c_str());
                return nullptr;
            }
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
                        microlensOffsets.push_back(toVec2(offset));
                    }
                }

                auto mljSurfaces = microlens["surfaces"];

                if (mljSurfaces.is_array() && mljSurfaces.size() > 0) {
                    for (auto jsurf : mljSurfaces) {
                        microlensData.push_back(toLensElementInterface(jsurf));
                    }
                }
                else {
                    Error("Error, microlens specification without a valid surface array in \"%s\".",
                        lensFile.c_str());
                    return nullptr;
                }
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
                               lensInterfaceData, microlensData, microlensDims, microlensOffsets, microlensSensorOffset, film, medium);
}

}  // namespace pbrt
