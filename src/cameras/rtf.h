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

#ifndef PBRT_CAMERAS_BlackBox_H
#define PBRT_CAMERAS_BlackBox_H

// cameras/BlackBox.h*
#include "pbrt.h"
#include "camera.h"
#include "film.h"

namespace pbrt {

// RTFCamera Declarations
class RTFCamera : public Camera {
  public:
    /*
    struct LensElementInterface {
        LensElementInterface() {}
        LensElementInterface(Float cRadius, Float aRadius,
            Float thickness, Float ior) :
            curvatureRadius(cRadius,cRadius),
            apertureRadius(aRadius ,aRadius),
            conicConstant((Float)0.0, (Float)0.0),
            transform(Transform()),
            thickness(thickness),
            eta(ior) {}
        Vector2f curvatureRadius;
        Vector2f apertureRadius;
        Vector2f conicConstant;
        Transform transform;
        Float thickness;
        Float eta;
        std::vector<Float> asphericCoefficients;
        Float zMin;
        Float zMax;
    };
    struct MicrolensData {
        std::vector<LensElementInterface> elementInterfaces;
        float offsetFromSensor;
        std::vector<Vector2f> offsets;
        Vector2i dimensions;
        // Non-physical term
        int simulationRadius;
    };
     */
    struct LensPolynomialTerm {
        LensPolynomialTerm() {}
        LensPolynomialTerm(std::string n, std::vector<Float> tr,
                        std::vector<Float> tu, std::vector<Float> tv,
                        std::vector<Float> coeff) :
                        name(n), termr(tr), termu(tu),
                        termv(tv), coeff(coeff) {}
        std::string name;
        std::vector<Float> termr;
        std::vector<Float> termu;
        std::vector<Float> termv;
        std::vector<Float> coeff;

    

    };

    struct RTFVignettingTerms {
        RTFVignettingTerms() {}
        RTFVignettingTerms(Float circlePlaneZ,int exitpupilIndex, std::vector<Float> pupilPos, std::vector<Float> pupilRadii, std::vector<Float> circleRadii,  std::vector<Float> circleSensitivities):
        circlePlaneZ(circlePlaneZ), exitpupilIndex(exitpupilIndex), pupilPos(pupilPos), pupilRadii(pupilRadii), circleRadii(circleRadii),circleSensitivities(circleSensitivities) {}
        Float circlePlaneZ;
        int exitpupilIndex; // Index that indicates main exit pupil
        std::vector<Float> pupilPos;
        std::vector<Float> pupilRadii;
        std::vector<Float> circleRadii;
        std::vector<Float> circleSensitivities;
    };

    /*
    // RTFCamera Public Methods
    RTFCamera(const AnimatedTransform &CameraToWorld, Float shutterOpen,
                    Float shutterClose, Float apertureDiameter, Float filmdistance,
                    Float focusDistance, bool simpleWeighting, bool noWeighting,
                    bool caFlag, const std::vector<RTFCamera::LensElementInterface> &lensData,
                    const std::vector<RTFCamera::LensElementInterface> &microlensData,
                    Vector2i microlensDims, const std::vector<Vector2f> & microlensOffsets,
                    float microlensSensorOffset, int microlensSimulationRadius, Film *film, const Medium *medium);
     */
    
    
    // RTFCamera Public Methods

    RTFCamera(const AnimatedTransform &CameraToWorld, Float shutterOpen, Float shutterClose, Float apertureDiameter, Float filmdistance, Float lensThickness, Float planeOffset, bool caFlag, Film *film, const Medium *medium,
        std::vector<std::map<std::string,RTFCamera::LensPolynomialTerm>> polynomialMaps, std::string bbmode,
        std::vector<RTFVignettingTerms>);
    
    Float GenerateRay(const CameraSample &sample, Ray *) const;

  private:
    // RTFCamera Private Declarations

    enum IntersectResult {MISS,CULLED_BY_APERTURE,HIT};
    const bool caFlag;
    const Float filmDistance;
    const Float planeOffset;
    const Float lensThickness;
    const Point2f exitPupilBounds;
    int pupilIndex;
    Ray ApplyPolynomial(const Ray &thisRay) const;
    Ray RotateRays(const Ray &thisRay, Float deg) const;
    Float PolynomialCal(const Ray &thisRay, std::string name, Vector2f &radiusRotation) const;
    Vector2f Pos2RadiusRotation(const Point3f pos) const;
    bool IsValidRay(const Ray &rCamera) const;
    bool IsValidRayCircles(const Ray &rotatedRay) const;
    bool TraceLensesFromFilm(const Ray &ray, Ray *rOut,
    const Transform CameraToLens) const;
    Point3f SampleExitPupil(const Point2f &pFilm, const Point2f &lensSample) const;
    std::map<std::string, RTFCamera::LensPolynomialTerm> poly;
    


    // Wavelength dependent RTF : vectorized. Each element corresponds to a given wavelength
    std::vector<Float> polyWavelengths_nm; // wavelengths read from file
    std::vector<std::map<std::string, RTFCamera::LensPolynomialTerm>> polynomialMaps; // Each element has corresponding wavelength
    std::vector<RTFCamera::RTFVignettingTerms> vignettingTerms;

    std::vector<Float> pupilPos;
    std::vector<Float> pupilRadii;
    std::vector<Float> circleRadii;
    std::vector<Float> circleSensitivities;
    Float circlePlaneZ; // Plane where circles were fitted
    std::string bbmode;
    
    /*
    // RTFCamera Private Data
    const bool simpleWeighting;
    const bool noWeighting;
    const bool caFlag;
    std::vector<LensElementInterface> elementInterfaces;
    std::vector<Bounds2f> exitPupilBounds;

    MicrolensData microlens;

    struct MicrolensElement {
        Point2f center;
        ConvexQuadf centeredBounds;
        Point2i index;
        Transform ComputeCameraToMicrolens() const;
    };

    // RTFCamera Private Methods
    Float LensRearZ() const { return elementInterfaces.back().thickness; }
    Float LensFrontZ() const {
        Float zSum = 0;
        for (const LensElementInterface &element : elementInterfaces)
            zSum += element.thickness;
        return zSum;
    }
    Float RearElementRadius() const {
        return elementInterfaces.back().apertureRadius.x;
    }
    bool TraceLensesFromFilm(const Ray &ray, const std::vector<LensElementInterface>& interfaces, Ray *rOut,
        const Transform CameraToLens, const ConvexQuadf& bounds) const;
    float TToBackLens(const Ray &ray, const std::vector<LensElementInterface>& interfaces,
        const Transform CameraToLens, const ConvexQuadf& bounds) const;

    static bool IntersectSphericalElement(Float radius, Float zCenter,
                                          const Ray &ray, Float *t,
                                          Normal3f *n);
    bool TraceLensesFromScene(const Ray &rCamera, Ray *rOut) const;
    void DrawLensSystem() const;
    void DrawRayPathFromFilm(const Ray &r, bool arrow,
                             bool toOpticalIntercept) const;
    void DrawRayPathFromScene(const Ray &r, bool arrow,
                              bool toOpticalIntercept) const;
    static void ComputeCardinalPoints(const Ray &rIn, const Ray &rOut, Float *p,
                                      Float *f);
    void ComputeThickLensApproximation(Float pz[2], Float f[2]) const;
    Float FocusThickLens(Float focusDistance);
    Float FocusBinarySearch(Float focusDistance);
    Float FocusDistance(Float filmDist);
    Bounds2f BoundExitPupil(Float pFilmX0, Float pFilmX1) const;
    void RenderExitPupil(Float sx, Float sy, const char *filename) const;
    Point3f SampleExitPupil(const Point2f &pFilm, const Point2f &lensSample,
                            Float *sampleBoundsArea) const;

    IntersectResult TraceElement(const LensElementInterface &element, const Ray& rLens, const Float& elementZ,
         Float& t, Normal3f& n, bool& isStop, const ConvexQuadf& bounds) const;

    void TestExitPupilBounds() const;

    Point2i MicrolensIndex(const Point2f& p) const;
    Point2f MicrolensCenterFromIndex(const Point2i& idx) const;
    MicrolensElement MicrolensElementFromIndex(const Point2i& idx) const;
    MicrolensElement ComputeMicrolensElement(const Ray & filmRay) const;

    bool TraceFullLensSystemFromFilm(const Ray & rIn, Ray * rOut) const;

    Point3f SampleMicrolensPupil(const Point2f &pFilm, const Point2f &lensSample,
        Float *sampleBoundsArea) const;


    bool HasMicrolens() const;
    */
};

RTFCamera *CreateRTFCamera(const ParamSet &params,
                            const AnimatedTransform &cam2world,
                            Film *film, const Medium *medium);

}  // namespace pbrt

#endif  // PBRT_CAMERAS_BlackBox_H
