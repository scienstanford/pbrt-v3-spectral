
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

// cameras/BlackBox.cpp*
#include "cameras/rtf.h"
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
#include <math.h>       /* pow */
#include "rng.h"
#include <gsl/gsl_roots.h> // For solving biconic surface intersections.
#include <gsl/gsl_errno.h>


using json = nlohmann::json;

namespace pbrt {

/*
STAT_PERCENT("Camera/Rays vignetted by lens system", vignettedRays, totalRays);
Float sgn(Float val) {
    return Float((Float(0) < val) - (val < Float(0)));
}
*/

// RTFCamera Method Definitions
// Pupil index: choose which pupil in the list of pupils counts as the exit pupil
RTFCamera::RTFCamera(const AnimatedTransform &CameraToWorld, Float shutterOpen, Float shutterClose, Float apertureDiameter, Float filmdistance, Float lensThickness, Float planeOffset, bool caFlag, Film *film, const Medium *medium,
        std::vector<std::map<std::string,RTFCamera::LensPolynomialTerm>> polynomialMaps, std::string bbmode,
        std::vector<RTFVignettingTerms>):
        Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
        caFlag(caFlag), filmDistance(filmdistance), planeOffset(planeOffset),lensThickness(lensThickness),polynomialMaps(polynomialMaps),bbmode(bbmode),vignettingTerms(vignettingTerms)
        {}

Vector2f RTFCamera::Pos2RadiusRotation(const Point3f pos) const{
    // Convert position to radius
    Float radius = std::sqrt(pos.x * pos.x + pos.y * pos.y);
    Float deg = std::atan2(pos.y,pos.x);
    
    // Determine whether to rotate 180 degrees % ATAN2 covers this already
    //if ( (pos.x < 0 && pos.y < 0) || (pos.x < 0 && pos.y > 0)) {
        //deg += Pi;
    //}
    
    Vector2f res = Vector2f(radius, Degrees(deg));
    
    return res;
}

Float RTFCamera::PolynomialCal(const Ray &thisRay, std::string name, Vector2f &radiusRotation) const{
    
    Point3f pos= thisRay.o;
    
//    pos = Point3f(0.4596, 0.3857, -2.0224);
    radiusRotation = Pos2RadiusRotation(pos);
    
    // Rotate rays
    Ray rotatedRay = RotateRays(thisRay, 90 - radiusRotation.y);
    Vector3f dir = Normalize(rotatedRay.d);
    

    LensPolynomialTerm thisOut = poly.at(name);
    Float res = 0;
    for (int i = 0; i < thisOut.termr.size(); i++) {
        res += (std::pow(radiusRotation.x * 1000, thisOut.termr[i]) * std::pow(dir.x, thisOut.termu[i]) * std::pow(dir.y, thisOut.termv[i])) * thisOut.coeff[i];
    }
    
    return res;
}

Ray RTFCamera::ApplyPolynomial(const Ray &thisRay) const{
    // Get radius and rotation
    Vector2f radiusRotation;
    
    // Main body:
    // Use 'poly' parameter - a map from name to the terms.
    // Optional keys are: 'outx', 'outy', 'outu', 'outv', 'outw'
    // Term contains: termr, termu, termv, coeff (input ray parameters)
    Float x = PolynomialCal(thisRay, "outx", radiusRotation) * 0.001f;
    Float y = PolynomialCal(thisRay, "outy", radiusRotation) * 0.001f;
    Float u = PolynomialCal(thisRay, "outu", radiusRotation);
    Float v = PolynomialCal(thisRay, "outv", radiusRotation);
    if (1 - u * u - v * v < 0) {
        Warning("Problemetic ray fitting.");
    }
    Float w = std::sqrt(std::abs(1 - u * u - v * v));
    
    // The outocming Ray starts at the output plane: 
//    Float outputPlane_z = -(lensThickness+filmDistance); // Lens thickness includes 2x planeOffsets
    Float outputPlane_z = -(filmDistance+lensThickness+planeOffset); /// Lens thickness: WIHOUT input plane
    // coordinate systems of lens and film are opposite -> that's why there is an inversion

    Ray rOut = Ray(Point3f(x, y, outputPlane_z), Vector3f(u, v, -w)); // Original
    
    // Rotate the output ray
    rOut = RotateRays(rOut, radiusRotation.y - 90);
    
    return rOut;
}

Ray RTFCamera::RotateRays(const Ray &thisRay, Float deg) const{
    Transform rot = Rotate(deg, Vector3f(0, 0, 1));
    Ray rRot = rot(thisRay);
    return rRot;
}

bool RTFCamera::IsValidRay(const Ray &rCamera) const{
    Vector3f dir = Normalize(rCamera.d);

    
    for (int i = 0; i < pupilPos.size(); i++) {
        // Calculate the length
        Float alpha = pupilPos[i] / dir.z;
//        std::cout << pupilRadii[i] << "\n";
    
        Point3f validBound = rCamera.o + alpha * dir; // TG EDIT TEMP
        if (std::sqrt(validBound.x * validBound.x + validBound.y * validBound.y) >= pupilRadii[i]) {
            return false;
        }
    } 

    return true;
}

bool RTFCamera::IsValidRayCircles(const Ray &rotatedRay) const{
    Vector3f dir = Normalize(rotatedRay.d);

      Float alpha = -circlePlaneZ / dir.z;  // minus sign because rotatedRay already has reversal of Z axi
     Point3f validBound = rotatedRay.o + alpha * dir; // TG EDIT TEMP
     Float offaxisDistance=rotatedRay.o.y; // in input plane
     Float xsquared = validBound.x * validBound.x;
     
     // CirclePlane instead of pupils
         for (int i = 0; i < circleRadii.size(); i++) {
        // Calculate the length

        Float distanceFromCenterY = validBound.y-offaxisDistance*circleSensitivities[i];
        
        if (std::sqrt(xsquared+ distanceFromCenterY*distanceFromCenterY) >= circleRadii[i]) {
            return false;
        }
    } 
    return true;
}

bool RTFCamera::TraceLensesFromFilm(const Ray &rCamera,
    Ray *rOut,
    const Transform CameraToLens = Scale(1, 1, -1)) const {
    
    /*if (!IsValidRay(rCamera)) {
        return false;
    }*/
    
    Ray rLens = CameraToLens(rCamera);
    
    // Since the polynomial is fitted only for 10 deg, we have a check. If not, reject the ray
    Vector2f radiusRotation = Pos2RadiusRotation(rLens.o);
    // Rotate rays
    Ray rotatedRay = RotateRays(rLens, 90 - radiusRotation.y);

//if(IsValidRay(rCamera)){
    //std::cout << IsValidRay(rCamera) << "," << IsValidRayCircles(rotatedRay) << "\n";
//    std::cout << (rCamera) << "," << (rotatedRay) << "\n";
//}


    if (!IsValidRayCircles(rotatedRay)) {
        return false;
    }


    Vector3f dir = Normalize(rotatedRay.d);
    //Float deg = Degrees(std::atan(std::sqrt(dir.x * dir.x + dir.y * dir.y) /std::abs(dir.z)));
    Float deg = Degrees(std::atan2(std::sqrt(dir.x * dir.x + dir.y * dir.y),std::abs(dir.z)));
    //if (deg > 10) {
    //        std::cout<< "Degree is larger than 10\n";
    //    return false;
    //}
    //
    
    
    // Apply polynomial fitting
    // Debug: Hardcode one light ray
    // rLens = Ray(Point3f(0, 0.0005333, -0.0020224), Vector3f(0.3060, -0.2223, 0.9257));
    // When use the test ray, the return should be: (0.006703, -0.002456, 0.2218, -0.3247)
    rLens = ApplyPolynomial(rLens);
    rLens.d = Normalize(rLens.d);
    // Rotate rays back to camera space
    if ((rOut != nullptr)) {
        *rOut = Inverse(CameraToLens)(rLens);
    }
    else{
        return false;
    }
    return true;
    
    /*
    // Test case of just flipping the image
    Transform flipTransform = Scale(-1, -1, 1);
    rLens.d = Vector3f(-rLens.d.x, -rLens.d.y, rLens.d.z);
    
    if (rOut != nullptr) {
        *rOut = Inverse(CameraToLens)(rLens);
    }
    return true;
     */
    
    
    /*
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

        // Update ray path for element interface interaction
        if (!isStop) {
            rLens.o = rLens(t);
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
     */

}
/*
float RTFCamera::TToBackLens(const Ray &rCamera,
    const std::vector<LensElementInterface>& interfaces,
    const Transform CameraToLens = Scale(1, 1, -1), const ConvexQuadf& bounds = ConvexQuadf()) const {
    // Transform _rCamera_ from camera to lens system space
    Ray rLens = CameraToLens(rCamera);

    const LensElementInterface &element = interfaces[interfaces.size() - 1];
    // Compute intersection of ray with lens element
    Float t;
    Normal3f n;
    bool isStop;
    Float elementZ = -element.thickness;
    IntersectResult result = TraceElement(element, rLens, elementZ, t, n, isStop, bounds);
    if (result != HIT)
        return std::numeric_limits<Float>::infinity();
    return t;
}
*/

/*
// Might be useful for bdpt
bool RTFCamera::TraceLensesFromScene(const Ray &rCamera,
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

        // Update ray path for from-scene element interface interaction
        if (!isStop) {
            rLens.o = rLens(t);
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
*/




/*
bool RTFCamera::TraceFullLensSystemFromFilm(const Ray& rIn, Ray* rOut) const {
    
    if (HasMicrolens()) {
        MicrolensElement centerElement = ComputeMicrolensElement(rIn);
        Float tMin = std::numeric_limits<Float>::infinity();
        int R = microlens.simulationRadius;
        Point2i cIdx = centerElement.index;
        MicrolensElement toTrace = centerElement;
        // Check to find the first microlens we intersect with
        // Could be sped up by only checking the directions of the projection of the ray
        for (int y = -R; y <= R; ++y) {
            for (int x = -R; x <= R; ++x) {
                const MicrolensElement el = MicrolensElementFromIndex(cIdx + Vector2i(x,y));
                float newT = TToBackLens(rIn, microlens.elementInterfaces, el.ComputeCameraToMicrolens(), el.centeredBounds);
                if (newT < tMin) {
                    tMin = newT;
                    toTrace = el;
                }
            }
        }
        if (tMin < std::numeric_limits<Float>::infinity()) {
            Ray rAfterMicrolens;
            if (rOut) rAfterMicrolens = *rOut;
            if (!TraceLensesFromFilm(rIn, microlens.elementInterfaces, &rAfterMicrolens, toTrace.ComputeCameraToMicrolens(), toTrace.centeredBounds)) {
                return false;
            }

            bool result = TraceLensesFromFilm(rAfterMicrolens, elementInterfaces, rOut);
            / *
            static int paths = 0;
            if (result && cIdx != toTrace.index) {
                ++paths;
                if (paths % 100 == 0) {
                    printf("Contributing non-center microlens path count: %d\n", paths);
                }
            }* /
            return result;
        } else {
            return false;
        }
    } else {
        return TraceLensesFromFilm(rIn, elementInterfaces, rOut);
    }
    
}
*/


Point3f RTFCamera::SampleExitPupil(const Point2f &pFilm, const Point2f &lensSample) const {
    // lensSample comes from 
    
    /*
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
     */
    
    // ZHENG: 
    // lensSample.X encodes  RADIUS ?
    // lensSample.Y encodes Degree?
    
    // TG: Lerp(t,x1,x2)  : (1-t)*x1+x2
    // So this proves that lensSample.X is a normalized parameter
        // TODO: refactor this to have mre meanintful names?

        // lensSample.x: gives a uniform distribution between 0  and 1

        // exitPupilBounds.x = lowerbound (normally zero)
        // exitPupilBounds.y = lowerbGeound (normally exit pupil radius)
        // Question: Why take square of exitpupilbounds?
        
        //lensSample : samples in box [0 1)^2
    // Radius
//     Float rSample = std::sqrt(Lerp(lensSample.x, exitPupilBounds.x * exitPupilBounds.x,
//                                             exitPupilBounds.y * exitPupilBounds.y));

// // TG: Sample A box
//     //Float xsample = Lerp(lensSample.x, -exitPupilBounds.y, exitPupilBounds.y);
//     //Float ysample = Lerp(lensSample.y, -exitPupilBounds.y, exitPupilBounds.y);
//     //std::cout << exitPupilBounds.y*exitPupilBounds.y;
//     // Degree
//         Float dSample = Lerp(lensSample.y, 0, 2 * Pi);
    
//     //Point2f pLens = exitPupilBounds.Lerp(lensSample);
//     Point2f pSample = Point2f(rSample * std::cos(dSample), rSample * std::sin(dSample));



    // Uniformly sample the exit pupil using PBRT function (TG)
    // Note ther eis also a function "UniformSampleDisk" but the PBRT book does not recommend
    // this because it distortes area (something related to the monte carlo stratified Sample )

    Float lensRadius=exitPupilBounds.y;
    Point2f pSample = lensRadius * ConcentricSampleDisk(lensSample);
//    Point2f pSample = lensRadius * UniformSampleDisk(lensSample);

    // pSample = Point2f(xsample,ysample);

    // Define input plane explicitly TG
    //  The pupil positins are defined relative to the inputplane_z.  



    Float inputPlane_z=filmDistance-planeOffset;

    // The Sample has been taken at the exit pupil, and now has to be projected back to the input plane 
    // This done using the Lerp (linear interpolation) function. First the interpolation factor (ratio) is calculated:
    Float ratio = inputPlane_z / (inputPlane_z + pupilPos[pupilIndex]);
   // Float ratio = (filmDistance) / ((filmDistance) + pupilPos[pupilIndex]); // original
   //  The interpolation is done:
    Point2f pLens = Lerp(ratio, pFilm, pSample);
     
    //Float tmp = Lerp(ratio, 0, filmDistance + pupilPos[pupilIndex]);
    //return Point3f(pLens.x, pLens.y, filmDistance); // Original

    // A point is made 
    return Point3f(pLens.x, pLens.y, inputPlane_z);
}

// This function is supposed to generate the output ray
Float RTFCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    
    // sample.pFilm: Point2f on the film
    // sample.pLens: Point2f on the lens : ()

    ProfilePhase prof(Prof::GenerateCameraRay);

    // ++totalRays;
    // Find point on film, _pFilm_, corresponding to _sample.pFilm_
    Point2f s(sample.pFilm.x / film->fullResolution.x,
        sample.pFilm.y / film->fullResolution.y);
        // Some kind of coordinate transformation
        //GetPhysicalExtent() returns the actual extent of the film in the scene. This information is specifically needed by the RealisticCamera. 
    Point2f pFilm2 = film->GetPhysicalExtent().Lerp(s);
    Point3f pFilm(-pFilm2.x, pFilm2.y, 0); //  ZHENG: Why reverse the direction?

    // We want to find the sample on the lens, this information is in sample.pLens. But how?
    // 
    Point3f pRear;
    pRear = SampleExitPupil(Point2f(pFilm.x, pFilm.y), sample.pLens);
    
    Ray rFilm = Ray(pFilm, pRear - pFilm, Infinity,
        Lerp(sample.time, shutterOpen, shutterClose));
    Ray rRear = Ray(pRear, pRear - pFilm, Infinity,
                    Lerp(sample.time, shutterOpen, shutterClose));
    

    // This is the core function of blackbox camera model
    if(!TraceLensesFromFilm(rRear, ray)) {
            // ++vignettedRays;
            return 0;
    }
    
    // Finish initialization of _RTFCamera_ ray
    // This is the output ray
    *ray = CameraToWorld(*ray);
    ray->d = Normalize(ray->d);
    ray->medium = medium;
    
    /*
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

    // Finish initialization of _RTFCamera_ ray
    *ray = CameraToWorld(*ray);
    ray->d = Normalize(ray->d);
    ray->medium = medium;

    // Return weighting for _RTFCamera_ ray
    if (HasMicrolens()) {
        // TODO: Proper weighting
        Float cosTheta = Normalize(rFilm.d).z;
        Float cos4Theta = (cosTheta * cosTheta) * (cosTheta * cosTheta);
        if (simpleWeighting) {
            // Normalize by eaxpecting that on average rays only make it through
            // one microlens to the outside world
            Float simDiameter = (microlens.simulationRadius*2.0 + 1.0);
            return cos4Theta * (simDiameter*simDiameter);
        }
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
     */
    
    Float cosTheta = Normalize(rFilm.d).z;
    Float cos4Theta = (cosTheta * cosTheta) * (cosTheta * cosTheta);
    Float pupilArea = pupilRadii[pupilIndex] * pupilRadii[pupilIndex] * Pi;
    return cos4Theta * pupilArea / (filmDistance * filmDistance);
}

RTFCamera *CreateRTFCamera(const ParamSet &params,
    const AnimatedTransform &cam2world,
    Film *film, const Medium *medium) {
    Float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    Float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
            shutterclose, shutteropen);
        std::swap(shutterclose, shutteropen);
    }
    
    Float apertureDiameter = params.FindOneFloat("aperturediameter", 1.0);

    // Add functionality to hard set the film distance
    Float filmDistance = params.FindOneFloat("filmdistance", 0);


    // Chromatic aberration flag
    bool caFlag = params.FindOneBool("chromaticAberrationEnabled", false);
    
    // Create bound
    // Point2f topLeft = params.FindOnePoint2f("topleftplane", Point2f(-0.0008, 0.0008));
    // Point2f bottomRight = params.FindOnePoint2f("bottomrightplane", Point2f(0.0008, -0.0008));
    // ZLY: Changed exit pupil bounds definition: Point2f now contains (0, radius^2)
    Point2f exitPupilBounds(0, INFINITY);
    
    int pupilIndex = -1;
    /*
    return new RTFCamera(cam2world, shutteropen, shutterclose,
                               apertureDiameter, filmDistance, focusDistance, simpleWeighting, noWeighting, caFlag,
                               lensInterfaceData, microlensData, microlensDims, microlensOffsets, microlensSensorOffset, microlensSimulationRadius, film, medium);
     */
    
    // Read in lens file
    std::string lensFile = params.FindOneFilename("lensfile", "");
    if (lensFile == "") {
        Error("No lens description file supplied!");
        return nullptr;
    }
    
    
    
    std::map<std::string, RTFCamera::LensPolynomialTerm> poly;

    std::vector<Float> polyWavelengths_nm; // wavelengths read from file
    std::vector<std::map<std::string, RTFCamera::LensPolynomialTerm>> polynomialMaps; // Each element has corresponding wavelength
    std::vector<RTFCamera::RTFVignettingTerms> vignettingTerms;
   
    // A pretty cool function
    auto endsWith = [](const std::string& str, const std::string& suffix) {
        return str.size() >= suffix.size() && 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
    };
    
    Float lensThickness = 0;
    Float planeOffset = 0;
    std::vector<Float> pupilPos;
    std::vector<Float> pupilRadii;
    std::vector<Float> circleRadii;
    std::vector<Float> circleSensitivities;


    Float circlePlaneZ;
    
    
    if (endsWith(lensFile, ".json")) {
        /*
         The format of polynomial lens would be:
         {
          "description": "equivalent lens poly",
          "name": "polynomial",
          "thickness": xxxx,
          "poly": [
           {
            "outputname": "outx",
            "termr":
            "termu":
            "termv":
            "coeff":
           },

          ]
         }
         
         Term and coeff will be stored in map.
         */
        // Polynomial lens case
        std::ifstream i(lensFile);
        json j;
        if (i && (i>>j)) {
            // Write name and description
            if (j["name"].is_string()) {
                LOG(INFO) << StringPrintf("Loading polynomial lens %s.\n",
                    j["name"].get<std::string>().c_str());
                    std::cout << j["name"].get<std::string>().c_str() << "\n";
            }
            if (j["description"].is_string()) {
                LOG(INFO) << StringPrintf("%s\n",
                    j["description"].get<std::string>().c_str());
            }
            
            if (j["thickness"].is_number()) {
                lensThickness = (Float) j["thickness"] * 0.001f;
            }
            if (j["planeoffset"].is_number()) {
                planeOffset = (Float) j["planeoffset"] * 0.001f; 
            }
            if (j["circlePlaneZ"].is_number()) {
                circlePlaneZ = (Float) j["circlePlaneZ"] * 0.001f; 
            }



            
            auto toTerms = [] (json jterms) {
                std::vector<Float> res;
                if (jterms.is_null()) {
                    return res;
                }
                
                for (int i = 0; i < jterms.size(); ++i) {
                    json thisnum = jterms[i];
                    if (!thisnum.is_number()) {
                        Error("Invalid polynomial in lens specification file must be an array of floats");
                    }
                    res.push_back((Float)thisnum);
                }
                return res;
            };

        auto toPolynomialStruct = [toTerms] (json jp) {
                RTFCamera::LensPolynomialTerm result;
                result.name = jp["outputname"].get<std::string>();
                result.termr = toTerms(jp["termr"]);
                result.termu = toTerms(jp["termu"]);
                result.termv = toTerms(jp["termv"]);
                result.coeff = toTerms(jp["coeff"]);
                return result;
            };

            auto toLensPolynomialTerms = [toTerms] (json jp) {
                RTFCamera::LensPolynomialTerm result;
                result.name = jp["outputname"].get<std::string>();
                result.termr = toTerms(jp["termr"]);
                result.termu = toTerms(jp["termu"]);
                result.termv = toTerms(jp["termv"]);
                result.coeff = toTerms(jp["coeff"]);
                return result;
            };
            

            auto toVignetTerms = [exitPupilBounds,pupilIndex,toTerms] (json sp) {
                RTFCamera::RTFVignettingTerms result; 
                Float mm_to_meter = 0.001f;

                result.circlePlaneZ = mm_to_meter*(Float)sp["circlePlaneZ"];
                result.pupilPos = toTerms(sp["pupilpos"]);
                result.pupilRadii = toTerms(sp["pupilradii"]);
                result.circleRadii = toTerms(sp["circleRadii"]);
                result.circleSensitivities = toTerms(sp["circleSensitivities"]);



                Float smallestBound  = INFINITY; //for initialization to find smallest pupil
                for (int i = 0; i < result.pupilPos.size(); i++) {
                        // Convert from mm to meter
                     

                    result.pupilPos[i] = mm_to_meter*result.pupilPos[i];
                    result.pupilRadii[i] = mm_to_meter*result.pupilRadii[i];
                    result.circleRadii[i] = mm_to_meter*result.circleRadii[i];
                    result.circleSensitivities[i] = mm_to_meter*result.circleSensitivities[i]; // sensitivies do not need to be converted because they are a ratio of distances
                    
                     // TG: in this implementation y is set to infinity, x to 0
                    // Goal: Find the smallest pupil, only go into the if statement if the current exit pupil bound is still
                    // larger than the pupilRadii[i]
                    
                    if (smallestBound > result.pupilRadii[i]) {
                     // Take the max of the two in case pupilRadii is even smaller than value of interest.
                     Float r = result.pupilRadii[i];
                      smallestBound = std::max(r, 0.0f);
                      result.exitpupilIndex   = i;
 
                   }
                    
                           
                }
                  return result;
            };

       
            // TG: so after this loop exitpupilBounds = Point2f  (0, pupilRadii[pupilIndex])


        
            // Loop over all wavelengths
            int wlIndex=0;
                auto polynomials = j["polynomials"];
                polyWavelengths_nm = std::vector<Float>(polynomials.size());
                  vignettingTerms = std::vector<RTFCamera::RTFVignettingTerms>(polynomials.size());

                polynomialMaps = std::vector<std::map<std::string, RTFCamera::LensPolynomialTerm>>(polynomials.size()); 
                if (polynomials.is_array() && polynomials.size() > 0) {
                    for (auto sp : polynomials) {
                    Float wavelength = (Float) sp["wavelength_nm"];   
                    polyWavelengths_nm[wlIndex]= (Float) wavelength;                 
                    

                    // Read vignetting terms
                    vignettingTerms[wlIndex] = toVignetTerms(sp);

                    // TODO: Determine exix pupil
                    // pupilIndx
                    // ExitpupilBounds
                     // Update exitPupil. Note pMin here is not topleft point but the radius.
          


                    // Read Polynomial Terms
                    auto jpoly = sp["poly"];
                    if (jpoly.is_array() && jpoly.size() > 0) {
                            for (auto jp : jpoly) {
                                auto curname = jp["outputname"].get<std::string>();
                                
                                polynomialMaps[wlIndex][curname] = toLensPolynomialTerms(jp);
                            }
                    } else {
                            Error("Error, invalid polynoial specification \"%s\".",
                                lensFile.c_str());
                            return nullptr;
                    }
                wlIndex=wlIndex+1;
                }   
            } else {
                Error("Error, invalid polynoial specification \"%s\".",
                    lensFile.c_str());
                return nullptr;
            }
    
            
            
        } else {
            Error("Error reading lens specification file \"%s\".",
                lensFile.c_str());
            return nullptr;
        }
    } else if (endsWith(lensFile, ".py")) {
        // Might use a neural network as equivalent lens model in the future
    } else {
        Error("Invalid format for lens specification file \"%s\".",
            lensFile.c_str());
        return nullptr;
    }
    
    std::string bbmode = params.FindOneFilename("bbmode", "polynomial");
    
    return new RTFCamera(cam2world, shutteropen, shutterclose,
                              apertureDiameter, filmDistance, lensThickness, planeOffset, caFlag, film, medium, polynomialMaps, bbmode, vignettingTerms);
}
 //   return new RTFCamera(cam2world, shutteropen, shutterclose,
                              //apertureDiameter, filmDistance, lensThickness, planeOffset, caFlag, film, medium, exitPupilBounds, poly, bbmode, pupilPos, pupilRadii, pupilIndex,circleRadii,circleSensitivities,circlePlaneZ);
//}

}  // namespace pbrt
