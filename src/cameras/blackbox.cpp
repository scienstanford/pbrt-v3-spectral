
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
#include "cameras/blackbox.h"
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

// BlackBoxCamera Method Definitions
BlackBoxCamera::BlackBoxCamera(const AnimatedTransform &CameraToWorld, Float shutterOpen,
    Float shutterClose, Float apertureDiameter, Float filmdistance, Float lensthickness, bool caFlag, Film *film, const Medium *medium, Bounds2f exitPupilBounds, std::map<std::string,BlackBoxCamera::LensPolynomialTerm> poly, std::string bbmode)
    : Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
      caFlag(caFlag), filmDistance(filmdistance), lensThickness(lensthickness), exitPupilBounds(exitPupilBounds),poly(poly), bbmode(bbmode) {
}
Vector2f BlackBoxCamera::Pos2RadiusRotation(const Point3f pos) const{
    // Convert position to radius
    Float radius = std::sqrt(pos.x * pos.x + pos.y + pos.y);
    Float deg = std::atan(pos.y/pos.x);
    
    // Determine whether to rotate 180 degrees
    if ( (pos.x < 0 && pos.y < 0) || (pos.x < 0 && pos.y > 0)) {
        deg += Pi/2;
    }
    
    Vector2f res = Vector2f(radius, deg);
    
    return res;
}

Float BlackBoxCamera::PolynomialCal(const Ray &thisRay, const std::string name) const{
    LensPolynomialTerm thisOut = poly[name];
    return 1.0;
}

Ray BlackBoxCamera::ApplyPolynomial(const Ray &thisRay) const{
    // Get radius and rotation
    Point3f pos= thisRay.o;
    Vector2f radiusRotation = Pos2RadiusRotation(pos);
    
    // Main body:
    // Use 'poly' parameter - a map from name to the terms.
    // Optional keys are: 'outx', 'outy', 'outu', 'outv', 'outw'
    // Term contains: termr, termu, termv, coeff (input ray parameters)
    Float x = PolynomialCal(thisRay, "outx");
    Float y = PolynomialCal(thisRay, "outy");
    Float u = PolynomialCal(thisRay, "outu");
    Float v = PolynomialCal(thisRay, "outv");
    Float w = PolynomialCal(thisRay, "outw");
    
    Ray rOut = Ray(Point3f(x, y, -lensThickness), Vector3f(u, v, w));
    
    // Rotate the output ray
    rOut = RotateRays(rOut, radiusRotation.y);
    
    return rOut;
}

Ray BlackBoxCamera::RotateRays(const Ray &thisRay, const Float deg) const{
    Transform rot = Rotate(deg, Vector3f(0, 0, 1));
    Ray rRot = rot(thisRay);
    return rRot;
}

bool BlackBoxCamera::TraceLensesFromFilm(const Ray &rCamera,
    Ray *rOut,
    const Transform CameraToLens = Scale(1, 1, -1)) const {
    
    Ray rLens = CameraToLens(rCamera);
        
    // Apply polynomial fitting
    *rOut = ApplyPolynomial(rLens);

    // Rotate rays back to camera space
    if (rOut != nullptr) {
        *rOut = Inverse(CameraToLens)(rLens);
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
float BlackBoxCamera::TToBackLens(const Ray &rCamera,
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
bool BlackBoxCamera::TraceLensesFromScene(const Ray &rCamera,
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
bool BlackBoxCamera::TraceFullLensSystemFromFilm(const Ray& rIn, Ray* rOut) const {
    
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


Point3f BlackBoxCamera::SampleExitPupil(const Point2f &lensSample) const {
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
    Point2f pLens = exitPupilBounds.Lerp(lensSample);
    return Point3f(pLens.x, pLens.y, filmDistance);
}

Float BlackBoxCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    
    ProfilePhase prof(Prof::GenerateCameraRay);
    // ++totalRays;
    // Find point on film, _pFilm_, corresponding to _sample.pFilm_
    Point2f s(sample.pFilm.x / film->fullResolution.x,
        sample.pFilm.y / film->fullResolution.y);
    Point2f pFilm2 = film->GetPhysicalExtent().Lerp(s);
    Point3f pFilm(-pFilm2.x, pFilm2.y, 0);

    Point3f pRear;
    pRear = SampleExitPupil(sample.pLens);

    Ray rFilm = Ray(pFilm, pRear - pFilm, Infinity,
        Lerp(sample.time, shutterOpen, shutterClose));
    Ray rRear = Ray(pFilm, pRear - pFilm, Infinity,
                    Lerp(sample.time, shutterOpen, shutterClose));
    
    // This is the core function of blackbox camera model
    if(!TraceLensesFromFilm(rRear, ray)) {
            // ++vignettedRays;
            return 0;
    }
    
    // Finish initialization of _BlackBoxCamera_ ray
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

    // Finish initialization of _BlackBoxCamera_ ray
    *ray = CameraToWorld(*ray);
    ray->d = Normalize(ray->d);
    ray->medium = medium;

    // Return weighting for _BlackBoxCamera_ ray
    if (HasMicrolens()) {
        // TODO: Proper weighting
        Float cosTheta = Normalize(rFilm.d).z;
        Float cos4Theta = (cosTheta * cosTheta) * (cosTheta * cosTheta);
        if (simpleWeighting) {
            // Normalize by expecting that on average rays only make it through
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
    return cos4Theta;
}

BlackBoxCamera *CreateBlackBoxCamera(const ParamSet &params,
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
    Float lensThickness = params.FindOneFloat("lensthickness", 0.0021224);
    // Chromatic aberration flag
    bool caFlag = params.FindOneBool("chromaticAberrationEnabled", false);
    
    // Create bound
    Point2f topLeft = params.FindOnePoint2f("topleftplane", Point2f(0, 0));
    Point2f bottomRight = params.FindOnePoint2f("bottomrightplane", Point2f(0, 0));
    Bounds2f exitPupilBounds(topLeft, bottomRight);
    /*
    return new BlackBoxCamera(cam2world, shutteropen, shutterclose,
                               apertureDiameter, filmDistance, focusDistance, simpleWeighting, noWeighting, caFlag,
                               lensInterfaceData, microlensData, microlensDims, microlensOffsets, microlensSensorOffset, microlensSimulationRadius, film, medium);
     */
    
    // Read in lens file
    std::string lensFile = params.FindOneFilename("lensfile", "");
    if (lensFile == "") {
        Error("No lens description file supplied!");
        return nullptr;
    }
    
    std::map<std::string, BlackBoxCamera::LensPolynomialTerm> poly;

    // A pretty cool function
    auto endsWith = [](const std::string& str, const std::string& suffix) {
        return str.size() >= suffix.size() && 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
    };
    
    if (endsWith(lensFile, ".json")) {
        /*
         The format of polynomial lens would be:
         {
          "description": "equivalent lens poly",
          "name": "polynomial",
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
            }
            if (j["description"].is_string()) {
                LOG(INFO) << StringPrintf("%s\n",
                    j["description"].get<std::string>().c_str());
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
            auto toLensPolynomialTerms = [toTerms] (json jp) {
                BlackBoxCamera::LensPolynomialTerm result;
                result.name = jp["outputname"].get<std::string>();
                result.termr = toTerms(jp["termr"]);
                result.termu = toTerms(jp["termu"]);
                result.termv = toTerms(jp["termv"]);
                result.coeff = toTerms(jp["coeff"]);
                return result;
            };
            
            // Parse polynomial term expressions
            auto jpoly = j["poly"];
            if (jpoly.is_array() && jpoly.size() > 0) {
                for (auto jp : jpoly) {
                    auto curname = jp["outputname"].get<std::string>();
                    poly[curname] = toLensPolynomialTerms(jp);
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
    
    return new BlackBoxCamera(cam2world, shutteropen, shutterclose,
                              apertureDiameter, filmDistance, lensThickness, caFlag, film, medium, exitPupilBounds, poly, bbmode);
}

}  // namespace pbrt
