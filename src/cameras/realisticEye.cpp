
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

// cameras/realistic.cpp*
#include "cameras/realisticEye.h"
#include "paramset.h"
#include "sampler.h"
#include "sampling.h"
#include "floatfile.h"
#include "imageio.h"
#include "reflection.h"
#include "stats.h"
#include "lowdiscrepancy.h"
#include <array>

#include <math.h>
#include <cmath>
#include <gsl/gsl_roots.h> // For solving biconic surface intersections.
#include <gsl/gsl_errno.h>

// NOTE: This code, alonge with realisticEye.h was adapted from pbrt-v2-spectral code by Trisha around Feb 2018.

namespace pbrt {
    
     //STAT_PERCENT("Camera/Rays vignetted by lens system", vignettedRays, totalRays);
    
    // -----------------------------------------
    // Needed for solving intersection of ray with biconic surface
    // -----------------------------------------
    struct biconic_params {
        float Rx, Ry, Cx, Cy;
        Ray ray;
    };
    double BiconicSag(double t, void *params){
        
        struct biconic_params *p;
        p = (struct biconic_params *)params;
        
        Point3f intersect = p->ray(t);
        float x,y,z;
        z = intersect.z;
        x = intersect.x;
        y = intersect.y;
        
        float f,g,g_term;
        f = (x*x)/p->Rx + (y*y)/p->Ry;
        g_term = 1 - (1+p->Cx)*(x*x)/(p->Rx*p->Rx) - (1+p->Cy)*(y*y)/(p->Ry*p->Ry);
        
        if(g_term < 0){
            // TODO: What to do here? Let's just make it a small number to ensure that f/g =/= z.
            g_term = 0.001;
        }
        
        g = 1 + sqrt(g_term);
        
        return z-f/g;
    }
    
    // -----------------------------------------
    // -----------------------------------------
    
    RealisticEye *CreateRealisticEye(const ParamSet &params,
                                     const AnimatedTransform &cam2world,
                                     Film *film, const Medium *medium) {
        
        // Needed for default camera class.
        Float shutteropen = params.FindOneFloat("shutteropen", 0.f);
        Float shutterclose = params.FindOneFloat("shutterclose", 1.f);
        if (shutterclose < shutteropen) {
            Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                    shutterclose, shutteropen);
            std::swap(shutterclose, shutteropen);
        }
        
        
        //The lens file should include  columns for:
        // [radiusX radiusY thickness mediumIndex semiDiameter conicConstantX conicConstantY].
        // The surfaces should be ordered from the world to the retina.
        // (E.g. The last surface in the text file is the retina.)
        std::string specfile = params.FindOneString("specfile", "");
        if (specfile == "") {
            // Try looking for "lensfile" instead of "specfile"
            specfile = params.FindOneString("lensfile", "");
            // If it's still empty, we have an issue.
            if(specfile == ""){
                Error( "No lens specification file supplied!\n" );
            }
        }
        
        // These are additional parameters we need to specify in the PBRT file.
        Float pupDiameter = params.FindOneFloat("pupilDiameter", 4.0); //mm
        Float retinaDistance = params.FindOneFloat("retinaDistance",16.32); //mm
        Float retinaRadius = params.FindOneFloat("retinaRadius",12); //mm
        Float retinaSemiDiam = params.FindOneFloat("retinaSemiDiam",4); //mm
        
        // Check for IORspectra slots
        Spectrum ior1 = params.FindOneSpectrum("ior1", 0);
        Spectrum ior2 = params.FindOneSpectrum("ior2", 0);
        Spectrum ior3 = params.FindOneSpectrum("ior3", 0);
        Spectrum ior4 = params.FindOneSpectrum("ior4", 0);
        Spectrum ior5 = params.FindOneSpectrum("ior5", 0);
        Spectrum ior6 = params.FindOneSpectrum("ior6", 0);
        
        // Put all the spectra into a single vector.
        std::vector<Spectrum> iorSpectra;
        iorSpectra.push_back(ior1);
        iorSpectra.push_back(ior2);
        iorSpectra.push_back(ior3);
        iorSpectra.push_back(ior4);
        iorSpectra.push_back(ior5);
        iorSpectra.push_back(ior6);
        
        // Flags
        bool flipRad = params.FindOneBool("flipLensRadius", 0.0);
        bool mmUnits = params.FindOneBool("mmUnits",1.0);
        bool diffractionEnabled = params.FindOneBool("diffractionEnabled", 0.0);
        
        // Weighting parameters for lens shading
        bool simpleWeighting = params.FindOneBool("simpleweighting", true);
        bool noWeighting = params.FindOneBool("noweighting", false); // Added by TL for depth maps.
        
        return new RealisticEye(cam2world,
                                shutteropen,
                                shutterclose,
                                simpleWeighting,
                                noWeighting,
                                film,
                                medium,
                                specfile,
                                pupDiameter,
                                retinaDistance,
                                retinaRadius,
                                retinaSemiDiam,
                                iorSpectra,
                                flipRad,
                                mmUnits,
                                diffractionEnabled);

    }
    
    // RealisticEye Method Definitions
    RealisticEye::RealisticEye(const AnimatedTransform &CameraToWorld,
                               Float shutterOpen,
                               Float shutterClose,
                               bool simpleWeighting,
                               bool noWeighting,
                               Film *film,
                               const Medium *medium,
                               std::string specfile,
                               Float pD,
                               Float rD,
                               Float rR,
                               Float rSD,
                               std::vector<Spectrum> iorS,
                               bool flipRad,
                               bool mmUnits,
                               bool diffEnab)
    : Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
    simpleWeighting(simpleWeighting), noWeighting(noWeighting) {
        
        // Scale units depending on the units of the scene
        if(mmUnits){
            lensScaling = 1;
        }else{
            lensScaling = 0.001;
        }
        
        pupilDiameter = pD*lensScaling;
        retinaDistance = rD*lensScaling;
        retinaRadius = rR*lensScaling;
        retinaSemiDiam = rSD*lensScaling;
        iorSpectra = iorS;
        diffractionEnabled = diffEnab;
        
        // -------------------------
        // --- Read in lens file ---
        // -------------------------
        
        /*
         Note on sign convention:
         For the realistic eye code, we have gone with the Zemax convention for the lens radius. A positive lens radius means the center of the spherical lens is in the positive direction (toward the scene) relative to the location of the lenes. This is different than the other camera classes where the convention is flipped. If necessary, one can turn on the [flipLensRadius] flag.
         */
        
        // Find the complete path for the specfile
        std::string lensFileName = AbsolutePath(ResolveFilename(specfile));
        
        std::vector<float> vals;
        
        // Check to see if there is valid input in the lens file.
        if (!ReadFloatFile(lensFileName.c_str(), &vals)) {
            Warning("Unable to read lens file!");
            return;
        }
        
        // The lens file should include  columns for [radiusX radiusY thickness mediumIndex semiDiameter conicConstantX conicConstantY].
        // Let's check then that the file is a multiple of 7 (not including the effective focal length at the top)
        if ((vals.size()-1) % 7 != 0)
        {
            Warning("Wrong number of float values in lens file! Did you forget to specify the focal length? Is this a lens file with biconic surfaces? Do you have a carriage return at the end of the data?");
            return;
        }
        
        effectiveFocalLength = vals[0]*lensScaling;   // Read the effective focal length.
        
        for (int i = 1; i < vals.size(); i+=7)
        {
            LensElementEye currentLensEl;
            currentLensEl.radiusX = vals[i]*lensScaling;
            currentLensEl.radiusY = vals[i+1]*lensScaling;
            currentLensEl.thickness = vals[i+2]*lensScaling;
            currentLensEl.mediumIndex = vals[i+3];
            currentLensEl.semiDiameter = vals[i+4]*lensScaling;
            currentLensEl.conicConstantX = vals[i+5];
            currentLensEl.conicConstantY = vals[i+6];
            
            // Note: Zemax and PBRT-spectral seem to have different conventions for what the positive and negative sign of the radius is. In Zemax, a positive radius means that the center of lens sphere is directed toward the positive Z-axis, and vice versa. In previous PBRT-spectral iterations, this was flipped. Here I've rewritten the lens tracing code to go with the Zemax convention, however to be backward compatible we might want to have this ability to flip the radii.
            if(flipRad){
                currentLensEl.radiusX = -1*currentLensEl.radiusX;
                currentLensEl.radiusY = -1*currentLensEl.radiusY;
                currentLensEl.conicConstantX = -1*currentLensEl.conicConstantX;
                currentLensEl.conicConstantY = -1*currentLensEl.conicConstantY;
                Warning("Flipping lens radius & conic convention.");
            }
            
            // A radius of zero in BOTH x and y directions indicates an aperture. We should be careful of this though, since sometimes we may want to define a flat surface...
            // If the surface is an aperture, we set it's size to be equal to the pupil diameter specified.
            if (currentLensEl.radiusX == 0 && currentLensEl.radiusY == 0 ){
                currentLensEl.semiDiameter = pupilDiameter/2;
            }
            else{
                // We have to do a semi-diameter check here. As we change accommodation, we also change the radius of curvature, and we don't want the semi-diameter to be bigger than the radius. This manifests as a square root of a negative number in equation 1 on Einighammer et al. 2009.
                // TODO: This check is sort of hack-y, is there a better mathematical way to do this?
                // Note: This calculation should be done in millimeters.
                float smallerR = std::min(currentLensEl.radiusX,currentLensEl .radiusY)*(1/lensScaling);
                float biggerK = std::max(currentLensEl.conicConstantX,currentLensEl.conicConstantY);
                if(currentLensEl.semiDiameter*(1/lensScaling)*currentLensEl.semiDiameter*(1/lensScaling)*(1+biggerK)/(smallerR*smallerR) > 1.0f ){
                    currentLensEl.semiDiameter = 0.95 * sqrt((smallerR*smallerR/(1+biggerK))); // 0.95 is to add some buffer zone, since rays act very strangely when they get too close to the edge of the conical surface.
                }
            }
            
            
            lensEls.push_back(currentLensEl);
        }
        
        // Check thickness of last element. It should be zero, since we use the "retina distance" parameter for this final "thickness."
        if(lensEls[lensEls.size()-1].thickness != 0){
            Error("Thickness of lens element closest to zero must be zero. Define thickness in 'retinaDistance' parameter instead.");
        }

        // To calculate the "film diagonal", we use the retina semi-diameter. The film diagonal is the diagonal of the rectangular image rendered out by PBRT, in real units. Since we restrict samples to a circular image, we can calculate the film diagonal to be the same as a square that circumscribes the circular image.
        retinaDiag = retinaSemiDiam*1.4142*2; // sqrt(2)*2
        
        // We are going to use our own error handling for gsl to prevent it from crashing the moment a ray doesn't intersect.
        gsl_set_error_handler_off();
        
        // Set up GSL random number generator for diffraction modeling
        const gsl_rng_type * T;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        
    }
    
    void RealisticEye::applySnellsLaw(Float n1, Float n2, Float lensRadius, Vector3f &normalVec, Ray * ray ) const
    {
        
        //                  |
        // scene ..... n2   |   n1 ..... sensor
        //                  |
        //                  |
        //                 lens
        
        
        Vector3f s1 = ray->d;
        if (lensRadius >0)
            normalVec = -normalVec;
        
        float radicand = 1 - (n1/n2) * (n1/n2) * Dot(Cross(normalVec, s1), Cross(normalVec, s1));
        if (radicand < 0){
            ray->d = Vector3f(0, 0,0);
            // Tlian: I put this return back in so the assert doesn't catch.
            return;   //reflection, no refraction - might want to change for lens flare!
        }
        
        Vector3f s2 = n1/n2 * (Cross(normalVec, Cross(-1 * normalVec, s1))) - normalVec * sqrt(radicand);
        ray->d = Normalize(s2);  //reassign the direction to the new direction
    }
    
    bool RealisticEye::IntersectLensElAspheric(const Ray &r, Float *tHit, LensElementEye currElement, Float zShift, Vector3f *n) const{
        
        // This code lets us find the intersection with an aspheric surface. At tHit, the ray will intersect the surface. Therefore:
        // If
        // (x,y,z) = ray_origin + thit * ray_direction
        // then:
        // z - u(x,y) = 0
        // where u(x,y) is the SAG of the surface, defined in Eq.1 of Einighammer et al. 2009 and in the Zemax help page under biconic surfaces.
        // We can use this fact to solve for thit. If the surface is not a sphere, this is a messy polynomial. So instead we use a numeric root-finding method (Van Wijingaarden-Dekker-Brent's Method.) This method is available in the GSL library.
        
        // DEBUG
        /*
         std::cout << "r.o = " << r.o.x << "," << r.o.y << "," << r.o.z << std::endl;
         std::cout << "r.d = " << r.d.x << "," << r.d.y << "," << r.d.z << std::endl;
         */
        
        // Move ray to object(lens) space.
        Ray objSpaceRay = r;
        objSpaceRay.o = objSpaceRay.o + Vector3f(0,0,zShift);
        
        int status;
        int iter = 0, max_iter = 100;
        const gsl_root_fsolver_type *T;
        gsl_root_fsolver *s;
        double root = 0;
        double x_lo = 0.0;
        double x_hi;
        if(currElement.thickness == 0){
            // Probably the surface closest to the retina.
            x_hi = retinaDistance*2;
        }else{
            x_hi = currElement.thickness*1.5; // thit will probably be less than this
        }
        gsl_function F;
        
        // DEBUG
        /*
         std::cout << "radiusX = " << currElement.radiusX << std::endl;
         std::cout << "radiusY = " << currElement.radiusY << std::endl;
         std::cout << "conicConstantX = " << currElement.conicConstantX << std::endl;
         std::cout << "conicConstantY = " << currElement.conicConstantY << std::endl;
         std::cout << "objSpaceRay.o = " << objSpaceRay.o.x << "," << objSpaceRay.o.y << "," << objSpaceRay.o.z << std::endl;
         std::cout << "objSpaceRay.d = " << objSpaceRay.d.x << "," << objSpaceRay.d.y << "," << objSpaceRay.d.z << std::endl;
         */
        
        struct biconic_params params = {currElement.radiusX,currElement.radiusY,currElement.conicConstantX,currElement.conicConstantY,objSpaceRay};
        F.function = &BiconicSag;
        F.params = &params;
        
        T = gsl_root_fsolver_brent;
        s = gsl_root_fsolver_alloc (T);
        
        status = gsl_root_fsolver_set (s, &F, x_lo, x_hi);
        if(status != 0){
            // Ray probably does not intersect. This might depend on the x_hi set above, i.e. if it's too small OR too large. TODO: Can we check this?
            gsl_root_fsolver_free (s);
            return false;
        }
        
        // DEBUG
        /*
         printf ("using %s method\n",
         gsl_root_fsolver_name (s));
         
         printf ("%5s [%9s, %9s] %9s %10s %9s\n",
         "iter", "lower", "upper", "root",
         "err", "err(est)");
         */
        
        do
        {
            iter++;
            gsl_root_fsolver_iterate (s);
            root = gsl_root_fsolver_root (s);
            x_lo = gsl_root_fsolver_x_lower (s);
            x_hi = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval (x_lo, x_hi,
                                             0, 0.0001);
            
            if (status == GSL_SUCCESS){
                
                // DEBUG
                /*
                 printf ("Converged:\n");
                 printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
                 iter, x_lo, x_hi,
                 root, root - root,
                 x_hi - x_lo);
                 */
                
                gsl_root_fsolver_free (s);
                
                // DEBUG
                /*
                 std::cout << "root = " << root << std::endl;
                 std::cout << "r.o = " << r.o.x << "," << r.o.y << "," << r.o.z << std::endl;
                 std::cout << "r.d = " << r.d.x << "," << r.d.y << "," << r.d.z << std::endl;
                 */
                
                *tHit = root;
                Point3f intersect = r(*tHit);
                
                // Check if intersection is within the semi-diameter of the lens, if not we return false (no intersect)
                // (If we don't do this here, we might get a complex normal which would crash the rendering.)
                if(intersect.x * intersect.x + intersect.y * intersect.y > (currElement.semiDiameter * currElement.semiDiameter)){
                    return false;
                }
                
                // Calculate normal at intersection
                // These equations are from Eq 2, 3, and 4 in Einighammer 2009
                float term1 = ((1+currElement.conicConstantX)*intersect.x*intersect.x)/(currElement.radiusX*currElement.radiusX);
                float term2 = ((1+currElement.conicConstantY)*intersect.y*intersect.y)/(currElement.radiusY*currElement.radiusY);
                
                float fprime_x = 2*intersect.x/currElement.radiusX;
                float gprime_x = (-1*(1+currElement.conicConstantX)*intersect.x)/(currElement.radiusX*currElement.radiusX*sqrt(1-term1-term2));
                
                float fprime_y = 2*intersect.y/currElement.radiusY;
                float gprime_y = (-1*(1+currElement.conicConstantY)*intersect.y)/(currElement.radiusY*currElement.radiusY*sqrt(1-term1-term2));
                
                float f = (intersect.x*intersect.x)/currElement.radiusX + (intersect.y*intersect.y)/currElement.radiusY;
                float g = 1+sqrt(1-term1-term2);
                
                float zprime_y = (fprime_y*g-gprime_y*f)/(g*g);
                float zprime_x = (fprime_x*g-gprime_x*f)/(g*g);
                
                Vector3f v_x = Vector3f(1,0,zprime_x);
                Vector3f v_y = Vector3f(0,1,zprime_y);
                
                *n = Normalize(Cross(v_x,v_y));
                *n = Faceforward(*n, -r.d);
                
                return true;
                
            }
            
            // For debugging
            /*
             printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
             iter, x_lo, x_hi,
             root, root - root,
             x_hi - x_lo);
             */
        }
        while (status == GSL_CONTINUE && iter < max_iter);
        
        gsl_root_fsolver_free (s);
        
        return false;
    }
    
    Float RealisticEye::GenerateRay(const CameraSample &sample, Ray *ray) const {
        
        
        //                         z=0         z= -filmDistance
        //                  |  |    ||           |
        //                  |  |    ||           |
        //  Scene <---------|--|----||---------> |
        //            +z    |  |    ||    -z     |
        //                  |  |    ||           |
        //                Lens Elements        Sensor
        //
        
        
        // Determine the size of the sensor in real world units (i.e. convert from pixels to millimeters).
        
        Point2i filmRes = film->fullResolution;
        float aspectRatio = (float)filmRes.x/(float)filmRes.y;
        float width = retinaDiag /sqrt((1.f + 1.f/(aspectRatio * aspectRatio)));
        float height = width/aspectRatio;
        
        Point3f startingPoint;
        
        startingPoint.x = -((sample.pFilm.x) - filmRes.x/2.f - .25)/(filmRes.y/2.f);
        startingPoint.y = ((sample.pFilm.y) - filmRes.y/2.f - .25)/(filmRes.y/2.f);
        
        // Convert starting point units to millimeters
        startingPoint.x = startingPoint.x * width/2.f;
        startingPoint.y = startingPoint.y * height/2.f;
        startingPoint.z = -retinaDistance;
        
        
        // Project sampled points onto the curved retina
        if (retinaRadius != 0)
        {
            // Right now the code only lets you curve the sensor toward the scene and not the other way around. See diagram:
            /*
             
             The distance between the zero point on the z-axis (i.e. the lens element closest to the sensor) and the dotted line will be equal to the "retinaDistance." The retina curvature is defined by the "retinaRadius" and it's height in the y and x direction is defined by the "retinaSemiDiam."
                      
         
                                        :
                                     |  :
                                      | :
                         | | |         |:
           scene <------ | | | <----   |:
                         | | |         |:
                     Lens System      | :
                                     |  :
                                        :
                                  retina
            <---- +z
         
             */
            
            // Limit sample points to a circle within the retina semi-diameter
            if((startingPoint.x*startingPoint.x + startingPoint.y*startingPoint.y) > (retinaSemiDiam*retinaSemiDiam)){
                return 0.f;
            }
            
            // Calculate the distance of a disc that fits inside the curvature of the retina.
            float zDiscDistance = -1*sqrt(retinaRadius*retinaRadius-retinaSemiDiam*retinaSemiDiam);
            
            // If we are within this radius, project each point out onto a sphere. There may be some issues here with even sampling, since this is a direct projection...
            float el = atan(startingPoint.x/zDiscDistance);
            float az = atan(startingPoint.y/zDiscDistance);
            
            // Convert spherical coordinates to cartesian coordinates (note: we switch up the x,y,z axis to match our conventions)
            float xc,yc,zc, rcoselev;
            xc = -1*retinaRadius*sin(el); // TODO: Confirm this flip?
            rcoselev = retinaRadius*cos(el);
            zc = -1*(rcoselev*cos(az)); // The -1 is to account for the curvature described above in the diagram
            yc = -1*rcoselev*sin(az); // TODO: Confirm this flip?
            
            zc = zc + -1*retinaDistance + retinaRadius; // Move the z coordinate out to correct retina distance
            
            startingPoint = Point3f(xc,yc,zc);
            
        }
        
        float lensU, lensV;
        Point2f lens = ConcentricSampleDisk(sample.pLens);
        lensU = lens.x;
        lensV = lens.y;
        
        //We need to shoot rays toward the disc that fits inside the curvature of first lens surface. Since we no longer have spherical elements, we have to calculate it as follows:
        // TODO: This is just a guess. Is there a correct way to do this? Should we look at the exitPupil code in PBRTv3?
        
        float lensP_semiDiam = lensEls[lensEls.size()-1].semiDiameter;
        float lensP_radius = lensEls[lensEls.size()-1].radiusX;
        
        // sgn(lensP_radius)
        // It's very rare for the first lens element to have a negative radius (spherical center toward sensor)...but just in case:
        float sgn_radius = (lensP_radius > 0) - (lensP_radius < 0);
        float discDistance = sgn_radius * BiconicZ(lensP_semiDiam, 0, lensEls[lensEls.size()-1]);
        
        // Scale the normalized lens coordinates by the size of the first lens element
        lensU *= lensP_semiDiam;
        lensV *= lensP_semiDiam;
        
        Point3f pointOnLens = Point3f(lensU, lensV, discDistance);   // We aim the ray at a flat disk, will that cause problems later?
        
        //float tempWavelength = ray->wavelength;
        //*ray = Ray(Point3f(0,0,0), Normalize(Vector3f(startingPoint)),INFINITY,0.f);
        ray->o = startingPoint;    //initialize ray origin
        ray->d = Normalize(pointOnLens - ray->o);
        //ray->wavelength = tempWavelength;  //so that wavelength information is retained.
        
        
        // DEBUG
        /*
         // Scene to retina
         startingPoint = Point(0,0.00081193431816,-16.319999973);
         pointOnLens = Point(0,1.544122632,0.21476753379);
         ray->o = startingPoint;
         ray->d = Normalize(pointOnLens - ray->o);
         ray->wavelength = 550;
         
         
         // Retina to scene
         startingPoint = Point(0,0,-16.3200);
         pointOnLens = Point(0,1.5294,0.2084);
         ray->o = startingPoint;
         ray->d = Normalize(pointOnLens - ray->o);
         ray->wavelength = 550;
         */
        
        // --------------------------------------------------------
        // --- Trace through the lens elements of the main lens ---
        // --------------------------------------------------------
        
        float lensDistance = 0; // How far we are from the "0" on the z-axis (see diagram above)
        
        for (int i = lensEls.size()-1; i>=0 ; i--)
        {
            
            ray->o = startingPoint;
            lensDistance += lensEls[i].thickness;
            
            // If the ray direction is zero, there is probably internal reflection going on somewhere. We will just terminate the ray here to avoid assert errors.
            if(ray->d == Vector3f(0,0,0)){
                return 0.f;
            }
            
            // DEBUG
            // ----
            /*
             std::cout << "\n" << std::endl;
             std::cout << "i = " << i << std::endl;
             std::cout << "start " << ray->o.x << " " << ray->o.y << " " << ray->o.z << std::endl;
             std::cout << "dir " << ray->d.x << " " << ray->d.y << " " << ray->d.z << std::endl;
             */
            // ----
            
            float tHit = 0;
            bool intersected = false;
            Vector3f normalVec(0,0,1);
            Point3f intersectPoint(0,0,0);
            
            if (lensEls[i].radiusX == 0 && lensEls[i].radiusY == 0)
            {
                // ---------------------
                // --- APERTURE CASE ---
                // ---------------------
                
                float tAperture = 0;
                if (i == lensEls.size()-1)
                    tAperture = retinaDistance/ray->d.z;   //special case for when aperture is the 1st element
                else
                    tAperture = (lensDistance - ray->o.z)/(ray->d.z);
                
                // Point where the ray intersects the aperture plane
                Point3f intersectPoint = (*ray)(tAperture);
                normalVec = Vector3f(0,0,1);
                
                // Check if ray makes it through the aperture
                if((intersectPoint.x * intersectPoint.x + intersectPoint.y * intersectPoint.y) > (lensEls[i].semiDiameter * lensEls[i].semiDiameter)){
                    return 0.f;
                }
                
                if(diffractionEnabled){
                    
                    // DEBUG: Check that direction did change
//                    std::cout << "ray->d = (" << ray->d.x << "," << ray->d.y << "," << ray->d.z << ")" << std::endl;
                    
                    // Adjust ray direction using HURB diffraction
                    Vector3f newDiffractedDir;
                    diffractHURB(intersectPoint, lensEls[i].semiDiameter,  ray->wavelength, ray->d, &newDiffractedDir);
                    ray->d = newDiffractedDir;
                    
                    // DEBUG: Check that direction did change
//                    std::cout << "ray->d = (" << ray->d.x << "," << ray->d.y << "," << ray->d.z << ")" << std::endl;
                    
                }
                
                startingPoint = intersectPoint;
                
                
                
            }
            else
            {
                // ----------------------------
                // --- REGULAR ELEMENT CASE ---
                // ----------------------------
                
                //---- Find intersection point ----
                
                // Since the surface is defined in object space with the edge at zero, so we need to move the ray accordingly from world space to object space. This is simply a translation according to the location of the surface.
                
                // z(x,y) is defined like this (flipped across the z axis for +r vs -r):
                //
                //   :   /
                //   :  /
                //   : /
                //   :|
                //   :|
                //   :|
                //   : \d_zema
                //   :  \
                //   :   \
                //  z=0
                //
                
                float zShift = -lensDistance;
                
                // Does the ray intersect the lens surface? If so, where?
                intersected = IntersectLensElAspheric(*ray, &tHit, lensEls[i], zShift, &normalVec);
                
                if (intersected)
                {
                    intersectPoint = (*ray)(tHit);
                    
                    //DEBUG
                    // std::cout << "Intersect point: " << intersectPoint.x << " " << intersectPoint.y << " " << intersectPoint.z << std::endl;
                    
                    
                    // ---- Apply Snell's Law ----
                    //                  |
                    // scene ..... n2   |   n1 ..... sensor
                    //           [i-1]  |  [i]
                    //                  |
                    //                 lens
                    //
                    
                    // The user can load IOR spectra for each ocular medium into ior1, ior2, etc. In the lens file, they can then specify with medium they would like to use. The number (X) corresponds to iorX.
                    
                    float n1,n2;
                    
                    n1 = lookUpIOR(lensEls[i].mediumIndex, *ray);
                    
                    // If we're at the lens surface closest to the scene, n2 should be air.
                    if (i-1 >= 0){
                        
                        n2 = lookUpIOR(lensEls[i-1].mediumIndex,*ray);
                        
                        // Trisha: If we're entering the aperture (sometimes we put n2 == 0 when that is the case) we skip the n2 == 0 and apply the next medium.
                        // (We put this in the current if statement so we can handle the unique case when the aperture is the first element.)
                        
                        //                 |    /
                        //                 |   |
                        //                 |  |
                        // <-----            |    <------
                        //  n2, [i-2]      |  |        n1, [i]
                        //                 |   |
                        //                 |    \
                        //
                        //           aperture, [i-1]
                        
                        if(n2 == 0)
                            n2 = lookUpIOR(lensEls[i-2].mediumIndex, *ray);
                        
                    }
                    else{
                        n2 = 1;
                    }
                    
                    // If n1 = n2 = 0, something is wrong!
                    if(n1 == 0 & n2 == 0){
                        Error("Index of refractions are set to zero. Something is wrong.");
                    }
                    applySnellsLaw( n1,  n2,  0, normalVec, ray );
                    // --- Update ray starting point ---
                    
                    startingPoint = intersectPoint;
                    
                }
                else
                {
                    return 0.f;
                }
            }
            
            
        }
        
        ray->o = startingPoint;
        
        // DEBUG
        
        // ----
        /*
         std::cout << "\n" << std::endl;
         std::cout << "final" << std::endl;
         std::cout << ray->o.x << " " << ray->o.y << " " << ray->o.z << std::endl;
         std::cout << ray->d.x << " " << ray->d.y << " " << ray->d.z << std::endl;
         */
        // ----
        
        *ray = CameraToWorld(*ray);
        ray->d = Normalize(ray->d);
        ray->medium = medium;
        
        // No weighting for now...we should add it in!
        return 1.f;
        
    }
    
    // Handy method to explicity solve for the z(x,y) at a given point (x,y),for the biconic SAG.
    float RealisticEye::BiconicZ(float x, float y, LensElementEye currElement) const{
        
        float f,g,g_term;
        float Rx,Ry,Cx,Cy;
        Rx = currElement.radiusX; Ry = currElement.radiusY;
        Cx = currElement.conicConstantX; Cy = currElement.conicConstantY;
        
        f = (x*x)/Rx + (y*y)/Ry;
        g_term = 1 - (1+Cx)*(x*x)/(Rx*Rx) - (1+Cy)*(y*y)/(Ry*Ry);
        
        if(g_term < 0){
            g_term = 0.001;
            Warning("Encountered a complex value when solving for the z(x,y) of the biconic.");
        }
        
        g = 1 + sqrtf(g_term);
        
        return f/g;
        
    }
    
    // Given the mediumIndex, load up the right spectra from the ones read in through ior1, ior2, etc. Then find the corresponding IOR for the given ray wavelength.
    float RealisticEye::lookUpIOR(int mediumIndex, const Ray &ray)const{
        
        float n;
        
            // Standard media
            // If spectral renderer is used, then ray.wavelength will be given a value. Otherwise, it will be a junk value. I'm not sure why it's not intialized to zero despite having wavelength = 0 in the constructer. Anyhow, we can error check by looking for values within the valid range.
            if((std::abs(ray.wavelength) > 400) & (std::abs(ray.wavelength) < 800)){
                iorSpectra[mediumIndex-1].GetValueAtWavelength(ray.wavelength,&n);
            }else{
                iorSpectra[mediumIndex-1].GetValueAtWavelength(550,&n);
            }
        
        return n;
    }
    
    void RealisticEye::diffractHURB(Point3f intersect, Float apertureRadius, const Float wavelength, const Vector3f oldDirection, Vector3f *newDirection) const {
        
//        std::cout << "wavelength = " << wavelength << std::endl;
        
        double dist2Int = sqrt(intersect.x*intersect.x + intersect.y*intersect.y);
        Vector3f dirS = Normalize(Vector3f(intersect.x, intersect.y, 0));
        Vector3f dirL = Normalize(Vector3f(-1*intersect.y, intersect.x, 0));
        Vector3f dirU = Vector3f(0,0,1); // Direction pointing normal to the aperture plane and toward the scene.
        
        double dist2EdgeS = apertureRadius - dist2Int;
        double dist2EdgeL = sqrt(apertureRadius*apertureRadius - dist2Int*dist2Int);
        
        // Calculate variance according to Freniere et al. 1999
        // If the scene is in meters, lensScaling = 0.001 and dist2Edge will be in meters.
        // if scene is in millimeters, lensScaling = 1 and dist2Edge will be in millimeters.
        double sigmaS = atan(1/(2 * dist2EdgeS * 2*Pi/(wavelength*10e-6*lensScaling) ));
        double sigmaL = atan(1/(2 * dist2EdgeL * 2*Pi/(wavelength*10e-6*lensScaling) ));
        
        // Sample from bivariate gaussian
        double initS = 0;
        double initL = 0;
        double *noiseS = &initS;
        double *noiseL = &initL;
        gsl_ran_bivariate_gaussian (r, sigmaS, sigmaL, 0, noiseS, noiseL);
        
        // DEBUG:
//        std::cout << "noiseS = " << *noiseS << std::endl;
//        std::cout << "noiseL = " << *noiseL << std::endl;
        
        // Decompose our original ray into dirS and dirL.
        double projS = Dot(oldDirection,dirS)/dirS.Length();
        double projL = Dot(oldDirection,dirL)/dirL.Length();
        double projU = Dot(oldDirection,dirU)/dirU.Length();
        
        /*
         We have now decomposed the original, incoming ray into three orthogonal
         directions: directionS, directionL, and directionU.
         directionS is the direction along the shortest distance to the aperture
         edge.
         directionL is the orthogonal direction to directionS in the plane of the
         aperture.
         directionU is the direction normal to the plane of the aperture, pointing
         toward the scene.
         To orient our azimuth and elevation directions, imagine that the
         S-U-plane forms the "ground plane." "Theta_x" in the Freniere paper is
         therefore the deviation in the azimuth and "Theta_y" is the deviation in
         the elevation.
         */
        
        // Calculate current azimuth and elevation angles
        double thetaA = atan(projS/projU); // Azimuth
        double thetaE = atan(projL/sqrt(projS*projS + projU*projU)); // Elevation
        
        // Deviate the angles
        thetaA = thetaA + *noiseS;
        thetaE = thetaE + *noiseL;
        
        // Recalculate the ray direction
        // Remember the ray direction is normalized, so it should have length = 1
        double newProjL = sin(thetaE);
        double newProjSU = cos(thetaE);
        double newProjS = newProjSU * sin(thetaA);
        double newProjU = newProjSU * cos(thetaA);
        
        // Add up the new projections to get a new direction
        *newDirection = Normalize(newProjS*dirS + newProjL*dirL + newProjU*dirU);
        
    }
    
}  // namespace pbrt
