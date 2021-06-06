
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


// media/homogeneous.cpp*
#include "media/uber.h"
#include "sampler.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "core/interpolation.h"

namespace pbrt {
 
UberPhaseFunction::UberPhaseFunction(std::vector<float> &data) : sigma_s(0.0)
{
    int numWavelenghts = data[0];
    //assert(data[1] == nAngularSamples, "VSF angular sampling mismatch");
    
    float *waves = (float *) malloc(sizeof(float) * numWavelenghts);
    float *vals = (float *) malloc(sizeof(float) * numWavelenghts);
    
    for (int j=0; j<nAngularSamples; j++)
    {
        phaseAngle[j] = (float) j / (float) (nAngularSamples - 1) * Pi;
        
        
        for (int k=0; k<numWavelenghts; k++)
        {
            waves[k] = data[j*2*numWavelenghts + 2*k + 2];
            vals[k] = data[j*2*numWavelenghts + 2*k + 3];
        }
        Spectrum cdf = SampledSpectrum::FromSampled(waves, vals, numWavelenghts);
        Spectrum pdf = SampledSpectrum::FromSampled(waves, vals, numWavelenghts);
        
        CDF.push_back(cdf);
        PDF.push_back(pdf);
        
        sigma_s += 2 * Pi * cdf * sin(Pi - phaseAngle[j]) * phaseAngle[1];
    }
    
    Spectrum normFact = Spectrum(sigma_s);
    for (int i=0; i<nSpectralSamples; i++) normFact[i] = normFact[i] == 0 ? 1.0f : normFact[i];
    
    // Normalize
    CDF[0] = CDF[0] * Pi * 2 * sin(Pi - phaseAngle[0]) * phaseAngle[1] / normFact;
    PDF[0] = PDF[0] / normFact;
    for (int j=1; j<nAngularSamples; j++)
    {
        CDF[j] = CDF[j-1] + CDF[j] * Pi * 2 * sin(Pi - phaseAngle[j]) * phaseAngle[1] / normFact;
        PDF[j] = PDF[j] / normFact * Pi * 2 * sin(Pi - phaseAngle[j]) * phaseAngle[1];
    }
    
    
    
    free(waves);
    free(vals);
}

// PhaseFunction Interface
UberPhaseFunction::~UberPhaseFunction() {
};

Spectrum UberPhaseFunction::getScatter()
{
    return sigma_s;
}

Spectrum UberPhaseFunction::p(const Vector3f &wo, const Vector3f &wi) const {

    ProfilePhase _(Prof::PhaseFuncEvaluation);
    
    Float dot = Dot(wo,wi);
    Float angle = acos(fmin(fmax(dot, -1.0), 1.0));
    
    int angleIndex = (int) (fmin(angle / Pi, 0.9999) * (nAngularSamples - 1));
    
    float angleFraction = angle - phaseAngle[angleIndex];
    float deltaAngle = phaseAngle[angleIndex+1] - phaseAngle[angleIndex];
    float weight = angleFraction / deltaAngle;
    
    return (PDF[angleIndex + 1] * weight + PDF[angleIndex] * (1-weight));
};

Float UberPhaseFunction::Sample_p(const Vector3f &wo, Vector3f *wi,
                       const Point2f &u, Float wavelength) const {
    
    ProfilePhase _(Prof::PhaseFuncSampling);
        
    Float cdf[nAngularSamples];
    for (int i=0; i<nAngularSamples; i++)
    {
        cdf[i] = CDF[i].GetValueAtWavelength(wavelength);
    }
    
    int i=1;
    while ((cdf[i] < u[0]) && (i < (nAngularSamples-1))) i++;
    Float w1 = (cdf[i] - u[0]) / (cdf[i] - cdf[i-1]);
    Float w2 = (u[0] - cdf[i-1]) / (cdf[i] - cdf[i-1]);
    Float sampledAngle = phaseAngle[i-1] * w1 + phaseAngle[i] * w2;
    
    Float cosTheta = cos(sampledAngle);
    Float sinTheta = sin(sampledAngle);
    
    Float phi = 2 * Pi * u[1];
    Vector3f v1, v2;
    CoordinateSystem(wo, &v1, &v2);
    *wi = SphericalDirection(sinTheta, cosTheta, phi, v1, v2, -wo);
    
    return p(-wo, *wi).GetValueAtWavelength(wavelength);
};


UberMedium* createUberMedium(std::string absFile, std::string vsfFile) {
    
    std::string filename = AbsolutePath(ResolveFilename(absFile));
    std::vector<Float> values;
    pbrt::ReadFloatFile(filename.c_str(), &values);
    
    float *lambda = (float *) malloc(sizeof(float) * values.size() / 2);
    float *vals   = (float *) malloc(sizeof(float) * values.size() / 2);
    for (int i=0; i< values.size() / 2; i++)
    {
        lambda[i] = values[2*i];
        vals[i] = values[2*i+1];
    }
    Spectrum sigma_a = SampledSpectrum::FromSampled(lambda, vals, values.size() / 2);
    
    free(lambda);
    free(vals);
    
    
    filename = AbsolutePath(ResolveFilename(vsfFile));
    values.clear();
    pbrt::ReadFloatFile(filename.c_str(), &values);
    
    UberPhaseFunction * ph = new UberPhaseFunction(values);
        
    return new UberMedium(sigma_a, ph->getScatter(), ph);
}

// HomogeneousMedium Method Definitions
Spectrum UberMedium::Tr(const Ray &ray, Sampler &sampler) const {
    ProfilePhase _(Prof::MediumTr);
    return Exp(-sigma_t * std::min(ray.tMax, MaxFloat) * ray.d.Length());
}

Spectrum UberMedium::Sample(const Ray &ray, Sampler &sampler,
                                   MemoryArena &arena,
                                   MediumInteraction *mi) const {
    ProfilePhase _(Prof::MediumSample);
    // Sample a channel and distance along the ray
    float sigma_t_wave = sigma_t.GetValueAtWavelength(ray.wavelength);
    
    Float dist = -std::log(1 - sampler.Get1D()) / sigma_t_wave;
    
    Float t = ray.tMax;
    bool sampledMedium = false;
    if (sigma_s.GetValueAtWavelength(ray.wavelength) > 0)
    {
        // Scattering is enabled
        t = std::min(dist / ray.d.Length(), ray.tMax);
        sampledMedium = t < ray.tMax;
    }
        
    if (sampledMedium)
    {
        *mi = MediumInteraction(ray(t), -ray.d, ray.time, this, this->ph);
    }
        
    // Compute the transmittance and sampling density
    Spectrum Tr = Exp(-sigma_t * std::min(t, MaxFloat) * ray.d.Length());
    
    Spectrum density = sampledMedium ? (sigma_t * Tr) : Tr;
    Float pdf = density.GetValueAtWavelength(ray.wavelength);    
    if (pdf == 0) {
        CHECK(Tr.IsBlack());
        pdf = 1;
    }

    return sampledMedium ? (Tr * sigma_s / pdf) : (Tr / pdf);
}

}  // namespace pbrt
