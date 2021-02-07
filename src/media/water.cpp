
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
#include "media/water.h"
#include "sampler.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "core/interpolation.h"

namespace pbrt {

static const int absSamples = 61;
const Float absWave[] = {400, 405, 410, 415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570, 575, 580, 585, 590, 595, 600, 605, 610, 615, 620, 625, 630, 635, 640, 645, 650, 655, 660, 665, 670, 675, 680, 685, 690, 695, 700};
const Float waterAbs[] = {0.035080, 0.033118, 0.030408, 0.028562, 0.026118, 0.024902, 0.023099, 0.021404, 0.019910, 0.018851, 0.017619, 0.017859, 0.018095, 0.018295, 0.018501, 0.018991, 0.019880, 0.020770, 0.021810, 0.023542, 0.025761, 0.029207, 0.032930, 0.037112, 0.040245, 0.042098, 0.044156, 0.046693, 0.049525, 0.052769, 0.056292, 0.061013, 0.065429, 0.070765, 0.076831, 0.085858, 0.100352, 0.121569, 0.148864, 0.180922, 0.221222, 0.243105, 0.257202, 0.267508, 0.277863, 0.285397, 0.292787, 0.299453, 0.306261, 0.314608, 0.325244, 0.348966, 0.374212, 0.393069, 0.407539, 0.422468, 0.441646, 0.470825, 0.505272, 0.557488, 0.617855};
const Float planktonAbs[] = {0.015500, 0.016200, 0.016900, 0.016950, 0.017000, 0.017400, 0.017800, 0.018100, 0.018400, 0.018100, 0.017800, 0.017950, 0.018100, 0.017600, 0.017100, 0.015850, 0.014600, 0.013850, 0.013100, 0.012600, 0.012100, 0.011450, 0.010800, 0.010250, 0.009700, 0.009250, 0.008800, 0.008300, 0.007800, 0.007100, 0.006400, 0.005800, 0.005200, 0.004900, 0.004600, 0.004700, 0.004800, 0.004850, 0.004900, 0.004500, 0.004100, 0.004150, 0.004200, 0.004550, 0.004900, 0.005400, 0.005900, 0.006000, 0.006100, 0.005750, 0.005400, 0.006500, 0.007600, 0.009500, 0.011400, 0.011250, 0.011100, 0.008650, 0.006200, 0.003900, 0.001600};

static const int numKopelevichAngles = 19;
const Float kopelevichAngles[] = {0.000000, 0.008727, 0.017453, 0.026180, 0.034907, 0.069813, 0.104720, 0.174533, 0.261799, 0.523599, 0.785398, 1.047198, 1.308997, 1.570796, 1.832596, 2.094395, 2.356194, 2.617994, 3.141593};
const Float kopelevichSmallParticles[] = {5.300000, 5.300000, 5.200000, 5.200000, 5.100000, 4.600000, 3.900000, 2.500000, 1.300000, 0.290000, 0.098000, 0.041000, 0.020000, 0.012000, 0.008600, 0.007400, 0.007400, 0.007500, 0.008100};
const Float kopelevichLargeParticles[] = {140.000000, 98.000000, 46.000000, 26.000000, 15.000000, 3.600000, 1.100000, 0.200000, 0.050000, 0.002800, 0.000620, 0.000380, 0.000200, 0.000063, 0.000044, 0.000029, 0.000020, 0.000020, 0.000070};

KopelevichPhaseFunction::KopelevichPhaseFunction(std::vector<float> &data) : cSmall(0), cLarge(0), sigma_s(0.0)
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
        PDF[j] = PDF[j] / normFact;
    }
    
    free(waves);
    free(vals);
    
    
}

KopelevichPhaseFunction::KopelevichPhaseFunction(Float cSmall, Float cLarge) :
    cSmall(cSmall), cLarge(cLarge), sigma_s(0.0) {
            
    Spectrum waterWaveWght;
    Spectrum smallWaveWght;
    Spectrum largeWaveWght;
        
    for (int j=0; j<Spectrum::nSamples; ++j)
    {
        // Compute average value of given SPD over $i$th sample's range
        Float lambda0 = Lerp(Float(j) / Float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
        Float lambda1 = Lerp(Float(j + 1) / Float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
        Float lambda = (lambda1 + lambda0) / 2;
        
        waterWaveWght[j] = pow(550 / lambda, 4.32);
        smallWaveWght[j] = pow(550 / lambda, 1.7);
        largeWaveWght[j] = pow(550 / lambda, 0.3);
    }
        
    for (int j=0; j<nAngularSamples; j++)
    {
        phaseAngle[j] = (float) j / (float) (nAngularSamples-1) * Pi;
        
        Spectrum cdf = 0.000093 * (1 + 0.835 * cos(phaseAngle[j]) * cos(phaseAngle[j])) * waterWaveWght;
        cdf += cSmall * CatmullRom(numKopelevichAngles, kopelevichAngles, kopelevichSmallParticles, phaseAngle[j]) * smallWaveWght;
        cdf += cLarge * CatmullRom(numKopelevichAngles, kopelevichAngles, kopelevichLargeParticles, phaseAngle[j]) * largeWaveWght;
        
        CDF.push_back(cdf);
        sigma_s += 2 * Pi * cdf * sin(phaseAngle[j]) * phaseAngle[1];
    }
        
    Spectrum normFact = Spectrum(sigma_s);
    for (int i=0; i<nSpectralSamples; i++) normFact[i] = normFact[i] == 0 ? 1.0f : normFact[i];
    
    // Normalize
    CDF[0] = CDF[0] * Pi * 2 * sin(phaseAngle[0]) * phaseAngle[1] / normFact;
    for (int j=1; j<nAngularSamples; j++)
    {
        CDF[j] = CDF[j-1] + CDF[j] * Pi * 2 * sin(phaseAngle[j]) * phaseAngle[1] / normFact;
    }
        
}

// PhaseFunction Interface
KopelevichPhaseFunction::~KopelevichPhaseFunction() {
};

Spectrum KopelevichPhaseFunction::getScatter()
{
    return sigma_s;
}

Spectrum KopelevichPhaseFunction::p(const Vector3f &wo, const Vector3f &wi) const {

    ProfilePhase _(Prof::PhaseFuncEvaluation);
    
    Float dot = Dot(wo,wi);
    Float angle = acos(fmin(fmax(dot, -1.0), 1.0)) / Pi;
    int angleInteger = (int) (fmin(angle, 0.9999) * (nAngularSamples - 1));
    float angleFloat = fmin(fmax(fmin(angle, 0.9999) * (nAngularSamples - 1) - angleInteger, 0.0f), 1.0f);
    return (PDF[angleInteger + 1] * angleFloat + PDF[angleInteger] * (1-angleFloat));
};

Float KopelevichPhaseFunction::Sample_p(const Vector3f &wo, Vector3f *wi,
                       const Point2f &u, Float wavelength) const {
    
    ProfilePhase _(Prof::PhaseFuncSampling);
        
    Float cdf[nAngularSamples];
    for (int i=0; i<nAngularSamples; i++)
        cdf[i] = CDF[i].GetValueAtWavelength(wavelength);
    
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
    
    return p(wo, *wi).GetValueAtWavelength(wavelength);
};

WaterMedium* createWaterMedium(Float cPlankton, Float aCDOM440, Float aNAP400, Float cSmall, Float cLarge) {
    
    // The total absorption spectrum = water absorpton + plankton absorption + CDOM absorption
    // + NAP absorption.
    Spectrum sigma_a = SampledSpectrum::FromSampled(absWave, waterAbs, absSamples);
    sigma_a += SampledSpectrum::FromSampled(absWave, planktonAbs, absSamples) * cPlankton; //aPlankton in mg/m3
    
    Float tmp[absSamples];
    for (int i=0; i<absSamples; i++)
        tmp[i] = aCDOM440 * exp(-0.014 * (absWave[i] - 440)) + aNAP400 * exp(-0.011 * (absWave[i] - 400));
    
    sigma_a += SampledSpectrum::FromSampled(absWave, tmp, absSamples);
    
    KopelevichPhaseFunction * ph = new KopelevichPhaseFunction(cSmall, cLarge);
        
    return new WaterMedium(sigma_a, ph->getScatter(), ph);
}

WaterMedium* createWaterMedium(std::string absFile, std::string vsfFile) {
    
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
    
    KopelevichPhaseFunction * ph = new KopelevichPhaseFunction(values);
        
    return new WaterMedium(sigma_a, ph->getScatter(), ph);
}

// HomogeneousMedium Method Definitions
Spectrum WaterMedium::Tr(const Ray &ray, Sampler &sampler) const {
    ProfilePhase _(Prof::MediumTr);
    return Exp(-sigma_t * std::min(ray.tMax, MaxFloat) * ray.d.Length());
}

Spectrum WaterMedium::Sample(const Ray &ray, Sampler &sampler,
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
    
    return sampledMedium ? (Tr * sigma_s) : (Tr);
}

}  // namespace pbrt
