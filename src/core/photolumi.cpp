//
//  photolumi.cpp
//  bsdftest
//
//  Created by ZhengLyu on 9/30/19.
//

#include "photolumi.h"
#include <algorithm>

namespace pbrt {

Float AveragePhotoLumiSamples(const Float *lambda, Float **vals, int n,
                              Float lambdaStart, Float lambdaEnd, int row) {
    if (row >= n) return 0.f;
    
    for (int i = 0; i < n - 1; ++i) CHECK_GT(lambda[i + 1], lambda[i]);
    CHECK_LT(lambdaStart, lambdaEnd);
    if (lambdaEnd <= lambda[0]) return vals[row][0];
    if (lambdaStart >= lambda[n - 1]) return vals[row][n - 1];
    if (n == 1) return vals[row][0];
    Float sum = 0;
    // Add contributions of constant segments before/after samples
    if (lambdaStart < lambda[0]) sum += vals[row][0] * (lambda[0] - lambdaStart);
    if (lambdaEnd > lambda[n - 1])
        sum += vals[row][n - 1] * (lambdaEnd - lambda[n - 1]);
    // Advance to first relevant wavelength segment
    int i = 0;
    while (lambdaStart > lambda[i + 1]) ++i;
    CHECK_LT(i + 1, n);
    
    // Loop over wavelength sample segments and add contributions
    auto interp = [lambda, vals, row](Float w, int i) {
        return Lerp((w - lambda[i]) / (lambda[i + 1] - lambda[i]), vals[row][i], vals[row][i + 1]);
    };
    for (; i + 1 < n && lambdaEnd >= lambda[i]; ++i) {
        Float segLambdaStart = std::max(lambdaStart, lambda[i]);
        Float segLambdaEnd = std::min(lambdaEnd, lambda[i + 1]);
        sum += 0.5 * (interp(segLambdaStart, i) + interp(segLambdaEnd, i)) *
                (segLambdaEnd - segLambdaStart);
    }
    return sum / (lambdaEnd - lambdaStart);
}

bool PhotoLumiSamplesSorted(const Float *lambda, Float **vals, int n) {
    for (int i = 0; i < n - 1; ++i)
        if (lambda[i] > lambda[i + 1]) return false;
    return true;
}

void SortPhotoLumiSamples(Float *lambda, Float **vals, int n) {
    std::vector<std::pair<Float, Float*>> sortVec;
    sortVec.reserve(n);
    for (int i = 0; i < n; ++i)
        sortVec.push_back(std::make_pair(lambda[i], vals[i]));
    std::sort(sortVec.begin(), sortVec.end());
    for (int i = 0; i < n; ++i) {
        lambda[i] = sortVec[i].first;
        vals[i] = sortVec[i].second;
    }
}

}
