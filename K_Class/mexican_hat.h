#ifndef MEXICAN_HAT_CWT_H
#define MEXICAN_HAT_CWT_H

#include <seiscomp/core/typedarray.h>
#include <vector>
#include <cmath>
#include <algorithm>

namespace Wavelet {

std::vector<double> MexicanHat(const Seiscomp::DoubleArray &data, size_t i1, size_t i2, double scale);
size_t findFirstProminentPeak(const std::vector<double> &cwtData, double prominenceThreshold);
size_t detectPWavePeak(const Seiscomp::DoubleArray &data, size_t i1, size_t i2, double noiseLevel, double scale);

}

#endif