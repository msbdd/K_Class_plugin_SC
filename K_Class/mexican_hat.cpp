#include "mexican_hat.h"

using namespace Seiscomp;
namespace Wavelet
{
// inline for possible speed boost; Mexican Hat Wavelet function
inline double
mexicanHat (double t, double scale)
{
	double s = scale;
	double ts = t / s;
	return (1 - ts * ts) * exp (-0.5 * ts * ts);
}

// Wavelet Transform using the Mexican Hat wavelet
std::vector<double>
MexicanHat (const DoubleArray &data, size_t i1, size_t i2, double scale)
{
	size_t dataSize = data.size ();
	if (i1 >= dataSize || i2 >= dataSize || i1 > i2)
	{
		throw std::out_of_range ("Invalid range specified for i1 and i2");
	}

	size_t n = i2 - i1 + 1;
	std::vector<double> cwt (n, 0.0);

	// Define half-window based on scale
	// The half window will be used to calulate the wavelet only in this window

	size_t halfWindow = static_cast<size_t> (std::ceil (3.0 * scale));

	for (size_t i = i1; i <= i2; ++i)
	{
		double sum = 0.0;

		// Convolve signal with the scaled wavelet
		for (size_t k = std::max (i1, i - halfWindow); k <= std::min (i2, i + halfWindow); ++k)
		{
			double t = static_cast<double> (k) - i;
			sum += data[k] * mexicanHat (t, scale);
		}

		cwt[i - i1] = sum / std::sqrt (scale); // Normalize by sqrt(scale)
	}

	return cwt;
}
// Find the first prominent peak
size_t
findFirstProminentPeak (
	const std::vector<double> &cwtData, double prominenceThreshold)
{
	size_t peakIndex = 0;
	double lastValue = cwtData[0];
	bool increasing = false;

	for (size_t i = 1; i < cwtData.size (); ++i)
	{
		if (cwtData[i] > lastValue)
		{
			increasing = true;
		}
		else if (increasing && (lastValue - cwtData[i] >= prominenceThreshold))
		{
			peakIndex = i - 1;
			break;
		}
		lastValue = cwtData[i];
	}
	return peakIndex;
}

// Main function to detect the P-wave peak
size_t
detectPWavePeak (
	const DoubleArray &data, size_t i1, size_t i2, double noiseLevel, double scale)
{
	double prominenceThreshold = 1.5 * noiseLevel;

	std::vector<double> cwtData = MexicanHat (data, i1, i2, scale);
	size_t peakIndex = findFirstProminentPeak (cwtData, prominenceThreshold);

	return i1 + peakIndex;
}
} // namespace Wavelet
