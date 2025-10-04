#define SEISCOMP_COMPONENT K_Class

#include "K-Class.h"
#include "mexican_hat.h"

#include <iostream>
#include <seiscomp/logging/log.h>
#include <seiscomp/math/geo.h>
#include <seiscomp/processing/amplitudes/ML.h>
#include <seiscomp/processing/magnitudeprocessor.h>

#include <boost/bind.hpp>

// Valid within 0-15 degrees
#define DELTA_MIN 0.
#define DELTA_MAX 15
#define DEPTH_MAX 180

#define MAG_TYPE "K_Class"

using namespace std;
using namespace Seiscomp;

// Using a custom void namespace to prevent name collisions and improve code
// organization

namespace
{

// Custom function to summ two amplitudes

Processing::AmplitudeProcessor::AmplitudeValue
summ (
	const Processing::AmplitudeProcessor::AmplitudeValue &v0,
	const Processing::AmplitudeProcessor::AmplitudeValue &v1)
{
	Processing::AmplitudeProcessor::AmplitudeValue v;
	v.value = (v0.value + v1.value);

	double l0 = v0.lowerUncertainty ? *v0.lowerUncertainty : 0.0;
	double u0 = v0.upperUncertainty ? *v0.upperUncertainty : 0.0;
	double l1 = v1.lowerUncertainty ? *v1.lowerUncertainty : 0.0;
	double u1 = v1.upperUncertainty ? *v1.upperUncertainty : 0.0;

	double l = 0, u = 0;

	v.lowerUncertainty = std::sqrt (l0 * l0 + l1 * l1);
	v.upperUncertainty = std::sqrt (u0 * u0 + u1 * u1);

	v.lowerUncertainty = l;
	v.upperUncertainty = u;

	return v;
}

// Custom function to find average time between two amplitudes picks (from SED
// MLh)

Processing::AmplitudeProcessor::AmplitudeTime
average (
	const Processing::AmplitudeProcessor::AmplitudeTime &t0,
	const Processing::AmplitudeProcessor::AmplitudeTime &t1)
{
	Processing::AmplitudeProcessor::AmplitudeTime t;
	t.reference = Core::Time ((double (t0.reference) + double (t1.reference)) * 0.5);

	// Compute lower and upper uncertainty
	Core::Time t0b = t0.reference + Core::TimeSpan (t0.begin);
	Core::Time t0e = t0.reference + Core::TimeSpan (t0.end);
	Core::Time t1b = t1.reference + Core::TimeSpan (t1.begin);
	Core::Time t1e = t1.reference + Core::TimeSpan (t1.end);

	Core::Time minTime = t.reference;
	Core::Time maxTime = t.reference;

	minTime = std::min (minTime, t0b);
	minTime = std::min (minTime, t0e);
	minTime = std::min (minTime, t1b);
	minTime = std::min (minTime, t1e);

	maxTime = std::max (maxTime, t0b);
	maxTime = std::max (maxTime, t0e);
	maxTime = std::max (maxTime, t1b);
	maxTime = std::max (maxTime, t1e);

	t.begin = (double)(minTime - t.reference);
	t.end = (double)(maxTime - t.reference);

	return t;
}

ADD_SC_PLUGIN ("K_Class magnitude", "Dmitry Sidorov-Biryukov", 0, 0, 1)

// We need to create a custom non abstract class for individual magnitude
// computation

class SimpleAmplitudeProcessor
	: public Processing::AbstractAmplitudeProcessor_ML
{
  public:
	double wavelet_scale;
	SimpleAmplitudeProcessor ()
		: Processing::AbstractAmplitudeProcessor_ML (MAG_TYPE)
	{
		setMinSNR (0);
		setMinDist (DELTA_MIN);
		setMaxDist (DELTA_MAX);
		setMaxDepth (DEPTH_MAX);
	}
	friend class AmplitudeProcessor_K_Class;

  private:
	bool
	computeAmplitude (
		const DoubleArray &data, size_t i1, size_t i2, size_t si1,
		size_t si2, double offset, AmplitudeIndex *dt,
		AmplitudeValue *amplitude, double *period, double *snr) override
	{
		SEISCOMP_DEBUG ("Custom computeAmplitude called");
		double maxAmplitude = 0.0;
		size_t amp_index = 0;
		if (usedComponent () == Vertical)
		{
			double noise = 1.0;
			if (_noiseAmplitude)
			{
				noise = *_noiseAmplitude;
			}
			SEISCOMP_DEBUG ("Vertical component noise is %f", noise);
			amp_index = Wavelet::detectPWavePeak (data, i1, i2, noise, wavelet_scale);
			SEISCOMP_DEBUG (
				"i1 = %u\nP wave max detected at %u\nWindow length is %u",
				i1, amp_index, i2);
			maxAmplitude = fabs (data[amp_index]) - offset;
		}
		if (usedComponent () == FirstHorizontal || usedComponent () == SecondHorizontal)
		{
			amp_index = find_absmax (
				data.size (), data.typedData (), si1, si2, offset);
			SEISCOMP_DEBUG (
				"i1 = %u\nS wave max detected at %u\nWindow length is %u",
				si1, amp_index, si2);
			maxAmplitude = fabs (data[amp_index] - offset);
		}

		if (maxAmplitude <= 0)
		{
			SEISCOMP_DEBUG (
				"Amplitude computation failed: No valid peak found");
			return false;
		}
		amplitude->value = maxAmplitude;
		if (_streamConfig[_usedComponent].gain != 0.0)
		{
			amplitude->value /= _streamConfig[_usedComponent].gain;
			amplitude->value *= 1E03;
			amplitude->value = std::abs (amplitude->value);
			dt->index = amp_index;
			*period = -1;
			*snr = -1;

			SEISCOMP_DEBUG (
				"Amplitude computed: value=%f, period=%f, snr=%f",
				amplitude->value, *period, *snr);

			return true;
		}
		else
		{
			return false;
		}
	}
};

class AmplitudeProcessor_K_Class : public Processing::AmplitudeProcessor
{
  public:
	AmplitudeProcessor_K_Class ()
		: Processing::AmplitudeProcessor (MAG_TYPE)
	{
		setUsedComponent (Any);

		_ampN.setUsedComponent (FirstHorizontal);
		_ampE.setUsedComponent (SecondHorizontal);
		_ampZ.setUsedComponent (Vertical);

		_ampN.setPublishFunction (boost::bind (
			&AmplitudeProcessor_K_Class::newAmplitude, this, _1, _2));
		_ampE.setPublishFunction (boost::bind (
			&AmplitudeProcessor_K_Class::newAmplitude, this, _1, _2));
		_ampZ.setPublishFunction (boost::bind (
			&AmplitudeProcessor_K_Class::newAmplitude, this, _1, _2));
	}

	bool
	setup (const Processing::Settings &settings) override
	{
		_ampN._type = _type;
		_ampE._type = _type;
		_ampZ._type = _type;

		_ampN.streamConfig (FirstHorizontalComponent) = streamConfig (FirstHorizontalComponent);
		_ampE.streamConfig (SecondHorizontalComponent) = streamConfig (SecondHorizontalComponent);
		_ampZ.streamConfig (VerticalComponent) = streamConfig (VerticalComponent);

		if (!Processing::AmplitudeProcessor::setup (settings))
		{
			return false;
		}

		_ampZ.wavelet_scale = settings.getDouble ("amplitudes.K_Class.scale");

		return true;
	}

	void
	setEnvironment (
		const DataModel::Origin *hypocenter,
		const DataModel::SensorLocation *receiver,
		const DataModel::Pick *pick) override
	{
		_ampE.setEnvironment (hypocenter, receiver, pick);
		_ampN.setEnvironment (hypocenter, receiver, pick);
		_ampZ.setEnvironment (hypocenter, receiver, pick);
	}

	void
	computeTimeWindow () override
	{
		_ampN.setConfig (config ());
		_ampE.setConfig (config ());
		_ampZ.setConfig (config ());

		_ampN.computeTimeWindow ();
		_ampE.computeTimeWindow ();
		_ampZ.computeTimeWindow ();

		setConfig (_ampE.config ());
		setTimeWindow (
			_ampE.timeWindow () | _ampN.timeWindow () | _ampZ.timeWindow ());
	}

	void
	setTrigger (const Core::Time &trigger) override
	{
		AmplitudeProcessor::setTrigger (trigger);
		_ampE.setTrigger (trigger);
		_ampN.setTrigger (trigger);
		_ampZ.setTrigger (trigger);
	}

	void
	reset () override
	{
		AmplitudeProcessor::reset ();
		_results[0] = _results[1] = _results[2] = Core::None;

		_ampE.reset ();
		_ampN.reset ();
		_ampZ.reset ();
	}

	bool
	feed (const Record *record) override
	{
		if (_ampE.isFinished () && _ampN.isFinished () && _ampZ.isFinished ())
		{
			return false;
		}
		if (status () > Finished)
		{
			return false;
		}

		if (record->channelCode () == _streamConfig[FirstHorizontalComponent].code ())
		{
			_ampN.feed (record);
		}
		else if (record->channelCode () == _streamConfig[SecondHorizontalComponent].code ())
		{
			_ampE.feed (record);
		}
		else if (record->channelCode () == _streamConfig[VerticalComponent].code ())
		{
			SEISCOMP_DEBUG (
				"Adding stream %s",
				_streamConfig[VerticalComponent].code ().c_str ());
			_ampZ.feed (record);
		}
		else
		{
			SEISCOMP_WARNING (
				"Record with channel %s does not match N, E, or Z channels",
				record->channelCode ().c_str ());
			return false;
		}
		return true;
	}

  private:
	bool
	computeAmplitude (
		const DoubleArray &data, size_t i1, size_t i2, size_t si1,
		size_t si2, double offset, AmplitudeIndex *dt,
		AmplitudeValue *amplitude, double *period, double *snr) override
	{
		return false;
	}
	void
	newAmplitude (
		const AmplitudeProcessor *proc,
		const AmplitudeProcessor::Result &res)
	{
		if (isFinished ())
		{
			return;
		}

		int idx = 0;
		if (proc == &_ampE)
		{
			idx = 0;
		}
		else if (proc == &_ampN)
		{
			idx = 1;
		}
		else if (proc == &_ampZ)
		{
			idx = 2;
		}

		_results[idx] = ComponentResult ();
		_results[idx]->value = res.amplitude;
		_results[idx]->time = res.time;

		if (_results[0] && _results[1] && _results[2])
		{
			setStatus (Finished, 100.);
			Result newRes;
			newRes.record = res.record;
			newRes.component = Any;

			if (_results[0]->value.value > _results[1]->value.value)
			{
				newRes.amplitude = summ (_results[0]->value, _results[2]->value);
				newRes.time = average (_results[0]->time, _results[2]->time);
			}
			else
			{
				newRes.amplitude = summ (_results[1]->value, _results[2]->value);
				newRes.time = average (_results[1]->time, _results[2]->time);
			}
			newRes.period = -1;
			newRes.snr = -1;

			emitAmplitude (newRes);
		}
	}

	struct ComponentResult
	{
		AmplitudeValue value;
		AmplitudeTime time;
	};

	mutable SimpleAmplitudeProcessor _ampE, _ampN, _ampZ;
	OPT (ComponentResult)
	_results[3];
};

class MagnitudeProcessor_K_Class : public Processing::MagnitudeProcessor
{

  public:
	double maxDepth;
	double dist;

  public:
	MagnitudeProcessor_K_Class ()
		: Processing::MagnitudeProcessor (MAG_TYPE), maxDepth (DEPTH_MAX) {}

	string
	amplitudeType () const override
	{
		return MAG_TYPE;
	}
	void
	setDefaults() override {
		_minimumDistanceDeg = DELTA_MIN;
		_maximumDistanceDeg = DELTA_MAX;
		_maximumDepthKm = DEPTH_MAX;
	}
	bool
	setup (const Processing::Settings &settings) override
	{
		return true;
	}

	Processing::MagnitudeProcessor::Status
	computeMagnitude (
		double amplitude, const std::string &unit, double period,
		double snr, double delta, double depth, const DataModel::Origin *,
		const DataModel::SensorLocation *, const DataModel::Amplitude *,
		const Locale *, double &value) override
	{
		if (delta < DELTA_MIN || delta > DELTA_MAX)
		{
			return DistanceOutOfRange;
		}

		if (depth > maxDepth)
		{
			return DepthOutOfRange;
		}

		return compute_K_Class (amplitude, delta, depth, &value);
	}

  private:
	MagnitudeProcessor::Status
	compute_K_Class (
		double amplitude, double delta, double depth, double *mag)
	{
		float epdistkm, hypdistkm, magcalc;

		if (amplitude <= 0.)
		{
			*mag = 0;
			return Error;
		}
		epdistkm = Math::Geo::deg2km (delta);
		hypdistkm = sqrt (epdistkm * epdistkm + depth * depth);

		magcalc = 2.94 + 1.935 * (log10 (amplitude) + 1.734 * log10 (hypdistkm));
		*mag = magcalc;

		return OK;
	}
};

REGISTER_AMPLITUDEPROCESSOR (AmplitudeProcessor_K_Class, MAG_TYPE);
REGISTER_MAGNITUDEPROCESSOR (MagnitudeProcessor_K_Class, MAG_TYPE);

} // namespace