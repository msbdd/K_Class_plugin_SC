#define SEISCOMP_COMPONENT K_Class
// Valid within 0-10 degrees from 0 to 80 km;
#define DELTA_MIN 0.
#define DELTA_MAX 10
#define DEPTH_MAX 80
#define MAG_TYPE "K_Class"

#include "K-Class.h"

#include <iostream>
#include <boost/bind/bind.hpp>

#include <seiscomp/logging/log.h>
#include <seiscomp/math/geo.h>
#include <seiscomp/processing/amplitudes/ML.h>
#include <seiscomp/processing/magnitudeprocessor.h>
#include <seiscomp/seismology/ttt.h>
#include <seiscomp/client/application.h>

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

	v.lowerUncertainty = std::sqrt (l0 * l0 + l1 * l1);
	v.upperUncertainty = std::sqrt (u0 * u0 + u1 * u1);

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

ADD_SC_PLUGIN ("K_Class magnitude", "Dmitry Sidorov-Biryukov", 0, 0, 3)

// We need to create a custom non abstract class for individual magnitude
// computation

class SimpleAmplitudeProcessor
	: public Processing::AbstractAmplitudeProcessor_ML
{
  public:
	SimpleAmplitudeProcessor ()
		: Processing::AbstractAmplitudeProcessor_ML (MAG_TYPE)
	{

	}
	void setTravelTimeTable(TravelTimeTableInterfacePtr ttt)
	{
        _ttt = ttt;
	}
	friend class AmplitudeProcessor_K_Class;

  private:
	bool _haveP;
	bool _haveS;
	Core::Time _sArrival;
	Core::Time _pArrival;
	TravelTimeTableInterfacePtr _ttt;

	void setEnvironment(const DataModel::Origin *hypocenter,
						const DataModel::SensorLocation *receiver,
						const DataModel::Pick *pick) 
		{
			AmplitudeProcessor::setEnvironment(hypocenter, receiver, pick);
			double dist, az, baz;
			double hypoLat, hypoLon, hypoDepth;
			double recvLat, recvLon;
			try {
				hypoLat   = _environment.hypocenter->latitude().value();
				hypoLon   = _environment.hypocenter->longitude().value();
				hypoDepth = _environment.hypocenter->depth().value();
				recvLat   = _environment.receiver->latitude();
				recvLon   = _environment.receiver->longitude();
			}
			catch (...) {
				setStatus(Error, 1);
				return;
			}
			Math::Geo::delazi_wgs84(hypoLat, hypoLon, recvLat, recvLon,
								&dist, &az, &baz);

			if ( dist < _config.minimumDistance || dist > _config.maximumDistance ) {
			setStatus(DistanceOutOfRange, dist);
			return;
			}

			// Horizontal → default behaviour
			if (usedComponent() == FirstHorizontal || usedComponent() == SecondHorizontal) {
				SEISCOMP_DEBUG("setEnvironment: Horizontal component → using default");
				return;
			}
			_haveP = _haveS = false;
			// Compute station elevation if present (optional safe)
			double recvElev = 0.0;
			try {
				recvElev = _environment.receiver->elevation();
			} catch (...) {}

			// Compute travel times
			TravelTimeList *ttlist =
				_ttt->compute(hypoLat, hypoLon, hypoDepth,
								recvLat, recvLon, recvElev, 1);

			if (!ttlist || ttlist->isEmpty()) {
				delete ttlist;
				setStatus(Error, 3);
				return;
			}

			// Fetch P and S
			const TravelTime *tp = getPhase(ttlist, "P");
			const TravelTime *ts = getPhase(ttlist, "S");
			if (tp) {
				_pArrival = _environment.hypocenter->time().value()
							+ Core::TimeSpan(tp->time);
				_haveP = true;
			}

			if (ts) {
				_sArrival = _environment.hypocenter->time().value()
							+ Core::TimeSpan(ts->time);
				_haveS = true;
			}

			delete ttlist;

			if (!_haveP || !_haveS) {
				setStatus(Error, 5); // cannot create accurate window
				return;
			}

		}

	void computeTimeWindow() override 
	{
		if (usedComponent() == FirstHorizontal || usedComponent() == SecondHorizontal) {
			SEISCOMP_DEBUG("computeTimeWindow: Horizontal component → using default");
			return AmplitudeProcessor::computeTimeWindow();
		}

		// Vertical component (Z): Custom window around P wave - from P to S arrivals - should be the P wave maximum
		SEISCOMP_DEBUG("computeTimeWindow: Vertical component → custom P-wave window");
		SEISCOMP_DEBUG("P wave arrival %s", _pArrival.toString("%F %T.%f00000"));
		SEISCOMP_DEBUG("S wave arrival %s", _sArrival.toString("%F %T.%f00000"));
		AmplitudeProcessor::setTimeWindow(Core::TimeWindow(_pArrival, _sArrival));
		}

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
				amp_index = find_absmax (data.size (), data.typedData (), 0, data.size (), offset);
				SEISCOMP_DEBUG ("Data size = %u",data.size ()); 
				SEISCOMP_DEBUG ("P wave max detected at %u",amp_index);
				maxAmplitude = fabs (data[amp_index]) - offset;
			}
			if (usedComponent () == FirstHorizontal || usedComponent () == SecondHorizontal)
			{
				amp_index = find_absmax (
					data.size (), data.typedData (), si1, si2, offset);
				SEISCOMP_DEBUG ("si1 = %u", si1);
				SEISCOMP_DEBUG ("S wave max detected at %u", amp_index);
				SEISCOMP_DEBUG ("si2 is %u", si2);
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
  	TravelTimeTableInterface *_ttt;
	AmplitudeProcessor_K_Class ()
		: Processing::AmplitudeProcessor (MAG_TYPE)
	{	
		setMinSNR (0);
		setMinDist (DELTA_MIN);
		setMaxDist (DELTA_MAX);
		setMaxDepth (DEPTH_MAX);
		setUsedComponent (Any);

		_ampN.setUsedComponent (FirstHorizontal);
		_ampE.setUsedComponent (SecondHorizontal);
		_ampZ.setUsedComponent (Vertical);

		_ampN.setPublishFunction (boost::bind (
			&AmplitudeProcessor_K_Class::newAmplitude, this, boost::placeholders::_1, boost::placeholders::_2));
		_ampE.setPublishFunction (boost::bind (
			&AmplitudeProcessor_K_Class::newAmplitude, this, boost::placeholders::_1, boost::placeholders::_2));
		_ampZ.setPublishFunction (boost::bind (
			&AmplitudeProcessor_K_Class::newAmplitude, this, boost::placeholders::_1, boost::placeholders::_2));
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
		// Create travel-time table interface
		// Model and interface are read from global config

		auto app = Seiscomp::Client::Application::Instance();
		std::string interface = app->configGetString("amplitudes.ttt.interface");
		std::string model = app->configGetString("amplitudes.ttt.model");
		_ttt = TravelTimeTableInterface::Create(interface.c_str());
		_ttt->setModel(model.c_str());
		_ampZ.setTravelTimeTable(_ttt);
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
		if (_ampZ.status() == DistanceOutOfRange) 
		{
            setStatus(DistanceOutOfRange, _ampZ.statusValue());
            return;
        }
		if (_ampZ.status() == DepthOutOfRange) 
		{
            setStatus(DepthOutOfRange, _ampZ.statusValue());
            return;
        }
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

	const AmplitudeProcessor *componentProcessor (Component comp) const override
	{
		switch (comp)
		{
		case FirstHorizontalComponent:
			return &_ampN;
		case SecondHorizontalComponent:
			return &_ampE;
		case VerticalComponent:
			return &_ampZ;
		default:
			break;
		}

		return nullptr;
	}

	const DoubleArray *processedData (Component comp) const override
	{
		switch (comp)
		{
		case FirstHorizontalComponent:
			return _ampN.processedData (comp);
		case SecondHorizontalComponent:
			return _ampE.processedData (comp);
		case VerticalComponent:
			return _ampZ.processedData (comp);
		default:
			break;
		}

		return nullptr;
	}

	int capabilities () const override
	{
		return _ampN.capabilities () | _ampE.capabilities () | _ampZ.capabilities ();
	}

	void reprocess (OPT (double) searchBegin, OPT (double) searchEnd) override
	{
		setStatus (WaitingForData, 0);
		_ampN.setConfig (config ());
		_ampE.setConfig (config ());
		_ampZ.setConfig (config ());

		_results[0] = _results[1] = _results[2] = Core::None;

		_ampN.reprocess (searchBegin, searchEnd);
		_ampE.reprocess (searchBegin, searchEnd);
		_ampZ.reprocess (searchBegin, searchEnd);

		if (!isFinished ())
		{
			if (_ampN.status () > Finished)
				setStatus (_ampN.status (), _ampN.statusValue ());
			else if (_ampE.status () > Finished)
				setStatus (_ampE.status (), _ampE.statusValue ());
			else if (_ampZ.status () > Finished)
				setStatus (_ampZ.status (), _ampZ.statusValue ());
		}
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
	double dist;
	double A, a1, a2, a3, a4;
	double b1, b2, b3, b4;
	double l1, l2, l3;

  public:
	MagnitudeProcessor_K_Class ()
		: Processing::MagnitudeProcessor (MAG_TYPE){}

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
		Processing::MagnitudeProcessor::setup(settings);
		try {
			l1 = settings.getDouble ("magnitudes.K_Class.l1");
		}
		catch ( ... ) {
			l1 = 75.0;
		}
		try {
			l2 = settings.getDouble ("magnitudes.K_Class.l2");
		}
		catch ( ... ) {
			l2 = 264.0;
		}
		try {
			l3 = settings.getDouble ("magnitudes.K_Class.l3");
		}
		catch ( ... ) {
			l3 = 800.0;
		}
		try {
			A = settings.getDouble ("magnitudes.K_Class.A");
		}
		catch ( ... ) {
			A = 1.84;
		}
		try {
			a1 = settings.getDouble ("magnitudes.K_Class.a1");
		}
		catch ( ... ) {
			a1 = 2.11;
		}
		try {
			a2 = settings.getDouble ("magnitudes.K_Class.a2");
		}
		catch ( ... ) {
			a2 = 1.1;
		}
		try {
			a3 = settings.getDouble ("magnitudes.K_Class.a3");
		}
		catch ( ... ) {
			a3 = 2.98;
		}
		try {
			a4 = settings.getDouble ("magnitudes.K_Class.a4");
		}
		catch ( ... ) {
			a4 = 0.0;
		}
		try {
			b1 = settings.getDouble ("magnitudes.K_Class.b1");
		}
		catch ( ... ) {
			b1 = 1.32;
		}
		try {
			b2 = settings.getDouble ("magnitudes.K_Class.b2");
		}
		catch ( ... ) {
			b2 = 3.21;
		}
		try {
			b3 = settings.getDouble ("magnitudes.K_Class.b3");
		}
		catch ( ... ) {
			b3 = -1.34;
		}
		try {
			b4 = settings.getDouble ("magnitudes.K_Class.b4");
		}
		catch ( ... ) {
			b4 = 8.0;
		}
		return true;
	}


	Processing::MagnitudeProcessor::Status
	computeMagnitude (
		double amplitude, const std::string &unit, double period,
		double snr, double delta, double depth, const DataModel::Origin *,
		const DataModel::SensorLocation *, const DataModel::Amplitude *,
		const Locale *, double &value) override
	{	
		SEISCOMP_DEBUG("Delta = %f", delta);
		if ((delta < DELTA_MIN)|| (delta > DELTA_MAX))
		{
			return DistanceOutOfRange;
		}

		if (depth > DEPTH_MAX)
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
		if (hypdistkm <= l1)
		{
			magcalc = A * (log10 (amplitude) + a1 * log10 (hypdistkm) + b1);
		}
		else if (hypdistkm > l1 && hypdistkm <= l2)
		{
			magcalc = A * (log10 (amplitude) + a2 * log10 (hypdistkm) + b2);
		}
		else if (hypdistkm > l2 && hypdistkm <= l3)
		{
			magcalc = A * (log10 (amplitude) + a3 * log10 (hypdistkm) + b3);
		}
		else
		{
			magcalc = A * (log10 (amplitude) + a4 * log10 (hypdistkm) + b4);

		}

			*mag = magcalc;

		return OK;
	}
};

REGISTER_AMPLITUDEPROCESSOR (AmplitudeProcessor_K_Class, MAG_TYPE);
REGISTER_MAGNITUDEPROCESSOR (MagnitudeProcessor_K_Class, MAG_TYPE);

} // namespace
