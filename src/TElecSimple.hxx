#ifndef TElecSimple_hxx_seen
#define TElecSimple_hxx_seen

#include <TEvent.hxx>
#include <TMCChannelId.hxx>

namespace CP {
    class TElecSimple;
};

/// This is a very simplistic electronics simulation.  It is not intended for
/// doing physic, but does capture enough of the behavior to develop software.
class CP::TElecSimple {
public:
    typedef std::vector<double> DoubleVector;
    TElecSimple();
    ~TElecSimple();

    void operator()(CP::TEvent& event);

private:
    /// The integration window for the trigger.
    double fIntegrationWindow;

    /// The energy threshold for a trigger.  This really should be in terms of
    /// photons, but this doesn't do a light collection simulation.
    double fThreshold;

    /// The time between triggers.
    double fIntegrationTime;

    /// The time between triggers.
    double fTriggerOffset;

    /// The start of the simulation window.
    double fStartIntegration;

    /// The end of the simulation window.  
    double fStopIntegration;

    /// The drift velocity.
    double fDriftVelocity;

    /// The diffusion coefficient;
    double fDiffusionCoeff;

    /// The recombination probability for electrons
    double fRecombination;

    /// The activation energy for Argon
    double fActivationEnergy;

    /// The electron life time.
    double fElectronLife;

    /// The wire noise level.
    double fNoiseSigma;

    /// The rise time for the amplifier
    double fAmplifierRise;

    /// The gain of the amplifier for the collection plane.  This must be
    /// matched to the range of the ADC.
    double fAmplifierCollectionGain;

    /// The gain of the amplifier for the induction planes.  This must be
    /// matched to the range of the ADC.
    double fAmplifierInductionGain;

    /// The gain of the amplifier for the ganged PMTs.  This must be
    /// matched to the range of the ADC.
    double fAmplifierPMTGain;

    /// The time step for each digitization bin.
    double fDigitStep;

    /// The threshold to start digitization.
    double fDigitThreshold;

    /// The pedestal
    double fDigitPedestal;

    /// The ADC range
    int fDigitRange;

    /// The amount of time to save before a threshold crossing.
    double fDigitPreTrigger;

    /// The amount of time to save after a threshold crossing.
    double fDigitPostTrigger;

    /// Add an elecSim "Header" to pass necessary constants to the
    /// calibration.  The header consists of several TRealDatum arrays indexed
    /// by the calibration type.  The index is defined as:
    ///
    /// * truth/eHeader/digitStep: The time step per bin.
    /// * truth/eHeader/gain:      The gain.
    /// * truth/eHeader/pedestal:  The pedestal.
    ///
    /// Each array is indexed by the plane [as defined in the MCChannelId, 0)
    /// X, 1) V, 2) U].  The light sensor calibrations are then saved in index
    /// 3.
    void AddElecSimHeader(CP::TEvent& ev);

    /// Find the trigger times for an event.  This is controlled by paramters
    /// in the elecSim.parameters.dat file.
    void GenerateTriggers(CP::TEvent& ev, DoubleVector& triggers);

    /// Fill a vector full of the charge arrival times for the PMTS.
    void LightSignal(CP::TEvent& ev, CP::TMCChannelId chan, DoubleVector& out);

    /// Fill a vector full of the charge arrival times for a particular wire.
    bool DriftCharge(CP::TEvent& ev, CP::TMCChannelId chan, DoubleVector& out);

    /// Add noise to a vector of charge arrival times.
    void AddNoise(CP::TMCChannelId channel, DoubleVector& out);
    
    /// Fill a vector with the amplified values for the charge.  
    void ShapeCharge(CP::TMCChannelId channel,
                     const DoubleVector& in, DoubleVector& out);

    /// Translate the shaped charge into digitized values.  This adds the
    /// digits to the event.
    void DigitizeCharge(CP::TEvent& ev, CP::TMCChannelId channel,
                        const DoubleVector& in);
};
#endif
