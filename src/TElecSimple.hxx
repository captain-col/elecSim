#ifndef TElecSimple_hxx_seen
#define TElecSimple_hxx_seen

#include <TEvent.hxx>
#include <TMCChannelId.hxx>

#include <vector>
#include <complex>

namespace CP {
    class TElecSimple;
};

class TVirtualFFT;

/// This is a very simple electronics simulation.  The main purpose of the
/// simple electronics simulation is to capture enough of the detector
/// response to develop software, but to predict the behavior.  It's a purely
/// emperical model with lots of tunable parameters, but given that it's a
/// "fit", it can do a rather quick simulation.  Most importantly, it's a good
/// test bed for understanding the interaction of the CAPTAIN reconstruction
/// and calibration to different electronics behaviors.
class CP::TElecSimple {
public:
    typedef std::vector<double> DoubleVector;
    typedef std::vector< std::complex<double> > ComplexVector;
    TElecSimple();
    ~TElecSimple();

    void operator()(CP::TEvent& event);

private:
    /// The start of the simulation window.  The FADC time steps are measured
    /// relative to this time.  The start of the simulation must always
    /// before the earliest possible trigger, or charge will not be properly
    /// collected.  In the detector, the start of the integration will be the
    /// trigger time minus a pre-trigger offset.
    double fStartSimulation;

    /// The end of the simulation window.  This needs to be after the start of
    /// the integration window plus the length of the integration window.
    double fStopSimulation;

    // Parameters controlling the "physics" of ionization and scintillation. 

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

    /// The wire noise level.  This is the amount of thermal noise generated
    /// in the wires (and is not associated with the amplifier).
    double fWireNoise;

    /// The charge induction factor for the induction wires.  This is assumed
    /// to be the same for the U and V planes since they have similar
    /// electrical properties.  A factor of 1.0 means that an average charge
    /// |e-| induced |e-| on the wire.
    double fWireInductionFactor;

    /// The detection efficiency for a photon.  This should be a function of
    /// position in the detector, but this is the *simple* electronics
    /// simulation!  This is expressed in terms of photons collected per MeV.
    double fPhotonCollection;

    /// The fraction of the light emitted in the short part of the
    /// scintillation
    double fShortFraction;

    /// The time constant for the short component of the scintillation.
    double fShortTime;

    /// The time constant for the long component of the scintillation.
    double fLongTime;

    /// The number of light sensors in the detector.
    int fLightSensorCount;

    // Parameters controlling the trigger simulation.

    /// The time window used to generate the trigger.  The current simulation
    /// trigger has two modes.  In the first mode, (a beam trigger), the
    /// trigger time is always zero.  In that case, the trigger window should
    /// be zero or negative.  In the second mode, the trigger window is
    /// positive, and the deposited energy in the window must be above
    /// fTriggerThreshold.
    double fTriggerWindow;

    /// The energy threshold for a trigger.  This really should be in terms of
    /// photons, but this doesn't do a light collection simulation.
    double fTriggerThreshold;

    /// The "dead" time between triggers.  The electronics aren't really dead,
    /// but this is the amount of time before we can accept a new trigger.
    double fTriggerDeadTime;

    /// The time between triggers.
    double fTriggerOffset;

    /// The amount of time to simulate before the trigger.
    double fPreTriggerTime;

    /// The amount of time to simulate after the trigger.
    double fPostTriggerTime;

    // Parameters controlling the amplifier simulations.

    /// The rise time for the TPC amplifier
    double fAmplifierRise;

    /// A flag for whether the amplification describes the pulse height, or
    /// area of a single sample pulse is described by the amplification.  If
    /// this is true, then the amplification is the ratio of areas for the
    /// input and the shaped, amplified pulse.  If this is false, then the
    /// amplification is the ratio pulse heights for the input and the shaped,
    /// amplified pulse.
    double fAmplifierConserveIntegral;

    /// The gain of the amplifier for the collection plane.  This must be
    /// matched to the range of the ADC.
    double fAmplifierCollectionGain;

    /// The gain of the amplifier for the induction planes.  This must be
    /// matched to the range of the ADC.
    double fAmplifierInductionGain;

    /// The gain of the amplifier for the ganged PMTs.  This must be
    /// matched to the range of the ADC.
    double fAmplifierPMTGain;

    /// The width of the PMT 1 pe peak.  The is the uncertainty on the gain.
    double fPMTPeak;

    // The digitization parameters.

    /// The digitization noise.  This is the total noise contributed by the
    /// electronics (after the shaping).  It is in units of ADC counts.
    double fDigitNoise;

    /// The time step for each digitization bin.
    double fDigitStep;
    
    /// The time step for the PMT digitiation bins.
    double fPMTStep;

    /// The amount of time before the trigger to start digitization.  This
    /// will specify the "zero" time bin for the event.
    double fDigitPreTriggerTime;

    /// The amount of time after the trigger to end digitization.  This
    /// will specify the total number of time bins in an event.
    double fDigitPostTriggerTime;

    /// The pedestal
    double fDigitPedestal;

    /// The slope of the ADC conversion in ADC/mV (i.e. ADC/(input signal)).
    double fDigitSlope;

    /// The ADC maximum value.
    int fDigitMaximum;

    /// The ADC minimum value.
    int fDigitMinimum;

    /// The threshold to start digitization.  If this is greater than zero,
    /// then there will be zero suppression applied to the digits.
    double fDigitThreshold;

    /// The amount of time to save before a threshold crossing.  This only
    /// applies when there is zero suppression.
    double fDigitPreThreshold;

    /// The amount of time to save after a threshold crossing.  This only
    /// applies when there is zero suppression.
    double fDigitPostThreshold;

    /// Add an elecSim "Header" to pass necessary constants to the
    /// calibration.  The header consists of several TRealDatum arrays indexed
    /// by the calibration type.  The index is defined as:
    ///
    /// * truth/elecSimple/digitStep: The time step per bin.
    /// * truth/elecSimple/gain:      The gain.
    /// * truth/elecSimple/pedestal:  The pedestal.
    /// * truth/elecSimple/slope:     The slope of the digitizer in adc/V.
    ///
    /// Each array is indexed by the plane [as defined in the MCChannelId, 0)
    /// X, 1) V, 2) U].  The light sensor calibrations are then saved in index
    /// 3.
    void AddElecSimHeader(CP::TEvent& ev);

    /// Find the trigger times for an event.  This is controlled by paramters
    /// in the elecSim.parameters.dat file.
    void GenerateTriggers(CP::TEvent& ev, DoubleVector& triggers);

    /// Fill a vector full of the photon arrival times for the PMTS.
    void LightSignal(CP::TEvent& ev, CP::TMCChannelId chan, DoubleVector& out);

    /// Build the digit for the PMT..  This adds the
    /// digits to the event.
    void DigitizeLight(CP::TEvent& ev, CP::TMCChannelId channel,
                       const DoubleVector& input,
                       const DoubleVector& triggers);

    /// Fill a vector full of the charge arrival times for a particular wire.
    bool DriftCharge(CP::TEvent& ev, CP::TMCChannelId chan, DoubleVector& out);

    /// Add the wire noise to a vector of charge arrival times.
    void AddWireNoise(CP::TMCChannelId channel, DoubleVector& out);
    
    /// The pulse shape from a delta function impulse.  This is not
    /// normalized.  It also gives the integral over the digit step.
    double PulseShaping(double t);
    
    /// The charge induced on the induction wires as a delta function of
    /// charge passes.  This is not normalized.
    double InducedCharge(double t);
    
    /// Fill a vector with the amplified values for the charge.  
    void ShapeCharge(CP::TMCChannelId channel,
                     const DoubleVector& in, DoubleVector& out);

    /// An FFT used to for the convolution.
    TVirtualFFT* fFFT;

    /// An FFT used to invert for the convolution
    TVirtualFFT* fInvertFFT;

    /// A buffer to hold the FFT of the delta function response.
    ComplexVector fResponseFFT;

    /// A buffer to hold the FFT of the induced charge response.  This is only
    /// used for induction wires.
    ComplexVector fInducedFFT;

    /// A buffer to hold the FFT of the capacitive coupling response.  This is
    /// only used for induction wires.
    ComplexVector fCurrentFFT;

    /// Translate the shaped wire charge into digitized values.  This adds the
    /// digits to the event.
    void DigitizeWires(CP::TEvent& ev, CP::TMCChannelId channel,
                        const DoubleVector& input,
                        const DoubleVector& triggers);

    /// Find the digits in the input range.  The start is the bin number in
    /// input to start looking for a new digit at.  The startBin and stopBin
    /// give the digitization range in the input.  The return value is a pair
    /// of integers giving the bins in the input to make into a digit.  This
    /// is where zero suppression will be implemented.
    std::pair<int,int> FindDigitRange(int start,
                                      int startBin,
                                      int stopBin,
                                      const DoubleVector& input);
};
#endif
