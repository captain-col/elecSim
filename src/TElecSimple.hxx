#ifndef TElecSimple_hxx_seen
#define TElecSimple_hxx_seen

#include <TEvent.hxx>
#include <TMCChannelId.hxx>
#include <TMCDigit.hxx>

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
    typedef std::vector<double> RealVector;
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

    /// The spectral noise RMS.  This is the noise that has a specific
    /// spectral dependence.  This is the noise associated with the
    /// amplifier/shaper.  The nominal value is 650 pe after shaping.
    double fSpectralNoise;

    /// The power value for the f^(-alpha) spectrum of the pink noise.
    double fSpectralAlpha;
    
    /// The low cut off for the pink noise.  It is specified in Hertz.  Below
    /// this index white noise is used.
    double fSpectralLowCut;

    /// The high cut off for the pink noise.  It is specified Hertz.  Above
    /// this index white noise is used.
    double fSpectralHighCut;
    
    /// The series resistance in the noise "short" impedance.  This must be
    /// between 0.0 and 1.0.
    double fSpectralResist;

    /// The half power frequency of the inductor in the parallel RLC part of
    /// the circuit.
    double fSpectralIndFreq;

    /// The internal resistance of the inductor in the parallel RLC part of
    /// the circuit.
    double fSpectralIndRes;

    /// An empirical factor modifying the shape of the impedance vs frequency
    /// for the inductor.
    double fSpectralIndPow;

    /// The half power frequency of the capacitor in the parallel RLC part of
    /// the circuit.
    double fSpectralCapFreq;

    /// The internal resistance of the capacitor in the parallel RLC part of
    /// the circuit.
    double fSpectralCapRes;
    
    /// An empirical factor modifying the shape of the impedance vs frequency
    /// for the capacitor.
    double fSpectralCapPow;

    /// The charge induction factor for the induction wires.  This is assumed
    /// to be the same for the U and V planes since they have similar
    /// electrical properties.  A factor of 1.0 means that an average charge
    /// |e-| induced |e-| on the wire.
    double fWireInductionFactor;

    /// The detection efficiency for a photon.  This should be a function of
    /// position in the detector, but this is the *simple* electronics
    /// simulation!  This is expressed in terms of photons collected per MeV.
    double fPhotonCollection;

    /// The dark current for the PMTs
    double fPMTDarkCurrent;
    
    /// The fraction of the light emitted in the short part of the
    /// scintillation
    double fShortFraction;

    /// The time constant for the short component of the scintillation.
    double fShortTime;

    /// The time constant for the long component of the scintillation.
    double fLongTime;

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

    /// The shape factor for the leading edge of the TPC pulse.
    double fAmplifierRiseShape;

    /// The shape factor for the trailing edge of the TPC pulse.
    double fAmplifierFallShape;
    
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
    
    /// The time step oversampling
    int fDigitOversample;
    
    /// The time step for the PMT digitiation bins.
    double fPMTStep;

    /// The amount of time before the trigger to start digitization.  This
    /// will specify the "zero" time bin for the event.
    double fDigitPreTriggerTime;

    /// The amount of time after the trigger to end digitization.  This
    /// will specify the total number of time bins in an event.
    double fDigitPostTriggerTime;

    /// The amount of time before the trigger to start digitization.  This
    /// will specify the "zero" time bin for the event.
    double fPMTPreTriggerTime;

    /// The amount of time after the trigger to end digitization.  This
    /// will specify the total number of time bins in an event.
    double fPMTPostTriggerTime;

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
    void GenerateTriggers(CP::TEvent& ev, RealVector& triggers);

    /// Fill a vector full of the photon arrival times for the PMTS.  The
    /// times are sorted from first to last.
    void LightSignal(CP::TEvent& ev, CP::TMCChannelId chan,
                     RealVector& out,
                     CP::TMCDigit::ContributorContainer& contrib,
                     CP::TMCDigit::InfoContainer& info);
    
    /// Build the digit for a PMT.  This adds the digits to the event.  The
    /// input vector contains the times of the photon arrival at the channel,
    /// The triggers vector is an input of the times that the trigger
    /// criteria was met.
    void DigitizeLight(CP::TEvent& ev, CP::TMCChannelId channel,
                       const RealVector& input,
                       const RealVector& triggers,
                       const CP::TMCDigit::ContributorContainer& contrib,
                       const CP::TMCDigit::InfoContainer& info);

    /// Fill a vector full of the charge arrival times for a particular wire.
    double DriftCharge(CP::TEvent& ev, CP::TMCChannelId chan,
                       RealVector& out,
                       CP::TMCDigit::ContributorContainer& contrib,
                       CP::TMCDigit::InfoContainer& info);

    /// Add the wire noise to a vector of charge arrival times.  This noise is
    /// added in the time domain, so it can be added directly to the drifted
    /// charge.  This is controlled by elecSim.simple.wire.noise.
    void AddWireNoise(CP::TMCChannelId channel, RealVector& out);
    
    /// Generate the background spectrum that can be added to the wire.  This
    /// is where in stocastic noise which has a specific frequency spectrum
    /// can be added as well as any backgrounds induced because of ground
    /// loops (etc).
    void GenerateBackgroundSpectrum(CP::TMCChannelId channel,
                                    ComplexVector& out);
    

    /// Generate the FFT of the electronics response.  This fills the
    /// necessary vectors for use later.
    void GenerateResponseFFT(std::size_t samples);
    
    /// Fill a vector with the shaped values for the charge.
    void ShapeCharge(CP::TMCChannelId channel,
                     const RealVector& in, const ComplexVector& bkg,
                     RealVector& out);

    /// Translate the shaped wire charge into digitized values.  This adds the
    /// digits to the event.  This step includes the amplification.
    void DigitizeWire(CP::TEvent& ev, CP::TMCChannelId channel,
                         const RealVector& input,
                         const RealVector& triggers,
                         const CP::TMCDigit::ContributorContainer& contrib,
                         const CP::TMCDigit::InfoContainer& info);

    /// Calculating the induced current on the wire needs the potential if the
    /// wire was held at 1*unit::volt and everything else was grounded, as
    /// well as the electron velocity as a function of position.  That's a
    /// fairly complicated calculation.  It's simplified to assume that the
    /// electron travels in a straight line with an impact parameter of
    /// corrected distance, and at constant velocity.  The potential is
    /// simplified to be 1/r - 1/r_max where rmax is the distance to the wire
    /// as the electron passes the grid plane.  This estimates the shape of
    /// (velocity)*(electric field), and is normalized to so that the integral
    /// from 0.0 to gDist is 1.0.  Notice the sign is set so that this starts
    /// out positive (so that it matchs the behavior of the electronics).
    double InducedShape(double dist, double impact);

    /// The pulse shape from a delta function impulse at t equal to zero.
    /// This is not normalized.  It takes three parameters: "t" is the time
    /// since the delta-function impulse, "w" is the window to average the
    /// shape over, and "samples" is the number of samples used to do the
    /// average.
    double PulseShaping(double t, double w, int samples=100);
    
    /// The charge induced on the induction wires as a delta function of
    /// charge passes.  This is normalized to the electron charge.
    double InducedCharge(bool isCollection, double impact,
                         double tSample, double sStep);

    /// Find the digits in the input range.  The start is the bin number in
    /// input to start looking for a new digit at.  The startBin and stopBin
    /// give the digitization range in the input.  The return value is a pair
    /// of integers giving the bins in the input to make into a digit.  This
    /// is where zero suppression will be implemented.
    std::pair<int,int> FindDigitRange(int start,
                                      int startBin,
                                      int stopBin,
                                      const RealVector& input);

    /// An FFT used to for the convolution.
    TVirtualFFT* fFFT;

    /// An FFT used to invert for the convolution
    TVirtualFFT* fInvertFFT;

    /// A buffer to hold the FFT of the delta function response.
    ComplexVector fResponseFFT;

};
#endif
