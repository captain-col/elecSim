#ifndef TElecSimple_hxx_seen
#define TElecSimple_hxx_seen

#include <TEvent.hxx>

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

    // The rise time for the amplifier
    double fAmplifierRise;

    // The gain of the amplifier for the collection plane.
    double fAmplifierCollectionGain;

    // The gain of the amplifier for the induction planes.
    double fAmplifierInductionGain;

    /// The time step for each digitization bin.
    double fDigitStep;

    /// The threshold to start digitization.
    double fDigitThreshold;

    /// The pedestal
    double fDigitPedestal;

    /// The ADC range
    double fDigitRange;

    /// The amount of time to save before a threshold crossing.
    double fDigitPreTrigger;

    /// The amount of time to save after a threshold crossing.
    double fDigitPostTrigger;

    /// Find the trigger times for an event.  This is controlled by paramters
    /// in the elecSim.parameters.dat file.
    void GenerateTriggers(CP::TEvent& ev, DoubleVector& triggers);

    /// Fill a vector full of the charge arrival times for a particular wire.
    bool CollectCharge(CP::TEvent& ev, int plane, int wire, DoubleVector& out);

    /// Add noise to a vector of charge arrival times.
    void AddNoise(int plane, DoubleVector& out);
    
    /// Fill a vector with the amplified values for the charge.  
    void ShapeCharge(int plane, int wire,
                     const DoubleVector& in, DoubleVector& out);

    /// Translate the shaped charge into digitized values.  This adds the
    /// digits to the event.
    void DigitizeCharge(CP::TEvent& ev, int plane, int wire, 
                        const DoubleVector& in);
};
#endif
