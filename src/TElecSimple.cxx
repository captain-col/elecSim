#include "TElecSimple.hxx"
#include "ComplexCauchy.hxx"

#include <TEvent.hxx>
#include <TG4HitSegment.hxx>
#include <TCaptLog.hxx>
#include <TRuntimeParameters.hxx>
#include <HEPUnits.hxx>
#include <HEPConstants.hxx>
#include <TGeomIdManager.hxx>
#include <TManager.hxx>
#include <TPulseMCDigit.hxx>
#include <TRealDatum.hxx>
#include <CaptGeomId.hxx>

#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TGeoBBox.h>

#include <TRandom.h>
#include <TVirtualFFT.h>

#include <TH1F.h>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>
#include <memory>
#include <string>
#include <iostream>
#include <fstream>

CP::TElecSimple::TElecSimple() {
    CaptLog("Starting the electronics simulation");

    // The integration window for the trigger.
    fTriggerWindow
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.trigger.window");

    // The energy threshold for a trigger.  This really should be in terms of
    // photons, but this doesn't do a light collection simulation.
    fTriggerThreshold
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.trigger.threshold");

    // The time between triggers.
    fTriggerDeadTime
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.trigger.deadTime");

    // The time between the trigger time being met and the trigger being
    // formed.
    fTriggerOffset
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.trigger.offset");

    // The amount of time to simulate before the trigger time.
    fPreTriggerTime
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.preTriggerTime");

    // The amount of time to simulate after the trigger time.
    fPostTriggerTime
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.postTriggerTime");

    // The rise time for the amplifier
    fAmplifierRise
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.amplifier.riseTime");

    // The shaping factor for the leading edge of the pulse.
    fAmplifierRiseShape
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.amplifier.riseShape");

    // The shaping factor for the trailing edge of the pulse.
    fAmplifierFallShape
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.amplifier.fallShape");

    // Set if the integral or the pulse height is conserved.  Based on lab
    // measurements, the pulse shaper conserves the peak voltage.
    fAmplifierConserveIntegral
        = CP::TRuntimeParameters::Get().GetParameterB(
            "elecSim.simple.amplifier.integral");

    // The gain of the amplifier for the collection plane.  This must be
    // matched to the range of the ADC.  The units are mV/fC
    fAmplifierCollectionGain
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.amplifier.gain.collection");
    fAmplifierCollectionGain *= unit::mV/unit::fC;

    // The gain of the amplifier for the induction planes.  The units are mV/fC
    fAmplifierInductionGain
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.amplifier.gain.induction");
    fAmplifierInductionGain *= unit::mV/unit::fC;

    // The gain of the amplifier for PMT.  This includes the actual PMT gain,
    // and any preamps.
    fAmplifierPMTGain
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.amplifier.gain.pmt");

    // The width of the 1 pe peak.
    fPMTPeak
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.amplifier.width");

    // The time step for each PMT digitization bin.
    fPMTStep
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.step.pmt");

    // The time step for each digitization bin.
    fDigitStep
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.step");

    // The time step for each digitization bin.
    fDigitOversample
        = CP::TRuntimeParameters::Get().GetParameterI(
            "elecSim.simple.digitization.oversample");
    if (fDigitOversample < 1) fDigitOversample = 1;

    // The amount of time to save before a threshold crossing.
    fDigitPreTriggerTime
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.preTriggerTime");
    if (fDigitPreTriggerTime < 0) fDigitPreTriggerTime = fPreTriggerTime;

    // The amount of time to save after a threshold crossing.
    fDigitPostTriggerTime
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.postTriggerTime");
    if (fDigitPostTriggerTime < 0) fDigitPostTriggerTime = fPostTriggerTime;

    // The amount of time to save before a threshold crossing.
    fPMTPreTriggerTime
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.preTriggerTime.pmt");
    if (fPMTPreTriggerTime < 0) fPMTPreTriggerTime = fPreTriggerTime;

    // The amount of time to save after a threshold crossing.
    fPMTPostTriggerTime
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.postTriggerTime.pmt");
    if (fPMTPostTriggerTime < 0) fPMTPostTriggerTime = fPostTriggerTime;

    // The pedestal.  This is then randomized per channel.
    fDigitPedestal
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.pedestal");

    // The slope of the digitizer in ADC/mV.
    fDigitSlope
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.slope");
    fDigitSlope *= 1.0/unit::mV;

    // The ADC maximum
    fDigitMaximum
        = CP::TRuntimeParameters::Get().GetParameterI(
            "elecSim.simple.digitization.maximum");

    // The ADC range
    fDigitMinimum
        = CP::TRuntimeParameters::Get().GetParameterI(
            "elecSim.simple.digitization.minimum");

    // The threshold to start digitizing a pulse.  This is specified in ADC
    // above pedestal.
    fDigitThreshold
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.threshold");

    // The noise introduced during digitization.  This is in terms of electron
    // equivalent noise.
    fDigitNoise
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.noise");

    // The amount of time to save before a threshold crossing.
    fDigitPreThreshold
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.preThreshold");

    // The amount of time to save after a threshold crossing.
    fDigitPostThreshold
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.postThreshold");

    // The recompination probability
    fRecombination
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.recombinationProb");

    // The ionization energy for argon.
    fActivationEnergy
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.activationEnergy");

    // The drift velocity
    fDriftVelocity
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.drift.velocity");
    fDriftVelocity = 1.6 *unit::mm/(1000*unit::ns);

    // The diffusion coefficient
    fDiffusionCoeff
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.drift.diffusion");

    // The electron lifetime.
    fElectronLife
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.drift.life");

    // The probability that a photon generate in LAr will make a
    // photo-electron.  This includes the PMT response, the light propagation,
    // and everything else.
    fPhotonCollection
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.light.collection");

    fPMTDarkCurrent = CP::TRuntimeParameters::Get().GetParameterD(
        "elecSim.simple.light.darkCurrent");

    fShortFraction
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.light.shortFraction");

    fShortTime
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.light.shortTime");

    fLongTime
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.light.longTime");

    fWireNoise
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.wire.noise");

    fSpectralNoise
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.spectral.noise");
    // The parameter has units of ADC, so turn it into voltage.
    fSpectralNoise /= fDigitSlope;
    
    fSpectralAlpha
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.spectral.alpha");
    
    fSpectralLowCut
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.spectral.lowCut");
    
    fSpectralHighCut
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.spectral.highCut");
    
    fSpectralResist
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.spectral.impedance.resistance");

    fSpectralIndFreq
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.spectral.impedance.inductor.frequency");

    fSpectralIndRes
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.spectral.impedance.inductor.resistance");

    fSpectralIndPow
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.spectral.impedance.inductor.auxPower");

    fSpectralCapFreq
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.spectral.impedance.capacitor.frequency");

    fSpectralCapRes
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.spectral.impedance.capacitor.resistance");
    
    fSpectralCapPow
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.spectral.impedance.capacitor.auxPower");

    fDiscreetScale
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.discreet.scale");
    
    if (CP::TRuntimeParameters::Get().HasParameter(
            "elecSim.simple.discreet.file")) {
        fDiscreetFile
            = CP::TRuntimeParameters::Get().GetParameterS(
                "elecSim.simple.discreet.file");
    }
        
    // The normalization factor for the charge induced on a wire.
    fWireInductionFactor
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.wire.induction");

    fFFT = NULL;
    fInvertFFT = NULL;
}

CP::TElecSimple::~TElecSimple() {}

void CP::TElecSimple::operator()(CP::TEvent& event) {
    CaptLog("Event " << event.GetContext());

    if (fDiscreetPeaks.empty() && !fDiscreetFile.empty()) {
        const char* root = std::getenv("ELECSIMROOT");
        if (root) {
            std::string name(root);
            name += "/parameters/";
            name += fDiscreetFile;
            OpenNoiseFile(name);
        }
        else {
            CaptError("ELECSIMROOT is undefined");
        }
    }
    
    if (event.GetContext().IsCAPTAIN()) {
        if (CP::TRuntimeParameters::Get().HasParameter(
                "elecSim.simple.drift.life.capt")) {
            fElectronLife
                = CP::TRuntimeParameters::Get().GetParameterD(
                    "elecSim.simple.drift.life.capt");
        }
    }
    else if (event.GetContext().IsMiniCAPTAIN()) {
        if (CP::TRuntimeParameters::Get().HasParameter(
                "elecSim.simple.drift.life.mini")) {
            fElectronLife
                = CP::TRuntimeParameters::Get().GetParameterD(
                    "elecSim.simple.drift.life.mini");
        }
    }
    
    // Check that the event has the truth hits.
    CP::THandle<CP::TDataVector> truthHits
        = event.Get<CP::TDataVector>("truth/g4Hits");

    if (!truthHits) {
        CaptError("No truth hits in event");
        return;
    }

    int segCount = 0;
    for (CP::TDataVector::iterator h = truthHits->begin();
         h != truthHits->end();
         ++h) {
        CP::THandle<CP::TG4HitContainer> g4Hits =
            (*h)->Get<CP::TG4HitContainer>(".");
        ++segCount;
    }
    if (segCount < 1) {
        CaptError("Event is missing the TG4HitContainer objects");
        return;
    }

    CP::TManager::Get().Geometry();

    // Add the elecSim header to truth.
    AddElecSimHeader(event);

    // Figure out when the triggers will be.  This is usually "just zero".
    RealVector triggerTimes;
    GenerateTriggers(event,triggerTimes);

    if (triggerTimes.empty()) {
        CaptLog("No trigger in the event");
        return;
    }

    // Find out the total integration period.  This is found based off of the
    // times of the first and last triggers.
    fStartSimulation = 100*unit::second;
    fStopSimulation = -100*unit::second;
    for (RealVector::iterator triggerTime = triggerTimes.begin();
         triggerTime != triggerTimes.end();
         ++triggerTime) {
        fStartSimulation = std::min(fStartSimulation,
                                    *triggerTime - fPreTriggerTime);
        fStopSimulation = std::max(fStopSimulation,
                                   *triggerTime + fPostTriggerTime);
    }

    // Simulate the PMT signals.
    for (int i = 0; i<100; ++i) {
        // Check to see if this PMT exists, quit if it doesn't.
        if (!CP::TManager::Get().GeomId().CdId(
                CP::GeomId::Captain::Photosensor(i))) break;
        TMCChannelId pmt(1,0,i);
        RealVector photonTimes;
        CP::TMCDigit::ContributorContainer contrib;
        CP::TMCDigit::InfoContainer info;
        LightSignal(event,pmt,photonTimes,contrib,info);
        DigitizeLight(event,pmt,photonTimes,triggerTimes,contrib,info);
    }

    // Calculate the charge bin step.  This can be changed to sub-sample the
    // charge distribution in each digitization step.
    double timeStep = fDigitStep/fDigitOversample;

    // Calculate the number of bins in the TPC samples.
    int chargeBins = (fStopSimulation - fStartSimulation)/timeStep;
    chargeBins = 2*(1+chargeBins/2);
    
    // Check if the stop time needs to be adjusted.
    fStopSimulation = fStartSimulation + chargeBins*timeStep;

    RealVector collectedCharge(chargeBins);
    ComplexVector frequencySpectrum(chargeBins);
    ComplexVector workSpace(chargeBins);
    RealVector measuredCharge(chargeBins);

    GenerateResponseFFT(chargeBins);

    // For each wire in the detector, figure out the signal.  The loop is done
    // this way so that we don't need to know how many planes and wires are
    // being simulated.
    int count = 0;
    for (int plane = 0; plane < 5; ++plane) {

        // Check to see if this plane exists, quit if it doesn't.
        if (!CP::TManager::Get().GeomId().CdId(
                CP::GeomId::Captain::Plane(plane))) break;
        for (int wire = 0; wire < 10000; ++wire) {
            // Check to see if the wire exists.  Quit if it doesn't
            if (!CP::TManager::Get().GeomId().CdId(
                    CP::GeomId::Captain::Wire(plane,wire))) break;
            TMCChannelId channel(0,plane,wire);
            if ((++count % 100) == 0) {
                CaptNamedInfo("channel",
                              "Channel " << channel);
            }
            CP::TMCDigit::ContributorContainer contrib;
            CP::TMCDigit::InfoContainer info;

            std::fill(collectedCharge.begin(), collectedCharge.end(), 0.0);
            std::fill(measuredCharge.begin(), measuredCharge.end(), 0.0);
            std::fill(frequencySpectrum.begin(), frequencySpectrum.end(),
                      std::complex<double>(0.0,0.0));
            std::fill(workSpace.begin(), workSpace.end(),
                      std::complex<double>(0.0,0.0));

            // Add the "time space" effect to the charge on the wire.  These
            // are applied before the shaping amplifier.
            DriftCharge(event,channel,collectedCharge,contrib,info);
            AddWireNoise(channel,collectedCharge);
            ShapeCharge(channel,collectedCharge,frequencySpectrum);

            /// Add the post-amplification effects.  These are applied using
            /// in "frequency space".
            GenerateBackgroundSpectrum(channel,frequencySpectrum,workSpace);
            GenerateDiscreetSpectrum(channel,frequencySpectrum,workSpace);

            InverseFFT(frequencySpectrum,measuredCharge);

            // Apply the digitization to the resulting charge.
            DigitizeWire(event,channel,measuredCharge,triggerTimes,
                         contrib,info);

            
        }
    }

    return;
}

void CP::TElecSimple::AddElecSimHeader(CP::TEvent& event) {
    CP::THandle<CP::TDataVector> truth = event.Get<CP::TDataVector>("truth");
    if (!truth) {
        CaptError("No truth information for this event)");
        return;
    }

    truth->AddDatum(new CP::TDataVector("elecSimple"));
    CP::THandle<CP::TDataVector> header
        = truth->Get<CP::TDataVector>("elecSimple");
    if (!header) {
        CaptError("Error adding the elecSimple header");
        return;
    }

    // Fill the digit steps.
    std::unique_ptr<CP::TRealDatum> digitStep(new CP::TRealDatum("digitStep"));
    digitStep->clear();
    digitStep->push_back(fDigitStep);  // X Plane
    digitStep->push_back(fDigitStep);  // V Plane
    digitStep->push_back(fDigitStep);  // U Plane
    digitStep->push_back(fPMTStep);  // PMT
    header->AddDatum(digitStep.release());

    // Fill the trigger offset.
    std::unique_ptr<CP::TRealDatum> triggerOffset(
        new CP::TRealDatum("triggerOffset"));
    triggerOffset->clear();
    triggerOffset->push_back(fDigitPreTriggerTime);  // X Plane
    triggerOffset->push_back(fDigitPreTriggerTime);  // V Plane
    triggerOffset->push_back(fDigitPreTriggerTime);  // U Plane
    triggerOffset->push_back(fPMTPreTriggerTime);  // PMT
    header->AddDatum(triggerOffset.release());

    // Fill the pedestals
    std::unique_ptr<CP::TRealDatum> pedestal(new CP::TRealDatum("pedestal"));
    pedestal->clear();
    pedestal->push_back(fDigitPedestal);  // X Plane
    pedestal->push_back(fDigitPedestal);  // V Plane
    pedestal->push_back(fDigitPedestal);  // U Plane
    pedestal->push_back(fDigitPedestal);  // PMT
    header->AddDatum(pedestal.release());

    // Fill amplifier gains (voltage/(input charge))
    std::unique_ptr<CP::TRealDatum> gain(new CP::TRealDatum("gain"));
    gain->clear();
    gain->push_back(fAmplifierCollectionGain); // X Plane
    gain->push_back(fAmplifierInductionGain);  // V Plane
    gain->push_back(fAmplifierInductionGain);  // U Plane
    gain->push_back(fAmplifierPMTGain);        // PMT
    header->AddDatum(gain.release());

    // Fill the digitizer slopes (ADC/mV).
    std::unique_ptr<CP::TRealDatum> slope(new CP::TRealDatum("slope"));
    slope->clear();
    slope->push_back(fDigitSlope);  // X Plane
    slope->push_back(fDigitSlope);  // V Plane
    slope->push_back(fDigitSlope);  // U Plane
    slope->push_back(1.0);          // PMT
    header->AddDatum(slope.release());

    /// Fill the shape times.
    std::unique_ptr<CP::TRealDatum> shapeTime(new CP::TRealDatum("shape"));
    shapeTime->clear();
    shapeTime->push_back(fAmplifierRise);    // X Plane
    shapeTime->push_back(fAmplifierRise);    // V Plane
    shapeTime->push_back(fAmplifierRise);    // U Plane
    shapeTime->push_back(1.0);               // PMT
    header->AddDatum(shapeTime.release());

    /// Fill the rising edge shape.
    std::unique_ptr<CP::TRealDatum> shapeRise(new CP::TRealDatum("shapeRise"));
    shapeRise->clear();
    shapeRise->push_back(fAmplifierRiseShape);    // X Plane
    shapeRise->push_back(fAmplifierRiseShape);    // V Plane
    shapeRise->push_back(fAmplifierRiseShape);    // U Plane
    shapeRise->push_back(1.0);                    // PMT
    header->AddDatum(shapeRise.release());

    /// Fill the falling edge shape.
    std::unique_ptr<CP::TRealDatum> shapeFall(new CP::TRealDatum("shapeFall"));
    shapeFall->clear();
    shapeFall->push_back(fAmplifierFallShape);    // X Plane
    shapeFall->push_back(fAmplifierFallShape);    // V Plane
    shapeFall->push_back(fAmplifierFallShape);    // U Plane
    shapeFall->push_back(1.0);                    // PMT
    header->AddDatum(shapeFall.release());

    /// Fill the drift velocity and electron lifetime.
    std::unique_ptr<CP::TRealDatum> argonState(new CP::TRealDatum("argon"));
    argonState->clear();
    argonState->push_back(fDriftVelocity);
    argonState->push_back(fElectronLife);
    header->AddDatum(argonState.release());

}

void CP::TElecSimple::FFT(const RealVector& in,ComplexVector& out) {

    if (in.size() != out.size()) {
        CaptError("Invalid input and output sizes");
        throw;
    }
    
    if (!fFFT) {
        int len = in.size();
        fFFT = TVirtualFFT::FFT(1,&len,"R2C K M");
        if (len != (int) in.size()) {
            CaptError("Invalid length for FFT");
            CaptError("     original length: " << in.size());
            CaptError("     allocated length: " << len);
        }
        CaptLog("FFT initialized with " << len << " elements");
    }

    for (std::size_t i = 0; i<in.size(); ++i) {
        fFFT->SetPoint(i,in[i]);
    }

    fFFT->Transform();

    double norm = 1.0/std::sqrt(1.0*out.size());
    for (std::size_t i = 0; i<out.size(); ++i) {
        double rl, im;
        fFFT->GetPointComplex(i,rl,im);
        out[i] = norm*std::complex<double>(rl,im);
    }
}

void CP::TElecSimple::InverseFFT(const ComplexVector& in, RealVector& out) {

    if (in.size() != out.size()) {
        CaptError("Invalid input and output sizes");
        throw;
    }
    
    if (!fInvertFFT) {
        int len = in.size();
        fInvertFFT = TVirtualFFT::FFT(1,&len,"C2R K M");
        if (len != (int) in.size()) {
            CaptError("Invalid length for FFT");
            CaptError("     original length: " << in.size());
            CaptError("     allocated length: " << len);
        }
        CaptLog("FFT initialized with " << len << " elements");
    }

    for (std::size_t i = 0; i<in.size(); ++i) {
        fInvertFFT->SetPoint(i,in[i].real(),in[i].imag());
    }

    fInvertFFT->Transform();
    
    double norm = 1.0/std::sqrt(1.0*out.size());
    for (std::size_t i=0; i<out.size(); ++i) {
        out[i] = norm*fInvertFFT->GetPointReal(i);
    }
}

void CP::TElecSimple::GenerateTriggers(CP::TEvent& event,
                                       RealVector& triggers) {
    triggers.clear();

    if (fTriggerWindow <= 0) {
        triggers.push_back(0.0);
        return;
    }

    // Check that the event has the truth hits.
    CP::THandle<CP::TDataVector> truthHits
        = event.Get<CP::TDataVector>("truth/g4Hits");

    typedef std::vector< std::pair<double, double> > ETMap;
    ETMap energyTime;
    for (CP::TDataVector::iterator h = truthHits->begin();
         h != truthHits->end();
         ++h) {
        CP::THandle<CP::TG4HitContainer> g4Hits =
            (*h)->Get<CP::TG4HitContainer>(".");
        for (CP::TG4HitContainer::const_iterator h = g4Hits->begin();
             h != g4Hits->end();
             ++h) {
            const CP::TG4HitSegment* seg
                = dynamic_cast<const CP::TG4HitSegment*>((*h));
            energyTime.push_back(
                std::make_pair(seg->GetStartT(), seg->GetEnergyDeposit()));
        }
    }
    std::sort(energyTime.begin(), energyTime.end());

    ETMap::iterator start = energyTime.begin();
    ETMap::iterator stop = energyTime.begin();
    double deadUntil = 0.0;
    while (stop != energyTime.end()) {
        // skip dead time.
        if (stop->first < deadUntil) {
            ++stop;
            start = stop;
            continue;
        }
        // Look for charge in a window.
        while ( (stop->first - start->first) > fTriggerWindow) ++start;
        double esum = 0.0;
        for (ETMap::iterator et = start; et != stop; ++et) {
            esum += et->second;
        }
        // If over threshold, then trigger.
        if (esum > fTriggerThreshold) {
            double trigger = stop->first + fTriggerOffset;
            CaptLog("Trigger at " << trigger/unit::ns << " ns" );
            triggers.push_back(trigger);
            deadUntil = stop->first + fTriggerDeadTime;
        }
        ++stop;
    }

}

void CP::TElecSimple::LightSignal(
    CP::TEvent& event,
    CP::TMCChannelId chan,
    RealVector& times,
    CP::TMCDigit::ContributorContainer& contrib,
    CP::TMCDigit::InfoContainer& info) {

    times.clear();

    // Check that the event has the truth hits.
    CP::THandle<CP::TDataVector> truthHits
        = event.Get<CP::TDataVector>("truth/g4Hits");

    int pmt = chan.GetNumber();

    CP::TManager::Get().Geometry(); // just in case...
    if (!CP::TManager::Get().GeomId().CdId(
            CP::GeomId::Captain::Photosensor(pmt))) {
        CaptError("Not a pmt " << chan << " " << pmt);
        return;
    }

    // Arrays for the transforms.
    double local[3];
    double master[3];
    double masterStart[3];
    double masterStop[3];

    int totalPhotons = 0.0;
    double meanPhotons = 0.0;
    for (CP::TDataVector::iterator h = truthHits->begin();
         h != truthHits->end();
         ++h) {
        CP::THandle<CP::TG4HitContainer> g4Hits =
            (*h)->Get<CP::TG4HitContainer>(".");
        for (CP::TG4HitContainer::const_iterator h = g4Hits->begin();
             h != g4Hits->end();
             ++h) {
            const CP::TG4HitSegment* seg
                = dynamic_cast<const CP::TG4HitSegment*>((*h));

            masterStart[0] = seg->GetStartX();
            masterStart[1] = seg->GetStartY();
            masterStart[2] = seg->GetStartZ();
            double startT = seg->GetStartT();
            masterStop[0] = seg->GetStopX();
            masterStop[1] = seg->GetStopY();
            masterStop[2] = seg->GetStopZ();

            for (int i=0; i<3; ++i) {
                master[i] = 0.5*masterStart[i] + 0.5*masterStop[i];
            }

            // The PMT is facing in the positive Z direction, and the origin
            // is at the center of the photocathode.
            gGeoManager->MasterToLocal(master,local);

            // Make sure the energy is in front of the PMT (no reflections are
            // done).
            if (local[2] < 0) continue;

            // Find the distance to the PMT;.
            double distance = std::sqrt(local[0]*local[0]
                                        +local[1]*local[1]
                                        +local[2]*local[2]);
            if (distance < 1*unit::mm) distance = 1*unit::mm;

            // Find the cosine to the normal of the PMT.
            double cangle = local[2]/distance;

            // Find the PMT area (this needs to be calculated from the
            // geometry).
            double pmtArea = -1;
            TGeoVolume *volume = gGeoManager->GetCurrentVolume();
            if (!volume) abort();
            do {
                TGeoTube *tube = dynamic_cast<TGeoTube*>(volume->GetShape());
                if (tube) {
                    double r = tube->GetRmax();
                    pmtArea = 2*unit::pi*r*r;
                    break;
                }
                TGeoBBox *box = dynamic_cast<TGeoBBox*>(volume->GetShape());
                if (box) {
                    double x = 2*tube->GetDX();
                    double y = 2*tube->GetDY();
                    pmtArea = x*y;
                    break;
                }
                abort();
            } while(false);

            // Find the solid angle for the pmt.  This uses the "small angle"
            // approximation and isn't right when the track is close to the
            // PMT.  To stop *serious* problems, this limits the covered area
            // to 2*pi (only half of total solid angle).
            double solidAngle = cangle*pmtArea/(4*unit::pi*distance*distance);
            if (solidAngle > 2*unit::pi) solidAngle = 2*unit::pi;

            // Estimate the mean number of photons generated at the point of
            // the energy deposition.
            double photons = fPhotonCollection*seg->GetEnergyDeposit();

            // Correct for solid angle
            photons *= solidAngle/(4*unit::pi);
            meanPhotons += photons;
            
            // Find the number of photons at a PMT.
            int nPhotons = gRandom->Poisson(photons);
            totalPhotons += nPhotons;
            
            if (nPhotons > 0) {
                contrib.push_back(seg);
                info.push_back(nPhotons);
            }

            for (int i=0; i< nPhotons; ++i) {
                double r = gRandom->Uniform();
                if (r < fShortFraction) {
                    double t = gRandom->Exp(fShortTime);
                    // This assumes that the UV index of refraction is 1.0
                    t += startT + distance/unit::c_light;
                    times.push_back(t);
                }
                else {
                    double t = gRandom->Exp(fLongTime);
                    // This assumes that the UV index of refraction is 1.0
                    t += startT + distance/unit::c_light;
                    times.push_back(t);
                }
            }
        }
    }

    // Add the effect of dark noise to the pmt.
    double deltaT = fStopSimulation - fStartSimulation;
    double meanNoise = deltaT * fPMTDarkCurrent;
    if (meanNoise > 0) {
        int noisePhotons = gRandom->Poisson(meanNoise);
        while (0 < noisePhotons--) {
            times.push_back(gRandom->Uniform(fStartSimulation,fStopSimulation));
        }
    }
    
    std::sort(times.begin(), times.end());

    return;
}

void CP::TElecSimple::DigitizeLight(
    CP::TEvent& ev, CP::TMCChannelId channel,
    const RealVector& input,
    const RealVector& triggers,
    const CP::TMCDigit::ContributorContainer& contrib,
    const CP::TMCDigit::InfoContainer& info) {

    // Get the digits container, and create it if it doesn't exist.
    CP::THandle<CP::TDigitContainer> digits
        = ev.Get<CP::TDigitContainer>("~/digits/pmt");
    if (!digits) {
        CP::THandle<CP::TDataVector> dv
            = ev.Get<CP::TDataVector>("~/digits");
        if (!dv) {
            CP::TDataVector* t = new CP::TDataVector("digits");
            ev.AddDatum(t);
            dv = ev.Get<CP::TDataVector>("~/digits");
        }
        CP::TDigitContainer* dg = new CP::TDigitContainer("pmt");
        dv->AddDatum(dg);
        digits = ev.Get<CP::TDigitContainer>("~/digits/pmt");
    }

    std::size_t pmtBins = (fStopSimulation-fStartSimulation)/fPMTStep;
    
    RealVector measuredCharge;
    if (measuredCharge.size() < pmtBins) {
        measuredCharge.resize(pmtBins);
    }

    // Add the signal for each photon.
    double signalWidth = 10.0*unit::ns;
    double sigNorm = 0.0;
    for (double v = 0.0; v<10.0*signalWidth; v+=fPMTStep) {
        double sig = v*std::exp(-1.0-v/signalWidth)/signalWidth;
        sigNorm += sig;
    }

    double totalSignal = 0.0;
    double gainSignal = 0.0;
    for (RealVector::const_iterator t = input.begin();
         t != input.end(); ++t) {
        double pulse = gRandom->Gaus(1.0,fPMTPeak);
        while (pulse < 0.1) pulse = gRandom->Gaus(1.0,fPMTPeak);
        double sigMax = 0.0;
        double dt = *t - fStartSimulation;
        for (double val = 0.0; val < 10.0*signalWidth; val += fPMTStep) {
            int bin = (dt + val)/fPMTStep;
            double v = dt + 2*val - bin*fPMTStep;
            double sig = v*std::exp(-1.0-v/signalWidth)/signalWidth;
            sig *= pulse/sigNorm;;
            totalSignal += sig;
            sig *= fAmplifierPMTGain;
            gainSignal += sig;
            sigMax = std::max(sigMax,sig);
            measuredCharge[bin] += sig;
        }
    }
    
    // double pedestal = fDigitPedestal + gRandom->Uniform(-0.5,0.5);
    double pedestal = fDigitPedestal + 0.5;

    for (RealVector::const_iterator trigger = triggers.begin();
         trigger != triggers.end(); ++trigger) {
        // This is the start time of the digitization window.
        double startTime = *trigger - fPMTPreTriggerTime;
        // This is the stop time of the digitization window.
        double stopTime = *trigger + fPMTPostTriggerTime;
        // This is the first bin in the digitization window
        int startBin = startTime/fPMTStep - fStartSimulation/fPMTStep;
        // This is the last bin in the digitization window
        int stopBin = startBin + (stopTime-startTime)/fPMTStep;
        if (startBin < 0) startBin = 0;
        if (stopBin > (int) measuredCharge.size()) stopBin = measuredCharge.size();
        // This is the first bin that might end up in a new digit.
        int lastStop = startBin;
        // The threshold to save a region of the FADC.  If this is negative,
        // then there isn't any zero suppression.
        double threshold = -5.0;
        // The number of bins to save before the threshold is crossed.
        int preThresholdBins = 1.0*unit::microsecond/fPMTStep;
        // The number of bins to save after the threshold is crossed.
        int postThresholdBins = 1.0*unit::microsecond/fPMTStep;
        // Find all of the possible digits in the digitization window.
        do {
            std::pair<int,int> digitRange = std::make_pair(lastStop,stopBin);

            if (threshold > 0) {
                // Find the first bin in the zero suppressed digit.  The tBin
                // variable is the bin where the threshold was crossed.
                int tBin;
                for (tBin=digitRange.first; tBin<digitRange.second; ++tBin) {
                    if (measuredCharge[tBin] > threshold) {
                        digitRange.first = std::max(digitRange.first,
                                                    tBin - preThresholdBins);
                        break;
                    }
                }
                // Check to make sure there is a non-zero digit.  If tBin is
                // at the end of the range, then there isn't a new digit.
                if (tBin == digitRange.second) {
                    digitRange.first = tBin;
                    break;
                }
                // Find the last bin in the zero suppressed digit.
                bool zeroRangeFound = false;
                do {
                    // Move to the first bin that is below threshold.
                    while (measuredCharge[tBin] > threshold
                           && tBin < digitRange.second) ++tBin;
                    // Check if there is a long section of zeros.  This allows
                    // room for the next preThresholdBins region.
                    zeroRangeFound = true;
                    for (int i = tBin;
                         i<tBin+postThresholdBins+preThresholdBins
                             && i<digitRange.second;
                         ++i) {
                        if (measuredCharge[i] > threshold) {
                            tBin = i;
                            zeroRangeFound = false;
                            break;
                        }
                    }
                } while (!zeroRangeFound);
                digitRange.second = std::min(digitRange.second,
                                             tBin + postThresholdBins);
            }

            CaptNamedInfo("Digitize", channel.AsString()
                          << " " << digitRange.first
                          << " " << digitRange.second);

            // Now copy to the output.  The new digit is between start and scan.
            double digitizedTotal = 0.0;
            CP::TPulseDigit::Vector adc;
            for (int bin = digitRange.first; bin < digitRange.second; ++bin) {
                double val = measuredCharge[bin];
                digitizedTotal += val;
                // Shift the baseline to the pedestal value.
                val += pedestal;
                int ival = val;
                ival = std::max(fDigitMinimum,std::min(ival,fDigitMaximum));
                adc.push_back(ival);
            }
            CP::TPulseMCDigit* digit
                = new TPulseMCDigit(channel,digitRange.first-startBin,
                                    adc,contrib,info);
            digits->push_back(digit);
            lastStop = digitRange.second;
        } while (lastStop < stopBin);
    }

}

// Calculating the induced current on the wire needs the potential if the
// wire was held at 1*unit::volt and everything else was grounded, as well
// as the electron velocity as a function of position.  That's a fairly
// complicated calculation.  It's simplified to assume that the electron
// travels in a straight line with an impact parameter of corrected
// distance, and at constant velocity.  The potential is simplified to be
// 1/r - 1/r_max where rmax is the distance to the wire as the electron
// passes the grid plane.  This estimates the shape of
// (velocity)*(electric field), and is normalized to so that the integral
// from 0.0 to gDist is 1.0.  Notice the sign is set so that this starts
// out positive (so that it matchs the behavior of the electronics).
double CP::TElecSimple::InducedShape(double dist, double impact) {
    // Distance from wire plane to grid.
    double gDist = 3.18*unit::mm;
    if (std::abs(dist) > gDist) return 0.0;
    // Magnitude of the electric field at a dist from the wire plane.
    double field = 1.0/std::sqrt(dist*dist + impact*impact)
        - 1.0/std::sqrt(gDist*gDist + impact*impact);
    // cosine of the angle between the velocity and the field vector
    double cosV = dist/sqrt(dist*dist + impact*impact);
    double current = - cosV*field;
    // Calculate the normalization (explicitly integrated)
    double a = sqrt(gDist*gDist+impact*impact);
    double b = std::log(gDist*gDist+impact*impact)
        -2.0*(std::log(impact)+1.0);
    double norm = 0.5*(a*b+2.0*impact)/a;
    return current/norm;
}

double CP::TElecSimple::PulseShaping(double tSample, double window,
                                     int samples) {
    double val = 0.0;
    double step = 1.0/samples;
    for (double t = tSample; t<tSample+window; t += step*window) {
        if (t<0.0) continue;
        double x = t/fAmplifierRise;
        if (x < 1.0) x = std::pow(x,fAmplifierRiseShape);
        else x = std::pow(x,fAmplifierFallShape);
        val += (x<40) ? x*std::exp(-x): 0.0;
    }
    return val;
}

/// The charge induced in a time bin by single electron interacting with a
/// wire.  This is just the current times the bin size.  For the collection
/// wires, this is simply the electron charge (i.e. 1.0) at tSample of zero
/// (possibly with small corrections for the impact distance).  For the
/// induction wires, the total time integral is 0.0.  The induction current is
/// calculated assuming the wire spacing is 3*mm and plane spacing is 3.18*mm.
/// The wireImpactDist is the impact distance assuming the electron path were
/// not distorted by the electric field (and ranges between 0*mm and 3*mm
double CP::TElecSimple::InducedCharge(bool isCollection,
                                      double wireImpactDist,
                                      double tSample,
                                      double tStep) {

    // The model assumes an infinite plane of wires between two "grid" planes
    // (which is applicable to both captain and minicaptain), so the induced
    // charge for an induction plain is symmetric around zero.  The induction
    // wire planes are assumed to be perfectly transparent (also a good
    // assumption), and the collection plane is perfectly opaque.

    // Estimate the distance to the wire at the time of the sample.  Don't
    // calculate an induced charge when the electron is past the grid plane.
    double dist = fDriftVelocity*tSample;
    if (dist > 3.18*unit::mm) return 0.0;
    if (dist < -3.18*unit::mm) return 0.0;

    // Make sure this electron actually passes in the vicinity of the wire.
    if (wireImpactDist>3.0*unit::mm) return 0.0;
    if (isCollection && wireImpactDist > 1.5*unit::mm) return 0.0;
    if (isCollection && tSample>0.0) return 0.0;

    // As the electron passes a wire, it induces current on the TWO closest
    // wires (the current on more distant wires is ignored).  As the electron
    // reachs the plane, the sum of the induced charge on the wires will be
    // equal to the charge of the electron.  The potential is adjusted so that
    // the closes drift passes within about 0.5 mm.
    double correctedImpact
        = 0.5*unit::mm + wireImpactDist*(2.0*unit::mm)/(3.0*unit::mm);

    // Estimate the fraction of the charge on the current wire.  The fraction
    // induced on this wire and the other wire sum to one.
    double fracWire = 1.0 - correctedImpact/3.0*unit::mm;
    double fracOther = correctedImpact/3.0*unit::mm;
#define SQUARE_AVERAGE
#ifdef SQUARE_AVERAGE
    // Crudely estimated to go as the square of the distance of closest
    // approach.  The EM calculation shouldn't be that hard, but I've not done
    // it.
    double charge = fracWire*fracWire/(fracWire*fracWire+fracOther*fracOther);
#else
    // Even more crudely estimated to go linearly as the distance of closest
    // approach.
    double charge = fracWire/(fracWire+fracOther);
#endif

    // Override the charge and corrected impact parameters for collection wires.
    if (isCollection) {
        correctedImpact = 0.01*unit::mm;
        charge = 1.0;
    }

    // Take care of the units.
    double norm = 1.0/fDriftVelocity;

    // Estimate the (normalized) charge in the current sample.
    double normCharge = InducedShape(dist,correctedImpact)*tStep/norm;

    // Return the induced charge in the current sample.
    return charge*normCharge;
}

double CP::TElecSimple::DriftCharge(
    CP::TEvent& event, CP::TMCChannelId channel,
    RealVector& out,
    CP::TMCDigit::ContributorContainer& contrib,
    CP::TMCDigit::InfoContainer& info) {

    // Get the drift hits.
    CP::THandle<CP::TG4HitContainer> g4Hits =
        event.Get<CP::TG4HitContainer>("truth/g4Hits/drift");
    if (!g4Hits) return false;

    if (channel.GetType() != 0) return false;
    int plane = channel.GetSequence();
    int wire = channel.GetNumber();

    // Move to the frame of the wire being simulated.
    CP::TManager::Get().Geometry(); // just in case...
    if (!CP::TManager::Get().GeomId().CdId(
            CP::GeomId::Captain::Wire(plane,wire))) return false;

    // Depending on the voltage applied to a plane, the wires will either
    // collect the drifting electrons, or have an induced current.  This is
    // true if the current wire will collect the charge.
    bool isCollection = false;
    if (plane == 0) isCollection = true;

    // Arrays for the transforms.
    double local[3];
    double master[3];
    double masterStart[3];
    double masterStop[3];

    // The fraction of the charge that is "collected".  This is normally 1 for
    // the collection wires and less for the induction wires.
    double collectionEfficiency = 1.0;

    // The zone around a wire where change is collected.
    double collectionZone = 3.0*unit::mm;
    if (isCollection) collectionZone = 1.5*unit::mm;

    // The maximum diffusion.
    double maxDiffusion = collectionZone + 3*unit::mm;

    // The weight to for each simulated electron.  If this is greater than
    // one, the electrons are under sampled so that the simulation can run
    // faster.
    double weight = 4.0;

    double startedElectrons = 0.0;
    double totalCharge = 0.0;
    // Check for every hit.
    for (CP::TG4HitContainer::const_iterator h = g4Hits->begin();
         h != g4Hits->end();
         ++h) {
        const CP::TG4HitSegment* seg
            = dynamic_cast<const CP::TG4HitSegment*>((*h));
        masterStart[0] = seg->GetStartX();
        masterStart[1] = seg->GetStartY();
        masterStart[2] = seg->GetStartZ();
        double startT = seg->GetStartT();
        masterStop[0] = seg->GetStopX();
        masterStop[1] = seg->GetStopY();
        masterStop[2] = seg->GetStopZ();
        double stopT = seg->GetStopT();

        for (int i=0; i<3; ++i) {
            master[i] = 0.5*masterStart[i] + 0.5*masterStop[i];
        }

        // Check if we should simulate this hit.  If it's too far from the
        // wire, just skip it.  The electrons drift from local positive Z
        // towards zero.  The length of the wire is along Y.
        gGeoManager->MasterToLocal(master,local);

        if (std::abs(local[0]) > maxDiffusion) {
            continue;
        }

        // This segment close, so find the mean number of quanta generated.
        double meanQuanta = seg->GetEnergyDeposit()/fActivationEnergy;

        // Check to see if this version of DETSIM is saving the energy
        // deposited as scintillation.  This handles the recombination.
        double nonIonizing = seg->GetSecondaryDeposit();

        // Estimate the mean number of electrons (first assuming we don't have
        // help from DETSIM, then overriding the value if DETSIM told us the
        // answer).
        double meanElectrons = (1-fRecombination)*meanQuanta;
        if (nonIonizing>0.0) {
            // There is help from DETSIM (by way of NEST), so use the
            // recombination fraction calculated by NEST for this particular
            // step.
            meanElectrons = (1.0-nonIonizing/meanQuanta/fActivationEnergy);
            meanElectrons *= meanQuanta;
        }

        // Fluctuate the number of electrons.
        int nElectrons = 0.5 + gRandom->Poisson(meanElectrons)/weight;

        startedElectrons += nElectrons;

        // Now simulate each electron...
        double timeStep = (fStopSimulation-fStartSimulation)/out.size();
        double segCharge = 0.0;
        for (int e = 0; e<nElectrons; ++e) {
            double v = gRandom->Uniform();
            // Find the location
            for (int i=0; i<3; ++i) {
                master[i] = v*masterStart[i] + (1.0-v)*masterStop[i];
            }
            double depositT = v*startT + (1.0 - v)*stopT;

            // The electrons drift from local positive Z towards zero.  The
            // length of the wire is along Y.
            gGeoManager->MasterToLocal(master,local);
            double driftDistance = local[2];

            // It's behind this drift plane, so assume no signal.  Not true,
            // but this is the "SIMPLE" electronics simulation.
            if (driftDistance < 0) continue;

            // Find the diffusion in position.
            double driftSigma
                = std::sqrt(2.0*fDiffusionCoeff*driftDistance/fDriftVelocity
                            + 0.001*unit::mm);

#ifdef FIX_DRIFT_DISTANCE
            // This fixes the drift and diffusion so that we don't have any
            // geometry effects while debugging.
            driftDistance = 1.0*unit::cm;
            driftSigma = 1E-7;
#endif

            // Find the drift time.  If necessary, do the diffusion time
            // spread.
            double driftTime = depositT + driftDistance/fDriftVelocity;
            if (driftSigma > 1E-6) {
                double timeSigma = driftSigma/fDriftVelocity;
                driftTime = gRandom->Gaus(driftTime,timeSigma);
            }

            // Remove electrons that don't survive to the wires.
            double driftLife = gRandom->Exp(fElectronLife);
            if (driftLife < driftTime) continue;

            // Find the distance to the wire.
            double wireDistance = local[0];

            if (driftSigma>1E-6) {
                wireDistance = gRandom->Gaus(wireDistance,driftSigma);
            }
            wireDistance = std::abs(wireDistance);

            // Charges outside of this zone are not collected at all.
            if (wireDistance > collectionZone) continue;

            // Add the electron to the collected charge.
            double deltaT = driftTime - fStartSimulation;
            int timeBin = deltaT/timeStep;
            if (timeBin < 0 || (int)out.size() <= timeBin) {
                CaptError("Drift out of time window " << driftTime/unit::ms
                          << " " << timeBin
                          << " " << out.size());
                continue;
            }

            // Add the induced charge to the output.
            if (isCollection) {
                // Handle the induced charge on the collection plane.
                int dBin = 3.18*unit::mm/fDriftVelocity/timeStep + 1;
                double integral = 0.0;
                double skewTime = wireDistance/fDriftVelocity;
                for (int i = timeBin-dBin;
                     i<timeBin+skewTime/timeStep+1;
                     ++i) {
                    if (i<0) continue;
                    if ((int) out.size() <= i) continue;
                    double charge = weight;
                    charge *= collectionEfficiency;
                    charge *= InducedCharge(isCollection,wireDistance,
                                            i*timeStep-deltaT-skewTime,
                                            timeStep);
                    out[i] += charge;
                    integral += charge;
                }
                totalCharge += integral;
                segCharge += integral;
            }
            else {
                // Handle the induced charge on the induction planes.
                int dBin = 3.18*unit::mm/fDriftVelocity/timeStep + 1;
                double integral = 0.0;
                for (int i = timeBin-dBin; i<=timeBin+dBin; ++i) {
                    if (i<0) continue;
                    if ((int) out.size() <= i) continue;
                    double charge = weight;
                    charge *= collectionEfficiency;
                    charge *= InducedCharge(isCollection,wireDistance,
                                            i*timeStep-deltaT,timeStep);
                    out[i] += charge;
                    integral += std::abs(charge);
                }
                totalCharge += integral/2.0;
                segCharge += integral/2.0;
            }
        }

        // This segment added charge to this wire, so add it to the
        // contributors.
        if (segCharge > 0) {
            contrib.push_back(seg);
            info.push_back(segCharge);
        }
    }

    return totalCharge;
}

void CP::TElecSimple::AddWireNoise(CP::TMCChannelId channel,
                                   RealVector& out) {
    // Add the wire nose.  This is before the electronics and covers "thermal"
    // noise from the wires.
    if (fWireNoise <= 0.1) return;

    // Translate the wire noise into electrons per bin and set the number of
    // samples that the noise is correlated over.  The number of correlated
    // samples is set by scaling, and should be parameterized so it can be fit
    // to the data.
    double timeStep = (fStopSimulation-fStartSimulation)/out.size();
    double sampleNoise = fWireNoise*std::sqrt(timeStep/(1000.0*unit::ns));

    for (RealVector::iterator o = out.begin(); o != out.end(); ++o) {
        (*o) += gRandom->Gaus(0,sampleNoise);
    }

}

namespace {
    double ColoredNoise(double freq,
                        double lowCut,
                        double highCut,
                        double alpha) {
        double mag = 1.0;
        double lowMag = 1.0;
        double highMag = 1.0;
        freq = std::abs(freq);
        lowMag = std::pow(lowCut,-alpha);
        highMag = std::pow(highCut,-alpha);
        if (freq < lowCut) {
            mag=lowMag/highMag;
        }
        else if (freq > highCut) {
            mag = highMag/highMag;
        }
        else {
            mag = std::pow(freq,-alpha)/highMag;
        }
        return mag;
    }

    /// Model the noise.  This started as a parallel rlc short for the noise,
    /// but the form has been modified to allow data to be fit.  This returns
    /// a value between 0.0 and 1.0 that is the equivalent to the RMS voltage
    /// between T1 and T2 when driven by a unit current source (e.g. the
    /// impedance between T1 and T2).
    ///
    ///      -------------------------T1
    ///                |
    ///                R
    ///                |
    ///           -----------
    ///           |    |    |
    ///           R'   L    C
    ///           |    |    |
    ///           |    r    r
    ///           |    |    |
    ///           -----------
    ///                |
    ///      -------------------------T2
    ///
    ///  All of the resistances can be different, and the inductance and
    ///  capitance are specified in units of frequency.  For this to actually
    ///  model the circuit show, indPow and capPow need to be 1.0.
    double ShortingImpedance(double freq,
                             double res,
                             double indFreq,
                             double indRes,
                             double indPow,
                             double capFreq,
                             double capRes,
                             double capPow) {
        freq = std::abs(freq);
        if (freq < 1.0*unit::hertz) freq = 1.0*unit::hertz;
        // Find impedance of the parallel rlc.
        double z = 1.0;
        z = 1.0
            /std::sqrt(
                1.0
                + std::pow(1.0/(indRes + std::pow(freq/indFreq,indPow))
                           - 1.0/(capRes + std::pow(capFreq/freq,capPow)),2.0));
        // Add in the series resistor
        if (res > 1.0) res = 1.0;
        if (res < 0.0) res = 0.0;
        z = res + (1.0-res)*z;
        return z;
    }
}

void CP::TElecSimple::GenerateBackgroundSpectrum(CP::TMCChannelId channel,
                                                 ComplexVector& inOut,
                                                 ComplexVector& work) {
    if (inOut.size() != work.size()) {
        CaptError("Wrong work vector size");
        throw;
    }

    // The sample size.
    double timeStep = (fStopSimulation-fStartSimulation)/work.size();

    // The nyquist frequency in hertz
    double nyquistFreq = (0.5/timeStep);

    // The frequency step per bin.
    double deltaFreq = 2.0*nyquistFreq/work.size();

    // Generate Gaussian noise.  This will be "sculpted" to have the right
    // power spectrum.
    for (std::size_t i=1; i<work.size()/2; ++i) {
        work[i] = std::complex<double>(gRandom->Gaus(),gRandom->Gaus());
        work[work.size()-i] = std::conj(work[i]);
    }
    work[0] = std::complex<double>(0.0,0.0);
    work[work.size()/2] = std::complex<double>(0.0,0.0);
    
    double impedMax = 0.0;
    double norm = 0.0;

    // Sculpt the power spectrum and find the peak.
    for (std::size_t i=1; i<work.size()/2; ++i) {
        double freq = deltaFreq*i;
        double mag = ColoredNoise(freq,
                                  fSpectralLowCut,
                                  fSpectralHighCut,
                                  fSpectralAlpha);
        
        double imp = ShortingImpedance(freq,
                                       fSpectralResist,
                                       fSpectralIndFreq,
                                       fSpectralIndRes,
                                       fSpectralIndPow,
                                       fSpectralCapFreq,
                                       fSpectralCapRes,
                                       fSpectralCapPow);
        
        if (imp > impedMax) {
            norm = mag*imp;
            impedMax = imp;
        }

        work[i] *= mag*imp;
        work[work.size()-i] = std::conj(work[i]);
    }

    // Adjust the normalization of the curve at the frequency with maximum
    // shorting impedance.
    if (norm > 0.0) norm = fSpectralNoise/norm;
    else {
        CaptError("Invalid normalization");
        throw;
    }

    for (std::size_t i=1; i<work.size()/2; ++i) {
        inOut[i] += norm*work[i];
        inOut[inOut.size()-i] = std::conj(inOut[i]); 
    }
}

void CP::TElecSimple::GenerateDiscreetSpectrum(CP::TMCChannelId channel,
                                               ComplexVector& inOut,
                                               ComplexVector& work) {

    // The sample size.
    double timeStep = (fStopSimulation-fStartSimulation)/work.size();

    // The nyquist frequency in hertz
    double nyquistFreq = (0.5/timeStep);

    // The frequency step per bin.
    double deltaFreq = 2.0*nyquistFreq/work.size();
    
    std::fill(work.begin(), work.end(), std::complex<double>(0.0,0.0));
    
    for (std::vector<CP::TElecSimple::DiscreetPeak>::iterator p
             = fDiscreetPeaks.begin();
         p != fDiscreetPeaks.end(); ++p) {
        double power = gRandom->Gaus(p->fPower,p->fPowerSigma);
        if (power<0.0) continue;
        power = p->fPower+p->fBackground;
        power = std::sqrt(power*power - p->fBackground*p->fBackground);
        power *= fDiscreetScale;
        power /= fDigitSlope;
        power /= std::sqrt(work.size()/fDigitOversample);
        double freq = p->fFrequency;
        std::size_t iFreq = freq/deltaFreq - 0.5;
        if (iFreq < 1 || iFreq > work.size()/2 - 1) {
            CaptError("Bad Freq " << freq << " " << deltaFreq
                      << " " << iFreq);
            continue;
        }
        double phase = gRandom->Uniform(0.0,2.0*3.14159);
#define SKIP_CAUCHY_PEAK true
        if (SKIP_CAUCHY_PEAK || p->fHalfWidth < fDigitOversample*deltaFreq) {
            // Peak is narrow, so inject a single frequency.
            work[iFreq] += std::polar(2.0*power,phase);
            work[work.size()-iFreq] = std::conj(work[iFreq]);
        }
#ifndef SKIP_CAUCHY_PEAK
        else {
            // Peak is wide, so inject a Cauchy.  This is skipped here since
            // the resulting model doesn't reproduce the data very well.
            double gamma = p->fHalfWidth/2.0;
            // Determine the normalization so that the power determines the
            // peak height.
            double norm = fDigitOversample;
            norm *= std::abs(ComplexCauchy(iFreq*deltaFreq,1.0,freq,gamma));
            norm = 1.0/norm;
            // Only add the Cauchy within 10 half widths of the peak.
            int halfWidths = 10;
            int bins = work.size();
            int iLow = iFreq - 1 - halfWidths*(p->fHalfWidth/deltaFreq);
            iLow = std:max(1,iLow);
            int iHigh
                = std::min(bins,int(iFreq + 1 + 10*(p->fHalfWidth/deltaFreq)));
            for (int i=iLow; i<iHigh; ++i) {
                double freq = i*deltaFreq;
                std::complex<double> v
                    = ComplexCauchy(freq,power,freq,gamma,phase);
                work[i] += norm*v;
                work[work.size()-i] = std::conj(work[i]);
            }
        }
#endif
    }
    
    for (std::size_t i=1; i<work.size()/2; ++i) {
        inOut[i] += work[i];
        inOut[inOut.size()-i] = std::conj(inOut[i]); 
    }
}

void CP::TElecSimple::GenerateResponseFFT(std::size_t samples) {
    // Find out the time step per simulated sample.
    double timeStep = (fStopSimulation-fStartSimulation)/samples;

    // Take the FFT of the response model.  This will allocate new FFTs if
    // it's required.
    if (samples != fResponseFFT.size()) {
        CaptLog("Create the response FFT");
        fResponseFFT.resize(samples);
        RealVector pulseShape(samples);
        // Find the normalization for the response
        double responseNorm = 0.0;
        for (std::size_t i = 0; i<samples; ++i) {
            pulseShape[i] = PulseShaping(i*timeStep, timeStep);
            if (fAmplifierConserveIntegral) {
                responseNorm += pulseShape[i];
            }
            else {
                responseNorm = std::max(responseNorm,pulseShape[i]);
            }
        }
        responseNorm /= std::sqrt(1.0*fResponseFFT.size());
        for (std::size_t i = 0; i<samples; ++i) pulseShape[i] /= responseNorm;
        FFT(pulseShape,fResponseFFT);
    }
}

void CP::TElecSimple::ShapeCharge(CP::TMCChannelId channel,
                                  const RealVector& in,
                                  ComplexVector& out) {
    if (out.size() != in.size()) {
        CaptError("Output vector does not match input size.");
        out.resize(in.size());
    }
    
    // Get the gain for this channel.
    double gain = fAmplifierCollectionGain;
    if (channel.GetType() == 0 && channel.GetSequence() != 0) {
        gain = fAmplifierInductionGain;
    }
        
    FFT(in,out);

    // Take the convolution in frequency space.
    for (std::size_t i=0; i<out.size(); ++i) {
        out[i] *= gain*fResponseFFT[i];
    }

}

void CP::TElecSimple::DigitizeWire(
    CP::TEvent& ev, CP::TMCChannelId channel,
    const RealVector& in,
    const RealVector& triggers,
    const CP::TMCDigit::ContributorContainer& contrib,
    const CP::TMCDigit::InfoContainer& info) {
    // Get the digits container, and create it if it doesn't exist.
    CP::THandle<CP::TDigitContainer> digits
        = ev.Get<CP::TDigitContainer>("~/digits/drift");
    if (!digits) {
        CP::THandle<CP::TDataVector> dv
            = ev.Get<CP::TDataVector>("~/digits");
        if (!dv) {
            CP::TDataVector* t = new CP::TDataVector("digits");
            ev.AddDatum(t);
            dv = ev.Get<CP::TDataVector>("~/digits");
        }
        CP::TDigitContainer* dg = new CP::TDigitContainer("drift");
        dv->AddDatum(dg);
        digits = ev.Get<CP::TDigitContainer>("~/digits/drift");
    }

    // Check for numeric problems
    for (RealVector::const_iterator i = in.begin(); i != in.end(); ++i) {
        if (!std::isfinite(*i)) {
            CaptError("INVALID NUMERIC SAMPLE");
        }
    }
    
    // Randomize the pedestal since it's going to be slightly different for
    // each channel.  This doesn't change the integral value of the pedestal
    // (it always rounds to the same value, but it does shift the probability
    // of the actual ADC mean.
    double pedestal = fDigitPedestal + gRandom->Uniform(-0.5, 0.5);
    
    // Calculate the time step for the input.
    double timeStep = (fStopSimulation-fStartSimulation)/in.size();

    for (RealVector::const_iterator trigger = triggers.begin();
         trigger != triggers.end(); ++trigger) {
        // This is the start time of the digitization window.
        double startTime = *trigger - fDigitPreTriggerTime;
        // This is the stop time of the digitization window.
        double stopTime = startTime
            + fDigitPreTriggerTime + fDigitPostTriggerTime;
        // This is the first bin in the digitization window
        int startBin = startTime/timeStep;
        if (startBin < 0) startBin = 0;
        // This is the last bin in the digitization window
        int stopBin = startBin + (stopTime-startTime)/timeStep;
        if (stopBin > (int) in.size()) stopBin = in.size();
        // This is the first bin that might end up in a new digit.
        int lastStop = startBin;
        // Find all of the possible digits in the digitization window.
        do {
            std::pair<int,int> digitRange
                = FindDigitRange(lastStop,startBin,stopBin,in);

            if (digitRange.first == digitRange.second) {
                CaptLog("No hit ");
                break;
            }

            // Now copy to the output.  The new digit is between start and scan.
            CP::TPulseDigit::Vector adc;
            int stride = (int) (fDigitStep/timeStep + 0.5);
            for (int bin = digitRange.first;
                 bin < digitRange.second-stride; bin += stride) {
                double val = 0;
                val = in[bin];
                val *= fDigitSlope;
                // Shift the baseline to the pedestal value.
                val += pedestal;
                if (fDigitNoise>0.01) val += gRandom->Gaus(0,fDigitNoise);
                int ival = val+0.5;
                ival = std::max(fDigitMinimum,std::min(ival,fDigitMaximum));
                if (ival<1) {
                    CaptWarn("Sample underflow " << adc.size() << " " << ival);
                    ival = 1;
                }
                else if (ival>4095) {
                    CaptError("Sample overflow " << adc.size() << " " << ival);
                    ival = 4095;
                }
                adc.push_back(ival);
            }

            CP::TPulseMCDigit* digit
                = new TPulseMCDigit(channel,
                                    (digitRange.first-startBin)/stride,
                                    adc,contrib,info);
            digits->push_back(digit);

            // Setup for the next (possible) digit.
            lastStop = digitRange.second;
        } while (lastStop < stopBin);

    }
}

std::pair<int, int>
CP::TElecSimple::FindDigitRange(int start,
                                int startBin,
                                int stopBin,
                                const RealVector& input) {
    if (fDigitThreshold < 0.0) return std::pair<int,int>(startBin,stopBin);

    // There is a threshold, so do the zero suppression.
    std::pair<int,int> digitRange = std::make_pair(start,stopBin);

    // Calculate the time step for the input.
    double timeStep = (fStopSimulation-fStartSimulation)/input.size();

    // Convert the threshold from ADC counts to mV (the input is the mV coming
    // out of the amplifiers.
    double threshold = fDigitThreshold/fDigitSlope;
    // The number of bins to save before the threshold is crossed.
    int preThresholdBins = fDigitPreThreshold/timeStep;
    // The number of bins to save after the threshold is crossed.
    int postThresholdBins = fDigitPostThreshold/timeStep;

    // Find the first bin in the zero suppressed digit.  The tBin
    // variable is the bin where the threshold was crossed.
    int tBin;
    for (tBin=digitRange.first; tBin<digitRange.second; ++tBin) {
        if (input[tBin] > threshold) {
            digitRange.first = std::max(digitRange.first,
                                        tBin - preThresholdBins);
            break;
        }
    }
    // Check to make sure there is a non-zero digit.  If tBin is
    // at the end of the range, then there isn't a new digit.
    if (tBin == stopBin) {
        digitRange.first = stopBin;
        digitRange.second = stopBin;
        return digitRange;
    }
    // Find the last bin in the zero suppressed digit.
    bool zeroRangeFound = false;
    do {
        // Move to the first bin that is below threshold.
        while (std::abs(input[tBin]) > threshold
               && tBin < digitRange.second) ++tBin;
        // Check if there is a long section of zeros.  This allows
        // room for the next preThresholdBins region.
        zeroRangeFound = true;
        for (int i = tBin;
             i<tBin+postThresholdBins+preThresholdBins
                 && i<digitRange.second;
             ++i) {
            if (std::abs(input[i]) > threshold) {
                tBin = i;
                zeroRangeFound = false;
                break;
            }
        }
    } while (!zeroRangeFound);
    digitRange.second = std::min(digitRange.second,
                                 tBin + postThresholdBins);

    return digitRange;
}

void CP::TElecSimple::OpenNoiseFile(std::string name) {
    fDiscreetPeaks.clear();
    if (name.empty()) return;

    std::ifstream input(name.c_str());
    if (!input.is_open()) {
        CaptError("Unable to read " << name);
        std::exit(1);
    }

    CaptLog("Reading noise from " << name);
    
    
    std::string line;
    while (std::getline(input,line)) {
        line = line.substr(0,line.find("#"));
        if (line.empty()) continue;
        std::istringstream parseLine(line);
        CP::TElecSimple::DiscreetPeak peak;
        parseLine >> peak.fIndex
                  >> peak.fFrequency
                  >> peak.fPower
                  >> peak.fPowerSigma
                  >> peak.fHalfWidth
                  >> peak.fBackground;
        peak.fFrequency *= unit::hertz;
        peak.fHalfWidth *= unit::hertz;
        fDiscreetPeaks.push_back(peak);
    }
}
