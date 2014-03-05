#include "TElecSimple.hxx"

#include <TEvent.hxx>
#include <TG4HitSegment.hxx>
#include <TCaptLog.hxx>
#include <TRuntimeParameters.hxx>
#include <HEPUnits.hxx>
#include <TGeomIdManager.hxx>
#include <TManager.hxx>
#include <TPulseMCDigit.hxx>
#include <TRealDatum.hxx>
#include <CaptGeomId.hxx>

#include <TGeoManager.h>
#include <TRandom.h>
#include <TVirtualFFT.h>

#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>
#include <memory>

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

    // The gain of the amplifier for the collection plane.  This must be
    // matched to the range of the ADC.
    fAmplifierCollectionGain
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.amplifier.gain.collection");
    fAmplifierCollectionGain *= unit::mV/unit::fC;

    // The gain of the amplifier for the induction planes.
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
            "elecSim.simple.digitization.pmtStep");

    // The time step for each digitization bin.
    fDigitStep
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.digitization.step");

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

    // The noise introduced during digitization
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

    fShortFraction 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.light.shortFraction");

    fShortTime 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.light.shortTime");

    fLongTime 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.light.longTime");

    fLightSensorCount 
        = CP::TRuntimeParameters::Get().GetParameterI(
            "elecSim.simple.light.sensors");

    // The wire noise level.
    double noise 
        = CP::TRuntimeParameters::Get().GetParameterD(
            "elecSim.simple.wire.noise");

    if (noise > 0) {
        fWireNoise = std::sqrt(noise*fDigitStep);
    }
    else {
        fWireNoise = 0.0;
    }
    CaptLog("Noise " << noise << " " << fWireNoise);

    fFFT = NULL;
    fInvertFFT = NULL;
}

CP::TElecSimple::~TElecSimple() {}

void CP::TElecSimple::operator()(CP::TEvent& event) {
    CaptLog("Event " << event.GetContext());

    // Check that the event has the truth hits.
    CP::THandle<CP::TDataVector> truthHits 
        = event.Get<CP::TDataVector>("truth/g4Hits");

    if (!truthHits) {
        CaptLog("No truth hits in event");
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
        CaptLog("Event is missing the TG4HitContainer objects");
        return;
    }

    CP::TManager::Get().Geometry();

    // Add the elecSim header to truth.
    AddElecSimHeader(event);

    // Figure out when the triggers will be.  This is usually "just zero".
    DoubleVector triggerTimes;
    GenerateTriggers(event,triggerTimes);

    if (triggerTimes.empty()) {
        CaptLog("No trigger in the event");
        return;
    }

    // Find out the total integration period.  This is found based off of the
    // times of the first and last triggers.
    fStartSimulation = 100*unit::second;
    fStopSimulation = -100*unit::second;
    for (DoubleVector::iterator triggerTime = triggerTimes.begin();
         triggerTime != triggerTimes.end();
         ++triggerTime) {
        fStartSimulation = std::min(fStartSimulation, 
                                     *triggerTime - fPreTriggerTime);
        fStopSimulation = std::max(fStopSimulation,
                                    *triggerTime + fPostTriggerTime);
    }
    
    // Simulate the (ganged) PMT signal.
    for (int i = 0; i<fLightSensorCount; ++i) {
        TMCChannelId pmt(1,0,i);
        DoubleVector photonTimes;
        LightSignal(event,pmt,photonTimes);
        DigitizeLight(event,pmt,photonTimes,triggerTimes);
    }

    int chargeBins = (fStopSimulation - fStartSimulation)/fDigitStep;
    chargeBins = 2*(1+chargeBins/2);
    DoubleVector collectedCharge(chargeBins);
    DoubleVector shapedCharge(chargeBins);

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
            if (!DriftCharge(event,channel,collectedCharge)) continue;
            AddWireNoise(channel,collectedCharge);
            ShapeCharge(channel,collectedCharge,shapedCharge);
            DigitizeWires(event,channel,shapedCharge,triggerTimes);
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
    std::auto_ptr<CP::TRealDatum> digitStep(new CP::TRealDatum("digitStep"));
    digitStep->clear();
    digitStep->push_back(fDigitStep);  // X Plane
    digitStep->push_back(fDigitStep);  // V Plane
    digitStep->push_back(fDigitStep);  // U Plane
    digitStep->push_back(fPMTStep);  // PMT
    header->AddDatum(digitStep.release());

    // Fill the pedestals
    std::auto_ptr<CP::TRealDatum> pedestal(new CP::TRealDatum("pedestal"));
    pedestal->clear();
    pedestal->push_back(fDigitPedestal);  // X Plane
    pedestal->push_back(fDigitPedestal);  // V Plane
    pedestal->push_back(fDigitPedestal);  // U Plane
    pedestal->push_back(fDigitPedestal);  // PMT
    header->AddDatum(pedestal.release());

    // Fill amplifier gains (voltage/(input charge))
    std::auto_ptr<CP::TRealDatum> gain(new CP::TRealDatum("gain"));
    gain->clear();
    gain->push_back(fAmplifierCollectionGain); // X Plane
    gain->push_back(fAmplifierInductionGain);  // V Plane
    gain->push_back(fAmplifierInductionGain);  // U Plane
    gain->push_back(fAmplifierPMTGain);        // PMT
    header->AddDatum(gain.release());

    // Fill the digitizer slopes (ADC/mV).
    std::auto_ptr<CP::TRealDatum> slope(new CP::TRealDatum("slope"));
    slope->clear();
    slope->push_back(fDigitSlope);  // X Plane
    slope->push_back(fDigitSlope);  // V Plane
    slope->push_back(fDigitSlope);  // U Plane
    slope->push_back(1.0);          // PMT
    header->AddDatum(slope.release());

    /// Fill the shape times.
    std::auto_ptr<CP::TRealDatum> shapeTime(new CP::TRealDatum("shape"));
    shapeTime->clear();
    shapeTime->push_back(fAmplifierRise);    // X Plane
    shapeTime->push_back(fAmplifierRise);    // V Plane
    shapeTime->push_back(fAmplifierRise);    // U Plane
    shapeTime->push_back(1.0);               // PMT
    header->AddDatum(shapeTime.release());

    /// Fill the drift velocity and electron lifetime.
    std::auto_ptr<CP::TRealDatum> argonState(new CP::TRealDatum("argon"));
    argonState->clear();
    argonState->push_back(fDriftVelocity);
    argonState->push_back(fElectronLife);
    header->AddDatum(argonState.release());

}

void CP::TElecSimple::GenerateTriggers(CP::TEvent& event,
                                       DoubleVector& triggers) {
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

void CP::TElecSimple::LightSignal(CP::TEvent& event,
                                  CP::TMCChannelId chan,
                                  DoubleVector& times) {

    times.clear();

    // Check that the event has the truth hits.
    CP::THandle<CP::TDataVector> truthHits 
        = event.Get<CP::TDataVector>("truth/g4Hits");

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
            double startT = seg->GetStartT();
            
            // Estimate the mean number of photons generated.
            double photons = fPhotonCollection*seg->GetEnergyDeposit();

            // Find the number of photons at a PMT.
            int nPhotons = gRandom->Poisson(photons);

            for (int i=0; i< nPhotons; ++i) {
                double r = gRandom->Uniform();
                if (r < fShortFraction) {
                    double t = gRandom->Exp(fShortTime);
                    t += startT;
                    t -= fStartSimulation;
                    times.push_back(t);
                }
                else {
                    double t = gRandom->Exp(fLongTime);
                    t += startT;
                    t -= fStartSimulation;
                    times.push_back(t);
                }
            }
        }
    }
}

void CP::TElecSimple::DigitizeLight(CP::TEvent& ev, CP::TMCChannelId channel,
                                    const DoubleVector& input,
                                    const DoubleVector& triggers) {

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

    int pmtBins = (fStopSimulation-fStartSimulation)/fPMTStep;
    DoubleVector shapedCharge(pmtBins);
    
    // Add the signal for each photon.
    double signalWidth = 10.0*unit::ns;
    for (DoubleVector::const_iterator t = input.begin();
         t != input.end(); ++t) {
        double pulse = gRandom->Gaus(1.0,fPMTPeak);
        while (pulse < 0.1) pulse = gRandom->Gaus(1.0,fPMTPeak);
        for (double val = 0.0; val < 10.0*signalWidth; val += fPMTStep) {
            int bin = (*t + val)/fPMTStep;
            double v = *t + 2*val - bin*fPMTStep;
            double sig = v*std::exp(-1.0-v/signalWidth)/signalWidth;
            shapedCharge[bin] += sig * pulse;
        }
    }

    double pedestal = fDigitPedestal;

    CP::TPulseDigit::Vector adc;
    for (DoubleVector::const_iterator trigger = triggers.begin();
         trigger != triggers.end(); ++trigger) {
        // This is the start time of the digitization window.
        double startTime = *trigger - fDigitPreTriggerTime;
        // This is the stop time of the digitization window.
        double stopTime = *trigger + fDigitPostTriggerTime;
        // This is the first bin in the digitization window
        int startBin = startTime/fPMTStep;
        if (startBin < 0) startBin = 0;
        // This is the last bin in the digitization window
        int stopBin = 1 + stopTime/fPMTStep;
        if (stopBin > (int) shapedCharge.size()) stopBin = shapedCharge.size();
        // Now copy to the output.  The new digit is between start and scan.
        CP::TPulseDigit::Vector adc;
        CP::TMCDigit::ContributorContainer contrib;
        for (int bin = startBin; bin < stopBin; ++bin) {
            double val = shapedCharge[bin]*fAmplifierPMTGain;
            // Shift the baseline to the pedestal value.
            val += pedestal;
            // Add the electronics noise.  This is from the electronics, and
            // therefore comes after the shaping.
            val += gRandom->Gaus(0.0,0.05*fAmplifierPMTGain);
            int ival = val;
            ival = std::max(fDigitMinimum,std::min(ival,fDigitMaximum));
            adc.push_back(ival);
        }
        CP::TPulseMCDigit* digit 
            = new TPulseMCDigit(channel,0,adc,contrib);
        digits->push_back(digit);
    }

}

bool CP::TElecSimple::DriftCharge(CP::TEvent& event,
                                  CP::TMCChannelId channel,
                                  DoubleVector& out) {

    for (DoubleVector::iterator t = out.begin(); t != out.end(); ++t) {
        *t = 0;
    }

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

    // Arrays for the transforms.
    double local[3];
    double master[3];
    double masterStart[3];
    double masterStop[3];

        
    // The weight to add for each electon.  This is setup so that the
    // effect of the weighting is "canceled" by the effect of the
    // digitization.  That means that the weight is related to 1/gain,
    // where the gain is the number of electrons per ADC digit.
#ifdef APPLY_WEIGHT
    double weight = std::max(fAmplifierCollectionGain, 
                             fAmplifierInductionGain);
    weight = 0.5*std::min(8.0,std::max(1.0,1/weight));
#else
    double weight = 1.0;
#endif

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

        // Define 15 mm as a long ways from the wire...
        if (std::abs(local[0]) > 3*unit::mm) {
            continue;
        }

        // This is close, so find the mean number of electrons generated.
        double electrons 
            = (1-fRecombination)*seg->GetEnergyDeposit()/fActivationEnergy;

        int nElectrons = 0.5 + gRandom->Gaus(electrons,sqrt(electrons))/weight;

        startedElectrons += nElectrons;

        // Now simulate each electron...
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
                            + 0.001);
            
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

            // BAD, BAD, BAD... I'm assuming the wire spacing is 3 mm.
            if (std::abs(wireDistance) > 1.5*unit::mm) continue;

            // Add the electron to the collected charge.
            double deltaT = driftTime - fStartSimulation;
            std::size_t timeBin = deltaT/fDigitStep;
            if (timeBin >= out.size()) {
                CaptError("Drift out of time window " << driftTime/unit::ms
                          << " " << timeBin
                          << " " << out.size());
                continue;
            }
            
            // BAD! BAD! BAD! I'm assuming that the signal doesn't vary
            // depending on how far the electron is from the wire.  This
            // matters for the induction wires.
            out[timeBin] += weight;
            totalCharge += weight;
        }
    }

    // A particle generates 1 electron per 34 eV of deposited energy, so a
    // typical MIP will generate ~5000 electrons per mm, so this is a very low
    // threshold.  If a wire sees less than 10, then it really didn't see
    // anything (ie ~350 eV of deposited energy).
    return (10 < totalCharge);
}

void CP::TElecSimple::AddWireNoise(CP::TMCChannelId channel, 
                                   DoubleVector& out) {
    if (fWireNoise <= 0.1) return;
    for (DoubleVector::iterator o = out.begin(); o != out.end(); ++o) {
        (*o) += gRandom->Gaus(0,fWireNoise);
    }
}

void CP::TElecSimple::ShapeCharge(CP::TMCChannelId channel,
                                  const DoubleVector& in,
                                  DoubleVector& out) {
    if (out.size() != in.size()) {
        CaptError("Output vector does not match input size.");
        out.resize(in.size());
    }

    bool bipolar = false;
    double gain = fAmplifierCollectionGain;
    if (channel.GetType() == 0 && channel.GetSequence() != 0) {
        bipolar = true;
        gain = fAmplifierInductionGain;
    }
    
    for (DoubleVector::iterator t = out.begin(); t != out.end(); ++t) {
        *t = 0;
    }

    // The normalization factor per transform;
    double norm = 1.0/std::sqrt(1.0*in.size());

    // Take the FFT of the response model.  This will allocate new FFTs if
    // it's required.
    if (in.size() != fResponseFFT.size()) {
        CaptLog("Initialize the FFT for convolutions");
        int len = in.size();
        if (fFFT) delete fFFT;
        fFFT = TVirtualFFT::FFT(1,&len,"R2C K");
        if (len != (int) in.size()) {
            CaptError("Invalid length for FFT");
            CaptError("     original length: " << in.size());
            CaptError("     allocated length: " << len);
        }
        len = in.size();
        if (fInvertFFT) delete fInvertFFT;
        CaptLog("Initialize the Inverted FFT for convolutions");
        fInvertFFT = TVirtualFFT::FFT(1,&len,"C2R K");
        if (len != (int) in.size()) {
            CaptError("Invalid length for inverse FFT");
            CaptError("     original length: " << in.size());
            CaptError("     allocated length: " << len);
        }
        CaptLog("Create the response FFT");
        fResponseFFT.resize(in.size());
        // Find the normalization for the response
        double responseNorm = 0.0;
        for (int i = 0; i<len; ++i) {
            double val = i*fDigitStep/fAmplifierRise;
            val = (val<40) ? fDigitStep*val*std::exp(-val)/fAmplifierRise: 0.0;
            responseNorm += val;
        }
        // Fill the response FFT.
        for (int i = 0; i<len; ++i) {
            double val = i*fDigitStep/fAmplifierRise;
            val = (val<40) ? fDigitStep*val*std::exp(-val)/fAmplifierRise: 0.0;
            fFFT->SetPoint(i,val/responseNorm);
        }
        fFFT->Transform();
        for (int i = 0; i<len; ++i) {
            double rl, im;
            fFFT->GetPointComplex(i,rl,im);
            fResponseFFT[i] = std::complex<double>(rl,im);
        }
        CaptLog("FFT initialized with " << len << " elements");
    }

    // Take the FFT of the input.
    double last = 0.0;
    for (std::size_t i = 0; i<in.size(); ++i) {
        double val = in[i];
        if (bipolar) val -= last;
        val *= gain;
        fFFT->SetPoint(i,val);
        last = in[i];
    }
    fFFT->Transform();

    // Take the convolution in frequency space and transform back to time.
    for (std::size_t i=0; i<in.size(); ++i) {
        double rl, im;
        fFFT->GetPointComplex(i,rl,im);
        std::complex<double> v(rl,im);
        v *= norm*fResponseFFT[i];
        fInvertFFT->SetPoint(i,v.real(),v.imag());
    }
    fInvertFFT->Transform();
    for (std::size_t i=0; i<out.size(); ++i) {
        out[i] = norm*fInvertFFT->GetPointReal(i);
    }

}

void CP::TElecSimple::DigitizeWires(CP::TEvent& ev, 
                                    CP::TMCChannelId chan,
                                    const DoubleVector& in,
                                    const DoubleVector& triggers) {
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

    // Randomize the pedestal since it's going to be slightly different for
    // each channel.  This doesn't change the integral value of the pedestal
    // (it always rounds to the same value, but it does shift the probability
    // of the actual ADC mean.
    double pedestal = fDigitPedestal + gRandom->Uniform(-0.5, 0.5);

    for (DoubleVector::const_iterator trigger = triggers.begin();
         trigger != triggers.end(); ++trigger) {
        // This is the start time of the digitization window.
        double startTime = *trigger - fDigitPreTriggerTime;
        // This is the stop time of the digitization window.
        double stopTime = *trigger + fDigitPostTriggerTime;
        // This is the first bin in the digitization window
        int startBin = startTime/fDigitStep;
        if (startBin < 0) startBin = 0;
        // This is the last bin in the digitization window
        int stopBin = 1 + stopTime/fDigitStep;
        if (stopBin > (int) in.size()) stopBin = in.size();
        // This is the first bin that might end up in a new digit.
        int lastStop = startBin;
        // Find all of the possible digits in the digitization window.
        do {
            std::pair<int,int> digitRange 
                = FindDigitRange(lastStop,startBin,stopBin,in);

            // Now copy to the output.  The new digit is between start and scan.
            CP::TPulseDigit::Vector adc;
            CP::TMCDigit::ContributorContainer contrib;
            for (int bin = digitRange.first; bin < digitRange.second; ++bin) {
                double val = in[bin]*fDigitSlope;
                // Shift the baseline to the pedestal value.
                val += pedestal;
                // Add the electronics noise.  This is from the electronics, and
                // therefore comes after the shaping.
                val += gRandom->Gaus(0.0,fDigitNoise);
                int ival = val;
                ival = std::max(fDigitMinimum,std::min(ival,fDigitMaximum));
                adc.push_back(ival);
            }
            CP::TPulseMCDigit* digit 
                = new TPulseMCDigit(chan,digitRange.first-startBin,adc,contrib);
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
                                const DoubleVector& input) {
    if (fDigitThreshold < 1) return std::pair<int,int>(startBin,stopBin);

    return std::pair<int,int>(startBin,stopBin);
}
