#include "TElecSimple.hxx"

#include <eventLoop.hxx>

#include <TROOT.h>

class TElecSimLoop: public CP::TEventLoopFunction {
public:
    TElecSimLoop() {
        fElecSim = NULL;
    }

    virtual ~TElecSimLoop() {};

    void Usage(void) {     }

    virtual bool SetOption(std::string option,std::string value="") {
        if (option=="peaks") {
            fPeaksFile = value;
            return true;
        }
        if (option=="noise") {
            fNoiseFile = value;
            return true;
        }
        return false;
    }

    bool operator () (CP::TEvent& event) {
        // Make sure the electronics simulated is created.
        if (!fElecSim) {
            fElecSim = new CP::TElecSimple();
            if (!fPeaksFile.empty()) {
                fElecSim->OpenPeaksFile(fPeaksFile);
            }
            if (!fNoiseFile.empty()) {
                fElecSim->OpenNoiseFile(fNoiseFile);
            }
        }

        // Run the simulation on the event.
        (*fElecSim)(event);

        // Save everything.
        return true;
    }

private:
    CP::TElecSimple* fElecSim;

    std::string fPeaksFile;
    std::string fNoiseFile;
};

int main(int argc, char **argv) {
    TElecSimLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
