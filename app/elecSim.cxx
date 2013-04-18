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
        return true;
    }

    bool operator () (CP::TEvent& event) {
        // Make sure the electronics simulated is created.
        if (!fElecSim) fElecSim = new CP::TElecSimple();

        // Run the simulation on the event.
        (*fElecSim)(event);

        // Save everything.
        return true;
    }

private:
    CP::TElecSimple* fElecSim;
};

int main(int argc, char **argv) {
    TElecSimLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
