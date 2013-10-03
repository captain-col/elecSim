#include "TPiecewiseConvolution.hxx"

#include <TCaptLog.hxx>

#include <tut.h>

#include <cmath>

namespace tut {
    struct baseTPiecewiseConvolution {
        // Run before each test.
        baseTPiecewiseConvolution() {
            fSignal.resize(6000);
            fResponse.resize(100);
            // Fill the signal.
            double peak = 0.5*fSignal.size();
            for (std::size_t i=0; i<fSignal.size(); ++i) {
                double value = 0.0;
                double bin = i;
                double delta = 0.5*(bin-peak);
                delta = delta*delta;
                if (delta < 20) value = std::exp(-delta);
                fSignal[i] = value;
            }
            // Fill the response
#ifdef DELTA_FUNCTION_TEST
            for (std::size_t i=0; i<fResponse.size()/2; ++i) {
                fResponse[i] = 0.0;
            }
            fResponse[0] = 1.0;
#else
            for (std::size_t i=0; i<fResponse.size()/2; ++i) {
                double delta = 0.1*i;
                delta = delta*delta;
                double value = 0.0;
                if (delta < 20) value = std::exp(-delta);
                if (i==0) {
                    fResponse[i] = value;
                }
                else { 
                    fResponse[i] = value;
                    fResponse[fResponse.size()-i] = value;
                }
            }
            double sum=0.0;
            for (std::size_t i=0; i<fResponse.size(); ++i) {
                sum += fResponse[i];
            }
            for (std::size_t i=0; i<fResponse.size(); ++i) {
                fResponse[i] /= sum;
            }
#endif
        }

        // Run after each test.
        ~baseTPiecewiseConvolution() {
        }

        std::vector<double> fSignal;
        std::vector<double> fResponse;
    };

    // Declare the test
    typedef 
    test_group<baseTPiecewiseConvolution>::object testTPiecewiseConvolution;
    
    test_group<baseTPiecewiseConvolution> 
    groupTPiecewiseConvolution("TPiecewiseConvolution");

    // Test P0D geometry identifiers.
    template<> template<>
    void testTPiecewiseConvolution::test<1> () {
        CaptLog("Run test 1");
        CP::TPiecewiseConvolution piecewise(fResponse);

        std::vector<double> convolution(fSignal.size());
        std::vector<double> deconvolution(fSignal.size());

        piecewise.Convolution(fSignal,convolution);
        piecewise.Deconvolution(convolution,deconvolution);
        
        for (std::size_t i=0; i<convolution.size(); ++i) {
            double delta = std::abs(deconvolution[i] - fSignal[i]);
            if (delta > 1e-4) {
                CaptError("Difference between signal and deconvolution "
                          << " " << i 
                          << " " << fSignal[i]
                          << " " << deconvolution[i]);
            }
            // ensure_distance(deconvolution[i], fSignal[i], 1e-4);
        }
    }
};

