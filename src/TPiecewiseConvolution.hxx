#ifndef TPiecewiseConvolution_hxx_Seen
#define TPiecewiseConvolution_hxx_Seen

namespace CP {
    class TPiecewiseConvolution;
}

class TVirtualFFT;

#include <vector>
#include <complex>

/// Do a convolution on a very long input sequence by breaking the sequence
/// into segments.  This is most useful when the input sequence is long
/// compared to the response function.  For instance, if the input sequence is
/// N measurements, and the response function is M measurements, then the
/// complexity for the full convolution is O(N log N), but the complexity for
/// the piecewise convolution is something like O(N log M). 
class CP::TPiecewiseConvolution {
public:

    /// Initialize the convolution.  This accepts the (time domain) response
    /// function.
    TPiecewiseConvolution(const std::vector<double>& response) {
        std::vector<double> dummyFilter;
        Initialize(response, dummyFilter);
    }

    /// Initialize the convolution.  This is used when there are two stages of
    /// the response.  An example of where this is used is when there is a
    /// response function for the wire, and a response function for the
    /// electronics.  An alternative way to think about this is that there
    /// will be a response function (what the wire does), and a filter (what
    /// the electronics do).
    TPiecewiseConvolution(const std::vector<double>& response, 
                          const std::vector<double>& filter) {
        Initialize(response, filter);
    }

    /// Do the convolution on a long sequence of data.  If the "convolution"
    /// boolean is true, then this does a convolution.  If the "convolution"
    /// boolean is false, this does a deconvolution.  There is a convenience
    /// function Deconvolution() that just calls this with convolution set to
    /// false, but it makes the users code more readable.
    void Convolution(const std::vector<double>& signal,
                     std::vector<double>& result,
                     bool convolution = true);

    /// Do the deconvolution of a long sequence of data.  See Convolution()
    /// for more information.
    void Deconvolution(const std::vector<double>& signal,
                       std::vector<double>& result) {
        Convolution(signal,result,false);
    }

private:
    void Initialize(const std::vector<double>& response, 
                    const std::vector<double>& filter); 

    /// The size of the input response and filter functions.  If both the
    /// response and filter are provided, then this will be the maximum of the
    /// two vector sizes.  The fResponseLength determines how much of the
    /// convolution is spoiled by wrap around.
    int fResponseLength;

    /// The complex FFT of the response (or response*filter).  This is the FFT
    /// of the full convolution function.  The vector is twice as long as the
    /// original response function.
    std::vector< std::complex<double> > fResponseFFT;

    /// The FFT that is used to do the piecewise deconvolution.  It's much
    /// shorter than the input signal length.
    TVirtualFFT *fFFT;

    /// The inverse FFT that is used to do the piecewise deconvolution.  It's
    /// much shorter than the input signal length.
    TVirtualFFT *fInvertFFT;
};
#endif
