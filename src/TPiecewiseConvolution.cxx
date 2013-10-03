#include "TPiecewiseConvolution.hxx"
#include "TCaptLog.hxx"

#include <TVirtualFFT.h>

void CP::TPiecewiseConvolution::Initialize(const std::vector<double>& response,
                                           const std::vector<double>& filter) {
    fResponseLength = std::max(response.size(), filter.size());
    
    int fftLength = 3*(1+fResponseLength);
    if (fftLength < 512) fftLength = 512;
    if (fftLength > 512) fftLength = 1024;
    if (fftLength > 1024) fftLength = 2048;

    fFFT = TVirtualFFT::FFT(1,&fftLength,"R2C M K");

    int ifftLength = fftLength;
    fInvertFFT = TVirtualFFT::FFT(1,&ifftLength,"C2R M K");

    // Make sure there wasn't a strange screwup creating the ffts.
    if (fftLength != ifftLength) {
        CaptError("Invalid FFT lengths for piecewise convolution");
    }

    // Set the size of the response work area.
    fResponseFFT.resize(fftLength);

    // Find the fft of the response and put it into the work area.
    CaptLog("Fill the response");
    for (std::size_t i=0; i<fResponseFFT.size(); ++i) {
        double value = 0.0;
        int j = fResponseFFT.size() - i;
        if (i < response.size()/2) {
            value = response[i];
        }
        else if (j < (int) response.size()/2) {
            value = response[response.size()-j];
        }
        fFFT->SetPoint(i,value);
    }
    fFFT->Transform();
    for (std::size_t i = 0; i<fResponseFFT.size(); ++i) {
        double rl, im;
        fFFT->GetPointComplex(i,rl,im);
        fResponseFFT[i] = std::complex<double>(rl,im);
    }

    // If there is a filter, then include it into the convolution function
    // FFT.
    if (filter.size() > 0) {
        CaptLog("Fill the filter");
        for (std::size_t i=0; i<fResponseFFT.size(); ++i) {
            double value = 0.0;
            int j = fResponseFFT.size() - i;
            if (i < filter.size()/2) {
                value = filter[i];
            }
            else if (j < (int) filter.size()/2) {
                value = filter[filter.size()-j];
            }
            fFFT->SetPoint(i,value);
        }
        fFFT->Transform();
        for (std::size_t i = 0; i<fResponseFFT.size(); ++i) {
            double rl, im;
            fFFT->GetPointComplex(i,rl,im);
            fResponseFFT[i] *= std::complex<double>(rl,im);
        }
    }
}

void CP::TPiecewiseConvolution::Convolution(const std::vector<double>& signal,
                                            std::vector<double>& result,
                                            bool convolution) {

    double norm = 1.0/fResponseFFT.size();

    std::size_t stride = fResponseFFT.size() - 2*fResponseLength;
    // The base is where to start convoluting the current chunk.
    for (std::size_t base = 0; base < signal.size(); 
         base += stride) {
        // Find the fft of the response and put it into the work area.  This
        // handles the wrap around and zeroing of the signal.
        for (std::size_t i = 0; i<fResponseFFT.size(); ++i) {
            std::size_t j = fResponseFFT.size() - i;
            double value = 0.0;
            if (i < stride + fResponseLength && i+base < signal.size()) {
                value = signal[base+i];
            }
            else if (j<base) {
                value = signal[base-j];
            }
            fFFT->SetPoint(i,value);
        }
        fFFT->Transform();

        // Take the convolution in frequency space and transform back to time.
        for (std::size_t i=0; i<fResponseFFT.size(); ++i) {
            double rl, im;
            fFFT->GetPointComplex(i,rl,im);
            std::complex<double> v(rl,im);
            if (convolution) {
                v *= fResponseFFT[i];
            }
            else {
                v /= fResponseFFT[i];
            }
            fInvertFFT->SetPoint(i,v.real(),v.imag());
        }
        fInvertFFT->Transform();

        for (std::size_t i=0; i<stride; ++i) {
            if (result.size() <= i+base) continue;
            result[i+base] = norm*fInvertFFT->GetPointReal(i);
        }
    }    
}
