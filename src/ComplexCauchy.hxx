#ifndef ComplexCauchy_hxx_Seen
#define ComplexCauchy_hxx_Seen

#include <complex>

/// Return the complex value of the Cauchy function at a frequency "freq".
/// The power sets peak normalization (be careful, the Breit-Wigner
/// (i.e. Lorentzian) can be extremely narrow in real cases which makes the
/// plotted height "surprising").  The peak position and half-width is set
/// using "peak" and "gamma" (the usual definitions for a Lorentzian).  The
/// overall phase of the function is set using "phase" which has a default
/// value of zero.  As expected, this is defined for both positive and
/// negative frequency, so if you are using it to fill an FFT array, be sure
/// to get the frequency definitions right.  The usual FFT array frequency
/// definition where the bins run from 0 to nsize-1 are:
///
///        0: DC
///        1: 1/nsize/sample_length
///        2: 2/nsize/sample_length
///        3: 3/nsize/sample_length
///        ...
///        nsize/2-1: (nsize/2 - 1)/nsize/sample_length
///        nsize/2: nyquistFreq (2/sample_length)
///        nsize/2+1: -(nsize-(nsize/2 + 1))/nsize/sample_length
///        ...
///        nsize-3: -3/nsize/sample_length
///        nsize-2: -2/nsize/sample_length
///        nsize-1: -1/nsize/sample_length
///
std::complex<double> ComplexCauchy(double freq,
                                   double power,
                                   double peak,
                                   double gamma,
                                   double phase=0.0) {
    std::complex<double> val = (peak*peak
                                - freq*freq
                                + gamma*gamma/4.0
                                + std::complex<double>(0.0,gamma*freq));
    val = (gamma/2.0 - std::complex<double>(0.0,freq))/val;
    return power*val*std::polar(1.0,phase);
}
#endif
