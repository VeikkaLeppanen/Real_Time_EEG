#include "dsp.h"

// Function to compute autocorrelations up to a given lag
Eigen::VectorXd computeAutocorrelation(const Eigen::VectorXd& data, int maxLag, const std::string& norm) {
    Eigen::VectorXd autocorrelations(maxLag + 1);
    double mean = data.mean();
    double variance = (data.array() - mean).square().sum() / data.size();

    for (int lag = 0; lag <= maxLag; ++lag) {
        if (lag == 0) {
            autocorrelations(lag) = variance;  // Autocorrelation at lag 0 is the variance
        } else {
            double autocorrelation = 0.0;
            for (int i = 0; i < data.size() - lag; ++i) {
                autocorrelation += (data(i) - mean) * (data(i + lag) - mean);
            }
            if (norm == "biased") {
                autocorrelations(lag) = autocorrelation / data.size();
            } else if (norm == "unbiased") {
                autocorrelations(lag) = autocorrelation / (data.size() - lag);
            } else {
                throw std::invalid_argument("Norm must be either 'biased' or 'unbiased'");
            }
        }
    }

    return autocorrelations;
}

// Perform FFT using FFTW
std::vector<std::complex<double>> performFFT(const std::vector<double>& data) {
    int N = data.size();
    std::vector<std::complex<double>> out(N);
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out_fftw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; ++i) {
        in[i][0] = data[i];
        in[i][1] = 0;
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out_fftw, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int i = 0; i < N; ++i) {
        out[i] = std::complex<double>(out_fftw[i][0], out_fftw[i][1]);
    }

    fftw_free(in);
    fftw_free(out_fftw);

    return out;
}

// Eigen version of FFT
std::vector<std::complex<double>> performFFT(const Eigen::VectorXd& data) {
    int N = data.size();
    std::vector<std::complex<double>> out(N);
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out_fftw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; ++i) {
        in[i][0] = data(i);
        in[i][1] = 0;
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out_fftw, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int i = 0; i < N; ++i) {
        out[i] = std::complex<double>(out_fftw[i][0], out_fftw[i][1]);
    }

    fftw_free(in);
    fftw_free(out_fftw);
    
    return out;
}   

// Perform inverse FFT using FFTW
std::vector<std::complex<double>> performIFFT(const std::vector<std::complex<double>>& data) {
    int N = data.size();
    std::vector<std::complex<double>> out(N);
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out_fftw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; ++i) {
        in[i][0] = data[i].real();
        in[i][1] = data[i].imag();
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out_fftw, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int i = 0; i < N; ++i) {
        out[i] = std::complex<double>(out_fftw[i][0] / N, out_fftw[i][1] / N); // Normalize by N
    }

    fftw_free(in);
    fftw_free(out_fftw);

    return out;
}

// Apply the Hilbert transform
std::vector<std::complex<double>> hilbertTransform(const std::vector<double>& signal) {
    size_t N = signal.size();
    auto fft_signal = performFFT(signal);

    // Apply the Hilbert transform in the frequency domain
    fft_signal[0] = 0; // DC component zeroed
    for (size_t k = 1; k < N / 2; ++k) {
        fft_signal[k] *= 2;
    }

    if (N % 2 == 0) {
        fft_signal[N / 2] = 0; // Nyquist frequency zeroed if needed
    }

    for (size_t k = N / 2 + 1; k < N; ++k) {
        fft_signal[k] = 0; // Zero the negative frequencies
    }

    // Perform inverse FFT
    return performIFFT(fft_signal);
}

// Eigen version of the Hilbert transform
std::vector<std::complex<double>> hilbertTransform(const Eigen::VectorXd& signal) {
    size_t N = signal.size();
    auto fft_signal = performFFT(signal);

    // Apply the Hilbert transform in the frequency domain
    fft_signal[0] = 0; // DC component zeroed
    for (size_t k = 1; k < N / 2; ++k) {
        fft_signal[k] *= 2;
    }

    if (N % 2 == 0) {
        fft_signal[N / 2] = 0; // Nyquist frequency zeroed if needed
    }

    for (size_t k = N / 2 + 1; k < N; ++k) {
        fft_signal[k] = 0; // Zero the negative frequencies
    }

    // Perform inverse FFT
    return performIFFT(fft_signal);
}

Eigen::VectorXd hamming(unsigned int N) {
    Eigen::VectorXd h(N);
    double a0 = 0.54, a1 = 0.46;
    for (unsigned int i = 0; i < N; ++i) {
        h(i) = a0 - a1 * std::cos(2 * M_PI * i / (N - 1));
    }
    return h;
}


Eigen::VectorXcd spectrum(const Eigen::VectorXd& x, const Eigen::VectorXd& W) {
    int N = x.size();
    Eigen::VectorXcd Pxx(N);

    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int i = 0; i < N; ++i) {
        in[i][0] = x(i) * W(i); // Real part
        in[i][1] = 0;           // Imaginary part
    }

    fftw_execute(p);

    double wc = W.sum();
    for (int i = 0; i < N; ++i) {
        Pxx(i) = std::complex<double>(out[i][0], out[i][1]) / wc;
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return Pxx;
}


Eigen::MatrixXcd specgram_cx(const Eigen::VectorXd& x, unsigned int Nfft, unsigned int Noverl) {
    int N = x.size();
    int D = Nfft - Noverl;
    int U = std::max(1, static_cast<int>((N - Nfft) / D + 1)); // Calculate number of columns
    Eigen::MatrixXcd Pw(Nfft, U);
    Eigen::VectorXd W = hamming(Nfft);

    for (int k = 0, m = 0; k <= N - Nfft; k += D, ++m) {
        Eigen::VectorXd xk = x.segment(k, Nfft);
        Pw.col(m) = spectrum(xk, W);
    }

    if (N <= Nfft) {
        Eigen::VectorXd W = hamming(N);
        Pw.resize(N, 1);
        Pw.col(0) = spectrum(x, W);
    }

    return Pw;
}

Eigen::MatrixXd specgram(const Eigen::VectorXd& x, unsigned int Nfft, unsigned int Noverl) {
    Eigen::MatrixXcd Pw = specgram_cx(x, Nfft, Noverl);
    return (Pw.array() * Pw.conjugate().array()).real();
}

Eigen::VectorXd pwelch(const Eigen::VectorXd& x, unsigned int Nfft, unsigned int Noverl, bool doubleSided) {
    Eigen::MatrixXd Pxx = specgram(x, Nfft, Noverl);

    if (doubleSided) {
        return Pxx.rowwise().mean();
    } else {
        return Pxx.rowwise().mean().segment(0, Pxx.rows() / 2 + 1);
    }
}

double calculateSNR(const Eigen::VectorXd& data, int overlap, int nfft, double fs, double target_freq, double bandwidth) {
    Eigen::VectorXd psd = pwelch(data.array(), nfft, overlap);

    double df = fs / nfft; // Frequency resolution
    int target_index = static_cast<int>(target_freq / df);
    int half_bandwidth = static_cast<int>(bandwidth / (2 * df));

    double signal_power = 0.0;
    for (int i = target_index - half_bandwidth; i <= target_index + half_bandwidth; ++i) {
        if (i >= 0 && i < psd.size()) {
            signal_power += psd(i);
        }
    }

    double noise_power = 0.0;
    for (int i = 0; i < psd.size(); ++i) {
        if (i < target_index - half_bandwidth || i > target_index + half_bandwidth) {
            noise_power += psd(i);
        }
    }
    noise_power /= (psd.size() - 2 * half_bandwidth);

    return 10 * std::log10(signal_power / noise_power); // SNR in dB
}

double ang_diff(double x, double y) {
    std::complex<double> result = std::exp(std::complex<double>(0, x)) / std::exp(std::complex<double>(0, y));
    return std::arg(result);
}

Eigen::VectorXd ang_diff(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Input vectors must have the same size.");
    }

    Eigen::VectorXd result(x.size());
    for (int i = 0; i < x.size(); ++i) {
        std::complex<double> num = std::exp(std::complex<double>(0, x[i]));
        std::complex<double> den = std::exp(std::complex<double>(0, y[i]));
        std::complex<double> quotient = num / den;
        result[i] = std::arg(quotient);
    }
    return result;
}