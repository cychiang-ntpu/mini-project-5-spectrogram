# Spectrogram Generation from Audio Files

This project demonstrates how to generate spectrograms from audio files using Python. It utilizes the `librosa` library for audio processing and `matplotlib` for visualization.

## Prerequisites

Make sure you have the following libraries installed:

- `numpy`
- `matplotlib`
- `librosa`

You can install them using pip:

\`\`\`bash
pip install numpy matplotlib librosa
\`\`\`

## Usage

1.  **Prepare your audio files:** Place your `.wav` files in a directory.
2.  **Run the script:** The script reads audio files, computes their Short-Time Fourier Transform (STFT), and displays the spectrogram.

## Code Explanation

### 1. Loading Audio

The audio file is loaded using `librosa.load`, which returns the audio time series and the sampling rate.

\`\`\`python
y, sr = librosa.load(audio_path)
\`\`\`

### 2. Short-Time Fourier Transform (STFT)

The STFT is computed to convert the audio signal from the time domain to the time-frequency domain.

The formula for the discrete STFT is:

$$
X(m, k) = \sum_{n=0}^{N-1} x(n + mH) w(n) e^{-j \frac{2\pi}{N} kn}
$$

Where:
- $x(n)$ is the input signal.
- $w(n)$ is the window function (e.g., Hann window).
- $N$ is the frame length (FFT size).
- $H$ is the hop length.
- $X(m, k)$ is the STFT coefficient at time frame $m$ and frequency bin $k$.

In the code:

\`\`\`python
D = librosa.stft(y)
\`\`\`

### 3. Log-Amplitude Spectrogram

Since the dynamic range of audio signals is very large, we usually convert the amplitude to a decibel (log) scale for better visualization.

$$
S_{dB} = 10 \log_{10}(|X(m, k)|^2)
$$

Or effectively:

$$
S_{dB} = 20 \log_{10}(|X(m, k)|)
$$

In the code:

\`\`\`python
S_db = librosa.amplitude_to_db(np.abs(D), ref=np.max)
\`\`\`

### 4. Visualization

We use `librosa.display.specshow` to plot the spectrogram.

\`\`\`python
import matplotlib.pyplot as plt
import librosa.display

plt.figure(figsize=(10, 5))
librosa.display.specshow(S_db, sr=sr, x_axis='time', y_axis='log')
plt.colorbar(format='%+2.0f dB')
plt.title('Spectrogram')
plt.show()
\`\`\`

## References

- [Librosa Documentation](https://librosa.org/doc/latest/index.html)
- [Short-time Fourier transform - Wikipedia](https://en.wikipedia.org/wiki/Short-time_Fourier_transform)
