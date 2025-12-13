---
title: 'mini project-5: Spectrogram'

---

# mini project-5: Spectrogram 
> 由 https://hackmd.io/@ntpu-ce-mmsp/mmsp-2025-mini-project-5 同步

## 1. Objectives
* Be familiar with Discrete Fourier Transform (DFT).
* Understand how a spectrogram is generated.
* Mastery over the C language.
* Start to use popular python library “Matplotlib.”

## 2. Background
Fig.1 shows waveforms of a compound vowel **‘ai’** and a syllable **‘ya.’** The compound vowel **‘ai’** and the syllable **‘ya’** have similar sound components **‘a’** and **‘i’**, but are in a different order over time. It can be observed that the waveform at the beginning of ‘ai’ and the one at the end of **‘ya’** are very similar; that is, they may have similar harmonic components, since they all are the sound **‘i’** (vice versa for the case of the end of **‘ai’** vs. the beginning of **‘ya.’** These two waveforms are different speech segment. But, if you take Fourier transform for **whole segments of ‘ai’ and ‘ya,’** they may have very similar frequency responses in Discrete-Time Fourier Transforms (DTFTs), indicating that we cannot distinguish the two sounds **‘ai’** and **‘ya’** if the analysis window is too long to make Fourier analysis lose **time resolution** because speech signals change over time.  

![image](https://hackmd.io/_uploads/HyVNc3sqee.png)


To catch the change in frequency components over time, as shown in Fig. 2, **Short-Time Fourier Analysis or Short-Time Fourier Transform (STFT)** is adopted to analyze a speech signal by decomposing the speech signal into a series of short segments, referred to as **analysis windows/frames**, and analyze each one independently. The waveform in each short segment (analysis window) could be stationary since waveforms in a shorter window look more periodic and have more similar harmonic components. Also, each analysis frame could catch a different waveform to observe signal changing over time. We, therefore, can analyze speech phonemes and their transition by spectrogram.


A spectrogram is a visual representation of a **spectrum sequence**, capturing spectrum changes with time. Spectrograms are sometimes called sonogram, spectral waterfalls, voiceprints, or voicegrams. In Digital Signal Processing (DSP), we can regard spectrogram, $X[m,k]$, as a function of time (frame index) and frequency bin, i.e.,

$$ X[m,k]=20\log_{10}|\sum_{n=0}^{N-1} x_m[n]e^{-j\frac{2\pi kn}{N}}\tag1$$

where $m$ represents frame index; $k$ represents frequency bin index. Note that $\sum_{n=0}^{N-1} x_m[n]e^{-j\frac{2\pi kn}{N}}$ in Eq. (1) is in the same form with Discrete Fourier Transform (DFT) and Discrete Fourier Series (DFS). $x_m[n]$ is a short-time signal with a length of $N$ samples, indicating taking DFT of a discrete time signal of length $N$. We usually call the parameter $N$ ‘FFT length’. Note that FFT stands for Fast Fourier Transform which is an efficient implantation of DFT. In mathematics, $x_m[n]$ is expressed by

$$x_m[n]=x[m \cdot M+n]\cdot w[n], \qquad 0\leq n < N \tag2$$

where $x[m\cdot M+n]$ denotes a time-shifted signal of original input speech,$x[n]$; $M$ represents frame interval (or usually "hop length") in sample; $w[n]$ is a window function, a mathematical function with zero-valued outside of some chosen interval (usually called a compact support). Therefore, we can regard $x_m[n]$ as $m$-th frame signal of length $N$. There are many types of window functions, such as rectangular window (Eq. 3), Hamming window (Eq. 4), and Hann window (Eq. 5).

<br>
<br>

![image](https://hackmd.io/_uploads/HJnw9hj9eg.png)
##### Fig.2: Overview of obtaining spectrogram. The process in sequence, includes framing, taking window, zero padding, and DFT or FFT.


<br>

$$\text{Rectangular window: } w[n] = \begin{cases} 1, & 0 \leq n < P \\ 0, & \text{otherwise} \end{cases}\tag3$$

$$\text{Hamming window: } w[n] = \begin{cases} 0.54-0.46cos(\frac{2\pi n}{P-1}), & 0 \leq n < P \\ 0, & \text{otherwise} \end{cases}\tag4$$

$$\text{Hann window: } w[n] = \begin{cases} 0.5-0.5cos(\frac{2\pi n}{P-1}), & 0 \leq n < P \\ 0, & \text{otherwise} \end{cases}\tag5$$

For analysis of speech, we usually choose Hamming window as the window function in Eq. 2. Note that in Eqs 3-5, a window function has zero-value outside of an interval $0\leq n < P$ . We therefore define $P$ as **analysis window size**, regarding how many sample points are excerpted from the analyzed signal $x[n]$ for an analysis window. Thus, the value $P$ must be equal or smaller than the FFT window length $N$. An analyzed signal is said to be zero-padded as $P<N$  because $x_m[n] = x[m \cdot M + n]\cdot w[n]$ is zero for $P\leq n <N$ . Since we only have non-zero samples for $0\leq n<P$, only $P$ points excerpted from $x[n]$ are informative for frequency analysis. As FFT window length $N$ is set to be greater than analysis window size $P$, the $N-P$ zero-padded samples in $x_m[n] = x[m \cdot M + n]\cdot w[n]$ just let the spectrum of m-th frame $X[m,k]$ has a smoother outline extrapolating more frequency bins.

## 3. Parameters of spectrogram by DFT
### 3.1. Analysis Window Types
**Why not using rectangular window for spectrogram?** Assume we have a pure sinusoidal signal to be analyzed by spectrogram, e.g.

$$x[n]=\cos(\omega_0 n)\tag6$$

and the corresponding Fourier transform is

$$X(e^{j\omega})=\sum_{k=-\infty}^{\infty}\pi \delta(\omega -\omega_0+2\pi k)+\pi \delta(\omega+\omega_0+2\pi k)\tag7$$

But that spectrogram analyzes the windowed signal $x_m[n] = x[m \cdot M + n]\cdot w[n]$ and this windowed signal is a distorted version of original pure sinusoidal signal $x[n]$. Also, let’s recall the convolution theorem, i.e.

$$\begin{align}
x_0[n] &= x[n] \cdot w[n]\overset{FT}{\longrightarrow} \; X_0(e^{j\omega}) = \frac{1}{2\pi} X(e^{j\omega}) \ast W(e^{j\omega}) \\[1em]
        &= \frac{1}{2\pi} \int_{\omega=-\pi}^{\pi} X(e^{j\theta}) W(e^{j(\omega - \theta)}) \, d\theta\tag8
\end{align}
$$

From Eq. (8), we know that the DTFT $X_0(e^{j\omega})$ depends on the DTFT $W(e^{j\omega})$ (DTFT of the window function $w[n]$). So, let’s check the $W(e^{j\omega})$ if $w[n]$ is a rectangular function of length $P$:

$$ \begin{align}
W(e^{j\omega}) &= \sum_{n=0}^{P-1} w[n] e^{-j\omega n} = \sum_{n=0}^{P-1} 1\cdot e^{-j\omega n} = \frac{1 - e^{-j\omega P}}{1 - e^{-j\omega}} \\[1em]
                &=\frac{e^{-j0.5\omega P} \left( e^{+j0.5\omega P} - e^{-j0.5\omega P} \right)}{e^{-j0.5\omega} \left( e^{+j0.5\omega } - e^{-j0.5\omega } \right)} \\[1em]
                &= e^{-j0.5\omega(P-1)} \frac{\sin(0.5\omega P)}{\sin(0.5\omega)} \tag9
\end{align}
$$

Apparently, $W(e^{j\omega})$ is in a form of a sinc function.<span style="color: red;">**So...Think about the rest by yourself. I talked too much!**</span>

### 3.2. The Releationship between Analysis Window Length $P$ and FFT Window Length $N$ 

Does $N$ or $P$ control "Analysis Bandwidth"? The bandwidth is:

$$bw=\frac{k}{N}f_s \text{ (unit: Hz)}$$

or

$$bw=\frac{k}{P}f_s \text{ (unit: Hz)}$$

where $f_s$ is sampling frequency. Please find the answer by yourself.

### 3.3. Frame Interval M
The parameter $M$ is the **frame interval** or usually **frame hop length** regarding rate of DFT analysis on time axis. Smaller $M$ means more detailed time resolution for spectrogram. Conventionally, $M$ is set to be equal or smaller than $N$, indicating each analysis frame $x_m[n]$ could have some overlapping samples to the neighboring analysis frames $(e.g. \quad x_{m-?}[n]...x_{m-1}[n],\quad x_{m+1}[n],...x_{m+?}[n]).$

### 3.4. Typical Settings
Typically, **analysis window length** (regarding $P$) is set to be 4ms~40ms, that is $P\in[32,320]$ for an 8kHZ sample-rate speech signal, or $P\in[64,640]$ for a 16kHz sample rate signal. FFT window length is typically set to be $N=2^{\lceil \log_2(P) \rceil}$. **Frame Hop Length** is usually set to be 5ms~20ms, that is $M\in[40,160]$ for an 8kHZ sample-rate speech signal, or $M\in[80,320]$ for a 16kHz sample rate signal.

## 4. Implementations
### 4.1. Generate two waveforms (*.wav) that satisfy the following specifications by `signal_gen.c`:
* Sample rates: 16,000 Hz or 8,000 Hz
* Bit depth: 16 bits
* Digital signal that represents
$$x(t) = \sum_{j=0}^{3} \sum_{i=0}^{9} \left\{ a_i \left[ u(t - 0.1i - j) - u(t - 0.1(i+1) - j) \right] s_j(t - 0.1i - j, f_i) \right\}
$$
    * $u(t)$ is the unit step function,
    * $\left\{a_i\right\}|_i = 0, 1, 2, \dots, 9 = \{ 100, 2000, 1000, 500, 250, 100, 2000, 1000, 500, 250 \}$ are amplitudes
    * $\left\{f_i\right\}|_i = 0, 1, 2, \dots, 9 = \{ 0,31.25,500,2000,4000,44,220,440,1760,3960 \}$ are frequencies in Hz;
    * $s_j(t,f_i)$ for $j$ = 0,1,2,3 are the waveforms generated by the function generator:
$$
s_j(t, f_i) = \begin{cases} 
\sin(2\pi f_i t), & \text{for } j = 0 \quad \text{(sine wave)} \\[1em]
f_i t - \lfloor f_i t \rfloor, & \text{for } j = 1 \quad \text{(sawtooth wave)} \\[1em]
\text{sgn}(\sin(2\pi f_i t)), & \text{for } j = 2 \quad \text{(square wave)} \\[1em]
2 \left| 2(f_i t - \lfloor f_i t + 0.5 \rfloor) \right| - 1, & \text{for } j = 3 \quad \text{(triangle wave)}
\end{cases}
$$
* Please name the two wave files: `s-8kHz.wav` and `s-16kHz.wav`. 





### 4.2. Generate 16 spectrograms $X[m,k]$ by `spectrogram.c`
The usage of `spectrogram.c` is:
```
spectrogram w_size w_type dft_size f_itv wav_in spec_out
```
where 
* `w_size`: analysis window size (unit: millisecond)
* `w_type`: a string to be “hamming” or “rectangular”
* `dft_size`: DFT/FFT window size (unit: millisecond)
* `f_itv`: frame interval (unit: millisecond)
* `wav_in`: nput WAVE file
* `spec_out`: output spectrigram data (ascii with 15 decimal places)

The 16 spectrograms are generated by the combinations of the following four window settings and the four waveform files:
### Four Window Settings:

* Setting 1:
    * Analysis window size = 32ms
    * Analysis window type = rectangular
    * DFT/FFT window size = 32ms
    * Frame interval = 10ms
* Setting 2:
    * Analysis window size = 32ms
    * Analysis window type = hamming
    * DFT/FFT window size = 32ms
    * Frame interval = 10ms
* Setting 3:
    * Analysis window size = 30ms
    * Analysis window type = rectangular
    * DFT/FFT window size = 32ms
    * Frame interval = 10ms
* Setting 4:
    * Analysis window size = 30ms
    * Analysis window type = hamming
    * DFT/FFT window size = 32ms
    * Frame interval = 10ms

### Four Waveform Files:
1. `aeueo-8kHz.wav` (download from https://github.com/cychiang-ntpu/mini-project-5-spectrogram/blob/main/aeueo-8Hz.wav): it contains 5 basic vowels with 4 tones in 8kHz.
2. `aeueo-16kHz.wav`(download from https://github.com/cychiang-ntpu/mini-project-5-spectrogram/blob/main/aeueo-16kHz.wav): it contain 5 basic vowels with 4 tones in 16kHz.
3. `s-8kHz.wav` as generated by `signal_gen.c`.
4. `s-16kHz.wav` as generated by `signal_gen.c`.


Therefore, 16 ascii files are saved: s-16k.{Set1~Set4}.txt, s-8k.{Set1~Set4}.txt, aeueo-16kHz.{Set1~Set4}.txt, and aeueo-8kHz.{Set1~Set4}.txt,

### 4.3. Generate the 16 PDFs of Spectrograms by `spectshow.py.`
The usage of “`spectshow.py`” is :
```
python3 spectshow.py in_wav in_txt out_pdf
```
where
* `in_wav`: input WAVE file
* `in_txt`: input text file that contains spectrogram data in ascii
* `out_pdf`: a pdf file that shows the figures of a waveform and the corresponding spectrogram.

You can use Matplotlib (https://matplotlib.org/) library to generate the figures. The pdf files opened should look like Figs. 3 and 4.

![image](https://hackmd.io/_uploads/rJgj53oqxg.png)
Fig. 3: The waveform and spectrogram for aeueo-8kHz.wav.

<br><br>
![image](https://hackmd.io/_uploads/rks6c3scle.png)
Fig. 4: The waveform and spectrogram for aeueo-16kHz.wav.

### 4.3. Write GitHub a Workflow
The workflow need to
1. compile and build source codes,
2. run programs that the generate waveforms,  the 16 spectrogram ascii files, and save the pdf files of the spectrograms, and
3. the artifacts include the generated waveforms for $x(t)$, the 16 spectrogram ascii txt files, and 16 pdf files that display waveforms and spectrogram.

### 4.4. Write a `README.md` File
The README.ms file consisted of the following contents:
1. Division of labor
2. Compare the results by Settings 1-4 and discuss the differences and their significance.
3. Calculate how many multiplications and additions are executed for each frame. 
4. Thoughts and Reflections

### 4.5. Push the Above Files to GitHub
Please provide the link of the repo in [lms3](https://lms3.ntpu.edu.tw/).

## 5. Important Dates
* submission due: YYYY-MM-DD=2025-12-21
* notification of the 1st evaluation: 2025-12-23
* resubmission due: 2025-12-28
* notification of the 2nd evaluation: 2026-1-2

## 6. Grading policy
* The assignments submitted for the first time before the submission due are scored from 50 to 100 points if the acompanying codes can be compiled successfully, or the assignments need to be resubmitted before the resubmission due.
* The assignments submitted after the submission due and before the resubmission due are scored from 40 to 70 points if the accompanying codes are compiled unsuccessfully in the 1st evaluation.
* Students are encouraged to resubmit their codes for the 2nd evaluation to get higher scores if they improve the codes submitted for the 1st evaluation.


## References:
[1]	Xuedong Huang, Alex Acero, Hsiao-Wuen Hon (2001). Spoken Language Processing: a guide to theory, algorithm, and system development, page 274-281. Prentice Hall

[2]	https://www.speech.kth.se/wavesurfer/# mini project-5: Spectrogram 
> 由 https://hackmd.io/@ntpu-ce-mmsp/mmsp-2025-mini-project-5 同步

## 1. Objectives
* Be familiar with Discrete Fourier Transform (DFT).
* Understand how a spectrogram is generated.
* Mastery over the C language.
* Start to use popular python library “Matplotlib.”

## 2. Background
Fig.1 shows waveforms of a compound vowel **‘ai’** and a syllable **‘ya.’** The compound vowel **‘ai’** and the syllable **‘ya’** have similar sound components **‘a’** and **‘i’**, but are in a different order over time. It can be observed that the waveform at the beginning of ‘ai’ and the one at the end of **‘ya’** are very similar; that is, they may have similar harmonic components, since they all are the sound **‘i’** (vice versa for the case of the end of **‘ai’** vs. the beginning of **‘ya.’** These two waveforms are different speech segment. But, if you take Fourier transform for **whole segments of ‘ai’ and ‘ya,’** they may have very similar frequency responses in Discrete-Time Fourier Transforms (DTFTs), indicating that we cannot distinguish the two sounds **‘ai’** and **‘ya’** if the analysis window is too long to make Fourier analysis lose **time resolution** because speech signals change over time.  

![image](https://hackmd.io/_uploads/HyVNc3sqee.png)


To catch the change in frequency components over time, as shown in Fig. 2, **Short-Time Fourier Analysis or Short-Time Fourier Transform (STFT)** is adopted to analyze a speech signal by decomposing the speech signal into a series of short segments, referred to as **analysis windows/frames**, and analyze each one independently. The waveform in each short segment (analysis window) could be stationary since waveforms in a shorter window look more periodic and have more similar harmonic components. Also, each analysis frame could catch a different waveform to observe signal changing over time. We, therefore, can analyze speech phonemes and their transition by spectrogram.


A spectrogram is a visual representation of a **spectrum sequence**, capturing spectrum changes with time. Spectrograms are sometimes called sonogram, spectral waterfalls, voiceprints, or voicegrams. In Digital Signal Processing (DSP), we can regard spectrogram, $X[m,k]$, as a function of time (frame index) and frequency bin, i.e.,

$$ X[m,k]=20\log_{10}|\sum_{n=0}^{N-1} x_m[n]e^{-j\frac{2\pi kn}{N}}\tag1$$

where $m$ represents frame index; $k$ represents frequency bin index. Note that $\sum_{n=0}^{N-1} x_m[n]e^{-j\frac{2\pi kn}{N}}$ in Eq. (1) is in the same form with Discrete Fourier Transform (DFT) and Discrete Fourier Series (DFS). $x_m[n]$ is a short-time signal with a length of $N$ samples, indicating taking DFT of a discrete time signal of length $N$. We usually call the parameter $N$ ‘FFT length’. Note that FFT stands for Fast Fourier Transform which is an efficient implantation of DFT. In mathematics, $x_m[n]$ is expressed by

$$x_m[n]=x[m \cdot M+n]\cdot w[n], \qquad 0\leq n < N \tag2$$

where $x[m\cdot M+n]$ denotes a time-shifted signal of original input speech,$x[n]$; $M$ represents frame interval (or usually "hop length") in sample; $w[n]$ is a window function, a mathematical function with zero-valued outside of some chosen interval (usually called a compact support). Therefore, we can regard $x_m[n]$ as $m$-th frame signal of length $N$. There are many types of window functions, such as rectangular window (Eq. 3), Hamming window (Eq. 4), and Hann window (Eq. 5).

<br>
<br>

![image](https://hackmd.io/_uploads/HJnw9hj9eg.png)
##### Fig.2: Overview of obtaining spectrogram. The process in sequence, includes framing, taking window, zero padding, and DFT or FFT.


<br>

$$\text{Rectangular window: } w[n] = \begin{cases} 1, & 0 \leq n < P \\ 0, & \text{otherwise} \end{cases}\tag3$$

$$\text{Hamming window: } w[n] = \begin{cases} 0.54-0.46cos(\frac{2\pi n}{P-1}), & 0 \leq n < P \\ 0, & \text{otherwise} \end{cases}\tag4$$

$$\text{Hann window: } w[n] = \begin{cases} 0.5-0.5cos(\frac{2\pi n}{P-1}), & 0 \leq n < P \\ 0, & \text{otherwise} \end{cases}\tag5$$

For analysis of speech, we usually choose Hamming window as the window function in Eq. 2. Note that in Eqs 3-5, a window function has zero-value outside of an interval $0\leq n < P$ . We therefore define $P$ as **analysis window size**, regarding how many sample points are excerpted from the analyzed signal $x[n]$ for an analysis window. Thus, the value $P$ must be equal or smaller than the FFT window length $N$. An analyzed signal is said to be zero-padded as $P<N$  because $x_m[n] = x[m \cdot M + n]\cdot w[n]$ is zero for $P\leq n <N$ . Since we only have non-zero samples for $0\leq n<P$, only $P$ points excerpted from $x[n]$ are informative for frequency analysis. As FFT window length $N$ is set to be greater than analysis window size $P$, the $N-P$ zero-padded samples in $x_m[n] = x[m \cdot M + n]\cdot w[n]$ just let the spectrum of m-th frame $X[m,k]$ has a smoother outline extrapolating more frequency bins.

## 3. Parameters of spectrogram by DFT
### 3.1. Analysis Window Types
**Why not using rectangular window for spectrogram?** Assume we have a pure sinusoidal signal to be analyzed by spectrogram, e.g.

$$x[n]=\cos(\omega_0 n)\tag6$$

and the corresponding Fourier transform is

$$X(e^{j\omega})=\sum_{k=-\infty}^{\infty}\pi \delta(\omega -\omega_0+2\pi k)+\pi \delta(\omega+\omega_0+2\pi k)\tag7$$

But that spectrogram analyzes the windowed signal $x_m[n] = x[m \cdot M + n]\cdot w[n]$ and this windowed signal is a distorted version of original pure sinusoidal signal $x[n]$. Also, let’s recall the convolution theorem, i.e.

$$\begin{align}
x_0[n] &= x[n] \cdot w[n]\overset{FT}{\longrightarrow} \; X_0(e^{j\omega}) = \frac{1}{2\pi} X(e^{j\omega}) \ast W(e^{j\omega}) \\[1em]
        &= \frac{1}{2\pi} \int_{\omega=-\pi}^{\pi} X(e^{j\theta}) W(e^{j(\omega - \theta)}) \, d\theta\tag8
\end{align}
$$

From Eq. (8), we know that the DTFT $X_0(e^{j\omega})$ depends on the DTFT $W(e^{j\omega})$ (DTFT of the window function $w[n]$). So, let’s check the $W(e^{j\omega})$ if $w[n]$ is a rectangular function of length $P$:

$$ \begin{align}
W(e^{j\omega}) &= \sum_{n=0}^{P-1} w[n] e^{-j\omega n} = \sum_{n=0}^{P-1} 1\cdot e^{-j\omega n} = \frac{1 - e^{-j\omega P}}{1 - e^{-j\omega}} \\[1em]
                &=\frac{e^{-j0.5\omega P} \left( e^{+j0.5\omega P} - e^{-j0.5\omega P} \right)}{e^{-j0.5\omega} \left( e^{+j0.5\omega } - e^{-j0.5\omega } \right)} \\[1em]
                &= e^{-j0.5\omega(P-1)} \frac{\sin(0.5\omega P)}{\sin(0.5\omega)} \tag9
\end{align}
$$

Apparently, $W(e^{j\omega})$ is in a form of a sinc function.<span style="color: red;">**So...Think about the rest by yourself. I talked too much!**</span>

### 3.2. The Releationship between Analysis Window Length $P$ and FFT Window Length $N$ 

Does $N$ or $P$ control "Analysis Bandwidth"? The bandwidth is:

$$bw=\frac{k}{N}f_s \text{ (unit: Hz)}$$

or

$$bw=\frac{k}{P}f_s \text{ (unit: Hz)}$$

where $f_s$ is sampling frequency. Please find the answer by yourself.

### 3.3. Frame Interval M
The parameter $M$ is the **frame interval** or usually **frame hop length** regarding rate of DFT analysis on time axis. Smaller $M$ means more detailed time resolution for spectrogram. Conventionally, $M$ is set to be equal or smaller than $N$, indicating each analysis frame $x_m[n]$ could have some overlapping samples to the neighboring analysis frames $(e.g. \quad x_{m-?}[n]...x_{m-1}[n],\quad x_{m+1}[n],...x_{m+?}[n]).$

### 3.4. Typical Settings
Typically, **analysis window length** (regarding $P$) is set to be 4ms~40ms, that is $P\in[32,320]$ for an 8kHZ sample-rate speech signal, or $P\in[64,640]$ for a 16kHz sample rate signal. FFT window length is typically set to be $N=2^{\lceil \log_2(P) \rceil}$. **Frame Hop Length** is usually set to be 5ms~20ms, that is $M\in[40,160]$ for an 8kHZ sample-rate speech signal, or $M\in[80,320]$ for a 16kHz sample rate signal.

## 4. Implementations
### 4.1. Generate two waveforms (*.wav) that satisfy the following specifications by `signal_gen.c`:
* Sample rates: 16,000 Hz or 8,000 Hz
* Bit depth: 16 bits
* Digital signal that represents
$$x(t) = \sum_{j=0}^{3} \sum_{i=0}^{9} \left\{ a_i \left[ u(t - 0.1i - j) - u(t - 0.1(i+1) - j) \right] s_j(t - 0.1i - j, f_i) \right\}
$$
    * $u(t)$ is the unit step function,
    * $\left\{a_i\right\}|_i = 0, 1, 2, \dots, 9 = \{ 100, 2000, 1000, 500, 250, 100, 2000, 1000, 500, 250 \}$ are amplitudes
    * $\left\{f_i\right\}|_i = 0, 1, 2, \dots, 9 = \{ 0,31.25,500,2000,4000,44,220,440,1760,3960 \}$ are frequencies in Hz;
    * $s_j(t,f_i)$ for $j$ = 0,1,2,3 are the waveforms generated by the function generator:
$$
s_j(t, f_i) = \begin{cases} 
\sin(2\pi f_i t), & \text{for } j = 0 \quad \text{(sine wave)} \\[1em]
f_i t - \lfloor f_i t \rfloor, & \text{for } j = 1 \quad \text{(sawtooth wave)} \\[1em]
\text{sgn}(\sin(2\pi f_i t)), & \text{for } j = 2 \quad \text{(square wave)} \\[1em]
2 \left| 2(f_i t - \lfloor f_i t + 0.5 \rfloor) \right| - 1, & \text{for } j = 3 \quad \text{(triangle wave)}
\end{cases}
$$
* Please name the two wave files: `s-8kHz.wav` and `s-16kHz.wav`. 





### 4.2. Generate 16 spectrograms $X[m,k]$ by `spectrogram.c`
The usage of `spectrogram.c` is:
```
spectrogram w_size w_type dft_size f_itv wav_in spec_out
```
where 
* `w_size`: analysis window size (unit: millisecond)
* `w_type`: a string to be “hamming” or “rectangular”
* `dft_size`: DFT/FFT window size (unit: millisecond)
* `f_itv`: frame interval (unit: millisecond)
* `wav_in`: nput WAVE file
* `spec_out`: output spectrigram data (ascii with 15 decimal places)

The 16 spectrograms are generated by the combinations of the following four window settings and the four waveform files:
### Four Window Settings:

* Setting 1:
    * Analysis window size = 32ms
    * Analysis window type = rectangular
    * DFT/FFT window size = 32ms
    * Frame interval = 10ms
* Setting 2:
    * Analysis window size = 32ms
    * Analysis window type = hamming
    * DFT/FFT window size = 32ms
    * Frame interval = 10ms
* Setting 3:
    * Analysis window size = 30ms
    * Analysis window type = rectangular
    * DFT/FFT window size = 32ms
    * Frame interval = 10ms
* Setting 4:
    * Analysis window size = 30ms
    * Analysis window type = hamming
    * DFT/FFT window size = 32ms
    * Frame interval = 10ms

### Four Waveform Files:
1. `aeueo-8kHz.wav` (download from https://github.com/cychiang-ntpu/mini-project-5-spectrogram/blob/main/aeueo-8Hz.wav): it contains 5 basic vowels with 4 tones in 8kHz.
2. `aeueo-16kHz.wav`(download from https://github.com/cychiang-ntpu/mini-project-5-spectrogram/blob/main/aeueo-16kHz.wav): it contain 5 basic vowels with 4 tones in 16kHz.
3. `s-8kHz.wav` as generated by `signal_gen.c`.
4. `s-16kHz.wav` as generated by `signal_gen.c`.


Therefore, 16 ascii files are saved: s-16k.{Set1~Set4}.txt, s-8k.{Set1~Set4}.txt, aeueo-16kHz.{Set1~Set4}.txt, and aeueo-8kHz.{Set1~Set4}.txt,

### 4.3. Generate the 16 PDFs of Spectrograms by `spectshow.py.`
The usage of “`spectshow.py`” is :
```
python3 spectshow.py in_wav in_txt out_pdf
```
where
* `in_wav`: input WAVE file
* `in_txt`: input text file that contains spectrogram data in ascii
* `out_pdf`: a pdf file that shows the figures of a waveform and the corresponding spectrogram.

You can use Matplotlib (https://matplotlib.org/) library to generate the figures. The pdf files opened should look like Figs. 3 and 4.

![image](https://hackmd.io/_uploads/rJgj53oqxg.png)
Fig. 3: The waveform and spectrogram for aeueo-8kHz.wav.

<br><br>
![image](https://hackmd.io/_uploads/rks6c3scle.png)
Fig. 4: The waveform and spectrogram for aeueo-16kHz.wav.

### 4.3. Write GitHub a Workflow
The workflow need to
1. compile and build source codes,
2. run programs that the generate waveforms,  the 16 spectrogram ascii files, and save the pdf files of the spectrograms, and
3. the artifacts include the generated waveforms for $x(t)$, the 16 spectrogram ascii txt files, and 16 pdf files that display waveforms and spectrogram.

### 4.4. Write a `README.md` File
The README.ms file consisted of the following contents:
1. Division of labor
2. Compare the results by Settings 1-4 and discuss the differences and their significance.
3. Calculate how many multiplications and additions are executed for each frame. 
4. Thoughts and Reflections

### 4.5. Push the Above Files to GitHub
Please provide the link of the repo in [lms3](https://lms3.ntpu.edu.tw/).

## 5. Important Dates
* submission due: YYYY-MM-DD=2025-12-21
* notification of the 1st evaluation: 2025-12-23
* resubmission due: 2025-12-28
* notification of the 2nd evaluation: 2026-1-2

## 6. Grading policy
* The assignments submitted for the first time before the submission due are scored from 50 to 100 points if the acompanying codes can be compiled successfully, or the assignments need to be resubmitted before the resubmission due.
* The assignments submitted after the submission due and before the resubmission due are scored from 40 to 70 points if the accompanying codes are compiled unsuccessfully in the 1st evaluation.
* Students are encouraged to resubmit their codes for the 2nd evaluation to get higher scores if they improve the codes submitted for the 1st evaluation.


## References:
[1]	Xuedong Huang, Alex Acero, Hsiao-Wuen Hon (2001). Spoken Language Processing: a guide to theory, algorithm, and system development, page 274-281. Prentice Hall

[2]	https://www.speech.kth.se/wavesurfer/