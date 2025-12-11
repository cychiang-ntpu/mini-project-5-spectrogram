# mini project-5: Spectrogram 

## 1. Objectives
(a) Be familiar with Discrete Fourier Transform (DFT).
(b) Understand how a spectrogram is generated.
\(c\) Mastery over the C language.
(d) Start to learn “shell script.”
(e) Start to use popular python library “Matplotlib.”

## 2. Background
Fig.1 shows waveforms of a compound vowel **‘ai’** and a syllable **‘ya.’** The compound vowel **‘ai’** and the syllable **‘ya’** have similar sound components **‘a’** and **‘i’**, but are in a different order over time. It can be observed that the waveform at the beginning of ‘ai’ and the one at the end of **‘ya’** are very similar; that is, they may have similar harmonic components, since they all are the sound **‘i’** (vice versa for the case of the end of **‘ai’** vs. the beginning of **‘ya.’** These two waveforms are different speech segment. But, if you take Fourier transform for **whole segments of ‘ai’ and ‘ya,’** they may have very similar frequency responses in Discrete-Time Fourier Transforms (DTFTs), indicating that we cannot distinguish the two sounds **‘ai’** and **‘ya’** if the analysis window is too long to make Fourier analysis lose **time resolution** because speech signals change over time.  

![image](https://hackmd.io/_uploads/HyVNc3sqee.png)


To catch the change in frequency components over time, as shown in Fig. 2, **Short-Time Fourier Analysis or Short-Time Fourier Transform (STFT)** is adopted to analyze a speech signal by decomposing the speech signal into a series of short segments, referred to as **analysis frames**, and analyze each one independently. The waveform in each short segment (analysis frame) could be stationary since waveforms in a shorter window look more periodic and have more similar harmonic components. Also, each analysis frame could catch a different waveform to observe signal changing over time. We, therefore, can analyze phonemes and their transition by spectrogram.


A spectrogram is a visual representation of a **spectrum sequence**, capturing spectrum changes with time. Spectrograms are sometimes called sonogram, spectral waterfalls, voiceprints, or voicegrams. In Digital Signal Processing (DSP), we can regard spectrogram, $X[m,k]$  , as a function of time (frame index) and frequency bin, i.e.

$$ X[m,k]=20log_{10}|\sum_{n=0}^{N-1} x_m[n]e^{-j\frac{2\pi kn}{N}}|               \;\qquad\qquad\qquad\qquad\qquad\;\;\;\;\;\;  (1) $$

where $m$ represents frame index; $k$ represents frequency bin index. Note that $$\sum_{n=0}^{N-1} x_m[n]e^{-j\frac{2\pi kn}{N}}$$  in Eq. 1 is in the same form with Discrete Fourier Transform (DFT) and Discrete Fourier Series (DFS). $x_m[n]$ is a short-time signal with a length of $N$ samples, indicating taking DFT of a discrete time signal of length $N$. We usually call the parameter $N$ ‘FFT length’. Note that FFT stands for Fast Fourier Transform which is an efficient implantation of DFT. In mathematics, $x_m[n]$ is expressed by

$$x_m[n]=x[m \cdot M+n]\cdot w[n], \qquad 0\leq n < N \qquad \qquad \qquad \qquad  (2)$$

where $x[m\cdot M+n]$ denotes a time-shifted signal of original input speech,$x[n]$; $M$ represents frame interval (or usually "hop length") in sample; $w[n]$ is a window function, a mathematical function with zero-valued outside of some chosen interval (usually called a compact support). Therefore, we can regard $x_m[n]$ as $m$-th frame signal of length $N$. There are many types of window functions, such as rectangular window (Eq. 3), Hamming window (Eq. 4), and Hann window (Eq. 5).

<br>
<br>

![image](https://hackmd.io/_uploads/HJnw9hj9eg.png)
##### Fig.2: Overview of obtaining spectrogram. The process in sequence, includes framing, taking window, zero padding, and DFT or FFT.


<br>

$$ \text{Rectangular window: } w[n] = \begin{cases} 1, & 0 \leq n < P \\ 0, & \text{otherwise} \end{cases} \qquad \qquad \qquad \qquad \qquad\,\,(3)$$

$$ \text{Hamming window: } w[n] = \begin{cases} 0.54-0.46cos(\frac{2\pi n}{P-1}), & 0 \leq n < P \\ 0, & \text{otherwise} \end{cases} \qquad \; \; \, \, (4)$$

$$ \text{Hann window: } w[n] = \begin{cases} 0.5-0.5cos(\frac{2\pi n}{P-1}), & 0 \leq n < P \\ 0, & \text{otherwise} \end{cases} \qquad \; \; \, \, \;\;\;\;\,\,\,\,(5)$$

For analysis of speech, we usually choose Hamming window as the window function in Eq. 2. Note that in Eqs 3-5, a window function has zero-value outside of an interval $0\leq n < P$ . We therefore define $P$ as **analysis window size**, regarding how many sample points are excerpted from the analyzed signal $x[n]$ for an analysis window. Thus, the value $P$ must be equal or smaller than the FFT window length $N$. An analyzed signal is said to be zero-padded as $P<N$  because $x_m[n] = x[m \cdot M + n]\cdot w[n]$ is zero for $P\leq n <N$ . Since we only have non-zero samples for $0\leq n<P$, only $P$ points excerpted from $x[n]$ are informative for frequency analysis. As FFT window length $N$ is set to be greater than analysis window size $P$, the $N-P$ zero-padded samples in $x_m[n] = x[m \cdot M + n]\cdot w[n]$ just let the spectrum of m-th frame $X[m,k]$ has a smoother outline extrapolating more frequency bins.

## 3. Parameters of spectrogram by DFT
### 3.1. Analysis Window Types
**Why not using rectangular window for spectrogram?** Assume we have a pure sinusoidal signal to be analyzed by spectrogram, e.g.

$$\hspace{10em}x[n]=\cos(\omega_0 n)\qquad\qquad\qquad\qquad\qquad\;\;\;\;\;\;\;\;(6)$$

and the corresponding Fourier transform is

$$X(e^{j\omega})=\sum_{k=-\infty}^{\infty}\pi \delta(\omega -\omega_0+2\pi k)+\pi \delta(\omega+\omega_0+2\pi k)\qquad\qquad\qquad\;\,(7)$$

But that spectrogram analyzes the windowed signal $x_m[n] = x[m \cdot M + n]\cdot w[n]$ and this windowed signal is a distorted version of original pure sinusoidal signal $x[n]$. Also, let’s recall the convolution theorem, i.e.

$$\begin{align}
x_0[n] &= x[n] \cdot w[n]\overset{FT}{\longrightarrow} \; X_0(e^{j\omega}) = \frac{1}{2\pi} X(e^{j\omega}) \ast W(e^{j\omega}) \\[1em]
        &= \frac{1}{2\pi} \int_{\omega=-\pi}^{\pi} X(e^{j\theta}) W(e^{j(\omega - \theta)}) \, d\theta \qquad\qquad\qquad\qquad\qquad\qquad\;\;\;\; (8)
\end{align}
$$

From Eq. (8), we know that the DTFT $X_0(e^{j\omega})$ depends on the DTFT $W(e^{j\omega})$ (DTFT of the window function $w[n]$). So, let’s check the $W(e^{j\omega})$ if $w[n]$ is a rectangular function of length $P$:

$$ \begin{align}
W(e^{j\omega}) &= \sum_{n=0}^{P-1} w[n] e^{-j\omega n} = \sum_{n=0}^{P-1} 1\cdot e^{-j\omega n} = \frac{1 - e^{-j\omega P}}{1 - e^{-j\omega}} \\[1em]
                &=\frac{e^{-j0.5\omega P} \left( e^{+j0.5\omega P} - e^{-j0.5\omega P} \right)}{e^{-j0.5\omega} \left( e^{+j0.5\omega } - e^{-j0.5\omega } \right)} \\[1em]
                &= e^{-j0.5\omega(P-1)} \frac{\sin(0.5\omega P)}{\sin(0.5\omega)} \qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad(9)
\end{align}
$$

Apparently, $W(e^{j\omega})$ is in a form of a sinc function.<span style="color: red;">**So...Think about the rest by yourself. I talked too much!**</span>

### 3.2. The Releationship between Analysis Window Length $P$ and FFT Window Length $N$ 

Does $N$ or $P$ control "Analysis Bandwidth"? The bandwidth is:

$$bw=\frac{k}{N}f_s \text{ (unit: Hz)}$$

or

$$bw=\frac{k}{P}f_s \text{ (unit: Hz)}$$

where $f_s$ is sampling frequency.

### 3.3. Frame Interval M
$\qquad$ The parameter $M$ is the **frame interval** or usually **frame hop length** regarding rate of DFT analysis on time axis. Smaller $M$ means more detailed time resolution for spectrogram. Conventionally, $M$ is set to be equal or smaller than $N$, indicating each analysis frame $x_m[n]$ could have some overlapping samples to the neighboring analysis frames $(e.g. \quad x_{m-?}[n]...x_{m-1}[n],\quad x_{m+1}[n],...x_{m+?}[n]).$

### 3.4. Typical Settings
Typically, **analysis window length** (regarding $P$) is set to be 4ms~40ms, that is $P\in[32,320]$ for an 8kHZ sample-rate speech signal, or $P\in[64,640]$ for a 16kHz sample rate signal. FFT window length is typically set to be $N=2^{\lceil \log_2(P) \rceil}$. **Frame Hop Length** is usually set to be 5ms~20ms, that is $M\in[40,160]$ for an 8kHZ sample-rate speech signal, or $M\in[80,320] for a 16kHz sample rate signal.

## 4. Implementations
### 4.1. Generate two waveforms (*.wav) that satisfy the following specifications:
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
* Please name the two wave files `s-{$smp_rate}.wav` for `{$smp_rate}={“8k”, “16k”}`. 

#### **$ii$.**	$\quad$**We also provide 2 wave files that contain 5 basic vowels with 4 tones in 8kHz and 16kHz sample rates. The two files are named aeueo-16kHz.wav and aeueo-8kHz.wav**

#### **$iii$.** 	**$\;\;\,$Save the spectrograms $X[m,k]$ for all the WAVE files with the following settings, to files in ascii:**

$\qquad\quad$Setting 1:
$\qquad\quad$Analysis window size = 32ms
$\qquad\quad$Analysis window type = rectangular
$\qquad\quad$DFT/FFT window size = 32ms
$\qquad\quad$Frame interval = 10ms

$\qquad\quad$Setting 2:
$\qquad\quad$Analysis window size = 32ms
$\qquad\quad$Analysis window type = hamming
$\qquad\quad$DFT/FFT window size = 32ms
$\qquad\quad$Frame interval = 10ms

$\qquad\quad$Setting 3:
$\qquad\quad$Analysis window size = 30ms
$\qquad\quad$Analysis window type = rectangular
$\qquad\quad$DFT/FFT window size = 32ms
$\qquad\quad$Frame interval = 10ms

$\qquad\quad$Setting 4:
$\qquad\quad$Analysis window size = 30ms
$\qquad\quad$Analysis window type = hamming
$\qquad\quad$DFT/FFT window size = 32ms
$\qquad\quad$Frame interval = 10ms

Therefore, 16 ascii files are saved: s-16k.{Set1~Set4}.txt, s-8k.{Set1~Set4}.txt, aeueo-16kHz.{Set1~Set4}.txt, and aeueo-8kHz.{Set1~Set4}.txt,

#### **$iv$.**	$\quad$ **Compare the results by Settings 1-4 and discuss the differences and their significance.**

#### **$v$.**	$\quad\;\,$ **Calculate how many multiplications and additions are executed for Settings 1-4.**

### B. Requirements

1) `signal_gen.c`: generate digital version of $x(t)$ with two wave files named `s-{$smp_rate}.wav` for {$smp_rate}={“8k”, “16k”}. 
2) `spectrogram.c`: generate the 16 spectrograms with the following usage:
`spectrogram w_size w_type dft_size f_itv wav_in spec_out`
arguments are:
`w_size`: analysis window size (unit: millisecond)
`w_type`: a string to be “hamming” or “rectangular”
`dft_size`: DFT/FFT window size (unit: millisecond)
`f_itv`: frame interval (unit: millisecond)
`wav_in`: nput WAVE file
`spec_out`: output spectrigram data (ascii)
5) `spectshow.py.`: save pdf files that display waveforms and spectrograms according to the files: `s-16k.{Set1~Set4}.txt`, `s-8k.{Set1~Set4}.txt`, `eueo-16kHz.{Set1~Set4}.txt`, and `aeueo-8kHz.{Set1~Set4}.txt`. You can use Matplotlib (https://matplotlib.org/) library to generate the figures. The usage of “`spectshow.py`” is :
`python3 spectshow.py in_wav in_txt out_pdf`
arguments are:
`in_wav`: input WAVE file
`in_txt`: input text file that contains spectrogram data in ascii
`out_pdf`: a pdf file that shows the figures of a waveform and the corresponding spectrogram.

The pdf files opened should look like Figs. 3 and 4.

![image](https://hackmd.io/_uploads/rJgj53oqxg.png)
Fig. 3: The waveform and spectrogram for aeueo-8kHz.wav.

<br><br>
![image](https://hackmd.io/_uploads/rks6c3scle.png)
Fig. 4: The waveform and spectrogram for aeueo-16kHz.wav.

<br><br>


在數位學苑只要繳交 GitHub 的 repository (repo) link 即可，記得權限先開成 private。repo 裡面請包含：
1. `README.md`: 簡單介紹這個 repo 在做什麼，作者有哪些人？如何分工？每個人貢獻哪些部分？佔多少比例？
2. `encoder.c`。
3. `decoder.c`。
4. `test_input_simple.txt`，裡面的文字內容是 ASCII 編碼的 `Do you regret study communication Engineering?`
5. `logger.c/h`，如果直接使用 https://github.com/cychiang-ntpu/logger-example 裡面的 codes。
6. `.github/workflows/c_build-simple.yml`：GitHub Actions 的工作流程設定檔（workflow file），用來設定開發者在 push code 至 remote (GitHub) 之後在 remote 上進行自動建置、自動執行及驗證、以及自動上傳建置產物。可參考 https://github.com/Bensonlllll/build_on_github_test 這個 repo 的 `feature/curl-download-and-test` branch。要包含以下動作：
      * 步驟 1: 取得 repo 中的程式碼
      * 步驟 2: 編譯/建置 C 語言程式 (encoder.c)
      * 步驟 3: 上傳建置成品 `encoder.exe`
      * 步驟 4: 編譯/建置 C 語言程式 (decoder.c)
      * 步驟 5: 上傳建置成品 `decoder.exe`
      * 步驟 6: 使用建置好的 encoder 對輸入 `test_input_simple.txt` 進行處理，輸出編碼 `test_encoded-simple.bin`、codebook `test_codebook-simple.csv`、以及 log `test_encoder-simple.log`。
      * 步驟 7: 上傳輸出編碼 `test_encoded-simple.bin`
      * 步驟 8: 上傳輸出 codebook `test_codebook-simple.csv`
      * 步驟 9: 上傳輸出 log `test_encoder-simple.log`
      * 步驟 10: 使用建置好的 decoder ，使用輸入 `test_encoded-simple.bin` 以及 `test_codebook-simple.csv` 進行處理，輸出解碼文字檔 `test_output-simple.txt` 以及 log `test_decoder-simple.log`。
      * 步驟 11: 上傳解碼文字檔 `test_output-simple.txt`
      * 步驟 12: 上傳 log `test_decoder-simple.log`
      * 步驟 13: 使用 `diff` 檢驗 `test_input_simple.txt` 以及 `test_output-simple.txt` 兩個檔案是否一樣？ 
7. `.github/workflows/c_build-complex.yml`：類似 `.github/workflows/c_build-simple.yml` 的做法，加入下載 cano.txt 文字檔的 step，也就是做 `curl -o test_input_complex.txt https://sherlock-holm.es/stories/plain-text/cano.txt`
      * 步驟 1: 取得 repo 中的程式碼
      * 步驟 2: 編譯/建置 C 語言程式 (encoder.c)
      * 步驟 3: 上傳建置成品 `encoder.exe`
      * 步驟 4: 編譯/建置 C 語言程式 (decoder.c)
      * 步驟 5: 上傳建置成品 `decoder.exe`
      * 步驟 6: 用 `curl` 指令抓取 https://sherlock-holm.es/stories/plain-text/cano.txt 文字檔，儲存成 `test_input_complex.txt`
      * 步驟 7: 使用建置好的 encoder 對輸入 `test_input_complex.txt` 進行處理，輸出編碼 `test_encoded-complex.bin`、codebook `test_codebook-complex.csv`、以及 log `test_encoder-complex.log`。
      * 步驟 8: 上傳輸出編碼 `test_encoded-complex.bin`
      * 步驟 9: 上傳輸出 codebook `test_codebook-complex.csv`
      * 步驟 10: 上傳輸出 log `test_encoder-complex.log`
      * 步驟 11: 使用建置好的 decoder ，使用輸入 `test_encoded-complex.bin` 以及 `test_codebook-complex.csv` 進行處理，輸出解碼文字檔 `test_output-complex.txt` 以及 log `test_decoder-complex.log`。
      * 步驟 12: 上傳解碼文字檔 `test_output-complex.txt`
      * 步驟 13: 上傳 log `test_decoder-complex.log`
      * 步驟 14: 使用 `diff` 檢驗 `test_input_complex.txt` 以及 `test_output-complex.txt` 兩個檔案是否一樣？ 



## 4. Important Dates
* submission due: YYYY-MM-DD=2025-12-01
* notification of the 1st evaluation: 2025-12-03
* resubmission due: 2025-12-08
* notification of the 2nd evaluation: 2025-12-12

## 5. Grading policy
* The assignments submitted for the first time before the submission due are scored from 50 to 100 points if the acompanying codes can be compiled successfully, or the assignments need to be resubmitted before the resubmission due.
* The assignments submitted after the submission due and before the resubmission due are scored from 40 to 70 points if the accompanying codes are compiled unsuccessfully in the 1st evaluation.
* Students are encouraged to resubmit their codes for the 2nd evaluation to get higher scores if they improve the codes submitted for the 1st evaluation.


## References:
[1]	Xuedong Huang, Alex Acero, Hsiao-Wuen Hon (2001). Spoken Language Processing: a guide to theory, algorithm, and system development, page 274-281. Prentice Hall

[2]	https://www.speech.kth.se/wavesurfer/