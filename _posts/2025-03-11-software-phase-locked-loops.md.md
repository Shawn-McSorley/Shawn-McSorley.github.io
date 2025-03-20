---
layout: post
title: "Software Phase Locked Loops"
categories: phase
tags: SDR
---
# Motivation
Often in my research I need to measure the phase of a signal. This signal will commonly be a radio-frequency (RF) tone, described as a time varying voltage $V(t)$,

$$V(t)=A\cos(2\pi f t +\alpha(t)).$$

Here, $A$ is the amplitude of the tone, $f$ is the frequency of the tone, $t$ is the measurement time, and $\alpha(t)$ is some arbitrary phase (which may depend on time). In this post, I want to run through some tools I've found useful for measuring $\alpha(t)$. The end-goal for this post will be a software phase-locked loop (PLL). In short, this is a control system that actively tracks the phase of a tone. The end-goal for this series will be implementing an all-digital PLL inside the field-programmable gate array (FPGA) of a Pluto software-defined radio (SDR). Rather than jumping straight to the FPGA, I want to take the time to describe some important signal processing theory. This can all be done in software, which is a bit easier to work with.
# Overview
Prior to reading this post, I highly recommend going through *PySDR: A Guide to SDR and DSP using Python* [1]. It goes into great detail on in-phase and quadrature (IQ) sampling, and signal processing on IQ samples. I'm also going to assume familiarity with Python, Fourier transforms, Z-transforms, power spectrums and power spectral densities. 

In this post, I want to cover the following theory:
1. How we can model noise in $V(t)$
2. Phasors, and what measuring phase actually mean
3. The Hilbert transform and its relationship to IQ data
4. Generating arbitrary noise from a power spectral density
5. How a phase-locked loop works (with IQ data)

After getting the theory out of the way, I'll start using Python to generate simulated data.

There's a bit to cover. We could just jump straight to code for a PLL, and of course you're welcome to do so (skip to the end). However, if you take the time to understand the tools I present here, you will be able to simulate and understand a multitude of experimental conditions.   

# Modelling Noise
There's two classes of noise that we need to consider. Firstly, there's additive noise. Additive noise is what we usually think about when we hear noise in a signal processing context. Simply put, it's random variations added on to the voltage of our signal. I like to model it as $\bar{A}(t)$, with our measured voltage being,

$$V(t)=A\cos(2\pi f t +\alpha(t))+\bar{A}(t).$$

Additive noise is often associated with electronic noise. I won't limit this generalisation to additive white Gaussian noise, but, more often than not, most additive noise sources I work with are white (i.e. shot noise on a photodetector).

Secondly, there's phase noise. I think of this as random variations in the phase delay of, or arrival time of, our signal. There are a great deal of phase noise sources. I recommend taking a quick read of Enrico's Chart for Phase Noise and Two-Sample Variances [2]. We can accommodate for phase noise fluctuations, $\bar{\epsilon}(t)$, in our voltage model as,

$$V(t)=A\cos(2\pi f t +\alpha(t)+\bar{\epsilon}(t))+\bar{A}(t).$$

It can be helpful to instead consider a phasor model of $V(t)$. I've drawn such a model in Figure. 1. 

![Alt text](/assets/phasor.png){: .fig }
<div style="text-align: center;">
Figure 1: Phasor diagram of the voltage model. $f$, the rate at which our phasor spins; $\bar{\phi}$ the phase that our PLL will measure; $A$, the amplitude (magnitude) of the measurement signal; $\bar{A}$, the additive amplitude noise in the measurement signal; $\bar{\epsilon}$, the phase noise of the signal.
</div>

Our phasor model is a complex exponential, 

$$\pmb{V}=\bar{A}e^{j{\bar{\phi}}}.$$

Here I've lumped the phase and amplitude terms into $\bar{A}$ and $\bar{\phi}$. We should note that this is no longer a voltage, rather it is a model for what the phase and amplitude of our signal are doing at any given time. The phasor model will hopefully make a lot more sense when I introduce the Hilbert transform. Visually, the phasor diagram tells us a few things. 

1. The phasor is rotating around the complex plane at a rate of $f$ Hz.
2. If we were to 'subtract' the rotating component of the phasor, we're left with a point on the complex plane that is, hopefully, slowly moving. This is the purple line. The phase we care about, $\alpha(t)$, is the angle of this purple line.
3. Random phase perturbations $\bar{\epsilon}$ cause the angle of the purple line (phasor) to fluctuate.
4. Amplitude fluctuations $\bar{A}$, couple into both the real and imaginary component of the phasor. This results in a circle of uncertainty around the phasor's point on the complex plane (orange circle). 
5. The phase that we will measure, $\bar{\phi}$, is the combined result of both the additive and phase noise. The measured phase is the angle of the dashed red line. 

Conceptually, these are the most crucial elements to understanding how we can measure the phase of our signal. To get a meaningful measurement, we need the ability to 'subtract' away the rotating component of the phasor. More intuitively, and practically, we can think of this as phase comparison. Consider if we had two rotating phasors,

$$\pmb{V_1}=A_1e^{j{\phi_1}},$$

$$\pmb{V_2}=A_1e^{j{\phi_2}},$$

$$\phi_1=2\pi f t,$$

$$\phi_2=2\pi f t + \alpha.$$

It is practically easier to compare the phase of two signals. We will call $\pmb{V_1}$ our local oscillator, which provides a stable phase reference to compare $\pmb{V_2}$ to. We can readily extract the phase of our signal as,

$$\alpha=\arg(\pmb{V_2V_1^*})=\arg(A_1A_2e^{j\alpha}).$$

The phasor model also tells us how noise will couple into our measurement. If we were to look at $\arg(\pmb{V_2V_1^*})$ in practice, we would not obtain $\alpha$. Rather we would obtain a combination of the phase we care about, $\alpha$, the phase noise, $\bar{\epsilon}$, and the additive noise, $\bar{A}$. To further describe how additive noise couples into our measurement, I would like to introduce the Hilbert transform and the Analytic signal.

# A Brief Review of IQ Sampling
I think it's most useful to mention IQ sampling here. Rather than sampling $V(t)$ directly, we can sample its I and Q components. Again, I highly recommend reading PySDR [1] for a more thorough description. For a signal $x(t)=A\cos(2\pi f t-\phi)$, we can, through some trigonometry, express it as,

$$ x(t)=I\cos(2\pi f t)+Q\sin(2\pi f t)$$

$$I=A\cos(\phi)$$

$$Q=A\sin(\phi)$$

Most SDRs will sample I and Q, some sample $V(t)$ directly. Soon, I will talk about going back and forth between the two. In Physics and Engineering, we like to use complex exponentials. They're slightly easier to work with mathematically, and, arguably also easier to work with numerically. What's nice about IQ data, is that it naturally describes a phasor. The magnitude of the IQ phasor, $A$, is obtained from,

$$A=\sqrt{I^2+Q^2},$$

and the phase, $\phi$,

$$\phi=\tan^{-1}(Q/I).$$

Thus, we can describe our IQ data as a phasor, $\pmb X$,

$$\pmb{X}=I+jQ.$$

$$A=|I+jQ|,$$

$$\phi=\angle(I+jQ)$$

What's nice about IQ data is that the rotating part of the phasor, herein called the carrier, is subtracted off for us in the analogue front end of the receiver. An important clarification is that the IQ phasor $\pmb X$ is used as it makes our life easy. We can readily analyse imaginary data in software. We can take its FFT to visualise what makes up our signal. We can frequency up-shift and down-shift simply by mixing the IQ data with a complex exponential. We can also readily extract phase and amplitude by taking the magnitude or argument of our IQ samples. It's important to remember that our data, $\pmb X=I+jQ$, represents the two independent quadratures of a physical signal, related back to $x(t)$ via,

$$ x(t)=I\cos(2\pi f t)+Q\sin(2\pi f t).$$


## Analytical Signals

With a good understanding of IQ sampling and the IQ phasor, we can introduce what's called an **analytic signal**. A good place to look for information is the SciPy reference [3] and [4]. To best understand this, let's look at a very important analytic signal, the complex exponential,

$$\hat{u}(t)=Ae^{j\phi(t)}=A\cos(\phi(t))+jA\sin{\phi(t)}.$$

The key properties of an analytic signal are:
- The imaginary component, $\sin(\phi)$ in this case, is the real component, $\cos(\phi)$ here, shifted by $90^\circ$ degrees.
- The analytic signal, $e^{j\phi}$ here, has no negative frequency components. This can readily be verified by adding the fourier transforms of $\cos(\phi)$ and $j\sin(\phi)$
- The imaginary component, $\sin(\phi)$, does not tell us anything new about our real signal $\cos(\theta)$

The reason we like the analytic signal is that we can obtain the phase and amplitude quite easily with,

$$A=|\hat{u}(t)|$$

$$\phi(t)=\angle \hat{u}(t)$$

The analytic signal is, for our purposes, a more technical term for phasor. The reason I introduce it, is because there's a useful tool called the Hilbert transform that can be used to generate the analytic signal from a real valued signal.

We can construct an analytic signal with,

$$\hat{u}(t)=x(t)+j\mathcal{H}\{x(t)\},$$

where $\mathcal{H}\\{ \dots \\}$ denotes the Hilbert transform. The Hilbert transform is, in and of itself, an integral transform. However, we can opt to relate it to the Fourier transform of $x(t)$. Essentially, the Hilbert transform phase shifts the signal $x(t)$. But, it actually applies a different phase shift to the positive and negative frequencies in the signal. The Fourier transform of $\mathcal{H}\\{x(t)\\}$, is,

$$\DeclareMathOperator{\sign}{sign}\mathcal{F}\{\mathcal{H}\{x\}\}=-j\sign(f)\mathcal{F}(x) $$

That is, the Hilbert transform rotates positive frequencies by $-90^\circ$ degrees, and negative frequencies by $+90^\circ$ degrees. Computationally, this is the easiest way to perform the Hilbert transform. You can take its FFT, and mix the positive frequency components by $-j$, and the negative frequency components by $+j$. Alternatively, you can take the FFT, zero out the negative frequencies, and double the positive frequencies. The inverse FFT would then provide the analytic signal directly.

The reason we would do this, is again convenience. It's easier to determine the amplitude and phase of the analytic signal, than it is to extract it from the real valued signal. Complex numbers make things easier. 

If we consider our original voltage model $V(t)=A\cos(2\pi f t + \alpha(t)+\bar{\epsilon}(t))+\bar{A}(t)$, we can readily construct its analytic signal (its phasor! shown in Figure 1) as,

$$\mathcal{H}\{V(t)\}=A\sin(2\pi f t+\alpha(t)+\bar\epsilon(t))+\mathcal{H}\{\bar{A}(t)\}$$

$$\hat{u}(t)=A\cos(2\pi f t + \alpha(t)+\bar{\epsilon}(t))+\bar{A}(t) + j[A\sin(2\pi f t+\alpha(t)+\bar\epsilon(t))+\mathcal{H}\{\bar{A}(t)\}]$$

I like this approach, as it readily highlights that the additive noise will couple into both the imaginary and real component of our phasor model. That is, our phasor is corrupted by $\bar A(t)+j\mathcal{H}\\{\bar A(t) \\}$. This will help us model the impact of additive noise on our phase measurement. Before we do some fun math, a word of caution. It's very easy to fall into the trap of thinking that IQ data and the analytic signal tell us the same thing (I've fallen into it myself). However, **in general**, they are different. The key distinction stems from a fundamental property of the analytic signal -- the imaginary component does not tell us anything new about our signal, it's dependent on the real signal we care about via the Hilbert transform. For IQ data however, we are looking at two independent quadratures, we cannot discern one from the other. Regardless, the analytic signal and IQ data that we would measure can be related to each other.

Here are some useful properties, derived using the analytic signal (hidden for clarity), that will be useful for our signal processing:
- - - -
#### Property 1: Relationship with $x(t)=I\cos(2\pi f t)+Q\sin(2\pi f t)$

The IQ components of our voltage $V(t)$, can be related to the analytic signal via,

$$I=\mathcal{Re}\{\hat{u}(t)e^{-j2\pi f t}\},$$

and,

$$Q=-\mathcal{Im}\{\hat{u}(t)e^{-j2\pi f t}\}.$$

Here, the analytic signal is down-converted by mixing with the carrier $e^{-j2\pi f t}$. This is useful if you want to generate IQ data from a simulated real signal.

<details>
  <summary>Click to expand</summary>

Consider a general signal to be transmitted,
$$ x(t)=A\cos(2\pi f t - \phi)$$
It can be expressed as,

$$ x(t)=I\cos(2\pi f t)+Q\sin(2\pi f t)$$
$$I=A\cos(\phi)$$
$$Q=A\sin(\phi)$$

The analytical signal is provided by,
$$\hat{u}(t)=A\cos(2\pi f t - \phi)+j\mathcal{H}\{x(t)\}=A\cos(2\pi f t - \phi) + jA\sin(2\pi f t-\phi)=Ae^{j(2\pi f t-\phi)}$$
Let's separate the analytical signal into,
$$ \hat{u}(t)=Ae^{-j\phi}e^{j2\pi f t}$$
We can also take $Ae^{-j\phi}=I_u+jQ_u$,
$$ \hat{u}(t)=(I_u+jQ_u)e^{j2\pi f t}$$
Now, recall that $x(t)=\mathcal{Re}\{\hat{u}(t)\}$,
$$ \hat{u}(t)=(I_u+jQ_u)e^{j2\pi f t}=(I_u+jQ_u)(\cos(2\pi f t)+j\sin(2\pi f t))=(I_u\cos(2\pi f t)-Q_u\sin(2\pi f t))+j(\dots)$$
From which we find that,
$$x(t)=I_u\cos(2\pi f t)-Q_u\sin(2\pi f t)$$
This tells us that the conventional IQ data approach isn't quite consistent. From comparison with,
$$x(t)=I\cos(2\pi f t)+Q\sin(2\pi f t),$$
we see that $I_u=I$ and $Q_u=-Q$.
In other words, to get IQ data, we actually want the conjugate of the down shifted analytical signal $(\hat{u}e^{-j2\pi f t})^*$. 

</details>
- - -
#### Property 2.1: Arbitrary Phase Modulation
We can generate a signal with arbitrary phase modulation , $\phi(t)$, by first generating the IQ samples,
$$ \pmb X = e^{j\phi(t)}=I+jQ.$$
We then have two equivalent options of converting this to a real signal:

$$x(t)=I\cos(2\pi f t)+Q\sin(2\pi f t)$$

$$x(t)=\mathcal{Re}\{\pmb X^* e^{j2\pi f t}\}$$

<details>
  <summary>Click to expand</summary>
We can verify $$x(t)=\mathcal{Re}\{\pmb X^* e^{j2\pi f t}\}$$ by direct evaluation.

$$x(t)=\mathcal{Re}\{\pmb X^* e^{j2\pi f t}\} = \mathcal{Re} \{ (I-jQ)(\cos(2\pi f t) + j\sin(2\pi f t))\}$$

$$=\mathcal{Re}\{I\cos(2\pi f t)-j^2Q\sin(2\pi f t) + j(\dots) \}= I\cos(2\pi f t)+Q\sin(2\pi f t)$$

</details>


---

#### Property 2.2: Arbitrary Amplitude Modulation
We can generate a signal with arbitrary amplitude modulate, $A(t)$, by generating IQ samples,

$$\pmb X = A(t)e^{j\phi(t)}.$$

If we want to instead model additive signal (or noise), $\bar{A}(t)$,

$$\pmb X =Ae^{j\phi(t)}+\bar{A}(t) - j\mathcal{H}\{\bar A (t)\}.$$

For this equation to make sense, $S_{\bar A}(f)$ should be defined from 0Hz to the sample rate $f_s$. The PSD would then be shifted up to the carrier when transmitting.

<details>
  <summary>Click to expand</summary>

Consider the signal,

$$ x(t)=A\cos(2\pi f t-\phi)+\bar{A}(t)$$

The analytical signal is, 

$$ \hat{u}(t)=(I_u+jQ_u)e^{j2\pi f t}+\bar{A}(t)+j\mathcal{H}\{\bar{A}\}$$

We can down convert this to figure out I and Q,

$$\hat{u}(t)e^{-j2\pi f t}=I_u+jQ_u+(\bar{A}(t)+j\mathcal{H}\{\bar{A}\})(\cos(2\pi f t)-j\sin(2\pi f t))$$

$$ = [I_u + \bar{A}(t)\cos(2\pi f t)+\mathcal{H}\{\bar{A}\}\sin(2\pi f t)] + j[Q_u+\mathcal{H}\{\bar{A}\}\cos(2\pi f t)-\bar{A}\sin(2\pi f t)]$$
And from property 1,

$$I=I_u + \bar{A}(t)\cos(2\pi f t)+\mathcal{H}\{\bar{A}\}\sin(2\pi f t)$$

$$Q=-(Q_u+\mathcal{H}\{\bar{A}\}\cos(2\pi f t)-\bar{A}\sin(2\pi f t))$$
Intuitively, our measurement will see additive noise around the carrier. In the IQ sampling process, the additive noise is frequency shifted down to 0Hz. This is why there are $\cos(2\pi f t)$ and $\sin(2\pi f t)$ terms. However, if we're directly generating noise from 0Hz to $f_s/2$, we can ignore the down-conversion process (it is implicitly assumed). This provides,

$$I=I_u + \bar{A}(t)$$

$$Q=-(Q_u+\mathcal{H}\{\bar{A}\})$$

From which it is clear that,

$$\pmb X =I+jQ=Ae^{j\phi(t)}+\bar{A}(t) - j\mathcal{H}\{\bar A (t)\}.$$

</details>

---
#### Property 3: Additive to Phase Noise Conversion

For an additive noise source with power spectral density (PSD) $S_\bar A(f)$, it is seen in our measurement as phase noise with a PSD $S_{\phi,\bar A}(f)$, given by,

$$ S_{\phi, \bar A}(f)=\frac{2}{A^2}S_\bar A(f),$$

for signal amplitude $A$. This essentially tells us how our signal-to-noise ratio limits our measurement. As our signal decreases, we see more additive noise in our phase measurement. This is a very important result!

<details>
  <summary>Click to expand</summary>
A simplified model of the noise sources is given by:

$$v(t)=A\cos(\omega t + \theta(t) + \bar{\epsilon}(t)) + \bar{A}(t)$$

where $v(t)$ is the voltage of the signal input to our phase measurement device, $A$ is the amplitude of the signal, $\omega$ is the frequency of the signal, $\theta$ is the deterministic phase of the signal we want to measure, $\bar{\epsilon}(t)$ is the phase noise, and $\bar{A}(t)$ is the amplitude noise. We can usually lump the phase we're interested in measuring, and its noise together.

A phase-locked loop (PLL) can be used to measure the phase of the signal. For the following, it is useful to think of the PLL as a device that outputs the phase $\phi(t)$, where,

$$ \phi(t) = \theta(t) + \bar{\epsilon}(t) + k \bar{A}(t). $$

Here, $\phi(t)$ is the total phase of the signal, $\theta(t)$ is the signal we want to measure, $\bar\epsilon$ is the phase noise, and $k\bar{A}(t)$ is the amplitude noise. Here, $k$, is just to show that the amplitude noise couples into our measurement. We will work out very shortly how it does so.

To avoid confusion, the phasor model can be separated into its time varying component and its phase component (ignoring additive noise):

$$ v(t) = \Re(Ae^{j\omega t}e^{j\phi})$$

and the ideal PLL measurement will be given by:

$$ v(t) = \text{Arg}(Ae^{j\omega t}e^{j\phi}\times e^{-j\omega t})=\phi$$

There are a few ways to figure out how $\bar A$ will couple into the phase measurement from the PLL. I'm going to run through a derivation that uses the Hilbert transform, as it is essentially a mathematical description of the phasor diagram above. The goal is to find how the noise output, $S_{PLL}(f)$, depends on the phase noise $S_{\bar\phi}(f)$ and the additive noise $S_{\bar A}(f)$, where $S(f)$ is the noise power spectral density.

While crude, this derivation will use:

$$S(f)=\frac{1}{\text{RBW}}|\mathcal{F}(v(t))|^2$$

where $\mathcal{F}$ is the Fourier transform, $v(t)$ is the signal of interest and $\text{RBW}$ is the resolution bandwidth of our FFT. I'm mixing up continuous and discrete here, but it helps get the end result. A more rigorous approach would be to consider the auto-correlation of our signals.

Considering our noise model,
$$v(t) = A \cos(\omega_i t + \theta_i + \bar{\epsilon}(t))+\bar{A}(t)$$
Take the Hilbert transform,
$$\mathscr{H}\{v(t)\}=A\sin(\omega_i t + \theta_i + \bar{\epsilon}(t))+\mathscr{H}\{\bar{A}(t)\}$$
Form the analytic signal,
$$u(t) = A \cos(\omega_i t + \theta_i + \bar{\epsilon}(t))+\bar{A}(t) + j[A\sin(\omega_i t + \theta_i + \bar{\epsilon}(t))+\mathscr{H}\{\bar{A}(t)\}]$$
The instantaneous phase $\bar{\phi}(t)$ is given by the argument of this analytic signal (which is the phase of the phasor drawn in Figure 2),
$$\bar{\phi}(t) = \tan^{-1} [\frac{A\sin(\omega_i t + \theta_i + \bar{\epsilon}(t))+\mathscr{H}\{\bar{A}(t)\}}{A \cos(\omega_i t + \theta_i + \bar{\epsilon}(t))+\bar{A}(t)}]$$
For small $x,y$ it can be shown using a taylor expansion of $x$ and $y$
that,
$$\tan^{-1}[\frac{\sin(\alpha)+x}{\cos(\alpha)+y}]\approx \alpha + x \cos(\alpha) - y \sin(\alpha)$$
from which it follows for small noise fluctuations $\bar{A}(t)$ that,
$$\bar{\phi}(t) = \omega_i t + \theta_i + \bar{\epsilon}(t) + \frac{\mathscr{H}\{\bar{A}(t)\}}{A}\cos(\omega_i t + \theta_i + \bar{\epsilon}(t)) - \frac{\bar{A}(t)}{A}\sin(\omega_i t + \theta_i + \bar{\epsilon}(t))$$
We can ignore the contribution from the frequency and static offset,
providing $$\begin{aligned}
    \bar{\phi(t)}=\bar{\epsilon}(t)+\frac{\mathscr{H}\{\bar{A}(t)\}}{A}\cos(\omega_i t + &\theta_i + \bar{\epsilon}(t)) - \frac{\bar{A}(t)}{A}\sin(\omega_i t + \theta_i + \bar{\epsilon}(t))\\
    \bar{\phi(t)} &= \bar{\epsilon}(t)+ A'(t)
\end{aligned}$$ From this it is clear that, $$\begin{aligned}
    S_{\phi,\epsilon}(f)=S_{\epsilon}(f)
\end{aligned}$$ which is expected. 

It is a bit more involved to work out the contribution from the additive noise. Let's just consider the additive noise now,

$$ \bar{\phi(t)}=\frac{\mathscr{H}\{\bar{A}(t)\}}{A}\cos(\omega_i t + \theta_i + \bar{\epsilon}(t)) - \frac{\bar{A}(t)}{A}\sin(\omega_i t + \theta_i + \bar{\epsilon}(t)) $$

Using the property of the Hilbert transform $\mathcal{F}(\mathcal{H}(x(t)))=-j \mathcal{F}(x(t))$,

$$ \mathcal{F}(\bar{\phi}(t)) = \mathcal{F}\{\frac{\mathcal{H}\{\bar{A}(t)\}}{A}\cos(\omega_i t + \theta_i + \bar{\epsilon}(t))\} + \mathcal{F}\{-\frac{\bar{A}(t)}{A}\sin(\omega_i t + \theta_i + \bar{\epsilon}(t))\} $$

Which when using the convolution property,

$$ \mathcal{F}(\bar{\phi}(t)) = \mathcal{F}\{\frac{\mathcal{H}\{\bar{A}(t)\}}{A}\} \ast \mathcal{F}\{\cos(\omega_i t + \theta_i + \bar{\epsilon}(t))\} + \mathcal{F}\{-\frac{\bar{A}(t)}{A}\} \ast \mathcal{F}\{\sin(\omega_i t + \theta_i + \bar{\epsilon}(t))\} $$
$$ \mathcal{F}(\bar{\phi}(t)) = -j\frac{\mathcal{F}\{\bar{A}(t)\}}{A} \ast \frac{1}{2}[\delta(f-f_i)+\delta(f+f_i)] + \mathcal{F}\{\frac{\bar{A}(t)}{A}\} \ast \frac{1j}{2}[-\delta(f-f_i)+\delta(f+f_i)] $$
Define $\pmb{A}(f)=\mathcal{F}\{\bar A(t)\}$,

$$ \mathcal{F}(\bar{\phi}(t)) = \frac{1}{2A}[-j\pmb{A}(f-f_i)-j\pmb{A}(f+f_i)-j\pmb{A}(f-f_i)+j\pmb{A}(f+f_i)]=\frac{-1j}{A}\pmb{A}(f-f_i) $$

Recalling that we can crudely relate the power spectral density to the magnitude of the Fourier transform, we can write the power spectral density of the phase noise as,

$$|\mathcal{F}(\bar \phi(t))|=\frac{1}{A^2} |\pmb A(f-f_i)| $$

$$ S_{\phi, A}(f) = \frac{1}{A^2} S_{\bar A}(f-f_i)$$

For our case, we only want to consider the single sided PSD, and we make the frequency shift by $f_i$ implicit. So we can write,

$$ S_{\phi, A}(f) = \frac{2}{A^2} S_{\bar A}(f).$$

</details>

- - - 

With this understanding of IQ data, the Analytic signal and the Hilbert transform, we're well equipped to simulate a variety of experimental conditions. Onto software!
# The rest is in progress!
# References
[1] “PySDR: A Guide to SDR and DSP using Python.”  [Online]. Available: [https://pysdr.org/#](https://pysdr.org/#)

[2] E. Rubiola and F. Vernotte, “The Companion of Enrico’s Chart for Phase Noise and Two-Sample Variances,” _IEEE Trans. Microwave Theory Techn._, vol. 71, no. 7, pp. 2996–3025, Jul. 2023, doi: [10.1109/TMTT.2023.3238267](https://doi.org/10.1109/TMTT.2023.3238267).

[3] “hilbert — SciPy v1.15.2 Manual.”  [Online]. Available: [https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert.html](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert.html)

[4] Mathuranathan, “Understanding Analytic Signal and Hilbert Transform,” GaussianWaves. [Online]. Available: [https://www.gaussianwaves.com/2017/04/analytic-signal-hilbert-transform-and-fft/](https://www.gaussianwaves.com/2017/04/analytic-signal-hilbert-transform-and-fft/)
