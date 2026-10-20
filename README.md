# Neuromorphic-Circuit-with-Resonant-tunneling-diode

RTD Neuromorphic Simulation Suite (MHz → THz)
An independent MATLAB reimplementation and extension of the resonant tunnelling diode (RTD) spiking-neuron model of Owen-Newns et al. (2026), plus a separate leaky integrate-and-fire (LIF) spiking neural network built from RTD-style neurons.
> **Status:** research/educational simulation. The experimental and GHz parameters follow the paper. The Near-THz and THz regimes are **projections** made by shrinking L and C, not measured devices.
![RTD neuromorphic SNN](figures/rtd_neuromorphic_snn.png)
---
Contents
Background
Neuromorphic computing and spiking neurons
Resonant tunnelling diodes
VCSELs and optical I/O
LIF neurons and SNNs
How the RTD neuron works
Model equations
Operating regimes
Repository layout
Quick start
Modules
Limitations
References
---
Background
Neuromorphic computing and spiking neurons
Neuromorphic circuits mimic the architecture and parallel processing principles of biological nervous systems using hardware rather than software. By building electronic counterparts to neurons and synapses, these circuits achieve ultra-low power consumption, event-driven asynchronous processing, and real-time adaptability.
Conventional computers separate memory from processing and run clocked, dense arithmetic. The brain instead computes with many simple units (neurons) that communicate through brief, sparse electrical pulses called spikes. Neuromorphic hardware tries to copy this principle:
Event-driven: a unit consumes energy mostly when it spikes.
Massively parallel: computation is spread over many simple neurons.
Temporal coding: information lives in when spikes occur, not only in how many.
Co-located memory and compute: state (membrane voltage, synaptic weights) is stored where it is processed.
A spiking neuron has three essential behaviours:
Behaviour	Meaning
Integration	Incoming inputs accumulate in an internal state
Threshold firing	When the state crosses a threshold, a spike is emitted
Refractory period	After a spike the neuron resets and cannot fire again immediately
Most neuromorphic chips emulate this with CMOS transistors and run at kHz–MHz rates. Photonic and high-speed electronic devices can reach GHz and potentially THz rates, which is the motivation for this project.
Resonant tunnelling diodes (RTDs)
A Resonant Tunneling Diode is an ultra-fast quantum-electronic semiconductor device composed of a double-barrier quantum well structure (such as InGaAs/AlAs or GaAs/AlGaAs). Electrons tunnel quantum-mechanically through thin potential barriers.When a forward bias is applied, the transmission coefficient peaks when the energy of incident electrons aligns with the discrete quantum-confined energy state inside the quantum well. As voltage increases further past this resonance point, alignment is lost, causing current to drop sharply. This creates a distinct Negative Differential Resistance (NDR) region in its current-voltage (I-V) characteristic, where increasing voltage yields decreasing current. This NDR property enables extremely high-frequency switching, self-sustained oscillations, and multi-stable logic states well into the Terahertz (THz) regime.
This produces an N-shaped current–voltage curve:
At low bias the well level is above the electron supply, so little current flows.
As bias rises, the level aligns with the emitter, so current rises to a peak.
Beyond the peak the level slips out of alignment and current falls even though voltage rises. This is negative differential resistance (NDR).
Past the valley, thermionic and other leakage currents make the current rise again.
The peak-to-valley current ratio (PVCR) measures how pronounced the NDR is.
Why NDR makes a neuron. Biased near the edge of the NDR region and connected to an inductive/capacitive circuit, the RTD is excitable: small perturbations decay away (sub-threshold), but a perturbation that crosses a threshold triggers a large, fast voltage excursion (a spike) before the circuit relaxes to rest. This is the same qualitative behaviour as the FitzHugh–Nagumo and Hodgkin–Huxley neuron models, and it comes from a single nanoscale device. Spike speed is set by the circuit's L and C, and the RTD itself can switch on picosecond to sub-picosecond timescales.
VCSELs and optical I/O
A VCSEL is a semiconductor laser whose laser beam emits perpendicularly from the top surface of the fabricated chip, unlike edge-emitting lasers. They feature distributed Bragg reflector (DBR) mirrors sandwiching an active quantum-well region. VCSELs offer low threshold currents, high modulation bandwidths, low power consumption, and easy 2D array integration, making them ideal for optical interconnects and optoelectronic neuromorphic computing.
A vertical-cavity surface-emitting laser (VCSEL) is a semiconductor laser that emits light perpendicular to the chip surface. It is cheap, compact, low-power and easy to arrange in 2D arrays.
In RTD-based photonic neurons:
Optical input: light falling on the RTD (via a photodetector or photosensitive RTD) shifts the current. In the model this is the `κ·S₀` term, where `S₀` is the optical input amplitude.
Optical output: the RTD's spike current modulates a laser, which emits a short optical pulse. In the code this is represented with a two-level laser rate-equation model (carrier number `N`, photon number `S`).
Why it matters: optical links between neurons avoid the RC delays, crosstalk and fan-out limits of metal wiring and allow all-to-all connectivity.
> The laser in this repo uses a generic rate-equation model with illustrative parameters. It is not a calibrated VCSEL model.
LIF neurons and spiking neural networks
A Spiking Neural Network (SNN) processes information using discrete temporal spikes rather than continuous analog values. The Leaky Integrate-and-Fire model is the most common abstraction for a biological neuron's membrane dynamics.
The standard continuous-time LIF membrane potential V(t) is governed by:
tau_m {dV(t)}{dt} = - (V(t) - V_{rest}) + R_m I_{in}(t)
where tau_m is the membrane time constant, V_{rest} is the resting potential, R_m is the membrane resistance, and I_{in}(t) is the incoming synaptic current. When V(t) crosses a predefined threshold V_{th}, the neuron fires an action potential (spike) and resets V(t) to V_{reset}.
A spiking neural network (SNN) connects layers of such neurons via weighted synapses. Each presynaptic spike produces a postsynaptic current with a double-exponential shape (fast rise `τ_r`, slower decay `τ_f`), and a neuron's input current is the weighted sum of its inputs' currents.
In `rtd_neuromorphic_snn.m` the plain leak is augmented with the RTD's own I–V characteristic as a nonlinear current term, giving an RTD-LIF neuron:
```
C_m · dV/dt = I_syn − V / R_leak − I_RTD(V)
```
The demo network is a 3-layer feed-forward SNN (8 input → 6 hidden → 4 output neurons) driven by Poisson spike trains.
---
How the RTD neuron works
In high-frequency optoelectronic neuromorphic architectures, RTDs and VCSELs are electrically coupled to physically realize the LIF neuron's integration and firing mechanism:Membrane Integration: Incoming electrical synaptic currents charge the capacitive node of the RTD circuit, mimicking the accumulation of charge across a biological cell membrane (the "integrate" phase).Thresholding via NDR: The nonlinear $I$-$V$ curve of the RTD—specifically its NDR region—acts as the firing threshold mechanism. When the integrated voltage reaches the peak voltage of the RTD's NDR region, the device undergoes a bifurcation-like rapid switching event.Spike Generation and Optical Output: This sharp electrical switching event discharges the circuit rapidly, producing a high-amplitude electrical transient that drives the coupled VCSEL. The VCSEL converts this electrical transient into a sub-picosecond or picosecond optical pulse (the "fire" phase).Leak and Reset: The circuit's inherent resistance provides the "leak" term, dissipating residual charge if input sub-threshold currents subside, returning the node to its baseline state.
```
 optical input S₀ ──►  RTD (NDR)  ──► spike in V,I ──► laser/VCSEL ──► optical output
                         ▲  │
 bias V₀ ──► R ──► L ────┘  C (parallel)
```
Bias the RTD at a voltage `V₀` just below the NDR threshold. The neuron sits at a stable resting point.
An optical input pulse injects extra current. If it is too small, `V` returns to rest (no spike).
If it is large enough to push the operating point past the peak, the RTD falls into the NDR region and the capacitor charges/discharges rapidly. `V` makes a large excursion. That is the spike.
The inductor `L` then pulls the current back, the system relaxes to rest, and during this recovery the neuron is refractory.
The spike current drives a laser, converting the electrical spike into an optical output that can excite other neurons.
Tuning knobs: the bias `V₀` sets the firing threshold; `L` and `C` set the spike duration and maximum rate; the optical amplitude sets how strongly the neuron is driven.
---
Model equations
RTD I–V (Schulman form):
```
f(V) = A · ln[(1 + e^((b − c + n₁V)/kT/q)) / (1 + e^((b − c − n₁V)/kT/q))]
         · [π/2 + atan((c − n₁V)/d)]
       + H · (e^(n₂V/kT/q) − 1)
```
Circuit dynamics (state `V`, `I`):
```
C · dV/dt = I − f(V) − κ·S₀(t)
L · dI/dt = V₀ − V − R·I
```
Laser rate equations (carriers `N`, photons `S`):
```
dN/dt = J_eff − (γ_m + γ_l + γ_nr)·N − γ_m·(N − N₀)·S
dS/dt = γ_m·(N − N₀)·S − S/τ_ph + γ_m·N
```
Timescales used by the code:
```
τ_spike = 2π·√(L·C)        f_max ≈ 1 / (3·τ_spike)
```
---
Operating regimes
Defined in `get_regime_params.m`. All regimes share the same I–V parameters and R = 10 Ω; only L and C change. Values below are computed from the formulas above.
Regime	L	C	τ_spike = 2π√(LC)	f_max = 1/(3τ)
Experiment	126 mH	2 pF	≈ 3.2 µs	≈ 0.11 MHz
GHz_Sim (paper Table 1)	126 nH	2 fF	≈ 100 ps	≈ 3.3 GHz
NearTHz (projected)	10 pH	0.05 fF	≈ 140 fs	≈ 2.4 THz
THz (projected)	1 pH	1 aF	≈ 6.3 fs	≈ 53 THz
> Note: the formula `2π√(LC)` is the LC *oscillation period*, a simple scaling estimate, not a full spike-width measurement. The numbers for the two projected regimes are far beyond what real devices can sustain (see Limitations); treat them as scaling trends.
---
Repository layout
```
rtd-thz-neuromorphic/
├── README.md
├── LICENSE
├── main_RTD_THz.m
├── src/
│   └── get_regime_params.m          
├── modules/
│   ├── RTD_IV_Model.m
│   ├── RTD_SingleNeuron.m
│   ├── RTD_EdgeDetection.m
│   ├── RTD_pSNN_ELM.m
│   ├── RTD_FadingMemory.m
│   └── RTD_RegimeComparison.m
├── standalone/                     
│   ├── photonic_snn_iris.m
│   ├── prl_node_simulation.m
│   ├── rtd_threshold_sweep.m
│   ├── rtd_neuromorphic_snn.m
│   └── rtd_neuron_schematic.m


---

## Quick start

Requirements: MATLAB (R2019b or newer recommended for `xline`/`yline`/`sgtitle`). No toolboxes are intended to be required; verify on your version.

```matlab
>> main_RTD_THz                 % runs all modules for REGIME = 'THz'
>> RTD_SingleNeuron(get_regime_params('GHz_Sim'))   % one module, one regime
>> rtd_neuromorphic_snn         % LIF SNN demo (saves figure to the current folder)
```
Change the regime by editing `REGIME` in `main_RTD_THz.m` (`'Experiment' | 'GHz_Sim' | 'NearTHz' | 'THz'`).
---
Modules
Module	What it shows
1. I–V model	Schulman I–V curve, NDR region, PVCR, firing-threshold map
2. Single neuron	Spike vs optical pulse, refractory behaviour, sub/supra-threshold excitability, RTD–laser output
3. Edge detection	Dual modulation (optical = y(t), bias = y(t−τ)) detects rising edges of a Mackey–Glass series
4. Photonic spiking ELM	Fixed random RTD-style spiking layer + trained linear readout classifies Iris
5. Fading memory	10 coupled neurons with delayed feedback; optical attenuation sets memory depth
6. Regime comparison	Spike timescales and design space across the four regimes
LIF SNN	8→6→4 network of RTD-LIF neurons with Poisson input, raster and firing-rate plots
---
Limitations
Please read before drawing conclusions from the results.
Projected regimes are not validated. NearTHz and THz simply shrink L and C. Real devices are limited by carrier transit time, quantum-well lifetime, contact resistance and parasitics, so sub-picosecond spiking at these values is not demonstrated here.
Simplified integration. Several modules use explicit Euler on a stiff system; results should be checked against `ode15s` or the Gauss–Legendre integrator before trusting spike shapes.
ELM module is simplified. It models the spiking layer as thresholded random projections rather than integrating the RTD ODE per sample, so accuracy is not directly comparable to the paper.
Laser model is generic. Parameters are illustrative, not fitted to a real VCSEL.
LIF SNN is phenomenological. It uses an empirical I–V curve and fixed positive weights with no learning; it demonstrates dynamics, not computation.
Independent work. This is not the authors' code, and the figures are not reproductions of the paper's data.
---
References
Owen-Newns et al., 2026, RTD-based photonic spiking neurons (the paper this project reimplements). Add the full citation/DOI here.
J. N. Schulman, H. J. De Los Santos, D. H. Choi, "Physics-based RTD current–voltage equation," IEEE Electron Device Lett., 1996.
R. Tsu and L. Esaki, "Tunneling in a finite superlattice," Appl. Phys. Lett., 1973.
W. Gerstner et al., Neuronal Dynamics, Cambridge University Press, 2014 (LIF and spiking models).
M. Mackey and L. Glass, "Oscillation and chaos in physiological control systems," Science, 1977.
