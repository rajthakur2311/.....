# BioHeat-Simulation
Coupled GPU-accelerated Monte Carlo (MCXLAB) and 3D finite-difference Pennes bio-heat solver to assess photothermal heating and tissue safety under continuous-wave LW-NIR (1300–2000 nm) optical interrogation. Evaluates steady-state thermal fields across broadband emitters, LEDs, and tunable lasers for non-invasive biophotonic pH sensing.


3D Bio-Thermal Modeling of Continuous-Wave LW-NIR Optical Interrogation in Biological TissueThis repository contains the numerical implementation of a coupled optical-thermal computational framework designed to evaluate localized tissue heating during continuous-wave (CW) long-wavelength near-infrared (LW-NIR, 1300–2000 nm) optical sensing.The framework couples GPU-accelerated Monte Carlo photon transport (MCXLAB) with a 3D finite-difference Pennes bio-heat solver (FDM/SOR) to evaluate the biological safety margins and thermal footprints of various optical interrogation sources targeting water vibrational overtone and combination absorption bands.

Key FeaturesCoupled Optical–Thermal Pipeline: 
Automatically maps 3D optical fluence rate distributions $\Phi(x,y,z)$ from GPU-accelerated Monte Carlo simulations directly into volumetric heat generation matrices $Q_{\text{laser}}(x,y,z)$.
3D Steady-State Pennes Bio-Heat Solver: Implements a second-order central finite-difference method solved via Successive Over-Relaxation (SOR) with a convergence residual threshold of $\le 10^{-6}\;^\circ\text{C}$.
Physiological Boundary Modeling: Incorporates active arterial blood perfusion thermoregulation ($\omega_b$), endogenous metabolic heating ($Q_m$), core baseline temperature ($T_b = 37.0\;^\circ\text{C}$), and a linearized Robin top boundary condition accounting for natural convection and thermal radiation ($h_{\text{eff}}$).
Multi-Source Spectral Discretization: Pre-configured routines for loss-compensated broadband and narrowband LW-NIR sources:Broadband Halogen Lamp (1300–2000 nm)Light-Emitting Diode (LED centered at 1450 nm)Thulium-Doped Fiber Amplifier ASE source (TDFA, 1800–2000 nm)Single-mode Tunable Laser Sources (TLS at 1450 nm and 1470 nm)
Validated Numerical Framework: Benchmarked against published literature results for CW laser irradiation in tissue domains.

Computational Workflow:

[Source Spectra & Loss Compensation (3.75 dB)]
                     │
                     ▼
  [GPU Monte Carlo Photon Transport (MCXLAB)]
   - 10^6 photons/band, Gaussian waist w0 = 0.5 mm
   - Tissue absorption µa(λ), scattering µs(λ), g = 0.95
                     │
                     ▼
  [3D Volumetric Heat Generation Matrix Q_laser(x,y,z)]
                     │
                     ▼
  [3D Finite-Difference Pennes Bio-Heat Solver]
   - Iterative SOR solver with Robin surface boundary
   - Conduction + Perfusion + Metabolic + Laser Source
                     │
                     ▼
  [3D Steady-State Temperature Field T(x,y,z)]
   - Axial decay profiles T(z), peak T_max, and safety margins

Dependencies & RequirementsMATLAB: 
R2023b or later (Parallel Computing Toolbox recommended).
Monte Carlo Package: MCXLAB (CUDA-enabled GPU acceleration).
GPU Hardware: NVIDIA CUDA-capable GPU (Compute Capability $\ge 6.0$).
