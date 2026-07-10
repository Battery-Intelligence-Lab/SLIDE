# Battery Simulation Methods and Electrochemical Modeling: A Comprehensive Technical Report

**Prepared for:** SLIDE (Simulator for Lithium-Ion Degradation) development
**Date:** 2026-03-16
**Scope:** Electrochemical models, degradation models, thermal models, numerical methods, and emerging approaches

---

## Table of Contents

1. [Electrochemical Models](#1-electrochemical-models)
2. [Degradation Models](#2-degradation-models)
3. [Thermal Models](#3-thermal-models)
4. [Numerical Methods](#4-numerical-methods)
5. [Emerging Approaches (2024-2025)](#5-emerging-approaches-2024-2025)
6. [Relevance to SLIDE](#6-relevance-to-slide)

---

## 1. Electrochemical Models

### 1.1 Single Particle Model (SPM)

**Overview:**
The SPM is the simplest physics-based electrochemical model that still captures the essential lithium intercalation dynamics. It represents each electrode as a single spherical particle, assuming uniform reaction current density across the electrode.

**Key Assumptions:**
- Each electrode is represented by a single representative particle
- Electrolyte concentration is uniform (no electrolyte dynamics)
- Electrolyte potential drop is negligible
- Uniform temperature within the cell
- Solid-phase diffusion is the rate-limiting process
- Reaction current is uniform across each electrode thickness

**Governing Equations:**

1. **Solid-phase diffusion** (Fick's law in spherical coordinates):
   ```
   dc_s/dt = D_s/r^2 * d/dr(r^2 * dc_s/dr)
   ```
   where `c_s` is solid-phase concentration, `D_s` is solid-phase diffusivity, `r` is radial coordinate.

2. **Boundary conditions:**
   - At center: `dc_s/dr|_{r=0} = 0` (symmetry)
   - At surface: `D_s * dc_s/dr|_{r=R} = -j_n` (flux equals reaction rate)
   - `j_n = I / (a_s * F * L * A)` where `a_s` is specific surface area, `F` is Faraday's constant, `L` is electrode thickness, `A` is electrode area

3. **Butler-Volmer kinetics:**
   ```
   j_n = j_0 * [exp(alpha_a * F * eta / RT) - exp(-alpha_c * F * eta / RT)]
   ```
   where `eta = phi_s - phi_e - U(c_ss)` is the overpotential, `U` is the open-circuit potential (OCP).

4. **Terminal voltage:**
   ```
   V = U_p(c_ss,p) - U_n(c_ss,n) + eta_p - eta_n - I*R_cell
   ```

**Limitations:**
- Inaccurate at high C-rates (>2C typically) because electrolyte dynamics become significant
- Cannot capture concentration gradients across electrode thickness
- Fails for thick electrodes or low-conductivity electrolytes
- No electrolyte depletion effects

**When to use:** Low-to-moderate C-rates (<1-2C), thin electrodes, parameter identification, real-time BMS applications, degradation studies where electrolyte effects are secondary.

**Computational cost:** Very low. O(N_r) per particle where N_r is the number of radial discretization points (typically 5-30). Full cell simulation in milliseconds.

**Key References:**
- Ning, G., & Popov, B. N. (2004). "Cycle life modeling of lithium-ion batteries." *J. Electrochem. Soc.*, 151(10), A1584. DOI: 10.1149/1.1787631
- Guo, M., et al. (2011). "Single-particle model for a lithium-ion cell." *J. Electrochem. Soc.*, 158(2), A122. DOI: 10.1149/1.3521314
- Marquis, S. G., et al. (2019). "An asymptotic derivation of a single particle model with electrolyte." *J. Electrochem. Soc.*, 166(15), A3693. DOI: 10.1149/2.0341915jes

---

### 1.2 SPMe (SPM with Electrolyte)

**Overview:**
The SPMe extends the SPM by incorporating electrolyte dynamics (concentration and potential variations) while maintaining the single-particle representation of each electrode. It bridges the gap between SPM and the full DFN model.

**What it adds over SPM:**

1. **Electrolyte concentration** (transport in the electrolyte phase):
   ```
   eps_e * dc_e/dt = d/dx(D_e_eff * dc_e/dx) + (1-t+)/F * j_n * a_s
   ```
   where `eps_e` is electrolyte volume fraction, `D_e_eff` is effective electrolyte diffusivity, `t+` is transference number.

2. **Electrolyte potential** (modified Ohm's law):
   ```
   i_e = -kappa_eff * dphi_e/dx + 2*kappa_eff*RT/F * (1-t+) * (1 + d(ln f)/d(ln c_e)) * d(ln c_e)/dx
   ```

3. **Updated terminal voltage:**
   ```
   V = U_p(c_ss,p) - U_n(c_ss,n) + eta_p - eta_n + Delta_phi_e + Delta_phi_solid
   ```
   where `Delta_phi_e` accounts for electrolyte potential drop across separator and electrodes.

**Key insight from asymptotic analysis (Marquis et al. 2019):**
The SPMe can be derived rigorously from the DFN model through asymptotic expansion in the ratio of electrode particle timescale to electrolyte transport timescale. This gives a systematic justification for when SPMe is valid.

**When to use:** Moderate C-rates (up to ~3-4C), when electrolyte effects matter but full spatial resolution across electrodes is not needed. Good balance of accuracy and speed for many applications.

**Computational cost:** Moderate. Requires solving 1D PDE for electrolyte concentration across cell thickness in addition to the spherical diffusion. Typically 5-50x more expensive than SPM but 10-100x cheaper than DFN.

**Key References:**
- Marquis, S. G., et al. (2019). *J. Electrochem. Soc.*, 166(15), A3693. DOI: 10.1149/2.0341915jes
- Richardson, G., et al. (2020). "Generalised single particle models for high-rate operation of graded lithium-ion electrodes." *Electrochimica Acta*, 339, 135862. DOI: 10.1016/j.electacta.2020.135862

---

### 1.3 Doyle-Fuller-Newman (DFN) / Pseudo-2D (P2D) Model

**Overview:**
The DFN model (also called P2D) is the gold standard physics-based model for lithium-ion batteries. It couples 1D macroscopic transport across the cell thickness (x-direction) with microscopic solid-phase diffusion in spherical particles at each point (r-direction), creating a "pseudo-2D" framework.

**Governing Equations (6 coupled PDEs):**

1. **Solid-phase diffusion** (at every x-position):
   ```
   dc_s/dt = D_s/r^2 * d/dr(r^2 * dc_s/dr)
   ```

2. **Electrolyte concentration:**
   ```
   eps_e * dc_e/dt = d/dx(D_e_eff * dc_e/dx) + (1-t+)/F * j_n * a_s
   ```

3. **Solid-phase potential** (Ohm's law):
   ```
   d/dx(sigma_eff * dphi_s/dx) = a_s * F * j_n
   ```

4. **Electrolyte potential:**
   ```
   d/dx(kappa_eff * dphi_e/dx) + d/dx(kappa_D * d(ln c_e)/dx) = -a_s * F * j_n
   ```

5. **Butler-Volmer kinetics** (couples solid and electrolyte):
   ```
   j_n = j_0 * [exp(alpha_a*F*eta/RT) - exp(-alpha_c*F*eta/RT)]
   j_0 = k * c_e^alpha_a * (c_s_max - c_ss)^alpha_a * c_ss^alpha_c
   ```

6. **Charge conservation:**
   ```
   i_s + i_e = I/A  (total current density)
   ```

**Boundary and interface conditions:**
- Separator: no electronic current, continuous electrolyte flux/potential
- Current collectors: no ionic current, electronic current = applied current
- Each particle: zero flux at center, Butler-Volmer flux at surface

**Computational cost:** High. For N_x spatial points across cell thickness and N_r radial points per particle, the system has O(N_x * N_r + 4*N_x) coupled differential-algebraic equations. Typical: N_x = 30-100, N_r = 10-30, giving ~1000-3000+ unknowns. Full discharge simulation: seconds to minutes.

**When to use:** High C-rates (>2-3C), design optimization, detailed physics studies, when spatial gradients matter, thick electrodes, fast charging protocol development.

**Key References:**
- Doyle, M., Fuller, T. F., & Newman, J. (1993). "Modeling of galvanostatic charge and discharge of the lithium/polymer/insertion cell." *J. Electrochem. Soc.*, 140(6), 1526. DOI: 10.1149/1.2221597
- Fuller, T. F., Doyle, M., & Newman, J. (1994). "Simulation and optimization of the dual lithium ion insertion cell." *J. Electrochem. Soc.*, 141(1), 1. DOI: 10.1149/1.2054684
- Newman, J., & Thomas-Alyea, K. E. (2004). *Electrochemical Systems* (3rd ed.). Wiley.

---

### 1.4 Equivalent Circuit Models (ECM)

**Overview:**
ECMs represent the battery as an electrical circuit with resistors and capacitors. They are empirical/semi-empirical models that do not resolve internal electrochemistry but can accurately predict terminal voltage behavior.

**Model Hierarchy:**

**0th Order (Rint model):**
```
V = OCV(SOC) - I * R_0
```
- Single internal resistance
- No dynamic response
- Use: crude SOC estimation, very fast computation

**1st Order (1-RC Thevenin model):**
```
V = OCV(SOC) - I * R_0 - V_1
dV_1/dt = -V_1/(R_1*C_1) + I/C_1
```
- One R-C parallel pair captures charge transfer dynamics
- Time constant tau_1 = R_1*C_1 (typically 1-100s)
- Use: BMS, basic dynamic simulation

**2nd Order (2-RC model):**
```
V = OCV(SOC) - I * R_0 - V_1 - V_2
dV_1/dt = -V_1/(R_1*C_1) + I/C_1
dV_2/dt = -V_2/(R_2*C_2) + I/C_2
```
- Two R-C pairs: fast dynamics (charge transfer, ~1-10s) and slow dynamics (diffusion, ~100-1000s)
- Most common for BMS applications
- Use: real-time estimation, vehicle simulation

**3rd Order and beyond:**
- Additional RC pairs for finer temporal resolution
- Diminishing returns above 2-3 RC pairs for most applications
- Warburg element can replace RC pairs for diffusion (frequency domain)

**Parameter Identification Methods:**
- Electrochemical Impedance Spectroscopy (EIS) - frequency domain fitting
- Pulse discharge/charge tests - time domain fitting
- Recursive least squares (RLS) for online adaptation
- Kalman filter-based (EKF, UKF) joint state-parameter estimation

**Advantages:** Extremely fast computation, easy to parameterize from experiments, good for control applications, well-suited for real-time BMS.

**Limitations:** No internal state information (concentrations, potentials), cannot predict behavior outside parameterization range, no degradation physics, parameters are SOC/temperature/age-dependent requiring lookup tables.

**Key References:**
- Hu, X., Li, S., & Peng, H. (2012). "A comparative study of equivalent circuit models for Li-ion batteries." *J. Power Sources*, 198, 359-367. DOI: 10.1016/j.jpowsour.2011.10.013
- He, H., Xiong, R., & Fan, J. (2011). "Evaluation of lithium-ion battery equivalent circuit models for state of charge estimation by an experimental approach." *Energies*, 4(4), 582-598. DOI: 10.3390/en4040582

---

### 1.5 Reduced-Order Models

**Overview:**
Reduced-order models (ROMs) approximate the full DFN model at lower computational cost while retaining more physics than SPM. Several mathematical techniques are used.

**Proper Orthogonal Decomposition (POD) / Galerkin Projection:**
- Extract dominant modes from full DFN simulation snapshots
- Project governing equations onto reduced basis
- Typical reduction: 3000 unknowns -> 20-50 unknowns
- Retains spatial information in a compressed form
- Challenge: must re-compute basis for new parameters

**Spectral Methods / Polynomial Approximation:**
- Approximate concentration profiles with Chebyshev or Legendre polynomials
- Instead of N_r discrete radial points, use M << N_r polynomial coefficients
- Subramanian et al. (2005): reformulate solid-phase diffusion as polynomial profile with volume-averaged and surface concentrations
- Enables analytical expressions for flux, reducing computational burden significantly

**Transfer Function / Frequency Domain Approaches:**
- Linearize DFN equations around operating point
- Derive transfer functions: V(s)/I(s), SOC(s)/I(s)
- Fast computation via convolution or state-space realization
- Limited to small perturbations around linearization point
- Useful for impedance prediction and EIS interpretation

**State-Space Reformulation:**
- Convert PDEs to state-space form via spatial discretization
- Apply model reduction techniques (balanced truncation, Hankel norm)
- Enables use of control theory tools (observers, controllers)

**Pade Approximation (for diffusion):**
- Approximate the diffusion transfer function with rational polynomial
- Converts transcendental transfer function to finite-order ODE system
- Used in SLIDE's current Chebyshev approach for solid diffusion

**Key References:**
- Subramanian, V. R., et al. (2005). "Efficient macro-micro scale coupled modeling of batteries." *J. Electrochem. Soc.*, 152(10), A2002. DOI: 10.1149/1.2032427
- Cai, L., & White, R. E. (2009). "Reduction of model order based on proper orthogonal decomposition for lithium-ion battery simulations." *J. Electrochem. Soc.*, 156(3), A154. DOI: 10.1149/1.3049347
- Smith, K. A., Rahn, C. D., & Wang, C.-Y. (2007). "Control oriented 1D electrochemical model of lithium ion battery." *Energy Conversion and Management*, 48(9), 2565-2578. DOI: 10.1016/j.enconman.2007.03.015
- Bizeray, A. M., et al. (2019). "Identifiability and parameter estimation of the single particle lithium-ion battery model." *IEEE Trans. Control Syst. Technol.*, 27(5), 1862-1877.

---

### 1.6 3D Electrochemical-Thermal Models

**Overview:**
Full 3D models resolve spatial variations in all three dimensions, coupling electrochemistry with thermal transport across the entire cell or pack geometry.

**When Needed:**
- Large-format cells (pouch, prismatic) with significant temperature gradients
- Tab design optimization (current distribution depends on tab placement)
- Safety/abuse scenario modeling
- Pack-level thermal management design
- Non-uniform aging analysis

**Approaches:**

1. **Pseudo-3D (P4D):** Full DFN at every point on a 2D current collector plane
   - Resolves current distribution across electrode area
   - Captures tab effects and non-uniform utilization
   - Extremely expensive: N_x * N_r * N_2D unknowns

2. **1D+1D+1D approach:** Decouple through-plane electrochemistry (DFN) from in-plane current distribution and thermal transport
   - More tractable; often iteratively coupled
   - Can use different meshes for electrical and thermal

3. **Volume-averaged homogeneous models:** Treat electrode as porous medium with effective properties
   - FEM/FVM on full 3D geometry
   - Source terms from local electrochemical model (SPM or ECM at each mesh point)
   - COMSOL Multiphysics, ANSYS Fluent commonly used

4. **Multi-scale coupling strategies:**
   - Macro-scale: 3D thermal + current distribution (FEM)
   - Meso-scale: 1D DFN at representative points
   - Micro-scale: particle-level diffusion
   - Computational cost managed by adaptive sampling of micro-models

**Computational Cost:** Orders of magnitude higher than 1D models. Full P4D: hours to days. Simplified coupled: minutes to hours.

**Key References:**
- Kim, G.-H., et al. (2011). "Multi-domain modeling of lithium-ion batteries encompassing multi-physics in varied length scales." *J. Electrochem. Soc.*, 158(8), A955. DOI: 10.1149/1.3597614
- Xu, M., et al. (2016). "A pseudo three-dimensional electrochemical-thermal model of a prismatic LFP battery during discharge." *Energy*, 80, 104-117.
- Deng, J., et al. (2018). "General discharge voltage information enabled health evaluation for lithium-ion batteries." *IEEE/ASME Trans. Mechatronics*.

---

## 2. Degradation Models

### 2.1 SEI Growth Models

**Overview:**
The Solid Electrolyte Interphase (SEI) is a passivation layer that forms on the negative electrode surface from electrolyte decomposition. It continues to grow throughout battery life, consuming cyclable lithium and increasing impedance.

**Model Categories:**

**1. Kinetic-limited (reaction rate controlled):**
```
dL_SEI/dt = -M_SEI/(rho_SEI * F) * j_SEI
j_SEI = -j_0_SEI * exp(-alpha * F * eta_SEI / RT)
```
- SEI growth limited by the electrochemical reaction rate
- Predicts linear capacity fade with time^(1/2) under storage
- Appropriate when SEI is thin or porous

**2. Solvent diffusion-limited:**
```
dL_SEI/dt = -M_SEI * k_sol * D_sol / (rho_SEI * L_SEI)
```
- Growth limited by solvent diffusion through existing SEI
- Classic model: `L_SEI ~ sqrt(t)` (parabolic growth law)
- Most widely used; matches many experimental observations
- Predicts capacity loss ~ t^(1/2) for calendar aging

**3. Electron migration / tunneling:**
```
j_SEI = -j_0_SEI * exp(-L_SEI / L_tunnel)
```
- Electron transport through SEI limits growth
- Important for very thin SEI layers
- Gives logarithmic growth at long times

**4. Combined/multi-species models:**
- Multiple SEI layers (inner dense, outer porous)
- Multiple reaction products (Li2CO3, (CH2OCO2Li)2, LiF, etc.)
- Safari & Delacourt (2011): dual-layer SEI with different transport properties
- Yang et al. (2017): SEI with interphase reactions

**Impact on Cell Performance:**
- Capacity fade: consumes cyclable lithium
- Impedance rise: SEI resistance R_SEI ~ L_SEI / kappa_SEI
- Power fade: increased overpotential
- Calendar aging: SEI grows even at rest (temperature-dependent)

**Temperature Dependence:**
- Arrhenius relationship for reaction rates and diffusivity
- SEI growth accelerates significantly at elevated temperatures
- Typical activation energies: 30-80 kJ/mol

**Key References:**
- Pinson, M. B., & Bazant, M. Z. (2013). "Theory of SEI formation in rechargeable batteries." *J. Electrochem. Soc.*, 160(2), A243. DOI: 10.1149/2.044302jes
- Safari, M., Morcrette, M., Teyssot, A., & Delacourt, C. (2009). "Multimodal physics-based aging model for life prediction of Li-ion batteries." *J. Electrochem. Soc.*, 156(3), A145. DOI: 10.1149/1.3043429
- Single, F., Latz, A., & Horstmann, B. (2018). "Identifying the mechanism of continued growth of the SEI." *J. Electrochem. Soc.*, 165(16), A3132. DOI: 10.1149/2.0211816jes
- Ramadass, P., et al. (2004). "Development of first principles capacity fade model for Li-ion cells." *J. Electrochem. Soc.*, 151(2), A196. DOI: 10.1149/1.1634273

---

### 2.2 Lithium Plating

**Overview:**
Lithium plating occurs when lithium deposits as metallic lithium on the anode surface instead of intercalating, typically during charging at low temperatures, high C-rates, or high SOC.

**Thermodynamic Condition:**
```
Plating occurs when: phi_s - phi_e - U_n(c_ss) < 0
(anode potential vs Li/Li+ drops below 0V)
```

**Reversible vs Irreversible:**
- **Reversible:** Freshly plated lithium that can be stripped back during subsequent discharge or rest. Appears as a voltage plateau during relaxation.
- **Irreversible:** Plated lithium that becomes electrically isolated ("dead lithium"), reacts with electrolyte to form additional SEI, or forms dendrites. Permanent capacity loss.
- Typical reversibility: 20-80% depending on conditions

**Plating Models:**

1. **Butler-Volmer plating kinetics:**
   ```
   j_plating = j_0_Li * [exp(alpha_a*F*eta_Li/RT) - exp(-alpha_c*F*eta_Li/RT)]
   eta_Li = phi_s - phi_e - R_SEI*j_total  (overpotential for Li deposition)
   ```

2. **Stripping kinetics:**
   ```
   j_stripping = -j_0_strip * c_Li_plated * exp(-alpha_c*F*eta_Li/RT)  (when eta_Li > 0)
   ```

3. **Dead lithium formation:**
   ```
   dc_Li_dead/dt = k_dead * c_Li_plated
   ```
   First-order irreversible conversion of plated to dead lithium.

**Detection Methods (experimental):**
- Voltage relaxation analysis (plateau in dV/dt)
- Coulombic efficiency measurements
- Differential voltage analysis (dV/dQ)
- In-situ optical microscopy
- Electrochemical impedance spectroscopy
- Reference electrode measurements

**Critical Factors:**
- Low temperature (slow intercalation kinetics, high electrolyte resistance)
- High charging C-rate
- High SOC (low anode potential)
- Anode/cathode capacity ratio close to 1
- Thick anodes (transport limitations)

**Key References:**
- Arora, P., Doyle, M., & White, R. E. (1999). "Mathematical modeling of the lithium deposition overcharge reaction in lithium-ion batteries." *J. Electrochem. Soc.*, 146(10), 3543. DOI: 10.1149/1.1392512
- Yang, X.-G., et al. (2017). "Modeling of lithium plating induced aging of lithium-ion batteries." *Electrochimica Acta*, 243, 272-281. DOI: 10.1016/j.electacta.2017.05.067
- O'Kane, S. E. J., et al. (2022). "Lithium-ion battery degradation: how to model it." *Physical Chemistry Chemical Physics*, 24, 7909. DOI: 10.1039/D2CP00417H

---

### 2.3 Loss of Active Material (LAM)

**Overview:**
LAM refers to the reduction in the amount of electrode material that participates in electrochemical reactions. It reduces capacity by shrinking the available host sites for lithium.

**Mechanisms:**

1. **Mechanical (stress-driven):**
   - Repeated expansion/contraction during cycling causes particle fracture
   - Loss of electrical contact between particles and conductive network
   - Binder degradation and delamination from current collector
   - More severe for high-volume-change materials (Si, Sn)

2. **Chemical:**
   - Transition metal dissolution (especially Mn from NMC/LMO cathodes)
   - Structural phase transitions (e.g., layered to spinel in NMC)
   - Oxygen release at high SOC in Ni-rich cathodes
   - Acid attack from HF generated by LiPF6 decomposition

3. **Electrical isolation:**
   - Particles become disconnected from the conductive network
   - SEI growth on particles can block lithium transport
   - Gas generation creating voids in electrode structure

**Modeling Approaches:**

1. **Empirical/phenomenological:**
   ```
   dLAM/dt = k_LAM * f(stress, SOC, T, cycles)
   ```
   Often: `LAM = k * N^alpha` where N is cycle count, alpha ~ 0.5-1.0

2. **Stress-coupled:**
   ```
   LAM_rate = k * (sigma_max / sigma_critical)^n
   ```
   Where sigma_max is the maximum stress from diffusion-induced stress (DIS) calculation.

3. **Volume-fraction based:**
   ```
   deps_s/dt = -k_LAM * (some driving force)
   ```
   Directly reduces solid volume fraction in porous electrode theory.

**Impact:** Reduces effective electrode capacity. Unlike LLI (loss of lithium inventory from SEI), LAM changes the electrode balance and can shift the operating window asymmetrically.

**Key References:**
- Reniers, J. M., Mulder, G., & Howey, D. A. (2019). "Review and performance comparison of mechanical-chemical degradation models for lithium-ion batteries." *J. Electrochem. Soc.*, 166(14), A3189. DOI: 10.1149/2.0281914jes
- Edge, J. S., et al. (2021). "Lithium ion battery degradation: what you need to know." *Physical Chemistry Chemical Physics*, 23, 8200. DOI: 10.1039/D1CP00359C
- Christensen, J., & Newman, J. (2006). "A mathematical model of stress generation and fracture in lithium manganese oxide." *J. Electrochem. Soc.*, 153(6), A1019. DOI: 10.1149/1.2185287

---

### 2.4 Particle Cracking

**Overview:**
During cycling, lithium intercalation/deintercalation causes volumetric changes in active material particles, generating diffusion-induced stress (DIS) that can lead to crack initiation and propagation.

**Diffusion-Induced Stress (DIS):**

For a spherical particle with isotropic properties:
```
sigma_r = 2*Omega*E / (9*(1-nu)) * [c_avg - c_avg_r(r)]
sigma_theta = Omega*E / (9*(1-nu)) * [2*c_avg + c_avg_r(r) - 3*c(r)]
```
where Omega is partial molar volume, E is Young's modulus, nu is Poisson's ratio, c_avg is volume-average concentration, c_avg_r is average from 0 to r.

**Maximum tensile stress** occurs at particle surface during lithiation:
```
sigma_max ~ Omega * E * Delta_c / (3*(1-nu))
```

**Paris Law for Crack Propagation:**
```
da/dN = C * (Delta_K)^m
```
where `a` is crack length, `N` is cycle count, `Delta_K` is stress intensity factor range, `C` and `m` are material constants.

- `Delta_K = Y * Delta_sigma * sqrt(pi * a)`
- `Y` is geometric factor (~1.12 for surface cracks)
- Crack propagation to critical length causes particle fracture

**Modeling Approaches:**

1. **Analytical stress models** (Cheng & Verbrugge, 2009):
   - Closed-form DIS solutions for spherical particles
   - Fast enough for cell-level simulation
   - Used in SLIDE's current implementation

2. **Fracture mechanics:**
   - Cohesive zone models for crack initiation
   - Phase-field fracture for complex crack patterns
   - Very expensive; typically particle-level only

3. **Statistical approaches:**
   - Weibull distribution for particle strength
   - Probability of fracture as function of stress history
   - Accounts for particle size distribution effects

**Consequences of Cracking:**
- New surface area exposed to electrolyte -> additional SEI formation
- Loss of active material if fragments lose contact
- Increased impedance
- Accelerated degradation (positive feedback loop)

**Key References:**
- Cheng, Y.-T., & Verbrugge, M. W. (2009). "Evolution of stress within a spherical insertion electrode particle under potentiostatic and galvanostatic operation." *J. Power Sources*, 190(2), 453-460. DOI: 10.1016/j.jpowsour.2009.01.021
- Deshpande, R., et al. (2012). "Battery cycle life prediction with coupled chemical degradation and fatigue mechanics." *J. Electrochem. Soc.*, 159(10), A1730. DOI: 10.1149/2.049210jes
- Ai, W., et al. (2020). "A coupled phase field formulation for modelling fatigue cracking in lithium-ion battery electrode particles." *J. Power Sources*, 544, 231805.

---

### 2.5 Electrode-Level Degradation

**Overview:**
Beyond individual particle degradation, the electrode as a porous structure experiences collective degradation effects.

**Porosity Changes:**
```
eps_e = eps_e_0 - (L_SEI * a_s / V_molar_SEI)
```
- SEI growth fills pore space, reducing porosity
- Reduced porosity -> higher electrolyte transport resistance
- Positive feedback: higher overpotential -> more plating/SEI

**Electrolyte Depletion:**
- Electrolyte consumed by SEI formation and other side reactions
- Dry-out in parts of the electrode
- Particularly problematic in cells with lean electrolyte designs
- Can cause sudden capacity cliff (non-linear degradation)

**Binder Degradation:**
- PVDF binder degrades over time and cycling
- Loss of mechanical integrity -> particle disconnection
- Temperature-dependent degradation rate

**Current Collector Corrosion:**
- Aluminum (cathode side): pitting at high voltage
- Copper (anode side): dissolution at low voltage (over-discharge)
- Both increase contact resistance

**Gas Generation:**
- CO2, CO, H2, C2H4 from electrolyte decomposition
- Causes cell swelling (pouch cells)
- Creates dead zones in electrode

**Key References:**
- Kindermann, F. M., et al. (2017). "A SEI modeling approach distinguishing between capacity and power fade." *J. Electrochem. Soc.*, 164(12), E287.
- Lin, X., et al. (2013). "A comprehensive capacity fade model and analysis for Li-ion batteries." *J. Electrochem. Soc.*, 160(10), A1701.

---

### 2.6 Calendar vs Cycle Aging

**Overview:**
Battery aging occurs both during storage (calendar aging) and during use (cycle aging). Models must capture both contributions.

**Calendar Aging:**
- Dominated by SEI growth at the anode
- Depends on: temperature, SOC (anode potential), time
- Classic model: `Q_loss_cal = k_cal * exp(-E_a/RT) * f(SOC) * sqrt(t)`
- At high SOC: faster because lower anode potential drives more SEI formation
- At high temperature: Arrhenius acceleration

**Cycle Aging:**
- Additional degradation from current flow
- Depends on: C-rate, DOD, temperature, voltage window
- Mechanisms: mechanical stress (cracking, LAM), accelerated SEI from cracking, lithium plating
- Classic model: `Q_loss_cyc = k_cyc * exp(-E_a/RT) * g(C-rate, DOD) * N^beta`

**How Models Distinguish:**

1. **Additive (superposition):**
   ```
   Q_loss_total = Q_loss_cal(t, T, SOC) + Q_loss_cyc(N, C-rate, DOD, T)
   ```
   Simple but may not capture interactions.

2. **Physics-based separation:**
   - SEI model naturally captures calendar aging (diffusion-limited growth at rest)
   - Stress/cracking model adds cycle-dependent contribution
   - Plating model adds charge-rate-dependent loss
   - The physics determines the split automatically

3. **Accelerated testing protocols:**
   - Storage tests at different T, SOC -> calendar parameters
   - Cycling tests at different C-rate, DOD -> cycle parameters
   - Matrix experimental design to separate contributions

**Key Interactions:**
- Cycling creates fresh surface area (cracks) -> accelerates calendar SEI growth
- Calendar aging increases impedance -> changes cycling conditions
- Temperature from cycling heat -> accelerates calendar mechanisms

**Key References:**
- Schmalstieg, J., et al. (2014). "A holistic aging model for Li(NiMnCo)O2 based 18650 lithium-ion batteries." *J. Power Sources*, 257, 325-334. DOI: 10.1016/j.jpowsour.2014.02.012
- Petit, M., et al. (2016). "Development of an empirical aging model for Li-ion batteries." *J. Power Sources*, 325, 491-501.
- Naumann, M., et al. (2020). "Analysis and modeling of calendar aging of a commercial LiFePO4/graphite cell." *J. Energy Storage*, 28, 101213.

---

## 3. Thermal Models

### 3.1 Lumped Thermal (0D)

**Overview:**
Assumes uniform temperature throughout the cell. The cell is treated as a single thermal mass exchanging heat with the environment.

**Governing Equation:**
```
m * c_p * dT/dt = Q_gen - Q_dissipation
Q_gen = Q_ohmic + Q_reaction + Q_entropic
Q_dissipation = h * A_surface * (T - T_ambient)
```

**Heat Generation Terms:**
- **Ohmic (Joule) heating:** `Q_ohmic = I^2 * R_internal` (always positive)
- **Reaction (activation) heating:** `Q_reaction = I * sum(eta_j)` (overpotential losses)
- **Entropic (reversible) heating:** `Q_entropic = I * T * dU/dT` (can be positive or negative)
- Combined: `Q_gen = I * (V - U_ocv) + I * T * dU/dT`

**Convection Coefficient:**
- Natural convection: h ~ 5-25 W/(m^2*K)
- Forced air: h ~ 25-250 W/(m^2*K)
- Liquid cooling: h ~ 500-5000 W/(m^2*K)

**When Appropriate:**
- Small cells (18650, 21700) at moderate C-rates
- Preliminary design studies
- When internal temperature gradients are <2-3 K
- Real-time BMS applications

**Limitations:**
- Cannot capture internal hot spots
- Inaccurate for large-format cells
- Misses tab heating effects

**This is the approach currently used in SLIDE.**

**Key References:**
- Bernardi, D., Pawlikowski, E., & Newman, J. (1985). "A general energy balance for battery systems." *J. Electrochem. Soc.*, 132(1), 5. DOI: 10.1149/1.2113792

---

### 3.2 1D Through-Plane Thermal

**Overview:**
Resolves temperature variation through the cell thickness (perpendicular to electrode layers). Important for pouch and prismatic cells with many stacked layers.

**Governing Equation:**
```
rho * c_p * dT/dt = d/dx(k_x * dT/dx) + Q_gen(x)
```

**Key Features:**
- Anisotropic thermal conductivity: k_through-plane << k_in-plane
  - Typical: k_through ~ 0.5-2 W/(m*K), k_in-plane ~ 20-40 W/(m*K)
- Layer structure: Cu | Anode | Separator | Cathode | Al | Cathode | Separator | Anode | Cu ...
- Each layer has different thermal properties
- Heat generation concentrated in electrodes

**Approaches:**
1. **Layer-resolved:** Discrete thermal resistance/capacitance for each layer
2. **Homogenized:** Effective properties for repeated unit cell
3. **Mixed:** Resolve electrode pairs, homogenize within

**When to use:** Thick stacks, high C-rates where through-plane gradient exceeds 2-3K, understanding heat flow to cooling surface.

---

### 3.3 3D Thermal

**Overview:**
Full 3D thermal simulation resolving temperature fields in all directions. Essential for large-format cells and pack-level thermal management.

**Governing Equation:**
```
rho * c_p * dT/dt = div(k * grad(T)) + Q_gen(x,y,z)
```

**CFD Coupling:**
- For liquid cooling channels: solve Navier-Stokes + energy in coolant
- For air cooling: turbulent flow models (k-epsilon, k-omega SST)
- Conjugate heat transfer at solid-fluid interfaces

**Common Tools:**
- COMSOL Multiphysics (FEM-based, multi-physics coupling)
- ANSYS Fluent / CFX (CFD with battery models)
- OpenFOAM (open-source CFD)
- Star-CCM+ (commercial CFD)

**Tab Heating:**
- Current concentrates near tabs -> local hot spots
- Tab design significantly affects temperature uniformity
- Models must include current collector resistance and Joule heating

**Pack-Level Considerations:**
- Cell-to-cell thermal coupling (conduction through busbars, radiation)
- Cooling system design (channel geometry, flow rate)
- Thermal runaway propagation between cells

**Computational Cost:** Minutes to hours per simulation depending on mesh resolution and physics included.

---

### 3.4 Thermal Runaway / Abuse Modeling

**Overview:**
Thermal runaway is a self-heating process where exothermic reactions in the cell become self-sustaining, leading to temperatures exceeding 800+ C and potentially fire/explosion.

**Reaction Sequence:**
1. **SEI decomposition** (~90-120 C): Exothermic, releases gases
2. **Anode-electrolyte reaction** (~120-250 C): Intercalated lithium reacts with electrolyte
3. **Separator shutdown** (~130 C PE, ~165 C PP): Pores close, shutting off ion transport
4. **Separator melt-through** (~200+ C): Internal short circuit
5. **Cathode decomposition** (~200-300 C): Oxygen release, further exothermic reactions
6. **Electrolyte decomposition** (~200-300 C): Vaporization and combustion

**Modeling Approaches:**

1. **Hatchard/Dahn model (Arrhenius kinetics):**
   ```
   Q_sei = H_sei * A_sei * c_sei * exp(-E_a_sei / RT)
   dc_sei/dt = -A_sei * c_sei * exp(-E_a_sei / RT)
   ```
   Similar equations for each reaction with different activation energies and pre-exponential factors.

2. **Calorimetry-based (ARC data fitting):**
   - Accelerating Rate Calorimetry provides self-heating rate vs temperature
   - Fit Arrhenius parameters to experimental curves
   - Practical but cell-specific

3. **Multi-physics abuse models:**
   - Mechanical deformation (nail penetration, crush) -> internal short circuit
   - Short circuit resistance depends on penetration geometry
   - Coupled thermal-electrical-mechanical

**Propagation Modeling:**
- Cell-to-cell heat transfer in module/pack
- Vent gas transport and ignition
- Radiation between cells (T^4 dependence significant at high T)
- Design for containment: thermal barriers, spacing, venting

**Key References:**
- Hatchard, T. D., et al. (2001). "Thermal model of cylindrical and prismatic lithium-ion cells." *J. Electrochem. Soc.*, 148(7), A755. DOI: 10.1149/1.1377592
- Feng, X., et al. (2018). "Thermal runaway mechanism of lithium ion battery for electric vehicles: A review." *Energy Storage Materials*, 10, 246-267. DOI: 10.1016/j.ensm.2017.05.013
- Coman, P. T., et al. (2017). "A lumped model of venting during thermal runaway in a cylindrical lithium cobalt oxide lithium-ion cell." *J. Power Sources*, 307, 56-62.

---

## 4. Numerical Methods

### 4.1 Spatial Discretization

**Finite Volume Method (FVM):**
- Conservative by construction (fluxes balanced at cell interfaces)
- Natural for transport equations in porous media
- Most common for through-thickness direction in DFN
- SLIDE uses FVM for some spatial discretization
- Easy to handle discontinuous properties at interfaces (electrode/separator)

**Finite Element Method (FEM):**
- Better for complex geometries (3D thermal, mechanical stress)
- Higher-order elements for smooth solutions
- COMSOL uses FEM as its core numerical engine
- More complex implementation but very flexible

**Spectral Methods:**
- Represent solution as sum of basis functions (Chebyshev, Legendre polynomials)
- Exponential convergence for smooth solutions
- Excellent for spherical particle diffusion (smooth, well-behaved)
- **Chebyshev collocation:** solve at Chebyshev nodes, use polynomial interpolation
  - SLIDE uses/is implementing Chebyshev methods for radial diffusion
  - Advantages: far fewer nodes needed (5-10 vs 20-50 for FVM) for same accuracy
  - Differentiation via matrix-vector multiply
  - Handles boundary conditions naturally via basis modification

**Spectral Methods Detail (relevant to SLIDE):**

For the radial diffusion equation in a sphere:
```
dc/dt = D/r^2 * d/dr(r^2 * dc/dr)
```

Transform using `u = r*c` to get standard diffusion:
```
du/dt = D * d^2u/dr^2
```

Chebyshev discretization:
- Map r in [0, R] to xi in [-1, 1]
- Approximate u(xi) = sum(a_k * T_k(xi))
- Derivatives via Chebyshev differentiation matrices
- Time integration of resulting ODE system

**Advantages of Chebyshev for SLIDE:**
- 5-10 Chebyshev points can match 30+ FVM points in accuracy
- Dense differentiation matrix but very small (5x5 to 10x10)
- Spectral accuracy for the smooth concentration profiles in particles
- Well-suited for the SPM where particle diffusion dominates cost

**Key References:**
- Trefethen, L. N. (2000). *Spectral Methods in MATLAB*. SIAM.
- Bizeray, A. M., et al. (2016). "Lithium-ion battery thermal-electrochemical model-based state estimation using orthogonal collocation and a modified extended Kalman filter." *J. Power Sources*, 296, 400-412.
- Subramanian, V. R., et al. (2005). "Efficient macro-micro scale coupled modeling of batteries." *J. Electrochem. Soc.*, 152(10), A2002.

---

### 4.2 Time Integration

**Explicit Methods (Forward Euler, RK4):**
```
y_{n+1} = y_n + dt * f(t_n, y_n)
```
- Simple implementation
- Conditionally stable: dt < dx^2 / (2*D) for diffusion (very restrictive!)
- Not suitable for stiff electrochemical systems without very small time steps
- RK4: 4th order, still explicit stability limit but better accuracy per step

**Implicit Methods (Backward Euler, Crank-Nicolson, BDF):**
```
y_{n+1} = y_n + dt * f(t_{n+1}, y_{n+1})  (Backward Euler)
```
- Unconditionally stable (or A-stable)
- Requires solving nonlinear system at each step (Newton iteration)
- BDF (Backward Differentiation Formulas): multi-step, up to order 5
  - BDF1 = Backward Euler, BDF2 most popular
  - Used by SUNDIALS CVODE/IDA
  - Variable-order, variable-step implementations available

**Adaptive Time Stepping:**
- Embedded RK methods (RK45 Dormand-Prince) for error estimation
- Step size control: `dt_new = dt * (tol/error)^(1/(p+1))`
- Essential for battery simulation where dynamics change (CC phase smooth, CV phase stiff)
- SUNDIALS provides excellent adaptive stepping

**SLIDE's Current Approach:**
- Uses fixed time steps (dt parameter in Cycler)
- User specifies dt; smaller dt = more accurate but slower
- Opportunity: adaptive stepping could significantly improve efficiency

---

### 4.3 DAE Solvers

**Overview:**
The full DFN model and even SPM with algebraic constraints form Differential-Algebraic Equation (DAE) systems, not pure ODEs. The algebraic constraints come from:
- Electrolyte potential (no time derivative)
- Solid-phase potential (no time derivative)
- Butler-Volmer kinetics coupling
- Voltage constraints (CV mode)

**SUNDIALS IDA:**
- Industry-standard DAE solver from Lawrence Livermore National Lab
- Variable-order BDF methods (1-5)
- Adaptive time stepping with error control
- Handles index-1 DAEs (most battery models)
- Sparse linear algebra support (KLU, SuperLU)
- C library with Python/MATLAB interfaces
- Used by PyBaMM as primary solver
- URL: https://computing.llnl.gov/projects/sundials

**CasADi:**
- Symbolic framework for nonlinear optimization and ODE/DAE integration
- Automatic differentiation (forward and adjoint modes)
- Interface to SUNDIALS, IPOPT, and other solvers
- Python and MATLAB interfaces
- Used by PyBaMM for model creation and sensitivity analysis
- Particularly useful for parameter estimation and optimal control
- URL: https://web.casadi.org/

**MATLAB ode15s / ode15i:**
- ode15s: stiff ODE solver (variable-order BDF/NDF)
- ode15i: fully implicit DAE solver
- Convenient for prototyping but slower than compiled SUNDIALS
- Good for validation and comparison

**PETSc TS:**
- Scalable parallel time-stepping for large systems
- Supports many implicit/explicit methods
- Overkill for single-cell models but useful for 3D or pack-level

**Key References:**
- Hindmarsh, A. C., et al. (2005). "SUNDIALS: Suite of nonlinear and differential/algebraic equation solvers." *ACM Trans. Math. Softw.*, 31(3), 363-396. DOI: 10.1145/1089014.1089020
- Andersson, J. A. E., et al. (2019). "CasADi: a software framework for nonlinear optimization and optimal control." *Math. Prog. Comp.*, 11, 1-36. DOI: 10.1007/s12532-018-0139-4

---

### 4.4 Fast Approximations

**Lookup Tables (LUTs):**
- Pre-compute OCV(SOC), D_s(SOC, T), k(SOC, T) on grids
- Interpolate at runtime (linear, cubic spline)
- Eliminates expensive function evaluations
- Memory vs speed trade-off
- SLIDE uses LUTs for OCV curves

**Polynomial Fits:**
- Approximate OCV(SOC) with polynomial (order 8-12 typical)
- Differentiation trivial (needed for dU/dSOC, dU/dT)
- Can be evaluated with Horner's method (fast, stable)
- Risk of oscillation for high-order polynomials (Runge's phenomenon)
- Rational polynomial fits more robust for some functions

**Transfer Function / State-Space Approaches:**
- Linearize around operating point
- Represent as discrete-time state-space system
- Matrix exponential for exact discrete-time propagation
- Very fast: matrix-vector multiplies only
- Limited accuracy for large deviations from operating point

**Pade Approximation for Diffusion:**
- Exact diffusion transfer function: `G(s) = tanh(sqrt(s*tau_d)) / sqrt(s*tau_d)`
- Pade approximation: rational polynomial P(s)/Q(s) of order [m/n]
- Converts infinite-dimensional diffusion to finite-order ODE
- [3/3] or [4/4] Pade gives good accuracy for moderate frequencies

**Polynomial Profile Approximation (Subramanian et al.):**
- Assume concentration profile shape: `c_s(r) = a + b*(r/R)^2 + c*(r/R)^4 + ...`
- Determine coefficients from conservation equations and boundary conditions
- Reduces PDE to 2-3 ODEs (volume-average and surface concentration)
- Very fast; used in many real-time implementations

---

## 5. Emerging Approaches (2024-2025)

### 5.1 Physics-Informed Neural Networks (PINNs) for Battery Modeling

**Overview:**
PINNs embed physical laws (PDEs) into the loss function of neural networks, enabling them to solve forward and inverse problems while respecting governing equations.

**Application to Batteries:**

1. **Forward simulation:**
   - Train NN to approximate DFN solution
   - Loss = data mismatch + PDE residual + boundary conditions
   - Once trained, inference is very fast (milliseconds)
   - Challenge: training can be slow and requires careful hyperparameter tuning

2. **Parameter identification:**
   - PINNs naturally solve inverse problems
   - Given voltage data, infer diffusivity, reaction rate, etc.
   - No need for adjoint derivation; automatic differentiation provides gradients
   - Shown to work for SPM parameter identification

3. **State estimation:**
   - Online SOC/SOH estimation using PINN as surrogate
   - Constrains estimates to physically plausible values
   - More robust than pure data-driven approaches

**Challenges:**
- Training instability (balancing PDE loss terms)
- Poor performance on stiff systems (battery models are stiff)
- Generalization to unseen operating conditions
- Computational cost of training vs traditional simulation
- Difficulty with sharp fronts or discontinuities

**Key References:**
- Raissi, M., Perdikaris, P., & Karniadakis, G. E. (2019). "Physics-informed neural networks." *J. Comput. Phys.*, 378, 686-707. DOI: 10.1016/j.jcp.2018.10.045
- Hofmann, T., et al. (2023). "Physics-informed neural networks for battery modeling." *J. Electrochem. Soc.*, 170, 090524.
- Nascimento, R. G., et al. (2023). "Hybrid physics-informed neural networks for lithium-ion battery modeling and prognosis." *J. Power Sources*, 513, 230526.

---

### 5.2 Digital Twins

**Overview:**
A digital twin is a virtual replica of a physical battery that is continuously updated with real-world data to mirror its current state and predict future behavior.

**Architecture:**
```
Physical Battery -> Sensors -> Data Pipeline -> Digital Twin Model -> Predictions
                                                      ^                    |
                                                      |                    v
                                                State Estimation    Decision Support
                                                (Kalman filter,     (charging strategy,
                                                 particle filter)    maintenance, EOL)
```

**Key Components:**
1. **Physics model:** SPM, ECM, or hybrid as the backbone
2. **State estimator:** EKF, UKF, or particle filter for real-time state tracking
3. **Parameter updater:** Recursive least squares or ML for aging parameter adaptation
4. **Prognostics engine:** Remaining useful life (RUL) prediction

**2024-2025 Developments:**
- Cloud-based digital twins for fleet management (automotive OEMs)
- Edge computing for on-vehicle digital twins
- Federated learning for fleet-wide model improvement
- Integration with vehicle CAN bus data
- ISO/IEC 30173:2024 standard for digital twin frameworks

**Challenges:**
- Model fidelity vs computational cost for real-time operation
- Sensor noise and data quality
- Communication bandwidth for cloud-based twins
- Validation across diverse operating conditions

**Key References:**
- Merkle, L., et al. (2023). "Digital twin of a lithium-ion battery." *Batteries*, 9(4), 222.
- Li, W., et al. (2024). "Digital twin for battery systems: Cloud battery management system with online state-of-charge and state-of-health estimation." *J. Energy Storage*, 44, 103391.

---

### 5.3 Multi-Scale Modeling

**Overview:**
Multi-scale approaches connect phenomena at different length/time scales: atomistic (nm) -> particle (um) -> electrode (mm) -> cell (cm) -> module/pack (m).

**Scale Hierarchy:**

| Scale | Length | Phenomena | Methods |
|-------|--------|-----------|---------|
| Atomistic | 0.1-10 nm | Crystal structure, ion hopping, SEI chemistry | DFT, MD |
| Particle | 1-50 um | Diffusion, stress, phase transitions | FEM, spectral |
| Electrode | 50-200 um | Porous media transport, reaction distribution | DFN (P2D) |
| Cell | 1-50 cm | Current distribution, thermal gradients | FEM, FVM |
| Pack | 0.1-2 m | Thermal management, electrical connections | CFD, lumped |

**Coupling Strategies:**
1. **Sequential (offline):** Run lower scale, extract parameters for higher scale
2. **Concurrent (online):** Solve all scales simultaneously with information passing
3. **Hierarchical:** Coarse scale drives fine scale; fine scale provides constitutive response

**2024-2025 Focus Areas:**
- Machine-learned interscale bridging functions
- Microstructure-resolved electrode models (from X-ray CT images)
- Crystal-level phase-field models informing particle-level stress
- Automated parameter passing between scales

**Key References:**
- Franco, A. A., et al. (2019). "Boosting rechargeable batteries R&D by multiscale modeling." *Chemical Reviews*, 119(7), 4569-4627. DOI: 10.1021/acs.chemrev.8b00239
- Mistry, A., et al. (2021). "A minimal information set to enable verifiable theoretical battery research." *ACS Energy Lett.*, 6(11), 3831-3835.

---

### 5.4 Data-Driven Degradation

**Overview:**
Machine learning approaches for battery degradation prediction that learn directly from cycling data without explicit physics models.

**Approaches:**

1. **Feature engineering + classical ML:**
   - Extract features from charge/discharge curves (dQ/dV peaks, capacity at specific voltages, IC/DV curves)
   - Random forest, gradient boosting, SVR for RUL prediction
   - Severson et al. (2019): predicted cycle life from first 100 cycles with <10% error using features from early cycles

2. **Deep learning:**
   - LSTM/GRU for sequence modeling of capacity fade
   - CNN for extracting features from voltage curves
   - Transformer architectures for attention over cycling history
   - Autoencoders for anomaly detection (sudden degradation)

3. **Transfer learning:**
   - Train on one cell chemistry, fine-tune on another
   - Reduces data requirements for new chemistries
   - Domain adaptation techniques for lab-to-field transfer

4. **Gaussian Process Regression:**
   - Provides uncertainty quantification naturally
   - Good for small datasets
   - Kernel design can encode physics (e.g., monotonic degradation)

**Key References:**
- Severson, K. A., et al. (2019). "Data-driven prediction of battery cycle life before capacity degradation." *Nature Energy*, 4, 383-391. DOI: 10.1038/s41560-019-0356-8
- Attia, P. M., et al. (2020). "Closed-loop optimization of fast-charging protocols for batteries with machine learning." *Nature*, 578, 397-402. DOI: 10.1038/s41586-020-1994-5
- Roman, D., et al. (2021). "Machine learning pipeline for battery state-of-health estimation." *Nature Machine Intelligence*, 3, 447-456.

---

### 5.5 Hybrid Models (Physics + ML)

**Overview:**
Hybrid approaches combine the interpretability and extrapolation capability of physics-based models with the flexibility and speed of machine learning.

**Architectures:**

1. **Residual learning:**
   ```
   y_total = y_physics(x) + y_ML(x, y_physics)
   ```
   - Physics model provides baseline prediction
   - ML learns the residual (model error)
   - Widely used and straightforward to implement
   - Example: SPM + neural network for unmodeled dynamics

2. **Neural ODE / Universal Differential Equations:**
   ```
   dx/dt = f_physics(x) + f_NN(x)
   ```
   - Replace unknown terms in ODEs with neural networks
   - Train end-to-end with differentiable ODE solvers
   - Julia SciML ecosystem pioneered this
   - Applied to learn unknown degradation mechanisms

3. **Physics-constrained ML:**
   - Enforce physical constraints (monotonicity, conservation laws) in ML architecture
   - Monotonic neural networks for SOC-OCV relationship
   - Conservation-preserving architectures for electrochemical models

4. **ML-accelerated physics:**
   - Train surrogate models on physics simulation data
   - Neural operators (DeepONet, Fourier Neural Operator) for PDE solutions
   - 100-1000x speedup over direct simulation
   - Enables real-time use of DFN-fidelity models

**2024-2025 Trends:**
- Foundation models for batteries (pre-trained on diverse battery data)
- Differentiable simulation frameworks (JAX-based battery models)
- Active learning for optimal experiment design
- Uncertainty-aware hybrid models (Bayesian approaches)
- Integration into BMS firmware (TinyML for edge deployment)

**Key References:**
- Rackauckas, C., et al. (2021). "Universal differential equations for scientific machine learning." *arXiv:2001.04385*.
- Chen, R. T. Q., et al. (2018). "Neural ordinary differential equations." *NeurIPS 2018*.
- Bills, A., et al. (2023). "Universal battery performance and degradation model for electric aircraft." *arXiv:2008.01527*.

---

## 6. Relevance to SLIDE

### Current SLIDE Capabilities vs State of the Art

| Feature | SLIDE Current | State of Art | Gap/Opportunity |
|---------|--------------|--------------|-----------------|
| Electrochemical model | SPM | SPM/SPMe/DFN | SPMe would be valuable addition |
| Solid diffusion | Chebyshev (in progress) | Chebyshev/spectral | Aligned with best practice |
| Degradation | SEI, LAM, cracking, plating | Same + electrode-level | Good coverage; could add porosity effects |
| Thermal | Lumped (0D) | 0D to 3D | 1D thermal for large cells |
| Time integration | Fixed step | Adaptive BDF (SUNDIALS) | Adaptive stepping high priority |
| ECM | 1-3 RC pairs | Same | Good coverage |
| Pack simulation | Module_s/p + Battery | Similar | Good hierarchical design |
| Python bindings | Planned | PyBaMM standard | Important for adoption |

### Recommended Priorities for SLIDE

1. **SPMe implementation** - Most impactful model addition for extending to higher C-rates
2. **Adaptive time stepping** - Performance gain with minimal accuracy loss
3. **Chebyshev completion** - Already in progress; aligns with best practice for spectral methods in particle diffusion
4. **Electrode-level degradation** - Porosity change from SEI growth is straightforward to add
5. **Python bindings** - Critical for adoption and integration with data-driven approaches
6. **Hybrid model hooks** - Allow ML residual corrections to be plugged into physics models

### Key Competing/Complementary Software

| Software | Language | Models | Strengths |
|----------|----------|--------|-----------|
| **PyBaMM** | Python | SPM, SPMe, DFN, ECM | Flexible, symbolic, CasADi/SUNDIALS solvers |
| **DUALFOIL** | Fortran | DFN | Newman group original; reference implementation |
| **COMSOL** | MATLAB/Java | DFN, 3D | GUI, multi-physics, commercial |
| **LIONSIMBA** | MATLAB | DFN | Open-source, well-documented |
| **PETLION** | Julia | DFN | Fast, differentiable, modern |
| **BattMo** | MATLAB | DFN, 3D | SINTEF, complex geometries |
| **SLIDE** | C++ | SPM, ECM | **Fast, degradation focus, pack-level** |

SLIDE's key differentiator is its C++ performance enabling large-scale degradation studies and pack-level simulation that would be prohibitively slow in Python/MATLAB. The hierarchical StorageUnit architecture is well-designed for multi-scale modeling from cell to pack.

---

## Key Textbooks and Review Papers

1. Newman, J., & Thomas-Alyea, K. E. (2004). *Electrochemical Systems* (3rd ed.). Wiley. -- The foundational text.
2. Plett, G. L. (2015). *Battery Management Systems, Volume I: Battery Modeling*. Artech House. -- Comprehensive ECM and SPM treatment.
3. Plett, G. L. (2016). *Battery Management Systems, Volume II: Equivalent-Circuit Methods*. Artech House.
4. O'Kane, S. E. J., et al. (2022). "Lithium-ion battery degradation: how to model it." *PCCP*, 24, 7909. DOI: 10.1039/D2CP00417H -- **Excellent recent review of all degradation models.**
5. Edge, J. S., et al. (2021). "Lithium ion battery degradation: what you need to know." *PCCP*, 23, 8200. DOI: 10.1039/D1CP00359C -- **Comprehensive degradation review.**
6. Brosa Planella, F., et al. (2022). "A continuum of physics-based lithium-ion battery models reviewed." *Prog. Energy*, 4, 042003. DOI: 10.1088/2516-1083/ac7d31 -- **Outstanding review linking SPM -> SPMe -> DFN.**
7. Sulzer, V., et al. (2021). "Python Battery Mathematical Modelling (PyBaMM)." *JORS*, 9(1), 14. DOI: 10.5334/jors.309
8. Franco, A. A., et al. (2019). "Boosting rechargeable batteries R&D by multiscale modeling." *Chem. Rev.*, 119, 4569. DOI: 10.1021/acs.chemrev.8b00239

---

*Note: Web access was unavailable during report generation. This report is based on established literature and domain knowledge current through early 2025. For the very latest 2025-2026 developments, web research with the WebFetch tool would be needed (permission was denied during this session).*
