# MonteCarloCollisionSim
This repository implements a Monte Carlo pipeline to estimate the collision probability between two orbiting objects (satellites / debris) using Conjunction Data Message (CDM) inputs and orbital propagation.
---

## Overview

As orbital density increases in Low Earth Orbit (LEO), probabilistic conjunction assessment becomes essential for mission safety. Nominal miss distance alone is insufficient to characterize collision risk — state uncertainty must be incorporated.

This repository implements a transparent and reproducible Monte Carlo pipeline that:

1. Parses Conjunction Data Messages (CDMs)
2. Samples from reported state covariance
3. Propagates sampled trajectories to time of closest approach (TCA)
4. Estimates collision probability from simulated outcomes

The goal is to demonstrate rigorous stochastic modeling applied to real-world aerospace safety problems.

---

## Methodology

### 1. CDM Parsing

Each CDM provides:

- Nominal state vectors (position & velocity)
- 6×6 covariance matrix
- Time of closest approach (TCA)
- Reference frame (e.g., ECI or RTN)
- Hard-body radius (if available)

The parser:

- Normalizes units
- Ensures consistent reference frames
- Constructs internal state objects for propagation

---

### 2. Uncertainty Modeling

We assume Gaussian state uncertainty:

x ~ N(x_nom, C)

Using Cholesky decomposition:

x_sample = x_nom + Lz  
where C = LLᵀ and z ~ N(0, I)

This produces Monte Carlo samples consistent with the reported covariance structure.

---

### 3. Orbital Propagation

Each sampled state is propagated to TCA using either:

- SGP4 (TLE-driven objects), or  
- Numerical integration (optional higher-fidelity models)

Propagation is vectorized and batch-processed for efficiency.

---

### 4. Collision Criterion

At TCA:

- Relative position: r_rel = r_A − r_B  
- Miss distance: d = ||r_rel||  
- Collision if: d ≤ R_A + R_B  

Estimated collision probability:

P_collision ≈ N_hits / N_samples

Confidence intervals are computed using binomial approximations.

---

## Example Usage

```bash
python -m src.jspoc_monte \
  --cdm data/example_cdm/case1.xml \
  --samples 100000 \
  --propagator sgp4 \
  --out results/run1.json
