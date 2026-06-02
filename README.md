# Quasi-1D Euler Solver for Convergent-Divergent Nozzle Flow

## Abstract

This project presents a quasi-1D C++ solver for convergent-divergent 
nozzle flow, developed from scratch. The goal is to solve the quasi-1D Euler 
equations with a geometric source term using the Rusanov flux scheme. 
The solver yields good agreement with isentropic theory: M = 1.0 at 
the throat and M = 3.3 at the exit, against an analytical value of 
M = 3.35 — a 1.5% error on the exit Mach number.

## Physics & Governing Equations
The solver resolves the quasi-1D Euler equations with a geometric source term,
in conservative vector form:

    ∂U/∂t + (1/A) · ∂(F·A)/∂x = S
where:

    U = (ρ, ρu, ρE) — conservative state vector
    F = (ρu, ρu²+p, (ρE+p)u) — flux vector
    S = (0, (p/A)·dA/dx, 0) — geometric source term
    A(x) = 1 + 2.2·(x − 1.5)² — nozzle cross-section (Anderson test case)

The analytical validation uses the isentropic area-Mach relation:

    A/A* = (1/M) · [(2/(γ+1)) · (1 + (γ−1)/2 · M²)]^((γ+1)/(2(γ−1)))

## Numerical Method

The Rusanov scheme estimates the interface flux as the average of left and 
right physical fluxes, stabilised by a dissipation term proportional to 
the maximum local wave speed. The time update reads:

    U_j^{n+1} = U_j^n − (Δt/Δx)·(1/A_j)·(F̂_{j+1}·A_{j+1} − F̂_j·A_j) + Δt·S_j
where F̂ denotes the Rusanov flux at cell interfaces. The time integration is 
explicit (forward Euler) with a CFL number of 0.5. At the subsonic inlet, 
total conditions are imposed (P₀ = 101 325 Pa, T₀ = 300 K, M = 0.1). 
At the supersonic outlet, all three characteristics exit the domain, so 
all variables are extrapolated from the interior (zero-gradient condition).

## Result and validation 

The steady-state Mach profile increases from M = 0.1 at the inlet to M ≈ 3.3 
at the exit, in good agreement with the isentropic area-Mach relation which predicts 
M = 3.35 (1.5% error). Total conditions are imposed at the subsonic inlet, while the 
supersonic outlet uses a zero-gradient extrapolation — the same approach as zeroGradient 
in OpenFOAM. Crucially, the sonic condition M = 1 at the throat is not imposed but emerges 
naturally from the solution, confirming the physical consistency of the solver.

![Validation](/results/validation.png)

## Build & Run

```bash
mkdir build && cd build
cmake ..
make
./Solver1D_TuyereLaval
```
Requires: C++17 compiler, CMake ≥ 3.16
## Project Structure

```
├── CMakeLists.txt
├── main.cpp
├── src/
│   ├── StateVector.h
│   ├── Cellule.h / .cpp
│   ├── Domaine.h / .cpp
│   └── NumericalSchemes.h / .cpp
├── scripts/
│   └── GifCreation.py
│   └── Validation.py
└── README.md
```

## Author & Context

Developed by Mathis TOPPAN, first-year engineering student at IMT Mines Albi (France). 
This project was entirely self-initiated as part of a personal training programme in computational 
fluid dynamics, with the long-term goal of understanding the internal mechanics of CFD solvers at the 
source-code level.

## Perspectives

- **HLLC flux scheme** — replace Rusanov for improved contact discontinuity resolution
- **MUSCL reconstruction** with minmod slope limiter — second-order spatial accuracy
- **Heat source term Q̇** in the energy equation, validated against Rayleigh flow theory — a stepping stone toward reacting flow simulations
- **2D structured Euler solver** with VTK export for ParaView visualisation (separate repository)

This solver is part of a broader self-initiated research programme 
in computational fluid dynamics, targeting a research internship at the von Karman Institute 
(VKI) on reacting flows.

## References

- Anderson, J.D. (1995). *Computational Fluid Dynamics: The Basics with Applications*. McGraw-Hill.
- Toro, E.F. (2009). *Riemann Solvers and Numerical Methods for Fluid Dynamics*. Springer.