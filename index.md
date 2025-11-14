# Introduction

Welcome to the hands-on training on adaptive mesh refinement methods applied to finite volume methods with the C++ library **samurai**. In a sequence of focused steps you will go from creating multi-resolution meshes and fields to implementing robust flux-based solvers for Burgers and the Euler equations, finishing with advanced shock problems (double Mach reflection) and custom positivity-preserving adaptation.

## Learning Path
1. Environment & tools: build, run, visualize (ParaView & Python).
2. Mesh & fields: multi-resolution meshes, interval-based loops, Gaussian initialization.
3. Naive Burgers: first FV update, discovering conservation pitfalls at level interfaces.
4. Flux mechanism: conservative / non-conservative fluxes, diffusion, vector extension.
5. Euler (part 1): Rusanov, HLL, HLLC, primitive ↔ conservative transforms, CFL control.
6. Euler (part 2): custom boundary conditions, double Mach reflection, safe prediction operator.

## You Will Be Able To
- Build conservative & non-conservative finite volume schemes rapidly.
- Use adaptive multi-resolution to concentrate resolution only where needed.
- Plug advanced Riemann solvers (Rusanov / HLL / HLLC) without rewriting infrastructure.
- Maintain physical validity (ρ, p > 0) under adaptation.
- Visualize and analyze results efficiently in 1D/2D.

## Key Strengths of samurai
- Compressed interval mesh representation (memory & speed).
- Dimension-agnostic flux definitions (same code 1D→2D→3D).
- Modular extension (new fluxes, BCs, operators) with minimal boilerplate.
- Multi-resolution error control independent of PDE specifics.

## Prerequisites
Basic C++ (templates, lambdas), finite volume concepts, CMake toolchain. ParaView and Python+matplotlib for visualization.

## Getting Ready
To start this hands-on, you have to install a proper development environment with the samurai library. You can follow the [Setup Instructions](00-setup.md) to install samurai via conda and verify the provided setup example. If you run into issues, reach out.

For each part, you have the solution files provided to help you check your implementation. The location of the solution is in the `solution` folder of the part name in the `practical_session/material` folder.
