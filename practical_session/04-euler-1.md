# Euler equations (part 1)

:::{note} Main Objectives
- Understand the Euler equations for compressible gas dynamics
- Implement three Riemann solvers: Rusanov, HLL, and HLLC
- Apply adaptive mesh refinement to capture shocks and contact discontinuities
- Visualize density, pressure, and velocity fields
:::

In this final part of the practical session, you will implement a complete compressible gas dynamics solver using the Euler equations. This builds upon everything you have learned: mesh creation, field manipulation, flux computation, and adaptive mesh refinement.

```{include} start_instructions.md
```

The Euler equations describe the motion of inviscid (non-viscous) compressible fluids. They are fundamental in aerodynamics, astrophysics, and shock physics. Unlike Burgers, the Euler system couples four conservation laws (mass, momentum in x and y, energy), requiring more sophisticated numerical schemes.

## Euler equations

The Euler equations in 2D read:

$$
\frac{\partial }{\partial t} \begin{pmatrix}
\rho \\\\ \rho u \\\\ \rho v \\\\ \rho E
\end{pmatrix}
+ \frac{\partial }{\partial x} \begin{pmatrix}
\rho u \\\\ \rho u^2 + p \\\\ \rho uv \\\\ u(\rho E + p)
\end{pmatrix}
+ \frac{\partial }{\partial y} \begin{pmatrix}
\rho v \\\\ \rho uv \\\\ \rho v^2 + p \\\\ v(\rho E + p)
\end{pmatrix}
= 0
$$

where:
- $\rho$ is the density (mass per unit volume)
- $u, v$ are the velocity components in x and y directions
- $\rho E$ is the total energy per unit volume (where $E$ is the specific total energy)
- $p$ is the pressure

where $E = e + \frac{1}{2}(u^2 + v^2)$ is the **specific** total energy (energy per unit mass), with $e$ the internal energy per unit mass.

The pressure is given by the equation of state for an ideal gas:

$$
p = (\gamma - 1) \rho e
$$

where $\gamma$ is the ratio of specific heats (typically 1.4 for air).

### Sound Speed

The sound speed (or speed of sound) is defined as:

$$
c = \sqrt{\frac{\gamma p}{\rho}}
$$

**Physical meaning:** This is the speed at which small pressure disturbances (acoustic waves) propagate through the fluid. It is crucial for:
- Determining the time step (CFL condition)
- Estimating wave speeds in Riemann solvers
- Identifying supersonic vs subsonic flow regimes (Mach number $M = |u|/c$)

## Explanation of the provided code

The Euler code becomes more complex due to the number of equations and the coupling between them. To save time, we provide you with starter code that you can find in the `euler` folder. We will now explain the different parts of the code that will help you in the next steps.

The configuration of the mesh and fields is defined in `config.hpp`.

### Naming convention of the variables

We have defined the numbering of the components of the vector field as follows:
- Component 0: $\rho$ (density)
- Component 1: $\rho E$ (total energy per unit volume)
- Component 2: $\rho u$ (momentum in x)
- Component 3: $\rho v$ (momentum in y)

Once you have created the vector field `conserved` containing the conservative variables, you can access each component using:

```cpp
// For cells
conserved[cell][EulerConsVar::rho]
conserved[cell][EulerConsVar::rhoE]
conserved[cell][EulerConsVar::mom(d)] // d = 0, 1 for x, y

// For intervals
conserved(EulerConsVar::rho, level, i, j)
conserved(EulerConsVar::rhoE, level, i, j)
conserved(EulerConsVar::mom(d), level, i, j)
```

### Conversion between conservative and primitive variables

Next, you need to be able to convert the primitive variables from the conservative ones and vice versa. The primitive variables are:
- $\rho$ (density)
- $u = \frac{(\rho u)}{\rho}$ (velocity in x-direction)
- $v = \frac{(\rho v)}{\rho}$ (velocity in y-direction)
- $p = (\gamma - 1) \left( \rho E - \frac{1}{2} \rho (u^2 + v^2) \right)$ (pressure)

where the notation $(\rho u)$ and $(\rho v)$ represent the momentum components.

The functions `cons2prim` and `prim2cons` defined in `variables.hpp` allow you to convert between these two sets of variables.

### Equation of state

The object `eos` defined in `eos.hpp` allows you to compute the pressure and sound speed from the primitive variables. Here is an example of its usage:

```cpp
auto p = eos::stiffened_gas::p(rho, e);
auto c = eos::stiffened_gas::c(rho, p);
auto e = eos::stiffened_gas::e(rho, p);
```

where `e` is the internal energy per unit mass.

### Time step computation

For explicit time integration, the CFL condition must be satisfied:

$$
\Delta t \leq \text{CFL} \cdot \frac{\Delta x}{\max(|u| + c, |v| + c)}
$$

where the maximum is taken over all cells. A default value is set to $\text{CFL} = 0.4$ in the main function. We provide the `get_max_lambda` function in `utils.hpp` that computes the maximum wave speed in the domain:

```cpp
auto get_max_lambda(const auto& u)
{
    static constexpr std::size_t dim = std::decay_t<decltype(u)>::dim;
    double res                       = 0.;

    const auto& mesh = u.mesh();

    samurai::for_each_cell(mesh,
                           [&](const auto& cell)
                           {
                               auto prim = cons2prim<dim>(u[cell]);

                               auto c = EOS::stiffened_gas::c(prim.rho, prim.p);
                               for (std::size_t d = 0; d < dim; ++d)
                               {
                                   res = std::max(std::abs(prim.v[d]) + c, res);
                               }
                           });
    return res;
}
```

## Numerical schemes

We propose implementing three different schemes for the Euler equations: Rusanov, HLL, and HLLC. Each of these schemes is described below.

### Rusanov scheme

The Rusanov scheme is a simple and robust approximate Riemann solver. The flux function for the Rusanov scheme can be defined as follows:

```math
\mathbf{F}_{\text{Rusanov}}(\mathbf{u}_L, \mathbf{u}_R) = \frac{1}{2} \left( \mathbf{F}(\mathbf{u}_L) + \mathbf{F}(\mathbf{u}_R) \right) - \frac{1}{2} s_{\max} (\mathbf{u}_R - \mathbf{u}_L)
```

where $s_{\max}$ is the maximum wave speed in the direction normal to the interface.

:::{note}$s_{\max}$ definition
**For the x-direction flux:**
$$
s_{\max} = \max(|u_L| + c_L, |u_R| + c_R)
$$

**For the y-direction flux:**
$$
s_{\max} = \max(|v_L| + c_L, |v_R| + c_R)
$$
:::

where $c_L = \sqrt{\gamma p_L / \rho_L}$ and $c_R = \sqrt{\gamma p_R / \rho_R}$ are the sound speeds on the left and right states.

### HLL scheme

The HLL (Harten-Lax-van Leer) scheme is another approximate Riemann solver that considers only the fastest left-going and right-going waves. The flux function for the HLL scheme is given by:

```math
\mathbf{F}_{\text{HLL}}(\mathbf{u}_L, \mathbf{u}_R) = \begin{cases}
\mathbf{F}(\mathbf{u}_L) & \text{if } s_L \geq 0 \\\\
\mathbf{F}(\mathbf{u}_R) & \text{if } s_R \leq 0 \\\\
\frac{s_R \mathbf{F}(\mathbf{u}_L) - s_L \mathbf{F}(\mathbf{u}_R) + s_L s_R (\mathbf{u}_R - \mathbf{u}_L)}{s_R - s_L} & \text{otherwise}
\end{cases}
```
where $s_L$ and $s_R$ are the estimated speeds of the left-going and right-going waves, respectively.

```math
s_L = \min(u_L - c_L, u_R - c_R)
```

```math
s_R = \max(u_L + c_L, u_R + c_R)
```

:::{note}
For the y-direction flux, replace $u$ with $v$ in the wave speed estimates:
$$
s_L = \min(v_L - c_L, v_R - c_R), \quad s_R = \max(v_L + c_L, v_R + c_R)
$$
:::

### HLLC scheme
The HLLC (Harten-Lax-van Leer-Contact) scheme is an extension of the HLL scheme that also captures the contact discontinuity. The flux function for the HLLC scheme is more complex and involves additional wave speed estimates.

#### Contact wave speed

The contact wave speed $s_M$ (also called star region velocity) is computed as:

```math
s_M = \frac{p_R - p_L + \rho_L u_L (s_L - u_L) - \rho_R u_R (s_R - u_R)}{\rho_L (s_L - u_L) - \rho_R (s_R - u_R)}
```

#### Star region states

Let's introduce the intermediate states $\mathbf{u}_L^*$ and $\mathbf{u}_R^*$ in the star region. We first compute the density in the star region:

```math
\rho_K^* = \rho_K \frac{s_K - u_K}{s_K - s_M}
```
where $K \in \{L, R\}$.

The velocity components in the star region remain the same as in the original state, except for the normal velocity which is set to $s_M$.

The fourth component (total energy per unit volume) in the star state is:

$$
(\rho E)_K^* = \rho_K^* \left[ E_K + (s_M - u_K)\left(s_M + \frac{p_K}{\rho_K(s_K - u_K)}\right) \right]
$$

where $E_K = e_K + \frac{1}{2}(u_K^2 + v_K^2)$ is the specific total energy.

Thus, for the x-direction flux, the star region states are given by:

```math
\mathbf{u}_K^* = \rho_K^* \begin{pmatrix}
1 \\\\
s_M \\\\ v_K \\\\
e_K + \frac{1}{2}(u_K^2 + v_K^2) + (s_M - u_K)\left(s_M + \frac{p_K}{\rho_K(s_K - u_K)}\right)
\end{pmatrix}
```

and for the y-direction flux, they are given by:

```math
\mathbf{u}_K^* = \rho_K^* \begin{pmatrix}
1 \\\\ u_K \\\\
s_M \\\\
e_K + \frac{1}{2}(u_K^2 + v_K^2) + (s_M - v_K)\left(s_M + \frac{p_K}{\rho_K(s_K - v_K)}\right)
\end{pmatrix}
```

#### Star region fluxes

The fluxes in the star region are computed as:

```math
\mathbf{F}_K^* = \mathbf{F}(\mathbf{u}_K) + s_K (\mathbf{u}_K^* - \mathbf{u}_K)
```

#### Final HLLC flux

The HLLC flux is then given by:

```math
\mathbf{F}_{\text{HLLC}}(\mathbf{u}_L, \mathbf{u}_R) = \begin{cases}
\mathbf{F}(\mathbf{u}_L) & \text{if } s_L \geq 0 \\
\mathbf{F}_L^* & \text{if } s_L < 0 \leq s_M \\
\mathbf{F}_R^* & \text{if } s_M < 0 < s_R \\
\mathbf{F}(\mathbf{u}_R) & \text{if } s_R \leq 0
\end{cases}
```

This four-wave structure allows the HLLC scheme to resolve contact discontinuities more accurately than the HLL scheme while maintaining robustness.

## Initial conditions

The physical domain is $[0,1] \times [0,1]$. The initial condition consists of a Riemann problem with four constant states separated by discontinuities:

$$
(\rho, u, v, p)(x,y,0) = \begin{cases}
(1.5, 0, 0, 1.5) & \text{if } x < 0.5 \text{ and } y < 0.5 \\\\
(0.5323, 1.206, 0, 0.3) & \text{if } x \geq 0.5 \text{ and } y < 0.5 \\\\
(0.5323, 0, 1.206, 0.3) & \text{if } x < 0.5 \text{ and } y \geq 0.5 \\\\
(0.138, 1.206, 1.206, 0.029) & \text{if } x \geq 0.5 \text{ and } y \geq 0.5
\end{cases}
$$

The boundary conditions are homogeneous Neumann conditions.

![Riemann Config 3 Setup](figures/riemann2d_config3_setup.svg)

The initial condition is provided in the `initial_condition.hpp` file.

## Exercises

To implement and simulate the Euler equations in 2D using samurai, you will need to follow the steps introduced in the previous sections. Here is a summary:

- Set the dimension of your problem
- Create a mesh
- Create a vector field with 4 components (density, x-momentum, y-momentum, energy)
- Define the initial condition
- Define the flux functions for the three schemes presented (Rusanov, HLL, HLLC)
- Create the time loop and apply the scheme at each time step
- Visualize the results using ParaView

If everything is working correctly, you can add adaptive mesh refinement using the multiresolution capabilities of samurai.

:::{note}
Go step by step and test each part of your code before moving to the next one. Start with a uniform mesh at level `8` for example and a simple scheme (like Rusanov) before implementing more complex schemes and adaptive refinement.
:::

```{exercise}
Implement the 2D Riemann problem (Configuration 3) with the following steps:
1. Create a uniform mesh at level 8 on the domain $[0, 1] \times [0, 1]$
2. Implement the initial condition with the four states shown in the figure
3. Implement the Rusanov scheme for the flux computation
4. Run the simulation until $t = 0.3$ with CFL = 0.4
5. Visualize the density field in ParaView

Once this is working, implement the HLL and HLLC schemes and compare the results.
```

## Conclusion

In this first part on the Euler equations, you have successfully implemented a complete solver for compressible gas dynamics. You learned to:

- Understand the **coupled system** of conservation laws (mass, momentum, energy)
- Work with **primitive variables** ($\rho$, $p$, $\mathbf{v}$) and **conservative variables** ($\rho$, $\rho E$, $\rho \mathbf{v}$)
- Implement three **Riemann solvers** with increasing sophistication:
  - **Rusanov**: Simple and robust, but dissipative
  - **HLL**: Better resolution of contact discontinuities
  - **HLLC**: Accurately captures all wave structures (shocks, contacts, rarefactions)
- Apply **adaptive mesh refinement** to efficiently resolve shocks and flow features
- Visualize complex multi-dimensional gas dynamics phenomena

In **Part 2**, you will tackle an even more challenging problem: the **double Mach reflection**. This classical benchmark will require you to implement custom boundary conditions (reflecting walls, inflow/outflow) and develop a robust prediction operator to maintain physical positivity during mesh adaptation. This next step will demonstrate how samurai's flexible framework can handle complex, realistic shock physics problems.
