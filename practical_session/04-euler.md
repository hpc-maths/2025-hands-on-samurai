# Euler equations

It's time to stand on your own two feet and implement the Euler equations using the flux mechanism we have just seen for the Burgers equation. The Euler equations in 2D read:

$$
\frac{\partial }{\partial t} \begin{pmatrix}
\rho \\\\ \rho u \\\\ \rho v \\\\ E
\end{pmatrix}
+ \frac{\partial }{\partial x} \begin{pmatrix}
\rho u \\\\ \rho u^2 + p \\\\ \rho uv \\\\ u(E + p)
\end{pmatrix}
+ \frac{\partial }{\partial y} \begin{pmatrix}
\rho v \\\\ \rho uv \\\\ \rho v^2 + p \\\\ v(E + p)
\end{pmatrix}
= 0
$$

where $\rho$ is the density, $u$ and $v$ are the velocity components in the x and y directions respectively, $E$ is the total energy per unit volume, and $p$ is the pressure given by the equation of state for an ideal gas:

$$
p = (\gamma - 1) \left( E - \frac{1}{2 } \rho (u^2 + v^2) \right)
$$

$\gamma$ is the ratio of specific heats (typically 1.4 for air).

We also need to introduce the sound speed (or speed of sound) is defined as:

```math
c = \sqrt{\frac{\gamma p}{\rho}}
```

It represents the speed at which pressure waves propagate through the fluid.

We propose to implement three different schemes for the Euler equations: Rusanov, HLL and HLLC. We describe in the following each of these schemes.

## Rusanov scheme

The Rusanov scheme is a simple and robust approximate Riemann solver. The flux function for the Rusanov scheme can be defined as follows:

```math
\mathbf{F}_{\text{Rusanov}}(\mathbf{u}_L, \mathbf{u}_R) = \frac{1}{2} \left( \mathbf{F}(\mathbf{u}_L) + \mathbf{F}(\mathbf{u}_R) \right) - \frac{1}{2} s_{\max} (\mathbf{u}_R - \mathbf{u}_L)
```

where $s_{\max}$ is the maximum wave speed, which can be estimated as:

```math
s_{\max} = \max(|u_L| + c_L, |u_R| + c_R)
```

## HLL scheme

The HLL (Harten-Lax-van Leer) scheme is another approximate Riemann solver that considers only the fastest left-going and right-going waves. The flux function for the HLL scheme is given by:

```math
\mathbf{F}_{\text{HLL}}(\mathbf{u}_L, \mathbf{u}_R) = \begin{cases}
\mathbf{F}(\mathbf{u}_L) & \text{if } s_L \geq 0 \\\\
\mathbf{F}(\mathbf{u}_R) & \text{if } s_R \leq 0 \\\\
\frac{s_R \mathbf{F}(\mathbf{u}_L) - s_L \mathbf{F}(\mathbf{u}_R) + s_L s_R (\mathbf{u}_R - \mathbf{u}_L)}{s_R - s_L} & \text{otherwise}
\end{cases}
```
where $s_L$ and $s_R$ are the estimated speeds of the left-going and right-going waves, respectively.

## HLLC scheme
The HLLC (Harten-Lax-van Leer-Contact) scheme is an extension of the HLL scheme that also captures the contact discontinuity. The flux function for the HLLC scheme is more complex and involves additional wave speed estimates.

### Wave speed estimates

First, we estimate the left and right wave speeds using the direct approach:

```math
s_L = \min(u_L - c_L, u_R - c_R)
```

```math
s_R = \max(u_L + c_L, u_R + c_R)
```

where $c_L = \sqrt{\gamma p_L / \rho_L}$ and $c_R = \sqrt{\gamma p_R / \rho_R}$ are the sound speeds on the left and right states.

### Contact wave speed

The contact wave speed $s_M$ (also called star region velocity) is computed as:

```math
s_M = \frac{p_R - p_L + \rho_L u_L (s_L - u_L) - \rho_R u_R (s_R - u_R)}{\rho_L (s_L - u_L) - \rho_R (s_R - u_R)}
```

### Star region states

The intermediate states $\mathbf{u}_L^*$ and $\mathbf{u}_R^*$ in the star region are computed as:

```math
\mathbf{u}_K^* = \rho_K \frac{s_K - u_K}{s_K - s_M} \begin{pmatrix}
1 \\
s_M \\
v_K \\
\frac{E_K}{\rho_K} + (s_M - u_K)\left(s_M + \frac{p_K}{\rho_K(s_K - u_K)}\right)
\end{pmatrix}
```

where $K \in \{L, R\}$.

### Star region fluxes

The fluxes in the star region are computed as:

```math
\mathbf{F}_K^* = \mathbf{F}(\mathbf{u}_K) + s_K (\mathbf{u}_K^* - \mathbf{u}_K)
```

### Final HLLC flux

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

## Exercise

In order to implement and simulate the Euler equations in 2D using samurai, you will need to follow the steps introduced in the previous sections. We summarize them here:

- Set the dimension of your problem
- Create a mesh
- Create a vector field with 4 components (density, momentum in x, momentum in y, energy)
- Define the initial condition
- Define the flux functions for the three schemes presented (Rusanov, HLL, HLLC)
- Make the time loop and apply the scheme at each time step
- Visualize the results using ParaView

If everything is working correctly, you can add the adaptive mesh refinement using the multiresolution capabilities of samurai.
