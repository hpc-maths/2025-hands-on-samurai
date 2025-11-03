# Flux mechanism in samurai

:::{note} Main Objectives
- Understand samurai's flux mechanism
- Write the solver for the inviscid Burgers equation using fluxes in ND
:::

In the previous practical session, we implemented a naive finite volume scheme for the inviscid Burgers equation using `for_each_interval`. We explained that this approach can be wrong when dealing with fluxes at the interfaces between different levels in a multi-resolution mesh. To address this issue, we introduced the concept of fluxes and how to handle them correctly using samurai's built-in flux mechanism. You will also see that the multi-dimensional case is handled in the same way.

A long documentation about the flux mechanism in samurai is available at: [finite volume schemes](https://hpc-math-samurai.readthedocs.io/en/latest/reference/finite_volume_schemes.html). We will give you a short glimpse of it here.

samurai provides three types of schemes to handle fluxes in finite volume schemes:

- the homogeneous linear scheme
- the heterogeneous linear scheme
- the nonlinear scheme

Following our example of the inviscid Burgers equation, we will focus on the nonlinear scheme in this practical session.

The definition of a nonlinear flux-based scheme requires defining a samurai::FluxDefinition<cfg>` where the config in our case is

```cpp
using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, stencil_size, field_t, field_t>;
```

We can observe that the scheme type is `NonLinear`. The stencil size defines the number of neighboring cells used to compute the flux at a given cell interface. In our case, we will use a stencil size of 2, which means that the flux at the interface between cells `i` and `i+1` will be computed using the values of `u` in cells `i` and `i+1`. The output field type and input field type are both `field_t`, which is the type of the solution field `u`. It happens that the output field doesn't have the same number of components as the input field as for example for the divergence operator where the input field is a vector and the output field is a scalar.

Now, the config is defined, we need to define the flux function for each Cartesian direction. In the scalar Burgers equation, the flux function is the same in each direction. Only the cells changed depending on the direction. To construct the flux in samurai, you have to build the following object:

```cpp
samurai::FluxDefinition<cfg> burgers_flux(
        [&](samurai::FluxValue<cfg>& flux, const samurai::StencilData<cfg>& data, const samurai::StencilValues<cfg>& u)
        {
            // Compute flux at the interface using values in u
        });

auto scheme = make_flux_based_scheme(burgers_flux);
scheme.set_name("upwind_flux");
```

Let us explain the arguments of the lambda function:

- `flux`: is the output flux at the interface between two cells. You have to fill this value using the values of `u`.
- `data`: contains information about the stencil, such as the indices of the cells in the stencil.
- `u`: contains the values of the solution field `u` in the stencil cells.

A figure illustrating the stencil of size 2 for a 1D flux is given below:

![](figures/cells.png)

`data` will contain an array `data.cells` with two `cell` instances representing `Cell L` and `Cell R` and an attribute `data.cell_length` containing the length of the cell. `u` will contain the values of the field in these cells which are given by `u[0]` for `uL` and `u[1]` for `uR`, respectively.

```{exercise}
Implement the inviscid Burgers equation using the flux mechanism in 1D.
```

```{exercise}
Set the dimension to 2D and see how the flux mechanism handles the multi-dimensional case.
```

So far we have only considered the scalar inviscid Burgers equation. We will now consider the viscous Burgers equation in 2D:

$$
\mathbf{u}_t + (\mathbf{u} \cdot \nabla)\mathbf{u} = \nu \Delta \mathbf{u}
$$

where $\mathbf{u} = (u,v)$ is the velocity vector and $\nu$ is the viscosity coefficient.

We have to modify the flux function to account for the vector nature of the solution. The flux in the x-direction is given by:

$$
\mathbf{F}(\mathbf{u}) = \displaystyle \begin{pmatrix}\frac{u^2}{2} \\\\ \frac{uv}{2}
\end{pmatrix}
$$

and in the y-direction by:

$$
\mathbf{G}(\mathbf{u}) = \begin{pmatrix}\frac{uv}{2} \\\\ \frac{v^2}{2}
\end{pmatrix}
$$

To handle the fact that the fluxes are not the same in each direction, we need to define in `samurai::FluxDefinition<cfg>` a lambda for each direction following this syntax:

```cpp
samurai::FluxDefinition<cfg> burgers_flux;

// Flux in x-direction
burgers_flux[0].cons_flux_function =
        [&](samurai::FluxValue<cfg>& flux, const samurai::StencilData<cfg>& data, const samurai::StencilValues<cfg>& u)
    {
        // Compute flux in x-direction
    };

// Flux in y-direction
burgers_flux[1].cons_flux_function =
        [&](samurai::FluxValue<cfg>& flux, const samurai::StencilData<cfg>& data, const samurai::StencilValues<cfg>& u)
    {
        // Compute flux in y-direction
    };
```

```{exercise}
Implement the viscous Burgers equation in 2D using the flux mechanism in samurai.
```

The diffusion operator can be implemented using the following configuration:

```cpp
using cfg = samurai::FluxConfig<samurai::SchemeType::LinearHomogeneous, stencil_size, field_t, field_t>;
```

The flux definition changes a little bit. Now, you have to define a lambda function with the following signature:

```cpp
samurai::FluxDefinition<cfg> diffusion_flux(
        [&](samurai::FluxStencilCoeffs<cfg>& coeffs, double h)
        {
            // Compute diffusion coefficients
        });
```

This is because samurai provides an explicit and implicit construction of the linear homogeneous and heterogeneous operators. In our case, we will use the explicit construction. The `coeffs` argument contains the coefficients of the diffusion operator at the interface between two cells. You have to fill these coefficients using the cell length `h`.

```{exercise}
Implement the diffusion operator using the flux mechanism in samurai.
````

```{caution} flux definition for diffusion
We recall that the finite volume flux for the diffusion operator is given at the interface by:
$$
F(u) = \nu \frac{u_R - u_L}{h}
$$
```

```{exercise}
Combine the convection and diffusion operators to solve the viscous Burgers equation in 2D using.
```


