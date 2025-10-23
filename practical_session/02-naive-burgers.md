# Naive Implementation of the Burgers Equation

:::{note} Main Objectives
- Improve your skills with samurai loop constructions
- Improve your skills with samurai field manipulations
:::

In this part of the practical session, we will implement a naive finite volume scheme to solve the one-dimensional inviscid Burgers equation using the samurai library. We will see in later parts how to improve this implementation. The Burgers equation is a fundamental partial differential equation that models various physical phenomena, including fluid dynamics and traffic flow. The inviscid Burgers equation is given by:

$$
\frac{\partial u}{\partial t} + \frac{1}{2} \frac{\partial u^2}{\partial x} = 0.
$$

## Finite Volume Scheme

We will use the simplest finite volume scheme to solve the Burgers equation. The domain is discretized into control volumes (cells), and the solution is approximated by its average value within each cell. The update formula for the cell average $u_i$ at time step $n+1$ is given by:

$$
u_i^{n+1} = u_i^n - \frac{\Delta t}{\Delta x} \left( F(u_{i+1/2}^n) - F(u_{i-1/2}^n) \right),
$$

where $\Delta t$ is the time step size, and $\Delta x$ is the cell size. We use a Forward Euler method here to discretize time for simplicity. However, more complex time schemes can be used, as we will show at the end of this practical session. $F(u_{i+1/2}^n)$ is the numerical flux function, which we will define using the upwind flux. We will assume that the velocity is always positive, leading to the following definition of the flux function:

$$
F(u_{i+1/2}) = \frac{1}{2} u_i^2
$$

```{exercise}
Take the provided code skeleton available in `material/02-naive-burgers/main.cpp` and complete the implementation of the naive finite volume scheme for the inviscid Burgers equation.

Ensure that you correctly implement the flux calculation by visualizing the evolution of the solution over time using the following command line:

````bash
python /path/to/read_mesh.py burgers_1d_ --field u --start 0 --end 340 --wait 10
````

```

## Adapt the mesh

It is time to adapt the mesh according the solution. We will use the multi-resolution capabilities of samurai to do so. We will not give you too much theoretical details about multi-resolution in this practical session, but you can refer to the [Thomas Bellotti thesis](https://hal.science/tel-04266822v1) for more information.

```{exercise}