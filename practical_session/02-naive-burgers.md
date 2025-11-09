# Naive implementation of the Burgers equation

:::{note} Main Objectives
- Develop skills in implementing finite volume schemes with samurai
- Learn to integrate multi-resolution mesh adaptation
- Recognize the challenges of flux calculation on non-uniform meshes
:::

In this part, we will implement a naive finite volume scheme to solve the one-dimensional inviscid Burgers equation. While this implementation is straightforward, we will see in later sections how to improve it to handle multi-resolution meshes correctly. We will use the Burgers equation as a demonstration. This equation is a fundamental partial differential equation that models various physical phenomena, including fluid dynamics and traffic flow. The inviscid Burgers equation is given by:

$$
\frac{\partial u}{\partial t} + \frac{1}{2} \frac{\partial u^2}{\partial x} = 0.
$$

## Finite volume scheme

We will use the simplest finite volume scheme to solve the Burgers equation. The domain is discretized into control volumes (cells), and the solution is approximated by its average value within each cell. The update formula for the cell average $u_i$ at time step $n+1$ is given by:

$$
u_i^{n+1} = u_i^n - \frac{\Delta t}{\Delta x} \left( F(u_{i+1/2}^n) - F(u_{i-1/2}^n) \right),
$$

where $u_i^n$ represents the cell average of $u$ in the cell $i$ at time step $n$, $\Delta t$ is the time step size, and $\Delta x$ is the cell size. We use a Forward Euler method here to discretize time for simplicity. However, more complex time schemes can be used, as we will show at the end of this practical session. $F(u_{i+1/2}^n)$ is the numerical flux function, which we will define using the upwind flux. For simplicity, we will assume that the velocity is always positive (upwind scheme). The flux at the right interface of cell $i$ is determined by the left state:

$$
F(u_{i+1/2}) = \frac{1}{2} u_i^2
$$

````{exercise}
A code skeleton is provided in `material/02-naive-burgers/main.cpp` to help you get started. Complete the implementation of the naive finite volume scheme for the inviscid Burgers equation.

Ensure that you have correctly implemented the flux calculation by visualizing the evolution of the solution over time using the following command:

```bash
python /path/to/read_mesh.py burgers_1d_ --field u --start 0 --end 340 --wait 10
```
````

## Mesh adaptation

It is time to adapt the mesh according to the solution. We will use the multi-resolution capabilities of samurai to do so. We will not go into the theoretical details of multi-resolution here, but you can refer to [Thomas Bellotti's thesis](https://hal.science/tel-04266822v1) for more information.

Adding the adaptation step performed by the multi-resolution framework of samurai is straightforward. You simply need to follow these steps:

- Include the multi-resolution adaptation header:

```cpp
#include <samurai/mr/adapt.hpp>
```

- Before the time loop, define the multi-resolution configuration:

```cpp
auto MRadaptation = samurai::make_MRAdapt(u);
auto mra_config   = samurai::mra_config();
````

- Inside the time loop, before updating the solution, call the adaptation function and resize the solution field:

```cpp
MRadaptation(mra_config);
unp1.resize();
```

```{note}
The process is always the same: first implement your solver on a fixed mesh, then add multi-resolution adaptation with just a few lines of code!
```

```{exercise}
Implement the multi-resolution adaptation step in your code.
```

:::{important}
Don't forget to change the `min_level` and `max_level` to allow for mesh adaptation. For example, you can set `min_level = 2` and `max_level = 8`.
:::

You can visualize the adapted mesh and the solution evolution over time using the same command as before. If you add the command-line option `--save-debug-fields` when running your program, samurai will save additional fields that can help you understand how the mesh is adapted over time, such as levels and coordinates.

If you want to visualize these additional fields, you can use the following command:

```bash
python /path/to/read_mesh.py burgers_1d_ --field levels u --start 0 --end 340 --wait 10
```

## Conclusion

In this part, you implemented a naive finite volume scheme for the inviscid Burgers equation using samurai. You learned how to set up the problem, implement the finite volume update, and adapt the mesh using samurai's multi-resolution capabilities. However, this naive implementation has several issues. Most notably, we do not compute the flux correctly at interfaces between cells of different refinement levels. In the next section, we will address this by implementing a flux calculation that properly accounts for the multi-resolution mesh structure.