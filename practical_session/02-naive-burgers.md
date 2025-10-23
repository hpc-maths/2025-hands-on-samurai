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

where $\Delta t$ is the time step size, and $\Delta x$ is the cell size. We use a Forward Euler method here to discretize time for simplicity. However, more complex time schemes can be used, as we will show at the end of this practical session. $F(u_{i+1/2}^n)$ is the numerical flux function, which we will define using the upwind flux. We will assume that the velocity is always positive, leading to the following flux function definition:

$$
F(u_{i+1/2}) = \frac{1}{2} u_i^2
$$

````{exercise}
Take the provided code skeleton available in `material/02-naive-burgers/main.cpp` and complete the implementation of the naive finite volume scheme for the inviscid Burgers equation.

Ensure that you have correctly implemented the flux calculation by visualizing the evolution of the solution over time using the following command:

```bash
python /path/to/read_mesh.py burgers_1d_ --field u --start 0 --end 340 --wait 10
```
````

## Adapt the Mesh

It is time to adapt the mesh according to the solution. We will use the multi-resolution capabilities of samurai to do so. We will not provide too many theoretical details about multi-resolution in this practical session, but you can refer to [Thomas Bellotti's thesis](https://hal.science/tel-04266822v1) for more information.

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
The process is always the same. This is why we like to say: first implement your solver for your favorite equation on a fixed mesh, and you will have multi-resolution adaptation for free!
```

```{exercise}
Implement the multi-resolution adaptation step in your code.
```

You can visualize the adapted mesh and the solution evolution over time using the same command as before. If you add the command-line option `--save-debug-fields` during the execution of your program, samurai will save additional fields that can help you understand how the mesh is adapted over time, such as levels and coordinates.

If you want to visualize these additional fields, you can use the following command:

```bash
python /path/to/read_mesh.py burgers_1d_ --field levels u --start 0 --end 340 --wait 10
```

## Conclusion

This part of the practical session introduced you to implementing a naive finite volume scheme for the inviscid Burgers equation using the samurai library. You learned how to set up the problem, implement the finite volume update, and adapt the mesh using samurai's multi-resolution capabilities. However, there are several issues with this naive implementation. One of the most significant is that we do not compute the flux at the interfaces between different levels correctly. We will address this issue in the next part of the practical session by implementing a more sophisticated flux calculation that accounts for the multi-resolution mesh structure.