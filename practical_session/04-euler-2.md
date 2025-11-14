# Euler equations (part 2)

:::{note} Main Objectives
- Implement custom boundary conditions for complex geometries
- Handle position-dependent and time-dependent boundary conditions
- Create a custom prediction operator to ensure physical positivity
- Solve the challenging double Mach reflection problem
:::

In this second part of the Euler equations practical session, we move from standard test cases to more complex and realistic scenarios. You will implement the double Mach reflection problem, a classical benchmark that tests the ability of numerical schemes to capture complex shock interactions. This problem requires:

You have to make the part 1 before starting this part.

## Adding a new case

We now want to implement a new test case: the double Mach reflection problem. This problem involves a Mach 10 shock wave reflecting off a solid wall at a 60-degree angle, creating complex shock interactions and flow features. The following figure illustrates the problem setup:

![Double Mach Reflection Setup](figures/double_mach_reflection_setup.svg)

This problem was introduced by Woodward and Colella (1984) as a challenging test case for high-resolution shock-capturing schemes.

:::{note} Reference
Woodward, P., & Colella, P. (1984). The numerical simulation of two-dimensional fluid flow with strong shocks. *Journal of Computational Physics*, 54(1), 115-173.
:::

The default boundary conditions implemented in samurai are not sufficient for this case, so you will need to create custom boundary condition classes. Samurai provides Dirichlet, Neumann, and periodic boundary conditions, but for this exercise, you will implement a new boundary condition that imposes specific values on the ghost cells.

Let's use the implementation of the Dirichlet boundary condition in samurai as an example to create your own custom boundary conditions.

```cpp
template <class Field>
struct DirichletImpl : public samurai::Bc<Field>
{
    INIT_BC(DirichletImpl, 2)

    apply_function_t get_apply_function(constant_stencil_size_t, const direction_t&) const override
    {
        return [](Field& u, const stencil_cells_t& cells, const value_t& dirichlet_value)
        {
            //      [0]   [1]
            //    |_____|.....|
            //     cell  ghost

            u[cells[1]] = 2 * dirichlet_value - u[cells[0]];
        };
    }
};

struct Dirichlet
{
    template <class Field>
    using impl_t = DirichletImpl<Field>;
};
```

:::{note}
- The number 2 in `INIT_BC(DirichletImpl, 2)` indicates that this boundary condition uses a stencil of size 2 (one cell inside the domain and one ghost cell).
- We omit the order parameter in the apply function for simplicity. You can find the full implementation of the Dirichlet boundary condition in samurai [here](https://github.com/hpc-maths/samurai/blob/master/include/samurai/bc/dirichlet.hpp).
:::

```{exercise}
Implement a boundary condition that imposes specific values on the ghost cells. Call this new boundary condition `Value`.
```

You now have all the ingredients to implement the double Mach reflection problem. Using samurai, you can select the boundary cells for each side of the domain as follows:

```cpp
const xt::xtensor_fixed<int, xt::xshape<dim>> bottom = {0, -1}; // Define the direction of the bottom boundary

// If you want to impose a value on all the bottom boundary cells
samurai::make_bc<Value>(u, rho, rhoE, momx, momy)
    ->on(bottom);

// If you want to impose a value depending on the position of the cell
samurai::make_bc<Value>(u,
                        [&](const auto& direction, const auto& cell_in, const auto& coord)
                        {
                            // Return an xt::xtensor_fixed<double, xt::xshape<dim + 2>> representing the value to impose (rho, rhoE, momx, momy)
                        })->on(bottom);
```

:::{note}
- `direction` is an array of integers of size `dim` that indicates the direction going out from `cell_in`.
- `cell_in` is of type `samurai::Cell` and gives the characteristics of the cell inside the domain that has a boundary face.
- `coord` is an array of doubles of size `dim` giving the center of the boundary face.
:::

The double Mach reflection problem requires careful implementation of boundary conditions on each side of the domain. Let's analyze what needs to be done for each boundary:

**Bottom boundary (y = 0):**
- For $x < x_0 = 2/3$: The reflecting wall is in contact with the post-shock state. Apply the post-shock state (left_state) on the ghost cells.
- For $x \geq x_0 = 2/3$: This is a reflecting wall. The density and energy should remain the same as in the interior cell, while the normal velocity component (v) should be reversed. You can achieve this by setting:
  ```cpp
  xt::xtensor_fixed<double, xt::xshape<dim + 2>>{
      u[cell][EulerConsVar::rho],
      u[cell][EulerConsVar::rhoE],
      u[cell][EulerConsVar::mom(0)],
      -u[cell][EulerConsVar::mom(1)]  // Note the minus sign for v-momentum
  }
  ```

**Top boundary (y = 1):**
The shock wave moves with time. At time $t$, the shock position along the top boundary is:
$$x_1(t) = x_0 + \frac{10 t}{\sin(60°)} + \frac{1}{\tan(60°)}$$

- For $x < x_1(t)$: The shock has passed; apply the post-shock state (left_state)
- For $x \geq x_1(t)$: The shock hasn't arrived yet; apply the pre-shock state (right_state)

**Left boundary (x = 0):**
This is an inflow boundary where the post-shock state enters the domain. Apply the post-shock state (left_state) on all ghost cells.

**Right boundary (x = 4):**
This is an outflow boundary. Use Neumann boundary conditions (zero gradient):
```cpp
samurai::make_bc<samurai::Neumann<1>>(u, 0, 0, 0, 0)->on(right);
```

:::{tip}
You can use the `cell.center(0)` and `cell.center(1)` methods to access the x and y coordinates of a cell's center, which is useful for implementing position-dependent boundary conditions.
:::

```{exercise}
Implement the boundary conditions required for the double Mach reflection problem.
```

The initial condition for the double Mach reflection problem divides the domain into two regions separated by the initial shock line. The shock starts at position $(x_0, 0) = (2/3, 0)$ and extends at a 60-degree angle.

**Initial shock geometry:**

The shock line is defined by the equation:
$$x < x_0 + \frac{y}{\tan(60°)}$$

- **Left of the shock line** (post-shock region): Apply the post-shock state with:
  - $\rho = 8.0$
  - $p = 116.5$
  - $\mathbf{v} = (8.25 \sin 60°, -8.25 \cos 60°)$

- **Right of the shock line** (pre-shock region): Apply the pre-shock state with:
  - $\rho = 1.4$
  - $p = 1.0$
  - $\mathbf{v} = (0, 0)$

:::{tip}
To implement this condition, you can iterate over each cell and check if the cell center $(x, y)$ satisfies the condition $x < x_0 + y/\tan(\alpha)$ where $\alpha = 60° = \pi/3$. Remember to convert the primitive variables $(\rho, p, \mathbf{v})$ to conservative variables $(\rho, \rho E, \rho u, \rho v)$ using the `prim2cons` function provided in the code.
:::

```{exercise}
Implement the initial condition for the double Mach reflection problem.
```

```{exercise}
Run the double Mach reflection simulation using the HLLC scheme with uniform mesh at level 8. Visualize the density field at different time steps to observe the shock reflections and interactions.
```

:::{tip}
For the double Mach reflection problem, typical simulation parameters are:
- Domain: $[0, 4] \times [0, 1]$
- Final time: $t = 0.2$
- CFL: 0.4
- Output interval: every 0.01 time units

You should observe the formation of a complex triple-point structure where the incident shock, reflected shock, and Mach stem meet.
:::

## Adaptation issues

If you try to run the double Mach reflection problem with adaptive mesh refinement, you may encounter issues where the density or pressure fields become negative in some regions due to how the multi-resolution algorithm works. The detail computation used for adaptation is based on wavelets, where we compute mean values over cells. This can lead to non-physical negative values for density or pressure during the refinement process when using the prediction operator provided by samurai. This phenomenon is only visible when the stencil size of the prediction operator is greater than or equal to 1 (the default is 1). If you use a prediction of order 0, you should not encounter this issue because you simply copy the value of the parent cell into the newly created child cells.

To avoid this issue, you can provide your own prediction operator that ensures positivity of density and pressure during the prediction step. The idea is to compute the new field values in the newly created child cells using the default prediction operator and then check if the density and pressure are positive. If not, you simply copy the field values from the parent cell into the child cells.

The prediction operator can be added during the creation of the multi-resolution object as follows:

```cpp
auto MRadaptation = samurai::make_MRAdapt(prediction_fn, u);
```

where `prediction_fn` is a function with the following definition:

```cpp
template <std::size_t dim, class TInterval>
class Euler_prediction_op : public samurai::field_operator_base<dim, TInterval>
{
  public:

    INIT_OPERATOR(Euler_prediction_op)

    inline void operator()(samurai::Dim<dim>, auto& dest, const auto& src) const
    {

    }
};

auto prediction_fn = [&](auto& new_field, const auto& old_field)
{
    return samurai::make_field_operator_function<Euler_prediction_op>(new_field, old_field);
};
```

:::{exercise}
Implement a custom prediction operator that ensures positivity of density and pressure during the adaptive mesh refinement process. Follow these steps:

**Step 1: Apply the default prediction**

Start by applying the default prediction of order `pred_order` to compute the values in the child cells:

```cpp
samurai::prediction<pred_order, true>(dest, src)(level, i, index);
```

**Step 2: Prepare for child cell inspection**

Create the indices for the child cells. Remember that when you refine a cell, each parent cell is divided into 4 children in 2D:

```cpp
auto i_f     = i << 1;      // Double the index in x-direction
i_f.step     = 2;           // Set step to 2
auto index_f = index << 1;  // Double the index in y-direction
```

**Step 3: Check for negative density**

For each of the 4 children, check if the density is negative. In 2D, the children are indexed by `(0,0)`, `(1,0)`, `(0,1)`, `(1,1)`:

```cpp
const auto mask_rho =
    (dest(EulerConsVar::rho, level + 1, i_f + 0, index_f + 0) < 0.0) ||
    (dest(EulerConsVar::rho, level + 1, i_f + 1, index_f + 0) < 0.0) ||
    (dest(EulerConsVar::rho, level + 1, i_f + 0, index_f + 1) < 0.0) ||
    (dest(EulerConsVar::rho, level + 1, i_f + 1, index_f + 1) < 0.0);
```

**Step 4: Compute pressure for each child**

For each child cell, compute the pressure from the conservative variables. Remember that:
- $e = E - \frac{1}{2}(u^2 + v^2)$ (internal energy from total energy)
- $p = (\gamma - 1) \rho e$ (equation of state)

Create an array to store the 4 pressure values and compute them:

```cpp
std::array<xt::xtensor<double, 1>, 4> pressure;
for (auto& p : pressure)
{
    p = xt::empty<double>({i.size()});
}

// For child (0,0)
auto rho_00 = dest(EulerConsVar::rho, level + 1, i_f + 0, index_f + 0);
auto e_00 = dest(EulerConsVar::rhoE, level + 1, i_f + 0, index_f + 0) / rho_00;
auto u_00 = dest(EulerConsVar::mom(0), level + 1, i_f + 0, index_f + 0) / rho_00;
auto v_00 = dest(EulerConsVar::mom(1), level + 1, i_f + 0, index_f + 0) / rho_00;
e_00 -= 0.5 * (u_00 * u_00 + v_00 * v_00);
pressure[0] = EOS::stiffened_gas::p(rho_00, e_00);

// Repeat for the other 3 children: (1,0), (0,1), (1,1)
// ...
```

**Step 5: Check for negative pressure**

```cpp
const auto mask_p = (pressure[0] < 0.0) || (pressure[1] < 0.0) ||
                    (pressure[2] < 0.0) || (pressure[3] < 0.0);
```

**Step 6: Apply prediction of order 0 on problematic cells**

For cells where either density or pressure is negative, copy the parent cell value to all children:

```cpp
samurai::apply_on_masked(mask_rho || mask_p,
                         [&](auto& ie)
                         {
                             // Copy parent to child (0,0)
                             xt::view(dest(level + 1, i_f + 0, index_f + 0), ie) =
                                 xt::view(src(level, i, index), ie);

                             // Copy parent to child (1,0)
                             xt::view(dest(level + 1, i_f + 1, index_f + 0), ie) =
                                 xt::view(src(level, i, index), ie);

                             // Copy parent to child (0,1)
                             xt::view(dest(level + 1, i_f + 0, index_f + 1), ie) =
                                 xt::view(src(level, i, index), ie);

                             // Copy parent to child (1,1)
                             xt::view(dest(level + 1, i_f + 1, index_f + 1), ie) =
                                 xt::view(src(level, i, index), ie);
                         });
```

This ensures that physical quantities remain positive during the mesh refinement process.
:::

```{exercise}
Run the double Mach reflection simulation using adaptive mesh refinement with your custom prediction operator. Verify that the simulation runs without issues related to negative density or pressure values.
```
