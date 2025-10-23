# First Steps

:::{note} Main Objectives
- Create a mesh using samurai
- Create and initialize a field on this mesh
:::

In this part of the practical session, we will create our first mesh using the samurai library. Then, we will create a field defined on this mesh and initialize it. We will start with a simple one-dimensional mesh and then extend our approach to two dimensions.

Recall that each step of this practical session builds upon the code from the previous step. If you are struggling but want to proceed to this step, you can find the complete code in the `material` folder under the corresponding step name. You can copy this code into your working directory to continue with the practical session.

## Creating a 1D Mesh

Several mesh types are defined in samurai, as described in [How-to: create a samurai mesh](https://hpc-math-samurai.readthedocs.io/en/latest/howto/mesh.html). You can also create your own mesh, but that is beyond the scope of this practical session.

In the following parts of the practical session, we will use adaptive mesh refinement techniques, particularly multi-resolution. Therefore, we will use a multi-resolution mesh for this first step. The class used to create such a mesh is `samurai::MRMesh`.

Let's start with a simple one-dimensional multi-resolution mesh. The first step is to create a box that defines the computational domain. In one dimension, a box is defined by its left and right boundaries.

````{exercise}
Create a one-dimensional box as described in [How-to: create a samurai domain using boxes](https://hpc-math-samurai.readthedocs.io/en/latest/howto/box.html) with the left boundary at -1 and the right boundary at 1.

Follow these steps:
- Include the header file needed to create a box
- Create the box using the appropriate constructor

To validate your code, print the box using `std::cout`:
```cpp
std::cout << box << std::endl;
```
````

You have created your first box using samurai! You can now use this box to create a one-dimensional multi-resolution mesh.

The construction of the multi-resolution mesh using the class `samurai::MRMesh` starts by defining the mesh configuration: the dimension, the stencil size of the numerical scheme, the number of cells for the graduation of the adapted mesh, and the stencil size of the prediction. Then, you can create the mesh using the box defined earlier and the minimum and maximum levels of the mesh.

````{exercise}
Create a one-dimensional multi-resolution mesh using the box you created in the previous exercise. Use the default configuration for the mesh.

For assistance, refer to the [How-to: create a samurai mesh](https://hpc-math-samurai.readthedocs.io/en/latest/howto/mesh.html).
````

Once again, you can validate your code by printing the mesh details using `std::cout`:

```cpp
std::cout << mesh << std::endl;
```

Take a moment to explore the output and understand the structure of the mesh you have created.

You have successfully created a one-dimensional multi-resolution mesh using samurai! In the next section, we will create a field defined on this mesh and initialize it.

## Creating and Initializing a 1D Scalar Field

In samurai, fields are defined on meshes using the `samurai::ScalarField` or `samurai::VectorField` classes. A field can represent various physical quantities, such as temperature, pressure, or velocity, depending on the problem being solved. It can be either scalar or vectorial. In this exercise, we will create a scalar field defined on the one-dimensional multi-resolution mesh we created earlier.

Two helper functions are available to create fields on meshes: `samurai::make_scalar_field` and `samurai::make_vector_field`. These functions simplify the process of creating fields by automatically handling the necessary configurations. For more details, refer to the [How-to: create a samurai field](https://hpc-math-samurai.readthedocs.io/en/latest/howto/field.html).

````{exercise}
Create a scalar field defined on the one-dimensional multi-resolution mesh you created earlier. Use the helper function `samurai::make_scalar_field` to create the field. The data type of the field should be `double` and its name should be `"u"`.

To validate your code, print the field details using `std::cout`:

```cpp
std::cout << field << std::endl;
```
````

You have successfully created a scalar field on the one-dimensional multi-resolution mesh! Now, let's initialize this field with a specific function.

We now want to initialize the field with a Gaussian function defined as:

$$
f(x) = \exp\left(-50 x^2\right)
$$

To do this, we need to loop over all the cells of the mesh, get the center of each cell, and set the field value at that cell to the Gaussian function value at the cell center.

samurai provides an easy way to loop over all the cells of the mesh using the `for_each_cell` function, as explained in the [How-to: loop over cells in a samurai mesh](https://hpc-math-samurai.readthedocs.io/en/latest/howto/loop.html).

We need to pause here to tell you that we love C++ lambdas! They are extremely useful for writing concise and readable code. If you are not familiar with them, we recommend taking a look at [Lambda expressions in C++](https://learn.microsoft.com/en-us/cpp/cpp/lambda-expressions-in-cpp?view=msvc-170). Don't worry though—we will guide you through the process.

A lambda is an anonymous function that can be defined inline. It can capture variables from the surrounding scope and be passed as an argument to other functions. The syntax of a lambda is as follows:

```cpp
[capture](parameters) -> return_type {
    // function body
}
```

The most important part for us is the capture list, which defines which variables from the surrounding scope are accessible inside the lambda. For example, to capture a variable `x` by reference, we can write:

```cpp
[&x](parameters) -> return_type {
    // function body
}
```

or to capture it by value, we can write:

```cpp
[x](parameters) -> return_type {
    // function body
}
```

You can also capture all variables by reference using `[&]` or by value using `[=]`.

The `for_each_cell` function takes the mesh and a lambda as arguments. The lambda is called for each cell of the mesh. The parameter of this lambda is an instance of `samurai::Cell`, which provides access to the cell's properties and methods. Therefore, if you want to access the field inside the lambda, you need to capture it by reference.

```{note}
Most of the time, we will capture all variables by reference using `[&]` to simplify the code.
```

Suppose we have a field named `field`. We can write:

```cpp
samurai::for_each_cell(mesh, [&](const auto& cell) {
    // Access the field and cell properties
});
```

`samurai::Cell` provides a method `center()` that returns the center of the cell as an `xt::xtensor`. In one dimension, this tensor has a single element, which is the x-coordinate of the cell center. To access this element, you can use the syntax `center()[0]` or `center(0)`. For more details, refer to the [`Cell` class documentation](https://hpc-math-samurai.readthedocs.io/en/latest/api/cell.html).

````{exercise}
Now, let's initialize the field with the Gaussian function. Do this by looping over all the cells of the mesh and setting the field value at each cell to the Gaussian function value at the cell center.
````

### Saving and Plotting the 1D Scalar Field

To visualize a samurai field, you can save it using the built-in samurai functions to export the data in formats compatible with visualization tools like ParaView or matplotlib. Refer to [How-to: save your samurai mesh and fields](https://hpc-math-samurai.readthedocs.io/en/latest/howto/save.html) and [How-to: plot samurai fields and meshes](https://hpc-math-samurai.readthedocs.io/en/latest/howto/plot.html) for detailed instructions on saving your field data and creating plots.

Since we are working with a one-dimensional field, we cannot use ParaView directly. However, we provide a Python script that allows you to visualize the data using matplotlib. This script is located in the samurai source directory, but we also include it in the `material/01-first-step` directory for your convenience. The script is named `read_mesh.py`.

Let's start by saving the mesh and the field to files that can be read by the Python script. You can use the `samurai::save` function to save both the mesh and the field.

````{exercise}
Write the code to save the mesh and the field to a file named `field_1d`.
````
Once you have saved the mesh and the field, you can use the provided Python script to visualize the field. Run the following command in your terminal:

```bash
python practical_session/material/01-first-step/read_mesh.py field_1d --field u
```

```{caution}
Make sure to replace `u` with the name you used when creating the field if it differs. The name is the string you provided when calling `samurai::make_scalar_field`.
```

You have successfully initialized and visualized the scalar field on the one-dimensional multi-resolution mesh with the Gaussian function! You can now proceed to the next part of the practical session, where we will extend our approach to two dimensions.

### A 2D Scalar Field

Using the knowledge you have gained so far, try to create a two-dimensional multi-resolution mesh, create a scalar field on this mesh, and initialize it with a two-dimensional Gaussian function defined as:

```math
u(x, y) = \exp\left(-50 (x^2 + y^2)\right)
```

At the end of this step, open ParaView to visualize the field you have created and initialized. You must open the xdmf file (not the h5 file) to see both the mesh and the field.

### Another Loop Approach

samurai provides an alternative way to loop over cells using its interval-based functionality. `samurai::for_each_cell` uses this functionality under the hood. This loop function is called `samurai::for_each_interval`. Here is how you can use it:

```cpp
samurai::for_each_interval(mesh, [&](std::size_t level, const auto& interval, const auto& index)
{
    // Access the field and interval properties
});
```

Let's explain the lambda parameters:
- `level`: the level of the cells in the current interval
- `interval`: an instance of `samurai::Interval` defining the cells in the x direction
- `index`: a container of size $dim - 1$ containing the indices of the cells in directions other than x

In the previous exercise, we used the accessor `field[cell]` to access the field value at the cell. In this loop, you can use the accessor `field(level, interval, index)` to access the field values at the interval.

```{caution}
You now have an entire interval of values rather than a single scalar value. We use xtensor to provide a NumPy-like array of values and perform lazy evaluation. For more details, refer to the [From NumPy to xtensor documentation](https://xtensor.readthedocs.io/en/latest/numpy.html).
```

```{note}
In the following exercises, we will use `i` instead of `interval` to simplify the notation.
```

````{exercise}
Rewrite the field initialization using the `samurai::for_each_interval` function instead of `samurai::for_each_cell`.

If you are familiar with NumPy, you can use the `xt::arange` function to create an array of cell centers for the interval. Then, compute the Gaussian function on this array and assign it to the field. Refer to [loop over intervals](https://hpc-math-samurai.readthedocs.io/en/latest/howto/loop.html#looping-over-intervals) to have an example.
````

### Conclusion

In this first step of the practical session, you have learned how to create one-dimensional and two-dimensional multi-resolution meshes using samurai. Recall that other mesh types exist in samurai, but they are beyond the scope of this session. You have also created scalar fields on these meshes and initialized them with Gaussian functions. Finally, you visualized the fields using matplotlib and ParaView.

Congratulations on completing this step! You are now ready to move on to the next part of the practical session, where we will explore how to implement a finite volume scheme to solve the Burgers equation.