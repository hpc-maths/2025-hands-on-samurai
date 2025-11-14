// Copyright 2025 the samurai team
// SPDX-License-Identifier:  BSD-3-Clause

#include <numbers>

#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>

#include "euler/config.hpp"
#include "euler/initial_condition.hpp"
#include "euler/utils.hpp"
#include "euler/variables.hpp"

int main(int argc, char* argv[])
{
    constexpr std::size_t dim = 2;
    using mesh_t              = config<dim>::mesh_t;

    samurai::initialize("Euler example", argc, argv);
    SAMURAI_PARSE(argc, argv);

    // Multiresolution parameters
    std::size_t min_level = 8;
    std::size_t max_level = 8;

    // Initialize the mesh
    auto box = ; // TO IMPLEMENT

    mesh_t mesh; // TO IMPLEMENT
    auto u = ;   // TO IMPLEMENT
    riemann_config_3(u);

    double Tf      = .25;
    double cfl     = 0.4;
    double t       = 0.;
    std::size_t nt = 0;
    double dx      = mesh.cell_length(max_level);

    samurai::save("euler_init", mesh, u);

    while (t != Tf)
    {
        double dt = cfl * dx / get_max_lambda(u);
        t += dt;

        if (t > Tf)
        {
            dt += Tf - t;
            t = Tf;
        }
        std::cout << fmt::format("iteration {}: t = {}, dt = {}", nt++, t, dt) << std::endl;
    }

    samurai::finalize();
    return 0;
}