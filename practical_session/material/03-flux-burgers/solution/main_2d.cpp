// Copyright 2025 the samurai team
// SPDX-License-Identifier:  BSD-3-Clause

#include <iostream>
#include <samurai/bc.hpp>
#include <samurai/box.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>
#include <samurai/schemes/fv.hpp>

template <typename field_t>
auto upwind_flux()
{
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, stencil_size, field_t, field_t>;

    samurai::FluxDefinition<cfg> burgers_flux(
        [&](samurai::FluxValue<cfg>& flux, const samurai::StencilData<cfg>&, const samurai::StencilValues<cfg>& u)
        {
            if (u[0] >= 0)
            {
                flux = u[0] * u[0];
            }
            else
            {
                flux = u[1] * u[1];
            }
        });

    auto scheme = make_flux_based_scheme(burgers_flux);
    return scheme;
}

int main(int argc, char* argv[])
{
    static constexpr std::size_t dim = 2;
    using config                     = samurai::MRConfig<dim>;
    using mesh_t                     = samurai::MRMesh<config>;

    samurai::initialize("samurai inviscid Burgers 2d using samurai flux", argc, argv);
    SAMURAI_PARSE(argc, argv);

    samurai::Box<double, dim> box({-2.0, -2.0}, {4.0, 4.0});

    std::size_t min_level = 8;
    std::size_t max_level = 8;
    mesh_t mesh{box, min_level, max_level};

    auto u = samurai::make_scalar_field<double>("u", mesh);
    samurai::make_bc<samurai::Dirichlet<1>>(u, 0);

    samurai::for_each_cell(mesh,
                           [&](const auto& cell)
                           {
                               auto x  = cell.center(0);
                               auto y  = cell.center(1);
                               u[cell] = std::exp(-20 * (x * x + y * y));
                           });

    double cfl     = 0.5;
    double min_dx  = mesh.cell_length(max_level);
    double dt      = cfl * min_dx;
    double t       = 0.0;
    double Tf      = 4;
    std::size_t nt = 0;

    auto unp1 = u; // u at the next time step

    auto MRadaptation = samurai::make_MRAdapt(u);
    auto mra_config   = samurai::mra_config();

    auto conv_flux = upwind_flux<decltype(u)>();

    while (t < Tf)
    {
        MRadaptation(mra_config);
        unp1.resize();

        t += dt;
        std::cout << "Time step " << nt << ", t = " << t << std::endl;

        unp1 = u - dt * conv_flux(u);

        samurai::swap(u, unp1);
        samurai::save("results", fmt::format("burgers_2d_{}", nt++), mesh, u);
    }
    samurai::finalize();
    return 0;
}
