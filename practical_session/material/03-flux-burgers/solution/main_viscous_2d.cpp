// Copyright 2025 the samurai team
// SPDX-License-Identifier:  BSD-3-Clause

// Copyright 2025 the samurai team
// SPDX-License-Identifier:  BSD-3-Clause

#include <iostream>
#include <numbers>

#include <samurai/bc.hpp>
#include <samurai/box.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>
#include <samurai/schemes/fv.hpp>

template <typename field_t>
auto upwind_conservative_flux()
{
    static constexpr std::size_t dim          = field_t::dim;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, stencil_size, field_t, field_t>;

    samurai::FluxDefinition<cfg> conv_flux;

    samurai::static_for<0, dim>::apply( // for each positive Cartesian direction 'd'
        [&](auto _d)
        {
            static constexpr std::size_t d = _d();

            conv_flux[d].flux_function =
                [&](samurai::FluxValuePair<cfg>& flux, const samurai::StencilData<cfg>&, const samurai::StencilValues<cfg>& u)
            {
                auto u_left  = u[0](d);
                auto u_right = u[1](d);

                if (u_left >= 0)
                {
                    flux[0] = u_left * u[0];
                    flux[1] = -flux[0];
                }
                else
                {
                    flux[0] = u_right * u[1];
                    flux[1] = -flux[0];
                }
            };
        });
    auto scheme = make_flux_based_scheme(conv_flux);
    return scheme;
}

template <typename field_t>
auto upwind_non_conservative_flux()
{
    static constexpr std::size_t dim          = field_t::dim;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, stencil_size, field_t, field_t>;

    samurai::FluxDefinition<cfg> conv_flux;

    samurai::static_for<0, dim>::apply( // for each positive Cartesian direction 'd'
        [&](auto _d)
        {
            static constexpr std::size_t d = _d();

            conv_flux[d].flux_function =
                [&](samurai::FluxValuePair<cfg>& flux, const samurai::StencilData<cfg>&, const samurai::StencilValues<cfg>& u)
            {
                auto u_left  = u[0](d);
                auto u_right = u[1](d);

                auto u_mean = 0.5 * (u_left + u_right);

                flux[0] = -u_mean * (std::copysign(1.0, u_mean) - 1.0) * 0.5 * (u[1] - u[0]);
                flux[1] = u_mean * (std::copysign(1.0, u_mean) + 1.0) * 0.5 * (u[1] - u[0]);
            };
        });
    auto scheme = samurai::make_flux_based_scheme(conv_flux);
    return scheme;
}

template <typename field_t>
auto diffusion_flux(double nu)
{
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::LinearHomogeneous, stencil_size, field_t, field_t>;

    samurai::FluxDefinition<cfg> diffusion(
        [nu](samurai::FluxStencilCoeffs<cfg>& coeffs, double h)
        {
            const std::size_t left  = 0;
            const std::size_t right = 1;

            coeffs[left].fill(0);
            coeffs[right].fill(0);
            for (std::size_t i = 0; i < field_t::n_comp; ++i)
            {
                coeffs[left](i, i)  = -1 / h;
                coeffs[right](i, i) = 1 / h;
            }
            coeffs *= nu;
        });
    auto scheme = samurai::make_flux_based_scheme(diffusion);
    return scheme;
}

auto compute_dt(const auto& u, double cfl, double nu)
{
    auto& mesh    = u.mesh();
    double min_dx = mesh.cell_length(mesh.max_level());

    // Compute the maximum velocity in the domain
    double u_max = 0.0;
    samurai::for_each_cell(mesh,
                           [&](const auto& cell)
                           {
                               // Norm of the velocity vector (u, v)
                               double u_magnitude = std::sqrt(u[cell][0] * u[cell][0] + u[cell][1] * u[cell][1]);
                               u_max              = std::max(u_max, u_magnitude);
                           });

    // Safety check if u_max is too small
    u_max = std::max(u_max, 1e-10);

    // CFL conditions
    // Convective constraint: dt ≤ CFL * dx / u_max
    double dt_conv = cfl * min_dx / u_max;

    // Diffusive constraint: dt ≤ CFL * dx² / (2*dim*ν)
    // In 2D, the factor 2*dim = 4
    double dt_diff = (nu > 1e-10) ? cfl * min_dx * min_dx / (4.0 * nu) : 1e10;

    // Take the minimum of the two constraints
    return std::min(dt_conv, dt_diff);
}

int main(int argc, char* argv[])
{
    static constexpr std::size_t dim = 2;
    using config                     = samurai::MRConfig<dim>;
    using mesh_t                     = samurai::MRMesh<config>;

    auto& app = samurai::initialize("samurai viscous vector Burgers 2d", argc, argv);

    auto pi = std::numbers::pi;
    samurai::Box<double, dim> box({0, 0}, {2 * pi, 2 * pi});
    double U0 = 1.0, k = 2.0;

    std::size_t min_level = 8;
    std::size_t max_level = 8;

    double Re_target       = 100.0;
    double L               = pi; // Characteristic length scale
    double nu              = U0 * L / Re_target;
    double cfl             = 0.45;
    double t               = 0.0;
    double Tf              = 5.0;
    std::size_t save_every = 10;
    std::size_t nt         = 0;

    app.add_option("--Re", Re_target, "Target Reynolds number")->capture_default_str()->group("Simulation parameters");
    app.add_option("--U0", U0, "Characteristic velocity")->capture_default_str()->group("Simulation parameters");
    app.add_option("--k", k, "Wavenumber")->capture_default_str()->group("Simulation parameters");
    app.add_option("--cfl", cfl, "CFL number")->capture_default_str()->group("Simulation parameters");
    app.add_option("--Tf", Tf, "Final time")->capture_default_str()->group("Simulation parameters");
    app.add_option("--save-every", save_every, "Save every n time steps")->capture_default_str()->group("Output");
    app.add_option("--min-level", min_level, "Minimum level of the multiresolution")->capture_default_str()->group("Multiresolution");
    app.add_option("--max-level", max_level, "Maximum level of the multiresolution")->capture_default_str()->group("Multiresolution");

    SAMURAI_PARSE(argc, argv);

    std::cout << "Reynolds number: Re = " << Re_target << std::endl;
    std::cout << "Viscosity: nu = " << nu << std::endl;

    mesh_t mesh{
        box,
        min_level,
        max_level,
        {true, true}
    };

    auto u = samurai::make_vector_field<double, 2>("u", mesh); // (u, v)

    samurai::for_each_cell(mesh,
                           [&](const auto& cell)
                           {
                               auto center = cell.center();
                               double x    = center[0];
                               double y    = center[1];

                               u[cell][0] = U0 * std::sin(k * x) * std::cos(k * y);
                               u[cell][1] = -U0 * std::cos(k * x) * std::sin(k * y);
                           });

    auto unp1 = u; // u at the next time step

    auto diff = diffusion_flux<decltype(u)>(nu);
    auto conv = upwind_conservative_flux<decltype(u)>();
    // auto conv = upwind_non_conservative_flux<decltype(u)>();

    auto MRadaptation = samurai::make_MRAdapt(u);
    auto mra_config   = samurai::mra_config();

    samurai::save("results", fmt::format("burgers_viscous_2d_{}", nt), mesh, u);
    nt++;

    while (t < Tf)
    {
        MRadaptation(mra_config);
        unp1.resize();

        double dt = compute_dt(u, cfl, nu);

        // Check that we do not exceed Tf
        if (t + dt > Tf)
        {
            dt = Tf - t;
        }

        t += dt;

        if (nt % save_every == 0)
        {
            std::cout << "Time step " << nt << ", t = " << t << ", dt = " << dt << std::endl;
        }

        unp1 = u + dt * (nu * diff(u) - conv(u));

        samurai::swap(u, unp1);

        if (nt % save_every == 0)
        {
            samurai::save("results", fmt::format("burgers_viscous_2d_{}", nt), mesh, u);
        }
        nt++;
    }

    std::cout << "Simulation finished at t = " << t << " after " << nt << " time steps." << std::endl;
    samurai::finalize();
    return 0;
}
