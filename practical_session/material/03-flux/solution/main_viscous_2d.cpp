#include <iostream>
#include <numbers>

#include <samurai/bc.hpp>
#include <samurai/box.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>
#include <samurai/schemes/fv.hpp>

template <typename field_t>
auto godunov_flux()
{
    static constexpr std::size_t dim          = field_t::dim;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, stencil_size, field_t, field_t>;

    auto f = [](auto u)
    {
        return 0.5 * u * u;
    };

    samurai::FluxDefinition<cfg> godunov_flux(
        [&](samurai::FluxValue<cfg>& flux, const samurai::StencilData<cfg>&, const samurai::StencilValues<cfg>& u)
        {
            auto u_left  = u[0];
            auto u_right = u[1];

            // Résolution du problème de Riemann pour Burgers: u_t + (0.5*u^2)_x = 0
            // Flux de Godunov = f(u*) où u* est la solution du problème de Riemann en x/t=0

            // Cas 1: u_left <= u_right (onde de raréfaction ou contact)
            if (u_left <= u_right)
            {
                if (u_left >= 0.0)
                {
                    flux = f(u_left); // Tout va vers la droite
                }
                else if (u_right <= 0.0)
                {
                    flux = f(u_right); // Tout va vers la gauche
                }
                else
                {
                    flux = 0.0; // 0 est dans la raréfaction, f(0) = 0
                }
            }
            else // Cas 2: u_left > u_right (choc)
            {
                // Le choc se propage à la vitesse s = (u_left + u_right)
                auto shock_speed = 0.5 * (u_left + u_right);
                if (shock_speed >= 0.0)
                {
                    flux = f(u_left); // Le choc va vers la droite, on prend u_left
                }
                else
                {
                    flux = f(u_right); // Le choc va vers la gauche, on prend u_right
                }
            }
        });

    auto scheme = make_flux_based_scheme(godunov_flux);
    scheme.set_name("godunov_flux");
    return scheme;
}

template <typename field_t>
auto burgers_convective_flux()
{
    static constexpr std::size_t dim          = field_t::dim;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, stencil_size, field_t, field_t>;

    auto f = [](auto u)
    {
        return 0.5 * u * u;
    };

    samurai::FluxDefinition<cfg> burgers_flux(

        [&](samurai::FluxValue<cfg>& flux, const samurai::StencilData<cfg>& data, const samurai::StencilValues<cfg>& u)
        {
            // Flux physique: f(u) = 0.5 * u^2
            auto f_left  = f(u[0]);
            auto f_right = f(u[1]);

            // Vitesse caractéristique maximale: |df/du| = |u|
            auto alpha = std::max(std::abs(u[0]), std::abs(u[1]));

            // Flux de Lax-Friedrichs (flux centré + diffusion numérique)
            // F_{i+1/2} = 0.5 * [f(u_L) + f(u_R)] - 0.5 * alpha * (u_R - u_L)
            flux = 0.5 * (f_left + f_right) - 0.5 * alpha * (u[1] - u[0]);
        });
    auto scheme = make_flux_based_scheme(burgers_flux);
    scheme.set_name("lax_friedrichs_flux");
    return scheme;
}

template <typename field_t>
auto diffusion_flux(double nu)
{
    static constexpr std::size_t dim          = field_t::dim;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::LinearHomogeneous, stencil_size, field_t, field_t>;

    samurai::FluxDefinition<cfg> diff_flux(
        [&](samurai::FluxStencilCoeffs<cfg>& coeffs, double h)
        {
            // Coefficients pour le schéma centré du second ordre
            coeffs[0] = -nu / h;
            coeffs[1] = nu / h;
        });

    auto scheme = make_flux_based_scheme(diff_flux);
    scheme.set_name("diffusion");
    return scheme;
}

int main(int argc, char* argv[])
{
    static constexpr std::size_t dim = 2;
    using config                     = samurai::MRConfig<dim>;
    using mesh_t                     = samurai::MRMesh<config>;

    samurai::initialize("samurai viscous Burgers 2d", argc, argv);
    SAMURAI_PARSE(argc, argv);

    auto pi = std::numbers::pi;
    // Domaine pour Taylor-Green: [-π, π] × [-π, π]
    samurai::Box<double, dim> box({-4 * pi, -4 * pi}, {4 * pi, 4 * pi});

    std::size_t min_level = 8;
    std::size_t max_level = 8;
    mesh_t mesh{
        box,
        min_level,
        max_level,
        {true, true}
    };

    auto u = samurai::make_scalar_field<double>("u", mesh);

    // Condition initiale: Taylor-Green vortex
    samurai::for_each_cell(mesh,
                           [&](const auto& cell)
                           {
                               auto x = cell.center(0);
                               auto y = cell.center(1);

                               // Vortex de Taylor-Green
                               double A = 2.0; // Amplitude
                               u[cell]  = A * std::sin(x) * std::cos(y);
                           });

    // Définir le nombre de Reynolds souhaité
    double Re_target = 100.0;             // Re modéré pour voir les structures
    double A         = 2.0;               // Amplitude
    double L         = 2 * pi;            // Longueur caractéristique
    double nu        = A * L / Re_target; // nu ≈ 0.126

    std::cout << "Reynolds number: Re = " << Re_target << std::endl;
    std::cout << "Viscosity: nu = " << nu << std::endl;

    double cfl = 0.25; // CFL réduit pour la stabilité
    double dx  = mesh.cell_length(max_level);

    // CFL pour Burgers visqueux: dt ≤ min(CFL*dx, 0.5*dx²/ν)
    double dt_conv = cfl * dx;           // Contrainte convective
    double dt_diff = 0.5 * dx * dx / nu; // Contrainte diffusive
    double dt      = std::min(dt_conv, dt_diff);

    std::cout << "dt_conv = " << dt_conv << ", dt_diff = " << dt_diff << ", dt = " << dt << std::endl;

    auto unp1 = u;
    double t  = 0.0;
    double Tf = 5.0; // Temps final

    std::size_t nt = 0;

    // auto conv_flux = burgers_convective_flux<decltype(u)>();
    auto conv_flux = godunov_flux<decltype(u)>();
    auto diff_flux = diffusion_flux<decltype(u)>(nu);

    // Sauvegarder l'état initial
    samurai::save(fmt::format("viscous_burgers_2d_{}", nt++), mesh, u);

    while (t < Tf)
    {
        t += dt;
        std::cout << "Time step " << nt << ", t = " << t << std::endl;

        // Mise à jour: u^{n+1} = u^n - dt*div(F_conv) + dt*div(F_diff)
        unp1 = u - dt * conv_flux(u) + dt * diff_flux(u);

        std::swap(u.array(), unp1.array());

        // Sauvegarder tous les 10 pas de temps
        if (nt % 10 == 0)
        {
            samurai::save(fmt::format("viscous_burgers_2d_{}", nt), mesh, u);
        }
        nt++;
    }

    samurai::finalize();
    return 0;
}