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

// Burgers vectoriel 2D: u_t + u*u_x + v*u_y = ν*Δu
//                       v_t + u*v_x + v*v_y = ν*Δv

template <typename field_t>
auto godunov_flux()
{
    static constexpr std::size_t dim          = field_t::dim;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, stencil_size, field_t, field_t>;

    samurai::FluxDefinition<cfg> godunov_flux;

    samurai::static_for<0, dim>::apply( // for each positive Cartesian direction 'd'
        [&](auto _d)
        {
            static constexpr std::size_t d = _d();

            auto f = [](auto u)
            {
                return u(d) * u;
            };

            godunov_flux[d].cons_flux_function =
                [&](samurai::FluxValue<cfg>& flux, const samurai::StencilData<cfg>&, const samurai::StencilValues<cfg>& u)
            {
                auto u_left  = u[0](d);
                auto u_right = u[1](d);

                // Résolution du problème de Riemann pour Burgers: u_t + (0.5*u^2)_x = 0
                // Flux de Godunov = f(u*) où u* est la solution du problème de Riemann en x/t=0

                // Cas 1: u_left <= u_right (onde de raréfaction ou contact)
                if (u_left <= u_right)
                {
                    if (u_left >= 0.0)
                    {
                        flux = f(u[0]); // Tout va vers la droite
                    }
                    else if (u_right <= 0.0)
                    {
                        flux = f(u[1]); // Tout va vers la gauche
                    }
                    // else
                    // {
                    //     flux = 0.0 * u[0]; // 0 est dans la raréfaction, f(0) = 0
                    // }
                }
                else // Cas 2: u_left > u_right (choc)
                {
                    // Le choc se propage à la vitesse s = (u_left + u_right)
                    auto shock_speed = 0.5 * (u_left + u_right);
                    if (shock_speed >= 0.0)
                    {
                        flux = f(u[0]); // Le choc va vers la droite, on prend u_left
                    }
                    else
                    {
                        flux = f(u[1]); // Le choc va vers la gauche, on prend u_right
                    }
                }
            };
        });
    auto scheme = make_flux_based_scheme(godunov_flux);
    scheme.set_name("godunov_flux");
    return scheme;
}

template <typename field_t>
auto burgers_vector_convective_flux()
{
    static constexpr std::size_t dim          = field_t::dim;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, stencil_size, field_t, field_t>;

    samurai::FluxDefinition<cfg> burgers_flux;

    samurai::static_for<0, dim>::apply(
        [&](auto integral_constant_d)
        {
            static constexpr std::size_t d = integral_constant_d();

            burgers_flux[d].cons_flux_function = [](auto& flux, const auto& data, const auto& coord)
            {
                auto& qL = coord[0];
                auto& qR = coord[1];

                // Flux pour composante u: f = u*u (direction x) ou g = u*v (direction y)
                // Flux pour composante v: f = u*v (direction x) ou g = v*v (direction y)
                if constexpr (d == 0) // Direction x
                {
                    auto f_u_L = qL(0) * qL(0); // u*u à gauche
                    auto f_u_R = qR(0) * qR(0); // u*u à droite
                    auto f_v_L = qL(0) * qL(1); // u*v à gauche
                    auto f_v_R = qR(0) * qR(1); // u*v à droite

                    auto alpha_u = std::max(std::abs(qL(0)), std::abs(qR(0)));
                    auto alpha_v = alpha_u; // Même vitesse caractéristique

                    flux(0) = 0.5 * (f_u_L + f_u_R) - 0.5 * alpha_u * (qR(0) - qL(0));
                    flux(1) = 0.5 * (f_v_L + f_v_R) - 0.5 * alpha_v * (qR(1) - qL(1));
                }
                else // Direction y
                {
                    auto g_u_L = qL(1) * qL(0); // v*u à gauche
                    auto g_u_R = qR(1) * qR(0); // v*u à droite
                    auto g_v_L = qL(1) * qL(1); // v*v à gauche
                    auto g_v_R = qR(1) * qR(1); // v*v à droite

                    auto alpha_u = std::max(std::abs(qL(1)), std::abs(qR(1)));
                    auto alpha_v = alpha_u;

                    flux(0) = 0.5 * (g_u_L + g_u_R) - 0.5 * alpha_u * (qR(0) - qL(0));
                    flux(1) = 0.5 * (g_v_L + g_v_R) - 0.5 * alpha_v * (qR(1) - qL(1));
                }
            };
        });

    auto scheme = make_flux_based_scheme(burgers_flux);
    scheme.set_name("burgers_vector_convective");
    return scheme;
}

template <typename field_t>
auto diffusion_flux(double nu)
{
    static constexpr std::size_t dim          = field_t::dim;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::LinearHomogeneous, stencil_size, field_t, field_t>;

    samurai::FluxDefinition<cfg> diff_flux;

    samurai::static_for<0, dim>::apply(
        [&](auto integral_constant_d)
        {
            static constexpr std::size_t d = integral_constant_d();

            diff_flux[d].cons_flux_function = [nu](samurai::FluxStencilCoeffs<cfg>& coeffs, double h)
            {
                // Coefficients pour le schéma centré du second ordre
                coeffs[0] = -1. / h;
                coeffs[1] = 1. / h;
            };
        });

    auto scheme = make_flux_based_scheme(diff_flux);
    scheme.set_name("diffusion_vector");
    return scheme;
}

int main(int argc, char* argv[])
{
    static constexpr std::size_t dim = 2;
    using config                     = samurai::MRConfig<dim>;
    using mesh_t                     = samurai::MRMesh<config>;

    samurai::initialize("samurai viscous vector Burgers 2d", argc, argv);
    SAMURAI_PARSE(argc, argv);

    // Box initiale: Tourbillon de Taylor-Green vectoriel
    auto pi = std::numbers::pi;
    // samurai::Box<double, dim> box({-2 * pi, -2 * pi}, {2 * pi, 2 * pi});
    // Box initiale: Dambreak 2D
    // samurai::Box<double, dim> box({-2.0, -2.0}, {2.0, 2.0});
    samurai::Box<double, dim> box({-1.0, -1.0}, {1.0, 1.0});

    // Double choc avec dissipation (N-wave)
    // samurai::Box<double, dim> box({-3.0, -3.0}, {3.0, 3.0});

    std::size_t min_level = 9;
    std::size_t max_level = 9;
    mesh_t mesh{
        box,
        min_level,
        max_level,
        {true, true}
    };

    auto u = samurai::make_vector_field<double, 2>("u", mesh); // (u, v)

    samurai::for_each_cell(mesh,
                           [&](auto& cell)
                           {
                               const double max = 1;
                               const double r   = 0.5;

                               double dist = 0;
                               for (std::size_t d = 0; d < dim; ++d)
                               {
                                   dist += std::pow(cell.center(d), 2);
                               }
                               dist = std::sqrt(dist);

                               double value = (dist <= r) ? (-max / r * dist + max) : 0;
                               u[cell]      = value;
                           });

    // Condition initiale: Tourbillon de Taylor-Green vectoriel
    // samurai::for_each_cell(mesh,
    //                        [&](const auto& cell)
    //                        {
    //                            auto x = cell.center(0);
    //                            auto y = cell.center(1);

    //                            double A = 1.0;
    //                            // Champ de vitesse qui tourne !
    //                            u[cell][0] = A * std::sin(x) * std::cos(y);  // u
    //                            u[cell][1] = -A * std::cos(x) * std::sin(y); // v
    //                        });

    double Re_target = 100.0;
    double U         = 1.0; // Vitesse caractéristique
    double L         = pi;  // Longueur caractéristique
    // double nu        = U * L / Re_target;

    // std::cout << "Reynolds number: Re = " << Re_target << std::endl;
    // std::cout << "Viscosity: nu = " << nu << std::endl;

    // // Dambreak 2D (rupture de barrage)
    // samurai::for_each_cell(mesh,
    //                        [&](const auto& cell)
    //                        {
    //                            auto x = cell.center(0);
    //                            auto y = cell.center(1);

    //                            // Carré de fluide qui se répand
    //                            u[cell] = (std::abs(x) < 0.5 && std::abs(y) < 0.5) ? 2.0 : 0.0;
    //                        });

    // double nu = 0.01; // Propagation avec coins

    // samurai::for_each_cell(mesh,
    //                        [&](const auto& cell)
    //                        {
    //                            auto x = cell.center(0);
    //                            auto y = cell.center(1);

    //                            // Positions des 4 gaussiennes
    //                            std::array<std::pair<double, double>, 4> centers = {
    //                                {{1.0, 1.0}, {-1.0, 1.0}, {1.0, -1.0}, {-1.0, -1.0}}
    //                            };

    //                            double total_amplitude = 0.0;

    //                            for (const auto& [cx, cy] : centers)
    //                            {
    //                                total_amplitude += std::exp(-10 * ((x - cx) * (x - cx) + (y - cy) * (y - cy)));
    //                            }

    //                            // Vitesse dirigée vers le centre (0,0)
    //                            double speed = 2.0;
    //                            double r     = std::sqrt(x * x + y * y) + 1e-10;

    //                            u[cell][0] = -speed * total_amplitude * x / r; // u vers le centre
    //                            u[cell][1] = -speed * total_amplitude * y / r; // v vers le centre
    //                        });

    double nu = 1e-12; // Interaction + fusion

    double cfl = 0.25;
    double dx  = mesh.cell_length(max_level);

    double dt_conv = dx;
    double dt_diff = cfl * (dx * dx) / nu;
    double dt      = 0.05 * std::min(dt_conv, dt_diff);

    std::cout << "dx = " << dx << std::endl;
    std::cout << "dt_conv = " << dt_conv << ", dt_diff = " << dt_diff << ", dt = " << dt << std::endl;

    auto unp1 = u;
    double t  = 0.0;
    double Tf = 5.0;

    std::size_t nt = 0;

    // auto conv = samurai::make_convection_upwind<decltype(u)>();
    auto diff = samurai::make_diffusion_order2<decltype(u)>();
    // auto conv_flux = burgers_vector_convective_flux<decltype(u)>();
    auto conv = godunov_flux<decltype(u)>();
    // auto diff = diffusion_flux<decltype(u)>(nu);

    auto MRadaptation = samurai::make_MRAdapt(u);
    auto mra_config   = samurai::mra_config();

    samurai::save(fmt::format("viscous_burgers_vector_2d_{}", nt), mesh, u);
    nt++;

    while (t < Tf)
    {
        MRadaptation(mra_config);
        unp1.resize();
        // Calculer la vitesse maximale dans le domaine
        double u_max = 0.0;
        samurai::for_each_cell(mesh,
                               [&](const auto& cell)
                               {
                                   // Norme du vecteur vitesse (u, v)
                                   double u_magnitude = std::sqrt(u[cell][0] * u[cell][0] + u[cell][1] * u[cell][1]);
                                   u_max              = std::max(u_max, u_magnitude);
                               });

        // Sécurité si u_max est trop petit
        u_max = std::max(u_max, 1e-10);

        double dx = mesh.cell_length(max_level);

        // CFL conditions
        double cfl_conv = 0.05; // Coefficient CFL pour la convection
        double cfl_diff = 0.05; // Coefficient CFL pour la diffusion

        // Contrainte convective: dt ≤ CFL * dx / u_max
        double dt_conv = cfl_conv * dx / u_max;

        // Contrainte diffusive: dt ≤ CFL * dx² / (2*dim*ν)
        // En 2D, le facteur 2*dim = 4
        double dt_diff = (nu > 1e-10) ? cfl_diff * dx * dx / (4.0 * nu) : 1e10;

        // Prendre le minimum des deux contraintes
        double dt = std::min(dt_conv, dt_diff);

        // Vérifier qu'on ne dépasse pas Tf
        if (t + dt > Tf)
        {
            dt = Tf - t;
        }

        t += dt;

        if (nt % 100 == 0)
        {
            std::cout << "Time step " << nt << ", t = " << t << std::endl;
        }

        // samurai::update_ghost_mr(u);
        // samurai::for_each_interval(
        //     mesh,
        //     [&](std::size_t level, const auto& i, const auto& index)
        //     {
        //         auto j = index[0];

        //         auto u_x    = xt::view(u(0, level, i, j), xt::all(), xt::newaxis());
        //         auto conv_x = xt::eval(
        //             xt::where(u_x > 0, u_x * (u(level, i, j) - u(level, i - 1, j)), u_x * (u(level, i + 1, j) - u(level, i, j))));
        //         auto u_y    = xt::view(u(1, level, i, j), xt::all(), xt::newaxis());
        //         auto conv_y = xt::eval(
        //             xt::where(u_y > 0, u_y * (u(level, i, j) - u(level, i, j - 1)), u_y * (u(level, i, j) - u(level, i, j + 1))));

        //         auto lap = (u(level, i + 1, j) - 2 * u(level, i, j) + u(level, i - 1, j)) / (dx * dx)
        //                  + (u(level, i, j + 1) - 2 * u(level, i, j) + u(level, i, j - 1)) / (dx * dx);

        //         unp1(level, i, j) = u(level, i, j) - dt * (conv_x + conv_y) + nu * dt * lap;
        //     });
        unp1 = u + dt * (-conv(u)); //+ nu * diff(u));
        // u1   = u - dt * conv(u) + nu * dt * diff(u);
        // u2   = 3. / 4 * u + 1. / 4 * (u1 - dt * conv(u1) + nu * dt * diff(u1));
        // unp1 = 1. / 3 * u + 2. / 3 * (u2 - dt * conv(u2) + nu * dt * diff(u2));

        std::swap(u.array(), unp1.array());

        if (nt % 100 == 0)
        {
            samurai::save(fmt::format("viscous_burgers_vector_2d_{}", nt), mesh, u);
        }
        nt++;
    }

    samurai::finalize();
    return 0;
}
