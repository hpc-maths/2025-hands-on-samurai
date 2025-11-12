#include <iostream>
#include <samurai/bc.hpp>
#include <samurai/box.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>

int main(int argc, char* argv[])
{
    static constexpr std::size_t dim = 1;
    using config                     = samurai::MRConfig<dim>;
    using mesh_t                     = samurai::MRMesh<config>;

    samurai::initialize("samurai naive inviscid Burgers 1d", argc, argv);
    SAMURAI_PARSE(argc, argv);

    samurai::Box<double, dim> box({-2.0}, {4.0});

    std::size_t min_level = 2;
    std::size_t max_level = 8;
    mesh_t mesh{box, min_level, max_level};

    auto u = samurai::make_scalar_field<double>("u", mesh);
    samurai::make_bc<samurai::Dirichlet<1>>(u, 0);
    samurai::for_each_cell(mesh,
                           [&](const auto& cell)
                           {
                               auto x  = cell.center(0);
                               u[cell] = std::exp(-20 * x * x);
                           });

    double cfl = 0.5;
    double dx  = mesh.cell_length(max_level);
    double dt  = cfl * dx;
    auto unp1  = u; // u at the next time step
    double t   = 0.0;
    double Tf  = 4;

    std::size_t nt = 0;

    auto MRadaptation = samurai::make_MRAdapt(u);
    auto mra_config   = samurai::mra_config();

    while (t < Tf)
    {
        MRadaptation(mra_config);
        unp1.resize();

        t += dt;
        std::cout << "Time step " << nt << ", t = " << t << std::endl;

        samurai::update_ghost_mr(u);
        samurai::for_each_interval(mesh,
                                   [&](std::size_t level, const auto& i, const auto&)
                                   {
                                       auto dx         = mesh.cell_length(level);
                                       auto flux_left  = 0.5 * u(level, i - 1) * u(level, i - 1);
                                       auto flux_right = 0.5 * u(level, i) * u(level, i);

                                       unp1(level, i) = u(level, i) - (dt / dx) * (flux_right - flux_left);
                                   });

        samurai::swap(u, unp1);
        samurai::save("results", fmt::format("burgers_1d_{}", nt++), mesh, u);
    }
    samurai::finalize();
    return 0;
}
