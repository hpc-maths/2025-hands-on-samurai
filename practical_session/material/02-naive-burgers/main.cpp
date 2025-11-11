#include <iostream>
#include <samurai/bc.hpp>
#include <samurai/box.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>
#include <samurai/update.hpp>

int main(int argc, char* argv[])
{
    static constexpr std::size_t dim = 1;
    using config                     = samurai::MRConfig<dim>;
    using mesh_t                     = samurai::MRMesh<config>;

    samurai::initialize("samurai naive inviscid Burgers", argc, argv);
    SAMURAI_PARSE(argc, argv);

    // Mesh initialization
    samurai::Box<double, dim> box({-2.0}, {4.0});

    std::size_t min_level = 8;
    std::size_t max_level = 8;
    mesh_t mesh{box, min_level, max_level};

    // Field initialization
    auto u = samurai::make_scalar_field<double>("u", mesh);
    samurai::make_bc<samurai::Dirichlet<1>>(u, 0);

    samurai::for_each_cell(mesh,
                           [&](const auto& cell)
                           {
                               auto x  = cell.center(0);
                               u[cell] = std::exp(-20 * x * x);
                           });

    // Time-stepping parameters
    double cfl = 0.5;
    double dx  = mesh.cell_length(max_level);
    double dt  = cfl * dx;
    auto unp1  = u; // u at the next time step
    double t   = 0.0;
    double Tf  = 4;

    // Time-stepping loop
    std::size_t nt = 0;
    while (t < Tf)
    {
        t += dt;
        std::cout << "Time step " << nt << ", t = " << t << std::endl;

        samurai::update_ghost_mr(u);

        samurai::for_each_interval(mesh,
                                   [&](std::size_t level, const auto& i, const auto&)
                                   {
                                       // PART TO IMPLEMENT
                                       // implement upwind fluxes here
                                   });

        samurai::swap(u, unp1);
        samurai::save("results", fmt::format("burgers_{}d_{}", dim, nt++), mesh, u);
    }
    samurai::finalize();
    return 0;
}
