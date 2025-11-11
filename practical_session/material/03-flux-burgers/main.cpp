#include <iostream>
#include <samurai/bc.hpp>
#include <samurai/box.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>
#include <samurai/schemes/fv.hpp>

//
// Flux definition for the inviscid (scalar) Burgers equation
//
template <typename field_t>
auto convective_flux()
{
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, stencil_size, field_t, field_t>;

    samurai::FluxDefinition<cfg> burgers_flux;

    burgers_flux[0].cons_flux_function =
        [&](samurai::FluxValue<cfg>& flux, const samurai::StencilData<cfg>& data, const samurai::StencilValues<cfg>& u)
    {
        // PART TO IMPLEMENT
    };
    auto scheme = make_flux_based_scheme(burgers_flux);
    scheme.set_name("convective_flux");
    return scheme;
}

int main(int argc, char* argv[])
{
    static constexpr std::size_t dim = 1;
    using config                     = samurai::MRConfig<dim>;
    using mesh_t                     = samurai::MRMesh<config>;

    samurai::initialize("samurai naive inviscid Burgers 1d", argc, argv);
    SAMURAI_PARSE(argc, argv);

    // Mesh initialization
    samurai::Box<double, dim> box({-2.0}, {4.0});

    std::size_t min_level = 2;
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
    double cfl    = 0.5;
    double min_dx = mesh.cell_length(max_level);
    double dt     = cfl * min_dx;
    double t      = 0.0;
    double Tf     = 4;

    // Scheme initialization
    auto conv_flux = convective_flux<decltype(u)>();
    auto unp1      = u; // u at the next time step

    // Time-stepping loop
    std::size_t nt = 0;
    while (t < Tf)
    {
        t += dt;
        std::cout << "Time step " << nt << ", t = " << t << std::endl;

        unp1 = un - dt * conv_flux(un);

        samurai::swap(u, unp1);
        samurai::save(fmt::format("burgers_1d_{}", nt++), mesh, u);
    }
    samurai::finalize();
    return 0;
}
