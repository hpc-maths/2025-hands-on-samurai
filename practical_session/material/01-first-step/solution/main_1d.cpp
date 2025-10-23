#include <iostream>
#include <samurai/box.hpp>
#include <samurai/field.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>

int main(int argc, char* argv[])
{
    static constexpr std::size_t dim = 1;
    using config                     = samurai::MRConfig<dim>;
    using mesh_t                     = samurai::MRMesh<config>;

    samurai::initialize("samurai first steps", argc, argv);
    SAMURAI_PARSE(argc, argv);

    samurai::Box<double, dim> box({-1.0}, {1.0});
    std::cout << box << std::endl;

    mesh_t mesh{box, 5, 5};
    std::cout << mesh << std::endl;

    auto u = samurai::make_scalar_field<double>("u", mesh);

    samurai::for_each_cell(mesh,
                           [&](const auto& cell)
                           {
                               auto x  = cell.center(0);
                               u[cell] = std::exp(-50 * x * x);
                           });
    std::cout << u << std::endl;

    samurai::save("field_1d", mesh, u);

    samurai::finalize();
    return 0;
}
