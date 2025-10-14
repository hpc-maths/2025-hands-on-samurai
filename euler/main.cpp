#include <samurai/algorithm/update.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/io/restart.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>

#include "euler/schemes.hpp"
#include "euler/variables.hpp"

double rhoL = 1.;
double pL   = 0.4;
double vL   = -2.;

double rhoR = 1.;
double pR   = 0.4;
double vR   = 2.;

void init(auto& u)
{
    static constexpr std::size_t dim = std::decay_t<decltype(u)>::dim;
    auto& mesh                       = u.mesh();

    u.resize();
    auto set_conserved = [](auto&& u, double rho, double p, double v)
    {
        u[EulerVariable::rho] = rho;
        double norm2          = 0.;
        for (std::size_t d = 0; d < dim; ++d)
        {
            u[EulerVariable::rhou + d] = rho * v;
            norm2 += v * v;
        }
        u[EulerVariable::rhoE] = rho * (EOS::e(rho, p) + 0.5 * norm2);
    };

    samurai::for_each_cell(mesh,
                           [&](auto& cell)
                           {
                               auto x = cell.center();

                               if (x[0] < 0.5)
                               {
                                   set_conserved(u[cell], rhoL, pL, vL);
                               }
                               else
                               {
                                   set_conserved(u[cell], rhoR, pR, vR);
                               }
                           });
}

auto get_max_lambda(const auto& u)
{
    static constexpr std::size_t dim = std::decay_t<decltype(u)>::dim;
    double res                       = 0.;

    const auto& mesh = u.mesh();

    samurai::for_each_cell(mesh,
                           [&](const auto& cell)
                           {
                               auto [rho, vel, e, p, c] = extract_primitive<dim>(u[cell]);
                               for (std::size_t d = 0; d < dim; ++d)
                               {
                                   res = std::max(std::abs(vel[d]) + c, res);
                               }
                           });
    return res;
}

// int main(int argc, char* argv[])
int main()
{
    constexpr std::size_t dim = 1;
    using Config              = samurai::MRConfig<dim>;

    // Simulation parameters
    xt::xtensor_fixed<double, xt::xshape<dim>> min_corner = {0.};
    xt::xtensor_fixed<double, xt::xshape<dim>> max_corner = {1.};

    // Multiresolution parameters
    std::size_t min_level = 10;
    std::size_t max_level = 10;

    // Initialize the mesh
    const samurai::Box<double, dim> box(min_corner, max_corner);

    samurai::MRMesh<Config> mesh;
    auto u = samurai::make_vector_field<double, 2 + dim>("euler", mesh);

    mesh = {box, min_level, max_level};
    init(u);

    const xt::xtensor_fixed<int, xt::xshape<1>> left  = {-1};
    const xt::xtensor_fixed<int, xt::xshape<1>> right = {1};

    samurai::make_bc<samurai::Dirichlet<1>>(u, rhoL, rhoL * (EOS::e(rhoL, pL) + 0.5 * vL * vL), rhoL * vL)->on(left);
    samurai::make_bc<samurai::Dirichlet<1>>(u, rhoR, rhoR * (EOS::e(rhoR, pR) + 0.5 * vR * vR), rhoR * vR)->on(right);

    auto hll     = make_euler_hll<decltype(u)>();
    auto rusanov = make_euler_rusanov<decltype(u)>();
    auto hllc    = make_euler_hllc<decltype(u)>();

    auto unp1 = samurai::make_vector_field<double, 2 + dim>("euler", mesh);

    double Tf  = .15;
    double cfl = 0.45;
    double dx  = mesh.cell_length(max_level);

    double t       = 0.;
    std::size_t nt = 0;

    samurai::save("euler_hll_init", mesh, u);

    while (t != Tf)
    {
        double dt = cfl * dx / get_max_lambda(u);
        t += dt;

        if (std::isnan(t))
        {
            std::cerr << "Error: Time became NaN, stopping simulation" << std::endl;
            break;
        }

        if (t > Tf)
        {
            dt += Tf - t;
            t = Tf;
        }

        std::cout << fmt::format("iteration {}: t = {}, dt = {}", nt++, t, dt) << std::endl;

        samurai::update_ghost_mr(u);
        unp1 = u - dt * hllc(u);
        std::swap(u.array(), unp1.array());
        samurai::save(fmt::format("euler_hll_{}", nt), mesh, u);
    }

    samurai::finalize();
    return 0;
}