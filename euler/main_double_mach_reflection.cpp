#include <numbers>

#include <samurai/algorithm/update.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/io/restart.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>

#include "euler/schemes.hpp"
#include "euler/utils.hpp"
#include "euler/variables.hpp"

double rhoL  = 8.;
double pL    = 116.5;
double alpha = std::numbers::pi / 3.;
double vxL   = 8.25 * std::sin(alpha);
double vyL   = -8.25 * std::cos(alpha);

double rhoR = 1.4;
double pR   = 1.;
double vxR  = 0.;
double vyR  = 0.;

double x0 = 2. / 3;

template <class Field>
struct ValueImpl : public samurai::Bc<Field>
{
    INIT_BC(ValueImpl, 2)

    stencil_t get_stencil(constant_stencil_size_t) const override
    {
        return samurai::line_stencil<dim, 0, 2>();
    }

    apply_function_t get_apply_function(constant_stencil_size_t, const direction_t&) const override
    {
        return [](Field& u, const stencil_cells_t& cells, const value_t& dirichlet_value)
        {
            u[cells[1]] = dirichlet_value;
        };
    }
};

struct Value
{
    template <class Field>
    using impl_t = ValueImpl<Field>;
};

void set_bc(auto& u, double& t)
{
    static constexpr std::size_t dim = std::decay_t<decltype(u)>::dim;

    const xt::xtensor_fixed<int, xt::xshape<dim>> bottom = {0, -1};
    samurai::make_bc<Value>(u,
                            [&](const auto&, const auto& cell, const auto&)
                            {
                                if (cell.center(0) < x0)
                                {
                                    return xt::xtensor_fixed<double, xt::xshape<dim + 2>>{
                                        rhoL,
                                        rhoL * (EOS::e(rhoL, pL) + 0.5 * (vxL * vxL + vyL * vyL)),
                                        rhoL * vxL,
                                        rhoL * vyL};
                                }
                                else
                                {
                                    return xt::xtensor_fixed<double, xt::xshape<dim + 2>>{u[cell][EulerVariable::rho],
                                                                                          u[cell][EulerVariable::rhoE],
                                                                                          u[cell][EulerVariable::rhou],
                                                                                          -u[cell][EulerVariable::rhou + 1]};
                                }
                            })
        ->on(bottom);

    const xt::xtensor_fixed<int, xt::xshape<dim>> top = {0, 1};
    samurai::make_bc<Value>(
        u,
        [&](const auto&, const auto& cell, const auto&)
        {
            double x1 = x0 + 10 * t / std::sin(alpha) + 1 / std::tan(alpha);
            if (cell.center(0) < x1)
            {
                return xt::xtensor_fixed<double, xt::xshape<dim + 2>>{rhoL,
                                                                      rhoL * (EOS::e(rhoL, pL) + 0.5 * (vxL * vxL + vyL * vyL)),
                                                                      rhoL * vxL,
                                                                      rhoL * vyL};
            }
            else
            {
                return xt::xtensor_fixed<double, xt::xshape<dim + 2>>{rhoR,
                                                                      rhoR * (EOS::e(rhoR, pR) + 0.5 * (vxR * vxR + vyR * vyR)),
                                                                      rhoR * vxR,
                                                                      rhoR * vyR};
            }
        })
        ->on(top);

    const xt::xtensor_fixed<int, xt::xshape<dim>> right = {1, 0};
    samurai::make_bc<samurai::Neumann<1>>(u, 0, 0, 0, 0)->on(right);

    const xt::xtensor_fixed<int, xt::xshape<dim>> left = {-1, 0};
    samurai::make_bc<Value>(u, rhoL, rhoL * (EOS::e(rhoL, pL) + 0.5 * (vxL * vxL + vyL * vyL)), rhoL * vxL, rhoL * vyL)->on(left);
}

void init(auto& u)
{
    auto& mesh = u.mesh();

    u.resize();
    auto set_conserved = [](auto&& u, double rho, double p, double vx, double vy)
    {
        u[EulerVariable::rho]      = rho;
        u[EulerVariable::rhou]     = rho * vx;
        u[EulerVariable::rhou + 1] = rho * vy;

        u[EulerVariable::rhoE] = rho * (EOS::e(rho, p) + 0.5 * (vx * vx + vy * vy));
    };

    samurai::for_each_cell(mesh,
                           [&](auto& cell)
                           {
                               auto x = cell.center();

                               if (x[0] < x0 + x[1] / std::tan(alpha))
                               {
                                   set_conserved(u[cell], rhoL, pL, vxL, vyL);
                               }
                               else
                               {
                                   set_conserved(u[cell], rhoR, pR, vxR, vyR);
                               }
                           });
}

int main(int argc, char* argv[])
{
    constexpr std::size_t dim = 2;
    using Config              = samurai::MRConfig<dim>;

    auto& app = samurai::initialize("Double mach reflection", argc, argv);

    // Simulation parameters
    xt::xtensor_fixed<double, xt::xshape<dim>> min_corner = {0., 0.};
    xt::xtensor_fixed<double, xt::xshape<dim>> max_corner = {4., 1.};

    // Multiresolution parameters
    std::size_t min_level = 8;
    std::size_t max_level = 8;

    double Tf  = .25;
    double cfl = 0.8;
    double t   = 0.;
    std::string restart_file;
    std::string scheme = "hllc";

    // Output parameters
    fs::path path        = fs::current_path();
    std::string filename = "double_mach_reflection";
    std::size_t nfiles   = 1;

    app.add_option("--min-corner", min_corner, "The min corner of the box")->capture_default_str()->group("Simulation parameters");
    app.add_option("--max-corner", max_corner, "The max corner of the box")->capture_default_str()->group("Simulation parameters");
    app.add_option("--cfl", cfl, "The CFL")->capture_default_str()->group("Simulation parameters");
    app.add_option("--Ti", t, "Initial time")->capture_default_str()->group("Simulation parameters");
    app.add_option("--Tf", Tf, "Final time")->capture_default_str()->group("Simulation parameters");
    app.add_option("--scheme", scheme, "Finite volume scheme")
        ->capture_default_str()
        ->check(CLI::IsMember({"rusanov", "hll", "hllc"}))
        ->group("Simulation parameters");
    app.add_option("--restart-file", restart_file, "Restart file")->capture_default_str()->group("Simulation parameters");
    app.add_option("--min-level", min_level, "Minimum level of the multiresolution")->capture_default_str()->group("Multiresolution");
    app.add_option("--max-level", max_level, "Maximum level of the multiresolution")->capture_default_str()->group("Multiresolution");
    app.add_option("--path", path, "Output path")->capture_default_str()->group("Output");
    app.add_option("--filename", filename, "File name prefix")->capture_default_str()->group("Output");
    app.add_option("--nfiles", nfiles, "Number of output files")->capture_default_str()->group("Output");

    SAMURAI_PARSE(argc, argv);

    // Initialize the mesh
    const samurai::Box<double, dim> box(min_corner, max_corner);

    samurai::MRMesh<Config> mesh;
    auto u = samurai::make_vector_field<double, 2 + dim>("euler", mesh);

    if (restart_file.empty())
    {
        mesh = {box, min_level, max_level};
        init(u);
    }
    else
    {
        samurai::load(restart_file, mesh, u);
    }
    set_bc(u, t);

    auto unp1 = samurai::make_vector_field<double, 2 + dim>("euler", mesh);

    double dx            = mesh.cell_length(max_level);
    const double dt_save = Tf / static_cast<double>(nfiles);
    std::size_t nsave    = 0;
    std::size_t nt       = 0;

    auto MRadaptation = samurai::make_MRAdapt(u);
    auto mra_config   = samurai::mra_config().epsilon(1e-5);
    MRadaptation(mra_config);

    samurai::save(fmt::format("{}_{}_init", filename, scheme), mesh, u);

    std::cout << "Using scheme: " << scheme << std::endl;
    auto fv_scheme = get_fv_scheme<decltype(u)>(scheme);

    while (t != Tf)
    {
        MRadaptation(mra_config);
        samurai::update_ghost_mr(u);

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

        unp1.resize();
        unp1 = u - dt * fv_scheme(u);

        std::swap(u.array(), unp1.array());

        if (t >= static_cast<double>(nsave + 1) * dt_save || t == Tf)
        {
            const std::string suffix = (nfiles != 1) ? fmt::format("_ite_{}", nsave++) : "";
            samurai::save(fmt::format("{}_{}{}", filename, scheme, suffix), mesh, u);
        }
    }

    samurai::finalize();
    return 0;
}