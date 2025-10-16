#include <numbers>

#include <samurai/algorithm/update.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/io/restart.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>

#include "euler/schemes.hpp"
#include "euler/utils.hpp"
#include "euler/variables.hpp"

double rhoL = 1.;
double pL   = 1;

double rhoR = 0.125;
double pR   = 0.1;

void init(auto& u, double theta)
{
    static constexpr std::size_t dim = std::decay_t<decltype(u)>::dim;
    auto& mesh                       = u.mesh();

    u.resize();
    auto set_conserved = [](auto&& u, double rho, double p)
    {
        u[EulerVariable::rho] = rho;
        for (std::size_t d = 0; d < dim; ++d)
        {
            u[EulerVariable::rhou + d] = 0;
        }
        u[EulerVariable::rhoE] = rho * EOS::e(rho, p);
    };

    double x0 = 0.5;
    double y0 = 0.5;

    samurai::for_each_cell(mesh,
                           [&](auto& cell)
                           {
                               auto x = cell.center();

                               const double x_theta = std::tan(theta) * (x[0] - x0);
                               const double y_theta = x[1] - y0;

                               const double p   = x_theta < y_theta ? pL : pR;
                               const double rho = x_theta < y_theta ? rhoL : rhoR;

                               set_conserved(u[cell], rho, p);
                           });
}

int main(int argc, char* argv[])
{
    constexpr std::size_t dim = 2;
    using Config              = samurai::MRConfig<dim>;

    auto& app = samurai::initialize("SOD shock tube 2D", argc, argv);

    // Simulation parameters
    xt::xtensor_fixed<double, xt::xshape<dim>> min_corner = {0., 0.};
    xt::xtensor_fixed<double, xt::xshape<dim>> max_corner = {1., 1.};

    double theta = std::numbers::pi / 4.;

    // Multiresolution parameters
    std::size_t min_level = 8;
    std::size_t max_level = 8;

    double Tf  = .25;
    double cfl = 0.45;
    double t   = 0.;
    std::string restart_file;
    std::string scheme = "hll";

    // Output parameters
    fs::path path        = fs::current_path();
    std::string filename = "sod_shock_2d";
    std::size_t nfiles   = 1;

    app.add_option("--min-corner", min_corner, "The min corner of the box")->capture_default_str()->group("Simulation parameters");
    app.add_option("--max-corner", max_corner, "The max corner of the box")->capture_default_str()->group("Simulation parameters");
    app.add_option("--cfl", cfl, "The CFL")->capture_default_str()->group("Simulation parameters");
    app.add_option("--theta", theta, "SOD angle")->capture_default_str()->group("Simulation parameters");
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
        init(u, theta);
    }
    else
    {
        samurai::load(restart_file, mesh, u);
    }

    samurai::make_bc<samurai::Neumann<1>>(u, 0, 0, 0, 0);

    auto unp1 = samurai::make_vector_field<double, 2 + dim>("euler", mesh);

    double dx            = mesh.cell_length(max_level);
    const double dt_save = Tf / static_cast<double>(nfiles);
    std::size_t nsave    = 1;
    std::size_t nt       = 0;

    samurai::save(fmt::format("{}_{}_init", filename, scheme), mesh, u);

    std::cout << "Using scheme: " << scheme << std::endl;
    auto fv_scheme = get_fv_scheme<decltype(u)>(scheme);

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