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

// Paramètres du problème de Sedov
double rho_ambient = 1.0;      // Densité ambiante
double p_ambient   = 1e-5;     // Pression ambiante (très petite)
double E_blast     = 0.244816; // Énergie de l'explosion
double r_blast     = 0.1;      // Rayon de la zone d'explosion

void init(auto& u, const xt::xtensor_fixed<double, xt::xshape<2>>& center)
{
    auto& mesh = u.mesh();

    u.resize();

    // Calculer le volume de la sphère d'explosion
    double V_blast = std::numbers::pi * r_blast * r_blast; // 2D: aire du disque

    samurai::for_each_cell(mesh,
                           [&](auto& cell)
                           {
                               auto coords = cell.center();
                               double dx   = coords[0] - center[0];
                               double dy   = coords[1] - center[1];
                               double r    = std::sqrt(dx * dx + dy * dy);

                               double rho = rho_ambient;
                               double p;
                               double vx = 0.;
                               double vy = 0.;

                               if (r < r_blast)
                               {
                                   // Zone d'explosion: énergie concentrée
                                   p = (EOS::gamma - 1.0) * E_blast / V_blast;
                               }
                               else
                               {
                                   // Zone ambiante
                                   p = p_ambient;
                               }

                               // Variables conservatives
                               u[cell][EulerVariable::rho]      = rho;
                               u[cell][EulerVariable::rhou]     = rho * vx;
                               u[cell][EulerVariable::rhou + 1] = rho * vy;
                               u[cell][EulerVariable::rhoE]     = rho * (EOS::e(rho, p) + 0.5 * (vx * vx + vy * vy));
                           });
}

int main(int argc, char* argv[])
{
    constexpr std::size_t dim = 2;
    using Config              = samurai::MRConfig<dim>;

    auto& app = samurai::initialize("Sedov blast wave", argc, argv);

    // Simulation parameters
    xt::xtensor_fixed<double, xt::xshape<dim>> min_corner = {-1., -1.};
    xt::xtensor_fixed<double, xt::xshape<dim>> max_corner = {1., 1.};
    xt::xtensor_fixed<double, xt::xshape<dim>> center     = {0., 0.};

    // Multiresolution parameters
    std::size_t min_level = 7;
    std::size_t max_level = 7;

    double Tf  = 1.0;
    double cfl = 0.45;
    double t   = 0.;
    std::string restart_file;
    std::string scheme = "hllc";

    // Output parameters
    fs::path path        = fs::current_path();
    std::string filename = "sedov_blast";
    std::size_t nfiles   = 20;

    app.add_option("--min-corner", min_corner, "The min corner of the box")->capture_default_str()->group("Simulation parameters");
    app.add_option("--max-corner", max_corner, "The max corner of the box")->capture_default_str()->group("Simulation parameters");
    app.add_option("--center", center, "Center of the blast")->capture_default_str()->group("Simulation parameters");
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

    // Sedov-specific parameters
    app.add_option("--rho-ambient", rho_ambient, "Ambient density")->capture_default_str()->group("Sedov parameters");
    app.add_option("--p-ambient", p_ambient, "Ambient pressure")->capture_default_str()->group("Sedov parameters");
    app.add_option("--E-blast", E_blast, "Blast energy")->capture_default_str()->group("Sedov parameters");
    app.add_option("--r-blast", r_blast, "Blast radius")->capture_default_str()->group("Sedov parameters");

    SAMURAI_PARSE(argc, argv);

    // Initialize the mesh
    const samurai::Box<double, dim> box(min_corner, max_corner);

    samurai::MRMesh<Config> mesh;
    auto u = samurai::make_vector_field<double, 2 + dim>("euler", mesh);

    if (restart_file.empty())
    {
        mesh = {box, min_level, max_level};
        init(u, center);
    }
    else
    {
        samurai::load(restart_file, mesh, u);
    }

    samurai::make_bc<samurai::Neumann<1>>(u, 0., 0., 0., 0.);

    auto unp1 = samurai::make_vector_field<double, 2 + dim>("euler", mesh);

    double dx            = mesh.cell_length(max_level);
    const double dt_save = Tf / static_cast<double>(nfiles);
    std::size_t nsave    = 1;
    std::size_t nt       = 0;

    samurai::save(fmt::format("{}_{}_init", filename, scheme), mesh, u);

    std::cout << "Using scheme: " << scheme << std::endl;
    std::cout << "Blast parameters: E = " << E_blast << ", r = " << r_blast << std::endl;
    std::cout << "Ambient: rho = " << rho_ambient << ", p = " << p_ambient << std::endl;

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
        std::cout << fmt::format("iteration {}: t = {:.6f}, dt = {:.6e}", nt++, t, dt) << std::endl;

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