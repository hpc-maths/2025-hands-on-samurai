#include <numbers>

#include <samurai/algorithm/update.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/io/restart.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>

#include "euler/config.hpp"
#include "euler/init/cases.hpp"
#include "euler/schemes.hpp"
#include "euler/utils.hpp"
#include "euler/variables.hpp"

template <class Field>
void init_sol(Field& u, const std::string& test_case_name)
{
    auto& registry  = test_case::TestCaseRegistry<Field>::instance();
    auto& test_case = registry.get(test_case_name);

    auto& mesh = u.mesh();
    u.resize();

    samurai::for_each_cell(mesh,
                           [&](auto& cell)
                           {
                               test_case.init(u, cell);
                           });
}

template <class Field>
void init_bc(Field& u, double& t, const std::string& test_case_name)
{
    auto& registry  = test_case::TestCaseRegistry<Field>::instance();
    auto& test_case = registry.get(test_case_name);
    test_case.bc(u, t);
}

template <class Field>
auto init_box(const std::string& test_case_name)
{
    auto& registry  = test_case::TestCaseRegistry<Field>::instance();
    auto& test_case = registry.get(test_case_name);
    return test_case.box();
}

int main(int argc, char* argv[])
{
    constexpr std::size_t dim = 2;
    using Config              = samurai::MRConfig<dim>;
    using mesh_t              = samurai::MRMesh<Config>;
    using field_t             = samurai::VectorField<mesh_t, double, dim + 2>;

    auto& app = samurai::initialize("Double mach reflection", argc, argv);

    // Multiresolution parameters
    std::size_t min_level = 8;
    std::size_t max_level = 8;

    double Tf  = .25;
    double cfl = 0.9;
    double t   = 0.;
    std::string restart_file;
    std::string scheme    = "hllc";
    std::string test_case = "double_mach_reflection";

    // Output parameters
    fs::path path      = fs::current_path();
    std::size_t nfiles = 1;

    auto available = test_case::TestCaseRegistry<field_t>::instance().available_test_cases();

    app.add_option("--cfl", cfl, "The CFL")->capture_default_str()->group("Simulation parameters");
    app.add_option("--Ti", t, "Initial time")->capture_default_str()->group("Simulation parameters");
    app.add_option("--Tf", Tf, "Final time")->capture_default_str()->group("Simulation parameters");
    app.add_option("--scheme", scheme, "Finite volume scheme")
        ->capture_default_str()
        ->check(CLI::IsMember({"rusanov", "hll", "hllc"}))
        ->group("Simulation parameters");
    app.add_option("--test-case", test_case, "Test case")->capture_default_str()->check(CLI::IsMember(available))->group("Simulation parameters");
    app.add_option("--restart-file", restart_file, "Restart file")->capture_default_str()->group("Simulation parameters");
    app.add_option("--min-level", min_level, "Minimum level of the multiresolution")->capture_default_str()->group("Multiresolution");
    app.add_option("--max-level", max_level, "Maximum level of the multiresolution")->capture_default_str()->group("Multiresolution");
    app.add_option("--path", path, "Output path")->capture_default_str()->group("Output");
    app.add_option("--nfiles", nfiles, "Number of output files")->capture_default_str()->group("Output");

    SAMURAI_PARSE(argc, argv);

    std::string filename = test_case;

    // Initialize the mesh
    auto box = init_box<field_t>(test_case);

    samurai::MRMesh<Config> mesh;
    auto u = samurai::make_vector_field<double, 2 + dim>("euler", mesh);

    if (restart_file.empty())
    {
        mesh = {box, min_level, max_level};
        init_sol(u, test_case);
    }
    else
    {
        samurai::load(restart_file, mesh, u);
    }
    init_bc(u, t, test_case);

    auto unp1 = samurai::make_vector_field<double, 2 + dim>("euler", mesh);

    double dx            = mesh.cell_length(max_level);
    const double dt_save = Tf / static_cast<double>(nfiles);
    std::size_t nsave    = 0;
    std::size_t nt       = 0;

    auto MRadaptation = samurai::make_MRAdapt(u);
    auto mra_config   = samurai::mra_config().epsilon(1e-3).relative_detail(true);
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