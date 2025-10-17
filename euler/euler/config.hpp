#pragma once

namespace config_2d
{
    static constexpr std::size_t dim = 2;
    using Config                     = samurai::MRConfig<dim>;
    using mesh_t                     = samurai::MRMesh<Config>;
    using field_t                    = samurai::VectorField<mesh_t, double, dim + 2>;
} // namespace config_2d