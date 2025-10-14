#pragma once

#include <samurai/schemes/fv.hpp>

#include "../variables.hpp"
#include "flux.hpp"

template <class Field>
auto make_euler_rusanov()
{
    static constexpr std::size_t dim          = Field::dim;
    static constexpr std::size_t field_size   = Field::n_comp;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, field_size, stencil_size, Field>;

    samurai::FluxDefinition<cfg> rusanov;

    samurai::static_for<0, dim>::apply( // for each positive Cartesian direction 'd'
        [&](auto _d)
        {
            static constexpr std::size_t d = _d();

            rusanov[d].cons_flux_function =
                [](samurai::FluxValue<cfg>& flux, const samurai::StencilData<cfg>& /*data*/, const samurai::StencilValues<cfg>& field)
            {
                static constexpr std::size_t left  = 0;
                static constexpr std::size_t right = 1;

                const auto& qL              = field[left];
                auto [rhoL, vL, eL, pL, cL] = extract_primitive<dim>(qL);

                const auto& qR              = field[right];
                auto [rhoR, vR, eR, pR, cR] = extract_primitive<dim>(qR);

                const auto lambda = std::max(std::abs(vL[d]) + cL, std::abs(vR[d]) + cR);

                flux = 0.5 * (compute_flux<dim, d>(rhoL, vL, eL, pL) + compute_flux<dim, d>(rhoR, vR, eR, pR) - lambda * (qR - qL));
            };
        });
    auto scheme = make_flux_based_scheme(rusanov);
    scheme.set_name("rusanov");

    return scheme;
}