#pragma once

#include <samurai/schemes/fv.hpp>

#include "../variables.hpp"
#include "flux.hpp"

template <class Field>
auto make_euler_hll()
{
    static constexpr std::size_t dim          = Field::dim;
    static constexpr std::size_t field_size   = Field::n_comp;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, field_size, stencil_size, Field>;

    samurai::FluxDefinition<cfg> hll;

    samurai::static_for<0, dim>::apply( // for each positive Cartesian direction 'd'
        [&](auto _d)
        {
            static constexpr std::size_t d = _d();

            hll[d].cons_flux_function =
                [](samurai::FluxValue<cfg>& flux, const samurai::StencilData<cfg>& /*data*/, const samurai::StencilValues<cfg>& field)
            {
                static constexpr std::size_t left  = 0;
                static constexpr std::size_t right = 1;

                const auto& qL = field[left];
                auto primL     = cons2prim<dim>(qL);

                const auto& qR = field[right];
                auto primR     = cons2prim<dim>(qR);

                double sL = std::min(primL.v[d] - primL.c, primR.v[d] - primR.c);
                double sR = std::max(primL.v[d] + primL.c, primR.v[d] + primR.c);

                if (sL >= 0)
                {
                    flux = compute_flux<d>(primL);
                }
                else if (sL < 0 && sR > 0)
                {
                    flux = (sR * compute_flux<d>(primL) - sL * compute_flux<d>(primR) + sL * sR * (qR - qL)) / (sR - sL);
                }
                else if (sR <= 0)
                {
                    flux = compute_flux<d>(primR);
                }
            };
        });
    auto scheme = make_flux_based_scheme(hll);
    scheme.set_name("hll");

    return scheme;
}