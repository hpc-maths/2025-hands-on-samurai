#pragma once

#include <samurai/schemes/fv.hpp>

#include "../variables.hpp"
#include "flux.hpp"

template <std::size_t d, std::size_t dim>
auto compute_star_state(const PrimState<dim>& prim, double s, double s_star)
{
    xt::xtensor_fixed<double, xt::xshape<dim + 2>> q_star;

    auto rho_star = prim.rho * (s - prim.v[d]) / (s - s_star);

    q_star[EulerConsVar::rho] = rho_star;

    for (std::size_t i = 0; i < dim; ++i)
    {
        q_star[EulerConsVar::rhou + i] = rho_star * prim.v[i];
    }
    q_star[EulerConsVar::rhou + d] = rho_star * s_star;
    q_star[EulerConsVar::rhoE]     = rho_star
                               * (prim.e + 0.5 * xt::sum(xt::square(prim.v))()
                                  + (s_star - prim.v[d]) * (s_star + prim.p / (prim.rho * (s - prim.v[d]))));

    return q_star;
}

template <class Field>
auto make_euler_hllc()
{
    static constexpr std::size_t dim          = Field::dim;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, stencil_size, Field, Field>;

    samurai::FluxDefinition<cfg> hllc;

    samurai::static_for<0, dim>::apply( // for each positive Cartesian direction 'd'
        [&](auto _d)
        {
            static constexpr std::size_t d = _d();

            hllc[d].cons_flux_function =
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
                double sM = (primL.rho * primL.v[d] * (sL - primL.v[d]) - primL.p - primR.rho * primR.v[d] * (sR - primR.v[d]) + primR.p)
                          / (primL.rho * (sL - primL.v[d]) - primR.rho * (sR - primR.v[d]));

                if (sL >= 0)
                {
                    flux = compute_flux<d>(primL);
                }
                else if (sL < 0 && sM >= 0)
                {
                    flux = compute_flux<d>(primL) + sL * (compute_star_state<d>(primL, sL, sM) - qL);
                }
                else if (sM < 0 && sR >= 0)
                {
                    flux = compute_flux<d>(primR) + sR * (compute_star_state<d>(primR, sR, sM) - qR);
                }
                else if (sR < 0)
                {
                    flux = compute_flux<d>(primR);
                }
            };
        });

    auto scheme = make_flux_based_scheme(hllc);
    scheme.set_name("hllc");

    return scheme;
}
