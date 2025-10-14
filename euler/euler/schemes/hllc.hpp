#pragma once

#include <samurai/schemes/fv.hpp>

#include "../variables.hpp"
#include "flux.hpp"

template <std::size_t dim, std::size_t d>
auto compute_star_state(const auto& q, double s, double s_star)
{
    using state_t = std::decay_t<decltype(q)>;

    state_t q_star;

    auto [rho, v, e, p, c] = extract_primitive<dim>(q);

    auto rho_star = rho * (s - v[d]) / (s - s_star);

    q_star[EulerVariable::rho] = rho_star;

    for (std::size_t i = 0; i < dim; ++i)
    {
        q_star[EulerVariable::rhou + i] = rho_star * v[i];
    }
    q_star[EulerVariable::rhou + d] = rho_star * s_star;
    q_star[EulerVariable::rhoE]     = rho_star * (e + 0.5 * xt::sum(xt::square(v))() + (s_star - v[d]) * (s_star + p / (rho * (s - v[d]))));

    return q_star;
}

template <class Field>
auto make_euler_hllc()
{
    static constexpr std::size_t dim          = Field::dim;
    static constexpr std::size_t field_size   = Field::n_comp;
    static constexpr std::size_t stencil_size = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, field_size, stencil_size, Field>;

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

                const auto& qL              = field[left];
                auto [rhoL, vL, eL, pL, cL] = extract_primitive<dim>(qL);

                const auto& qR              = field[right];
                auto [rhoR, vR, eR, pR, cR] = extract_primitive<dim>(qR);

                double sL = std::min(vL[d] - cL, vR[d] - cR);
                double sR = std::max(vL[d] + cL, vR[d] + cR);
                double sM = (rhoL * vL[d] * (sL - vL[d]) - pL - rhoR * vR[d] * (sR - vR[d]) + pR)
                          / (rhoL * (sL - vL[d]) - rhoR * (sR - vR[d]));

                if (sL >= 0)
                {
                    flux = compute_flux<dim, d>(rhoL, vL, eL, pL);
                }
                else if (sL < 0 && sM >= 0)
                {
                    flux = compute_flux<dim, d>(rhoL, vL, eL, pL) + sL * (compute_star_state<dim, d>(qL, sL, sM) - qL);
                }
                else if (sM < 0 && sR >= 0)
                {
                    flux = compute_flux<dim, d>(rhoR, vR, eR, pR) + sR * (compute_star_state<dim, d>(qR, sR, sM) - qR);
                }
                else if (sR < 0)
                {
                    flux = compute_flux<dim, d>(rhoR, vR, eR, pR);
                }
            };
        });

    auto scheme = make_flux_based_scheme(hllc);
    scheme.set_name("hllc");

    return scheme;
}
