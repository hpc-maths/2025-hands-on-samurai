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

                const auto& qL              = field[left];
                auto [rhoL, vL, eL, pL, cL] = extract_primitive<dim>(qL);

                const auto& qR              = field[right];
                auto [rhoR, vR, eR, pR, cR] = extract_primitive<dim>(qR);

                double sL = std::min(vL[d] - cL, vR[d] - cR);
                double sR = std::max(vL[d] + cL, vR[d] + cR);

                if (sL >= 0)
                {
                    flux = compute_flux<dim, d>(rhoL, vL, eL, pL);
                }
                else if (sL < 0 && sR > 0)
                {
                    // std::cout << "sL = " << sL << ", sR = " << sR << std::endl;
                    // std::cout << "left flux: " << compute_flux<dim, d>(rhoL, vL, eL, pL) << std::endl;
                    // std::cout << "right flux: " << compute_flux<dim, d>(rhoR, vR, eR, pR) << std::endl;
                    // std::cout << "qL: " << qL << std::endl;
                    // std::cout << "qR: " << qR << std::endl;
                    flux = (sR * compute_flux<dim, d>(rhoL, vL, eL, pL) - sL * compute_flux<dim, d>(rhoR, vR, eR, pR) + sL * sR * (qR - qL))
                         / (sR - sL);
                    // std::cout << "HLL Computed flux: " << flux << std::endl;
                }
                else if (sR <= 0)
                {
                    flux = compute_flux<dim, d>(rhoR, vR, eR, pR);
                }
                // std::cout << "Left state: rho = " << rhoL << ", v = " << vL[d] << ", p = " << pL << ", e = " << eL << ", c = " << cL
                //           << std::endl;
                // std::cout << "right state: rho = " << rhoR << ", v = " << vR[d] << ", p = " << pR << ", e = " << eR << ", c = " << cR
                //           << std::endl;

                // std::cout << "HLL Computed flux: " << flux << std::endl;
            };
        });
    auto scheme = make_flux_based_scheme(hll);
    scheme.set_name("hll");

    return scheme;
}