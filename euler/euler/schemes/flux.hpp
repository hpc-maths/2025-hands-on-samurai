#pragma once

#include "../variables.hpp"

template <std::size_t d, std::size_t dim>
auto compute_flux(const PrimState<dim>& prim)
{
    auto flux = xt::xtensor_fixed<double, xt::xshape<dim + 2>>{};

    flux[EulerConsVar::rho]  = prim.rho * prim.v[d];
    auto e                   = EOS::stiffened_gas::e(prim.rho, prim.p);
    flux[EulerConsVar::rhoE] = (prim.rho * e + prim.p) * prim.v[d];
    for (std::size_t i = 0; i < dim; ++i)
    {
        flux[EulerConsVar::rhou + i] = prim.rho * prim.v[i] * prim.v[d];
        flux[EulerConsVar::rhoE] += 0.5 * prim.rho * prim.v[i] * prim.v[i] * prim.v[d];
    }
    flux[EulerConsVar::rhou + d] += prim.p;
    return flux;
}