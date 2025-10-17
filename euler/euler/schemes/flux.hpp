#pragma once

#include "../variables.hpp"

template <std::size_t d, std::size_t dim>
auto compute_flux(const PrimState<dim>& prim)
{
    auto flux = xt::xtensor_fixed<double, xt::xshape<dim + 2>>{};

    flux[EulerConsVar::rho]  = prim.rho * prim.v[d];
    flux[EulerConsVar::rhoE] = (prim.rho * (prim.e + 0.5 * xt::sum(xt::square(prim.v))()) + prim.p) * prim.v[d];
    for (std::size_t i = 0; i < dim; ++i)
    {
        flux[EulerConsVar::rhou + i] = prim.rho * prim.v[i] * prim.v[d];
    }
    flux[EulerConsVar::rhou + d] += prim.p;
    return flux;
}