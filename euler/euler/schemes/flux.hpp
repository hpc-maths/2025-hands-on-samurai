#pragma once

#include "../variables.hpp"

template <std::size_t dim, std::size_t d>
auto compute_flux(double rho, const auto& v, double e, double p)
{
    auto flux = xt::xtensor_fixed<double, xt::xshape<dim + 2>>{};

    flux[EulerVariable::rho]  = rho * v[d];
    flux[EulerVariable::rhoE] = (rho * (e + 0.5 * xt::sum(xt::square(v))()) + p) * v[d];
    for (std::size_t i = 0; i < dim; ++i)
    {
        flux[EulerVariable::rhou + i] = rho * v[i] * v[d];
    }
    flux[EulerVariable::rhou + d] += p;
    return flux;
}