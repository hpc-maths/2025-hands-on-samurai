#pragma once

#include "eos.hpp"

namespace EulerConsVar
{
    static constexpr std::size_t rho  = 0;
    static constexpr std::size_t rhoE = 1;
    static constexpr std::size_t rhou = 2;
}

template <std::size_t dim>
struct PrimState
{
    double rho;
    double p;
    xt::xtensor_fixed<double, xt::xshape<dim>> v;
};

template <std::size_t dim>
auto cons2prim(const xt::xtensor_fixed<double, xt::xshape<dim + 2>>& conserved)
{
    PrimState<dim> primitives;

    primitives.rho = conserved[EulerConsVar::rho];
    auto e         = conserved[EulerConsVar::rhoE] / conserved[EulerConsVar::rho];
    for (std::size_t d = 0; d < dim; ++d)
    {
        primitives.v[d] = conserved[EulerConsVar::rhou + d] / conserved[EulerConsVar::rho];
        e -= 0.5 * (primitives.v[d] * primitives.v[d]);
    }
    primitives.p = EOS::stiffened_gas::p(primitives.rho, e);
    return primitives;
}

template <std::size_t dim>
auto prim2cons(const PrimState<dim>& primitives)
{
    xt::xtensor_fixed<double, xt::xshape<dim + 2>> conserved;

    conserved[EulerConsVar::rho]  = primitives.rho;
    auto e                        = EOS::stiffened_gas::e(primitives.rho, primitives.p);
    conserved[EulerConsVar::rhoE] = e * conserved[EulerConsVar::rho];
    for (std::size_t d = 0; d < dim; ++d)
    {
        conserved[EulerConsVar::rhou + d] = primitives.v[d] * conserved[EulerConsVar::rho];
        conserved[EulerConsVar::rhoE] += 0.5 * primitives.v[d] * primitives.v[d] * conserved[EulerConsVar::rho];
    }
    return conserved;
}