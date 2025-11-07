#pragma once

#include "eos.hpp"

template <std::size_t Dim, std::size_t NSpecies = 0>
struct EulerLayout
{
    static constexpr std::size_t rho  = 0;
    static constexpr std::size_t rhoE = 1;

    static constexpr std::size_t mom(std::size_t d)
    {
        assert(d < Dim);
        return 2 + d;
    } // d=0..Dim-1

    static constexpr std::size_t species(std::size_t i)
    {
        assert(i < NSpecies);
        return 2 + Dim + i;
    } // i=0..NSpecies-1

    static constexpr std::size_t first_mom     = 2;
    static constexpr std::size_t first_species = 2 + Dim;
    static constexpr std::size_t size          = 2 + Dim + NSpecies;
};

template <std::size_t Dim, std::size_t NSpecies = 0>
struct PrimState;

template <std::size_t Dim>
struct PrimState<Dim, 0>
{
    double rho;
    double p;
    xt::xtensor_fixed<double, xt::xshape<Dim>> v;
};

template <std::size_t Dim, std::size_t NSpecies>
    requires(NSpecies > 0)
struct PrimState<Dim, NSpecies>
{
    double rho;
    double p;
    xt::xtensor_fixed<double, xt::xshape<Dim>> v;
    xt::xtensor_fixed<double, xt::xshape<NSpecies>> Y;
};

template <std::size_t Dim, std::size_t NSpecies = 0>
auto cons2prim(const xt::xtensor_fixed<double, xt::xshape<EulerLayout<Dim, NSpecies>::size>>& conserved)
{
    using EulerConsVar = EulerLayout<Dim, NSpecies>;

    PrimState<Dim, NSpecies> primitives;

    primitives.rho = conserved[EulerConsVar::rho];
    auto e         = conserved[EulerConsVar::rhoE] / conserved[EulerConsVar::rho];
    for (std::size_t d = 0; d < Dim; ++d)
    {
        primitives.v[d] = conserved[EulerConsVar::mom(d)] / conserved[EulerConsVar::rho];
        e -= 0.5 * (primitives.v[d] * primitives.v[d]);
    }
    primitives.p = EOS::stiffened_gas::p(primitives.rho, e);
    return primitives;
}

template <std::size_t Dim, std::size_t NSpecies>
auto prim2cons(const PrimState<Dim, NSpecies>& primitives)
{
    using EulerConsVar = EulerLayout<Dim, NSpecies>;

    xt::xtensor_fixed<double, xt::xshape<EulerConsVar::size>> conserved;

    conserved[EulerConsVar::rho]  = primitives.rho;
    auto e                        = EOS::stiffened_gas::e(primitives.rho, primitives.p);
    conserved[EulerConsVar::rhoE] = e * conserved[EulerConsVar::rho];
    for (std::size_t d = 0; d < Dim; ++d)
    {
        conserved[EulerConsVar::mom(d)] = primitives.v[d] * conserved[EulerConsVar::rho];
        conserved[EulerConsVar::rhoE] += 0.5 * primitives.v[d] * primitives.v[d] * conserved[EulerConsVar::rho];
    }
    return conserved;
}