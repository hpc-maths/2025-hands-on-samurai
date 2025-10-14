#pragma once

#include "eos.hpp"

namespace EulerVariable
{
    static constexpr std::size_t rho  = 0;
    static constexpr std::size_t rhoE = 1;
    static constexpr std::size_t rhou = 2;
}

template <std::size_t dim>
auto extract_primitive(const auto& q)
{
    auto rho = q[EulerVariable::rho];
    auto v   = xt::xtensor_fixed<double, xt::xshape<dim>>{};
    for (std::size_t i = 0; i < dim; ++i)
    {
        v[i] = q[EulerVariable::rhou + i] / rho;
    }
    auto e = q[EulerVariable::rhoE] / rho - 0.5 * xt::sum(xt::square(v))();
    auto p = EOS::p(rho, e);
    auto c = EOS::c(rho, p);
    return std::make_tuple(rho, v, e, p, c);
}
