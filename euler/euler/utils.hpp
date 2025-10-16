#pragma once

#include <type_traits>

#include <samurai/algorithm.hpp>

#include "variables.hpp"

auto get_max_lambda(const auto& u)
{
    static constexpr std::size_t dim = std::decay_t<decltype(u)>::dim;
    double res                       = 0.;

    const auto& mesh = u.mesh();

    samurai::for_each_cell(mesh,
                           [&](const auto& cell)
                           {
                               auto [rho, vel, e, p, c] = extract_primitive<dim>(u[cell]);
                               double sum               = 0;
                               for (std::size_t d = 0; d < dim; ++d)
                               {
                                   sum += std::abs(vel[d]) + c;
                                   //    res = std::max(std::abs(vel[d]) + c, res);
                               }
                               res = std::max(sum, res);
                           });
    return res;
}