// Copyright 2025 SAMURAI TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
//

#pragma once

#include <samurai/numeric/prediction.hpp>
#include <samurai/operators_base.hpp>

#include "variables.hpp"

template <std::size_t dim, class TInterval>
class Euler_prediction_op : public samurai::field_operator_base<dim, TInterval>
{
  public:

    INIT_OPERATOR(Euler_prediction_op)

    inline void operator()(samurai::Dim<2>, auto& dest, const auto& src) const
    {
        using EulerConsVar = EulerLayout<dim>;
        using field_t      = std::decay_t<decltype(src)>;

        /*--- Compute unchanged prediction ---*/
        constexpr std::size_t pred_order = field_t::mesh_t::config::prediction_order;

        auto ii = i << 1;
        ii.step = 2;

        auto jj = j << 1;

        auto qs_i  = samurai::Qs_i<pred_order>(src, level, i, j);
        auto qs_j  = samurai::Qs_j<pred_order>(src, level, i, j);
        auto qs_ij = samurai::Qs_ij<pred_order>(src, level, i, j);

        dest(level + 1, ii, jj)         = src(level, i, j) + qs_i + qs_j - qs_ij;
        dest(level + 1, ii + 1, jj)     = src(level, i, j) - qs_i + qs_j + qs_ij;
        dest(level + 1, ii, jj + 1)     = src(level, i, j) + qs_i - qs_j + qs_ij;
        dest(level + 1, ii + 1, jj + 1) = src(level, i, j) - qs_i - qs_j - qs_ij;

        const auto mask_rho = (dest(EulerConsVar::rho, level + 1, ii, jj) < 0.0 || dest(EulerConsVar::rho, level + 1, ii + 1, jj) < 0.0
                               || dest(EulerConsVar::rho, level + 1, ii, jj + 1) < 0.0
                               || dest(EulerConsVar::rho, level + 1, ii + 1, jj + 1) < 0.0);

        samurai::apply_on_masked(mask_rho,
                                 [&](auto& ie)
                                 {
                                     xt::view(dest(level + 1, ii, jj), ie)         = xt::view(src(level, i, j), ie);
                                     xt::view(dest(level + 1, ii + 1, jj), ie)     = xt::view(src(level, i, j), ie);
                                     xt::view(dest(level + 1, ii, jj + 1), ie)     = xt::view(src(level, i, j), ie);
                                     xt::view(dest(level + 1, ii + 1, jj + 1), ie) = xt::view(src(level, i, j), ie);
                                 });

        samurai::static_nested_loop<dim - 1, 0, 2>(
            [&](const auto& stencil)
            {
                {
                    auto rho = dest(EulerConsVar::rho, level + 1, ii, jj + stencil[0]);
                    auto e   = xt::eval(dest(EulerConsVar::rhoE, level + 1, ii, jj + stencil[0])
                                      / dest(EulerConsVar::rho, level + 1, ii, jj + stencil[0]));

                    for (std::size_t d = 0; d < dim; ++d)
                    {
                        auto v_d = dest(EulerConsVar::mom(d), level + 1, ii, jj + stencil[0])
                                 / dest(EulerConsVar::rho, level + 1, ii, jj + stencil[0]);
                        e -= 0.5 * v_d * v_d;
                    }
                    auto p = EOS::stiffened_gas::p(rho, e);

                    samurai::apply_on_masked(p < 0,
                                             [&](auto& ie)
                                             {
                                                 xt::view(dest(level + 1, ii, jj), ie)         = xt::view(src(level, i, j), ie);
                                                 xt::view(dest(level + 1, ii + 1, jj), ie)     = xt::view(src(level, i, j), ie);
                                                 xt::view(dest(level + 1, ii, jj + 1), ie)     = xt::view(src(level, i, j), ie);
                                                 xt::view(dest(level + 1, ii + 1, jj + 1), ie) = xt::view(src(level, i, j), ie);
                                             });
                }
                {
                    auto rho = dest(EulerConsVar::rho, level + 1, ii + 1, jj + stencil[0]);
                    auto e   = xt::eval(dest(EulerConsVar::rhoE, level + 1, ii + 1, jj + stencil[0])
                                      / dest(EulerConsVar::rho, level + 1, ii + 1, jj + stencil[0]));

                    for (std::size_t d = 0; d < dim; ++d)
                    {
                        auto v_d = dest(EulerConsVar::mom(d), level + 1, ii + 1, jj + stencil[0])
                                 / dest(EulerConsVar::rho, level + 1, ii + 1, jj + stencil[0]);
                        e -= 0.5 * v_d * v_d;
                    }
                    auto p = EOS::stiffened_gas::p(rho, e);

                    samurai::apply_on_masked(p < 0,
                                             [&](auto& ie)
                                             {
                                                 xt::view(dest(level + 1, ii, jj), ie)         = xt::view(src(level, i, j), ie);
                                                 xt::view(dest(level + 1, ii + 1, jj), ie)     = xt::view(src(level, i, j), ie);
                                                 xt::view(dest(level + 1, ii, jj + 1), ie)     = xt::view(src(level, i, j), ie);
                                                 xt::view(dest(level + 1, ii + 1, jj + 1), ie) = xt::view(src(level, i, j), ie);
                                             });
                }
            });
    }
};
