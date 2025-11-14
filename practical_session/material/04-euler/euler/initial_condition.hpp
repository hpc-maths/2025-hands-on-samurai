// Copyright 2025 the samurai team
// SPDX-License-Identifier:  BSD-3-Clause

#pragma once

#include <numbers>

#include "variables.hpp"

auto riemann_config_3(auto& u)
{
    double x0 = 0.5;
    double y0 = 0.5;

    PrimState<2> quad1_state{
        1.5,
        1.5,
        xt::xtensor_fixed<double, xt::xshape<2>>{0., 0.}
    };

    PrimState<2> quad2_state{
        0.5323,
        0.3,
        xt::xtensor_fixed<double, xt::xshape<2>>{1.206, 0.}
    };

    PrimState<2> quad3_state{
        0.138,
        0.29,
        xt::xtensor_fixed<double, xt::xshape<2>>{1.206, 1.206}
    };

    PrimState<2> quad4_state{
        0.5323,
        0.3,
        xt::xtensor_fixed<double, xt::xshape<2>>{0, 1.206}
    };

    samurai::for_each_cell(u.mesh(),
                           [&](auto& cell)
                           {
                               auto x = cell.center();

                               if (x[0] >= x0 && x[1] >= y0)
                               {
                                   u[cell] = prim2cons<2>(quad1_state);
                               }
                               else if (x[0] < x0 && x[1] >= y0)
                               {
                                   u[cell] = prim2cons<2>(quad2_state);
                               }
                               else if (x[0] < x0 && x[1] < y0)
                               {
                                   u[cell] = prim2cons<2>(quad3_state);
                               }
                               else // (x[0] >= x0 && x[1] < y0)
                               {
                                   u[cell] = prim2cons<2>(quad4_state);
                               }
                           });
}
