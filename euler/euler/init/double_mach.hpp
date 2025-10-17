#pragma once

#include <samurai/bc.hpp>
#include <samurai/box.hpp>

#include "../variables.hpp"
#include "registry.hpp"

namespace test_case::double_mach_reflection
{
    double alpha = std::numbers::pi / 3.;
    double x0    = 2. / 3;

    PrimState<2> left_state{
        8.,
        EOS::e(8., 116.5),
        116.5,
        EOS::c(8., 116.5),
        xt::xtensor_fixed<double, xt::xshape<2>>{8.25 * std::sin(alpha), -8.25 * std::cos(alpha)}
    };

    PrimState<2> right_state{
        1.4,
        EOS::e(1.4, 1.),
        1.,
        EOS::c(1.4, 1.),
        xt::xtensor_fixed<double, xt::xshape<2>>{0., 0.}
    };

    auto init_fn = [](auto& u, auto& cell)
    {
        auto x = cell.center();

        if (x[0] < x0 + x[1] / std::tan(alpha))
        {
            u[cell] = prim2cons<2>(left_state);
        }
        else
        {
            u[cell] = prim2cons<2>(right_state);
        }
    };

    template <class Field>
    struct ValueImpl : public samurai::Bc<Field>
    {
        INIT_BC(ValueImpl, 2)

        stencil_t get_stencil(constant_stencil_size_t) const override
        {
            return samurai::line_stencil<dim, 0, 2>();
        }

        apply_function_t get_apply_function(constant_stencil_size_t, const direction_t&) const override
        {
            return [](Field& u, const stencil_cells_t& cells, const value_t& dirichlet_value)
            {
                u[cells[1]] = dirichlet_value;
            };
        }
    };

    struct Value
    {
        template <class Field>
        using impl_t = ValueImpl<Field>;
    };

    void bc_fn(auto& u, double& t)
    {
        static constexpr std::size_t dim = std::decay_t<decltype(u)>::dim;

        const xt::xtensor_fixed<int, xt::xshape<dim>> bottom = {0, -1};
        samurai::make_bc<Value>(u,
                                [&](const auto&, const auto& cell, const auto&)
                                {
                                    if (cell.center(0) < x0)
                                    {
                                        return prim2cons(left_state);
                                    }
                                    else
                                    {
                                        return xt::xtensor_fixed<double, xt::xshape<dim + 2>>{u[cell][EulerConsVar::rho],
                                                                                              u[cell][EulerConsVar::rhoE],
                                                                                              u[cell][EulerConsVar::rhou],
                                                                                              -u[cell][EulerConsVar::rhou + 1]};
                                    }
                                })
            ->on(bottom);

        const xt::xtensor_fixed<int, xt::xshape<dim>> top = {0, 1};
        samurai::make_bc<Value>(u,
                                [&](const auto&, const auto& cell, const auto&)
                                {
                                    double x1 = x0 + 10 * t / std::sin(alpha) + 1 / std::tan(alpha);
                                    if (cell.center(0) < x1)
                                    {
                                        return prim2cons(left_state);
                                    }
                                    else
                                    {
                                        return prim2cons(right_state);
                                    }
                                })
            ->on(top);

        const xt::xtensor_fixed<int, xt::xshape<dim>> right = {1, 0};
        samurai::make_bc<samurai::Neumann<1>>(u, 0, 0, 0, 0)->on(right);

        const xt::xtensor_fixed<int, xt::xshape<dim>> left = {-1, 0};
        samurai::make_bc<Value>(
            u,
            left_state.rho,
            left_state.rho * (left_state.e + 0.5 * (left_state.v[0] * left_state.v[0] + left_state.v[1] * left_state.v[1])),
            left_state.rho * left_state.v[0],
            left_state.rho * left_state.v[1])
            ->on(left);
    }

    template <std::size_t dim>
    auto box_fn()
    {
        xt::xtensor_fixed<double, xt::xshape<dim>> min_corner = {0., 0.};
        xt::xtensor_fixed<double, xt::xshape<dim>> max_corner = {4., 1.};

        return samurai::Box<double, dim>(min_corner, max_corner);
    }
}

REGISTER_TEST_CASE(double_mach_reflection,
                   test_case::double_mach_reflection::box_fn,
                   test_case::double_mach_reflection::init_fn,
                   test_case::double_mach_reflection::bc_fn)
