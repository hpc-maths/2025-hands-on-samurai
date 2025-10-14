#pragma once

struct EOS
{
    static constexpr double gamma  = 1.4;
    static constexpr double pi_inf = 0.;
    static constexpr double q_inf  = 0.;

    static double p(double rho, double e)
    {
        return (gamma - 1.0) * rho * (e - q_inf) - gamma * pi_inf;
    }

    static double c(double rho, double p)
    {
        return std::sqrt(gamma * (p + pi_inf) / rho);
    }

    static double e(double rho, double p)
    {
        return (p + gamma * pi_inf) / ((gamma - 1.0) * rho) + q_inf;
    }
};