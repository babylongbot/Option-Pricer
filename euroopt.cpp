#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <random>

struct optionParams
{
    double S;
    double K;
    double r;
    double v;
    double T;    
};

double montecarloCallPrice( const int nSim, const optionParams& p)
{
    double drift = (p.r - 0.5 * p.v * p.v)* p.T;
    double diffusion = p.v *sqrt(p.T);

    double payoffSum = 0.0;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);

    for(int i=0; i<nSim; i++)
    {
        double Z = distribution(generator);
        double SForward = p.S * exp(drift+ diffusion * Z);

        payoffSum += std::max(SForward - p.K, 0.0);
    }

    return (payoffSum / nSim) * exp(-p.r * p.T);
}

double montecarloPutPrice( const int nSim, const optionParams& p)
{
    double drift = (p.r - 0.5 * p.v * p.v)* p.T;
    double diffusion = p.v *sqrt(p.T);

    double payoffSum = 0.0;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);

    for(int i=0; i<nSim; i++)
    {
        double Z = distribution(generator);
        double SForward = p.S * exp(drift+ diffusion * Z);

        payoffSum += std::max(p.K - SForward, 0.0);
    }

    return (payoffSum / nSim) * exp(-p.r * p.T);
}

int main()
{
    int nSim = 100000;

    optionParams params;
    params.S = 25048.65;
    params.K = 25050;
    params.r = 0.0648;
    params.v = 0.129;
    params.T = 0.178;

    double montecarloCallPrice_ = montecarloCallPrice(nSim, params);
    double montecarloPutPrice_ = montecarloPutPrice(nSim, params);

    std::cout << "\n Monte Carlo call price: " << montecarloCallPrice_;
    std::cout << "\n Monte Carlo put price: " << montecarloPutPrice_;

    std::cout << std::endl;
    return 0;
}
