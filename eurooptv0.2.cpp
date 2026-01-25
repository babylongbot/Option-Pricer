#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <random>
#include "cxxopts.hpp"

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

int runpricing(double S, double K, double r, double v, double T, int nSim)
{
    optionParams p{S,K,r,v,T};
    auto call = montecarloCallPrice(nSim, p);
    auto put = montecarloPutPrice(nSim, p);
    std::cout << "Call: " << call <<  " Put: " << put << std::endl;
    return 0;
}


int main(int argc, char** argv)
{
    cxxopts::Options options("Option_pricer", "Monte Carlo option pricing");
    options.add_options()
        ("s,spot", "Spot price", cxxopts::value<double>()->default_value("25048.65"))
        ("K,strike", "Strike price", cxxopts::value<double>()->default_value("25050.00"))
        ("r", "Risk-free rate", cxxopts::value<double>()->default_value("0.0648"))
        ("v,iv", "Implied vol", cxxopts::value<double>()->default_value("0.129"))
        ("t,days", "Days to expiry", cxxopts::value<double>()->default_value("65"))
        ("n,sims", "MC simulations", cxxopts::value<int>()->default_value("100000"))
        ("h,help", "Print help");

    auto args = options.parse(argc, argv);
    if (args.count("help")) { std::cout << options.help() << std::endl; return 0; }

    runpricing(args["spot"].as<double>(), args["strike"].as<double>(),
    args["r"].as<double>(), args["v"].as<double>(), args["days"].as<double>()/365,
    args["sims"].as<int>());
       
    int nSim = 100000;

    optionParams params;
    params.S = 25048.65;
    params.K = 25050;
    params.r = 0.0648;
    params.v = 0.129;
    params.T = 0.178;

    double montecarloCallPrice_ = montecarloCallPrice(nSim, params);
    double montecarloPutPrice_ = montecarloPutPrice(nSim, params);

    //std::cout << "\n Monte Carlo call price: " << montecarloCallPrice_;
    //std::cout << "\n Monte Carlo put price: " << montecarloPutPrice_;

    std::cout << std::endl;
    return 0;
}
