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

double calcPrice(double S, double K, double r, double v, double T, int nSim, bool isCall = true)
{
    optionParams p{S,K,r,v,T};    
    return isCall ? montecarloCallPrice(nSim, p) : montecarloPutPrice(nSim, p);
}

double calcDelta(double S, double K, double r, double v, double T, int nSim, bool isCall = true)
{
    double h = 0.1; //bump size
    double price_up = calcPrice(S+h, K, r, v, T, nSim);
    double price_down = calcPrice(S-h, K, r, v, T, nSim);
    return (price_up - price_down) / (2 * h);
}

double calcGamma(double S, double K, double r, double v, double T, int nSim,  bool isCall = true)
{
    double h = 0.1; //bump size
    double price_up = calcPrice(S+h, K, r, v, T, nSim);
    double price_down = calcPrice(S-h, K, r, v, T, nSim);
    double price_base = calcPrice(S, K, r, v, T, nSim);
    return (price_up + price_down - 2 * price_base) / (h * h);
}

double calcVega(double S, double K, double r, double v, double T, int nSim,  bool isCall = true)
{
    double h = 0.001; //bump size
    double price_up = calcPrice(S+h, K, r, v, T, nSim);
    double price_down = calcPrice(S-h, K, r, v, T, nSim);
    return (price_up - price_down) / (2 * h) * 100;
}

double calcTheta(double S, double K, double r, double v, double T, int nSim,  bool isCall = true)
{
    double h = 1.0/365.0;
    double price_up = calcPrice(S+h, K, r, v, T, nSim);
    double price_down = calcPrice(S-h, K, r, v, T, nSim);
    return (price_up - price_down) / h ;
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
        ("h,help", "Print help") 
        ("delta", "Compute Delta")
        ("gamma", "Compute Gamma") 
        ("vega", "Compute Vega")
        ("theta", "Compute Theta")
        ("greeks", "Compute all Greeks");

    auto args = options.parse(argc, argv);
    if (args.count("help")) { std::cout << options.help() << std::endl; return 0; }

    double S = args["spot"].as<double>();
    double K = args["strike"].as<double>();
    double r = args["r"].as<double>();
    double v = args["v"].as<double>();
    double T = args["days"].as<double>() / 365.0;
    int nSim = args["sims"].as<int>();

    double callPrice = calcPrice(S, K, r, v, T, nSim, true);
    double PutPrice = calcPrice(S, K, r, v, T, nSim, false);
    std::cout << "Call: " << callPrice << std::endl;
    std::cout << "Put: " << PutPrice << std::endl;

    if (args.count("delta")) { std::cout << "Delta: " << calcDelta(S, K, r, v, T, nSim) << std::endl;}
    if (args.count("gamma")) { std::cout << "Gamma: " << calcGamma(S, K, r, v, T, nSim) << std::endl;}
    if (args.count("vega")) { std::cout << "Vega: " << calcVega(S, K, r, v, T, nSim) << std::endl;}
    if (args.count("theta")) { std::cout << "Theta: " << calcTheta(S, K, r, v, T, nSim) << std::endl;}
    if (args.count("greeks")) 
    {
        std::cout << "Delta: " << calcDelta(S, K, r, v, T, nSim) << std::endl;
        std::cout << "Gamma: " << calcGamma(S, K, r, v, T, nSim) << std::endl;
        std::cout << "Vega: " << calcVega(S, K, r, v, T, nSim) << std::endl;
        std::cout << "Theta: " << calcTheta(S, K, r, v, T, nSim) << std::endl;
    }
       
    return 0;
}
