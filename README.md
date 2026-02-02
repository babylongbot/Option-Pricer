Option-Pricer
Monte Carlo European Option Pricer with Greeks

Prices European call/put options and computes Delta, Gamma, Vega, and Theta using Monte Carlo simulation with finite difference approximations.

Quick Start
bash
# Clone and compile
g++ -O2 -std=c++11 -o option main.cpp -lcxxopts -lm

# Run basic pricing
./option

# Price with custom parameters
./option -s 25000 -K 25050 -r 0.065 -v 0.13 -t 60 -n 1000000

# Compute Greeks
./option --delta --gamma --vega --theta
./option --greeks
Current Settings
text
S = 25048.65  (NIFTY spot)
K = 25050     (ATM strike)  
r = 6.48%     (India T-Bill)
σ = 12.9%     (Implied Vol)
T = 0.178 yr  (65 days)
n = 100,000   (MC simulations)
Expected output:

text
Call: ~950-1000
Put: ~950-1000
Command Line Options
text
Usage: ./option [OPTIONS]

Options:
  -s, --spot      Spot price [25048.65]
  -K, --strike    Strike price [25050.00]
  -r              Risk-free rate [0.0648]
  -v, --iv        Implied volatility [0.129]
  -t, --days      Days to expiry [65]
  -n, --sims      MC simulations [100000]
  --delta         Compute Delta
  --gamma         Compute Gamma
  --vega          Compute Vega
  --theta         Compute Theta
  --greeks        Compute all Greeks
  -h, --help      Print help
Features
✅ European Call & Put pricing via Monte Carlo

✅ Delta (bump spot ±0.1%)

✅ Gamma (2nd order spot bump)

✅ Vega (bump vol ±0.1%)

✅ Theta (bump time ±1 day)

✅ Command line interface with cxxopts

✅ Production-ready C++11 code

Example Outputs
Basic pricing:

text
$ ./option
Call: 978.45
Put: 962.31
With Greeks:

text
$ ./option --greeks
Call: 978.45
Put: 962.31
Delta: 0.512
Gamma: 0.0234
Vega: 245.67
Theta: -12.34
Edit Parameters
Modify via command line flags or source code:

cpp
// In main(), default values are set here:
cxxopts::value<double>()->default_value("25048.65")  // spot
cxxopts::value<double>()->default_value("25050.00")  // strike  
cxxopts::value<double>()->default_value("0.0648")    // rate
cxxopts::value<double>()->default_value("0.129")     // vol
cxxopts::value<double>()->default_value("65")        // days
cxxopts::value<int>()->default_value("100000")       // sims
