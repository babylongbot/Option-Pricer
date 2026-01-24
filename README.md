# Option-Pricer
Monte Carlo Option Pricer
=======================

Prices European call and put options using Monte Carlo simulation.

QUICK START
-----------
1. Save code as euroopt.cpp
2. Compile: g++ -O2 -std=c++11 -o option euroopt.cpp -lm
3. Run: ./option

CURRENT SETTINGS
----------------
S = 25048.65  (NIFTY spot)
K = 25050     (ATM strike) 
r = 6.48%     (India T-Bill)
σ = 12.9%     (Implied Vol)
T = 0.178 yr  (65 days)

Expected output:
Monte Carlo call price: ~950-1000
Monte Carlo put price: ~950-1000

EDIT PARAMETERS
---------------
Change values in main() function:
params.S = YOUR_SPOT;
params.K = YOUR_STRIKE;
params.r = RISK_FREE_RATE;
params.v = IMPLIED_VOL;
params.T = DAYS_TO_EXPIRY / 365.0;
