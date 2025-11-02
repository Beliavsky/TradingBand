# TradingBand
Do simulations to study optimal trading bands in a trading system to maximize Sharpe ratio after transaction costs. `compile.bat` compiles the Fortran program. Run the Python program with `python xtrading_band.py`. Output of the Fortran program compiled by gfortran is

```
parameters
 days_year   252.00000000000000     
n=5, T=200, sims=100, mu= 0.001000, sigma= 0.100000, rho= 0.300000, c= 0.005000, gamma= 3.000000, seed=7
delta_base (cube-root) =   0.61979809

sharpe-oriented calibration
best k =  0.400000
best delta =   0.24791924
best after-cost sharpe =  0.263153

pre-cost results
strategy        pre_mean      pre_vol       pre_sharpe  
100*delta           0.001152      0.071496      0.245690
10*delta            0.001152      0.071496      0.245690
3*delta             0.001152      0.071496      0.245690
delta               0.001185      0.069613      0.265969
delta/3             0.001253      0.067235      0.296122
delta/10            0.001287      0.066665      0.306681
delta/100           0.001297      0.066564      0.310145
full_rebal (0)      0.001301      0.066561      0.311162
buy_and_hold        0.001152      0.071496      0.245690

after-cost results
strategy        after_mean    after_vol     after_sharpe
100*delta           0.001152      0.071496      0.245690
10*delta            0.001152      0.071496      0.245690
3*delta             0.001152      0.071496      0.245690
delta               0.001173      0.069608      0.263153
delta/3             0.001207      0.067230      0.285391
delta/10            0.001174      0.066663      0.279709
delta/100           0.001045      0.066565      0.249784
full_rebal (0)      0.001002      0.066560      0.239648
buy_and_hold        0.001152      0.071496      0.245690

trade statistics (turnover and trading frequency)
strategy        avg_turnover  trade_freq  
100*delta           0.000000      0.000000
10*delta            0.000000      0.000000
3*delta             0.000000      0.000000
delta               0.002460      0.071400
delta/3             0.009072      0.369800
delta/10            0.022665      0.790150
delta/100           0.050683      0.999850
full_rebal (0)      0.060013      1.000000
buy_and_hold        0.000000      0.000000
```
