import argparse
import numpy as np
import pandas as pd

from trading_band import (
    cube_root_band,
    calibrate_k_for_sharpe,
    compare_strategies,
)

def main():
    p = argparse.ArgumentParser(description="simulate no-trade band strategies with proportional costs")
    p.add_argument("--n", type=int, default=5, help="number of assets")
    p.add_argument("--T", type=int, default=200, help="periods per simulation path")
    p.add_argument("--sims", type=int, default=100, help="number of monte carlo paths")
    p.add_argument("--mu", type=float, default=0.001, help="per-period mean return")
    p.add_argument("--sigma", type=float, default=0.1, help="per-period volatility")
    p.add_argument("--rho", type=float, default=0.30, help="pairwise correlation")
    p.add_argument("--cost", type=float, default=0.005, help="proportional transaction cost c")
    p.add_argument("--gamma", type=float, default=3.0, help="risk aversion for mv objective")
    p.add_argument("--seed", type=int, default=7, help="random seed")
    p.add_argument("--calibrate_k", action="store_true", help="calibrate k for the cube-root band to maximize after-cost sharpe")
    p.add_argument("--kmin", type=float, default=0.4, help="min k in calibration grid")
    p.add_argument("--kmax", type=float, default=1.8, help="max k in calibration grid")
    p.add_argument("--ksteps", type=int, default=15, help="number of grid steps for k")
    p.add_argument("--save_csv_prefix", type=str, default="", help="if non-empty, save tables to <prefix>_*.csv")
    args = p.parse_args()

    n = args.n
    T = args.T
    sims = args.sims
    mu = args.mu
    sigma = args.sigma
    rho = args.rho
    c = args.cost
    gamma = args.gamma
    seed = args.seed

    # 1) cube-root base band
    delta_base = cube_root_band(mu, sigma, rho, c, gamma)

    print("parameters")
    print(f"n={n}, T={T}, sims={sims}, mu={mu:.6f}, sigma={sigma:.6f}, rho={rho:.6f}, c={c:.6f}, gamma={gamma:.6f}, seed={seed}")
    print(f"delta_base (cube-root) = {delta_base:.8f}")

    # 2) optionally calibrate k for after-cost sharpe
    if args.calibrate_k:
        k_grid = np.linspace(args.kmin, args.kmax, args.ksteps)
        out = calibrate_k_for_sharpe(n, T, sims, mu, sigma, rho, c, gamma, k_grid=k_grid, seed=seed)
        print("\nsharpe-oriented calibration")
        print(f"best k = {out['k']:.6f}")
        print(f"best delta = {out['delta']:.8f}")
        print(f"best after-cost sharpe = {out['best_sharpe']:.6f}")
        if args.save_csv_prefix:
            out["grid"].to_csv(f"{args.save_csv_prefix}_k_grid.csv", index=False)
        delta_ref = out["delta"]
    else:
        delta_ref = delta_base

    # 3) compare strategies: pre-/post-cost plus turnover stats
    pre_tab, post_tab, trade_tab = compare_strategies(
        n, T, sims, mu, sigma, rho, c, gamma, delta_ref=delta_ref, seed=seed
    )

    # formatting
    def fmt(df):
        z = df.copy()
        for col in z.columns:
            z[col] = z[col].astype(float)
        return z

    pre_fmt = fmt(pre_tab)
    post_fmt = fmt(post_tab)
    trade_fmt = fmt(trade_tab)

    print("\npre-cost results")
    print(pre_fmt.to_string(float_format=lambda x: f"{x:.6f}"))

    print("\nafter-cost results")
    print(post_fmt.to_string(float_format=lambda x: f"{x:.6f}"))

    print("\ntrade statistics (turnover and trading frequency)")
    print(trade_fmt.to_string(float_format=lambda x: f"{x:.6f}"))

    if args.save_csv_prefix:
        pre_fmt.to_csv(f"{args.save_csv_prefix}_pre_cost.csv")
        post_fmt.to_csv(f"{args.save_csv_prefix}_after_cost.csv")
        trade_fmt.to_csv(f"{args.save_csv_prefix}_trade_stats.csv")
        print(f"\nsaved: {args.save_csv_prefix}_pre_cost.csv, {args.save_csv_prefix}_after_cost.csv, {args.save_csv_prefix}_trade_stats.csv")

if __name__ == "__main__":
    main()
