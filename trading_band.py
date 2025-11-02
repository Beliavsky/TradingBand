import numpy as np
import pandas as pd

# symmetric mean-variance portfolio with proportional costs and band rebalancing

def make_cov(n, sigma, rho):
    """equi-correlation covariance."""
    cov = np.full((n, n), rho * sigma**2)
    np.fill_diagonal(cov, sigma**2)
    return cov

def project_to_band(w, delta):
    """project weights to [1/n - delta, 1/n + delta] with sum=1. delta=0 -> target."""
    n = w.size
    tgt = 1.0 / n
    if delta == 0.0:
        return np.full(n, tgt)

    lo, hi = tgt - delta, tgt + delta
    w_new = w.copy()

    over = w_new > hi
    under = w_new < lo
    w_new[over] = hi
    w_new[under] = lo

    rem = ~(over | under)
    need = 1.0 - w_new[~rem].sum()
    cur = w[rem].sum()

    if rem.any():
        if cur > 0:
            w_new[rem] = w[rem] * (need / cur)
        else:
            w_new[rem] = need / rem.sum()

    # guards
    w_new = np.maximum(w_new, 0.0)
    s = w_new.sum()
    if s <= 0:
        return np.full(n, tgt)
    return w_new / s

def simulate_strategy(n, T, sims, mu, sigma, rho, cost_rate, band_width, seed=None, tol=1e-12):
    """
    simulate a band policy.
    returns:
      pre  = (pre_mean, pre_vol, pre_sharpe)
      post = (after_mean, after_vol, after_sharpe)
      trade= (avg_turnover, trade_freq)
    avg_turnover = mean_t(traded_notional_t / W_pre_t); trade_freq = mean_t(1{turnover>0}).
    tol is a relative turnover tolerance to suppress fp artifacts.
    """
    rng = np.random.default_rng(seed) if seed is not None else np.random.default_rng()

    cov = make_cov(n, sigma, rho)
    mu_vec = np.full(n, mu)

    pre_means, pre_vols, pre_srs = [], [], []
    post_means, post_vols, post_srs = [], [], []
    avg_turnovers, trade_freqs = [], []

    for _ in range(sims):
        v = np.ones(n) / n  # holdings; wealth=1
        W = 1.0
        pre_rets, post_rets = [], []
        per_turn_list, per_flag_list = [], []

        for _t in range(T):
            # returns and mark-to-market
            r = rng.multivariate_normal(mu_vec, cov)
            v *= (1.0 + r)
            W_pre = v.sum()
            w = v / W_pre

            # pre-cost return
            pre_rets.append((W_pre - W) / W)

            # decide post-trade weights
            if band_width is None:
                w_new = w  # buy-and-hold
            else:
                w_new = project_to_band(w, band_width)

            v_post = W_pre * w_new

            # traded notional and relative turnover
            traded_notional = float(np.abs(v_post - v).sum())
            per_turn = traded_notional / max(W_pre, 1e-16)

            # tolerance snap for fp noise
            if per_turn <= tol:
                per_turn = 0.0
                traded_notional = 0.0
                fee = 0.0
                traded_flag = 0.0
            else:
                fee = cost_rate * traded_notional
                traded_flag = 1.0

            per_turn_list.append(per_turn)
            per_flag_list.append(traded_flag)

            # wealth after cost
            W_post = v_post.sum() - fee
            post_rets.append((W_post - W) / W)

            # scale positions to reflect fee taken from wealth
            ssum = v_post.sum()
            v = v_post * (W_post / ssum) if ssum > 0.0 else v_post
            W = W_post

        pre_rets = np.asarray(pre_rets)
        post_rets = np.asarray(post_rets)
        per_turn_arr = np.asarray(per_turn_list)
        per_flag_arr = np.asarray(per_flag_list)

        pre_m = pre_rets.mean()
        pre_s = pre_rets.std(ddof=1)
        post_m = post_rets.mean()
        post_s = post_rets.std(ddof=1)

        pre_sr = pre_m / pre_s if pre_s > 0 else np.nan
        post_sr = post_m / post_s if post_s > 0 else np.nan

        pre_means.append(pre_m); pre_vols.append(pre_s); pre_srs.append(pre_sr)
        post_means.append(post_m); post_vols.append(post_s); post_srs.append(post_sr)

        avg_turnovers.append(per_turn_arr.mean())
        trade_freqs.append(per_flag_arr.mean())

    pre = (float(np.mean(pre_means)), float(np.mean(pre_vols)), float(np.nanmean(pre_srs)))
    post = (float(np.mean(post_means)), float(np.mean(post_vols)), float(np.nanmean(post_srs)))
    trade = (float(np.mean(avg_turnovers)), float(np.mean(trade_freqs)))
    return pre, post, trade

def cube_root_band(mu, sigma, rho, cost_rate, gamma):
    """cube-root base: delta_base = (c / (gamma * sigma^2 * (1 - rho)))^(1/3)."""
    var_pull = sigma**2 * max(1.0 - rho, 1e-12)
    base = (cost_rate / (gamma * var_pull)) ** (1.0 / 3.0)
    return base

def calibrate_k_for_sharpe(n, T, sims, mu, sigma, rho, cost_rate, gamma, k_grid=None, seed=None, tol=1e-12):
    """grid search k in delta = k*base maximizing after-cost Sharpe; returns dict with grid df."""
    if k_grid is None:
        k_grid = np.linspace(0.2, 2.0, 19)

    base = cube_root_band(mu, sigma, rho, cost_rate, gamma)

    rows = []
    best = (-np.inf, None, None)
    for k in k_grid:
        delta = k * base
        _, post, trade = simulate_strategy(n, T, sims, mu, sigma, rho, cost_rate, delta, seed=seed, tol=tol)
        post_mean, post_vol, post_sr = post
        avg_turn, trade_freq = trade
        rows.append((k, delta, post_mean, post_vol, post_sr, avg_turn, trade_freq))
        if post_sr > best[0]:
            best = (post_sr, k, delta)

    grid = pd.DataFrame(
        rows,
        columns=["k", "delta", "after_mean", "after_vol", "after_sharpe", "avg_turnover", "trade_freq"]
    )
    return {"best_sharpe": best[0], "k": best[1], "delta": best[2], "base": base, "grid": grid}

def compare_strategies(n, T, sims, mu, sigma, rho, cost_rate, gamma, delta_ref=None, seed=None, tol=1e-12):
    """
    build three dataframes:
      pre_tab:   pre_mean, pre_vol, pre_sharpe
      post_tab:  after_mean, after_vol, after_sharpe
      trade_tab: avg_turnover, trade_freq
    strategies: wide(3*delta_ref), normal(delta_ref), narrow(delta_ref/3), full_rebal(0), buy_and_hold(None)
    """
    if delta_ref is None:
        delta_ref = cube_root_band(mu, sigma, rho, cost_rate, gamma)

    strategies = {
        "100*delta"   : 100.0 * delta_ref,
        "10*delta"    :  10.0 * delta_ref,
        "wide (3*delta)"   : 3.0 * delta_ref,
        "normal (delta)"   : delta_ref,
        "narrow (delta/3)" : delta_ref/3.0,
        "0.1*delta"        : 0.1 * delta_ref,
        "0.01*delta"        : 0.01 * delta_ref,
        "full_rebal (0)"   : 0.0,
        "buy_and_hold"     : None,
    }

    pre_rows, post_rows, trade_rows = [], [], []
    for name, bw in strategies.items():
        pre, post, trade = simulate_strategy(n, T, sims, mu, sigma, rho, cost_rate, bw, seed=seed, tol=tol)
        pre_rows.append((name, *pre))
        post_rows.append((name, *post))
        trade_rows.append((name, *trade))

    pre_tab = pd.DataFrame(pre_rows,  columns=["strategy", "pre_mean", "pre_vol", "pre_sharpe"]).set_index("strategy")
    post_tab = pd.DataFrame(post_rows, columns=["strategy", "after_mean", "after_vol", "after_sharpe"]).set_index("strategy")
    trade_tab = pd.DataFrame(trade_rows, columns=["strategy", "avg_turnover", "trade_freq"]).set_index("strategy")
    return pre_tab, post_tab, trade_tab
