module trading_band_mod
  !!* Core routines for simulating proportional-cost trading bands.
  !  Provides covariance construction, band projection, Monte Carlo simulation,
  !  calibration utilities, formatted printing, and CSV/file helpers.
  use iso_fortran_env, only: real64, int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_quiet_nan, ieee_value
  implicit none

  real(real64), parameter :: zero = 0.0_real64
  real(real64), parameter :: one = 1.0_real64
  real(real64), parameter :: tiny_val = 1.0e-16_real64
  real(real64), parameter :: pi_val = acos(-one)
  real(real64), parameter :: days_year = 252.0_real64
  integer, parameter :: stat_count = 3
  integer, parameter :: trade_stat_count = 2
  integer, parameter :: nstrat = 9
  integer, parameter :: strategy_name_len = 32

  public :: cube_root_band, calibrate_k_for_sharpe, compare_strategies
  public :: print_table, write_metrics_csv, simulate_strategy, days_year
  public :: stat_count, trade_stat_count, nstrat, strategy_name_len

contains

  subroutine seed_rng(seed)
    !!* Seed the intrinsic pseudo-random generator deterministically.
    !  Spreads a single integer seed across the state vector using an LCG.
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    integer(int64) :: state
    integer(int64), parameter :: modulus = 2147483647_int64

    call random_seed(size = n)
    allocate(put(n))

    state = int(seed, int64)
    if (state <= 0_int64) state = 1_int64

    do i = 1, n
      state = mod(48271_int64 * state, modulus)
      put(i) = int(mod(state, modulus))
    end do

    call random_seed(put = put)
    deallocate(put)
  end subroutine seed_rng

  subroutine make_cov(n, sigma, rho, cov)
    !!* Build an equi-correlation covariance matrix.
    !  Diagonal entries equal sigma**2; off-diagonals equal rho * sigma**2.
    integer, intent(in) :: n
    real(real64), intent(in) :: sigma, rho
    real(real64), intent(out) :: cov(n, n)
    integer :: i, j
    real(real64) :: sig2

    sig2 = sigma * sigma
    do i = 1, n
      do j = 1, n
        cov(i, j) = rho * sig2
      end do
      cov(i, i) = sig2
    end do
  end subroutine make_cov

  subroutine project_to_band(weights, delta, weights_new)
    !!* Project weights into a symmetric band around the equally weighted target.
    !  Enforces bounds [1/n - delta, 1/n + delta] and re-normalises to unity.
    real(real64), intent(in) :: weights(:)
    real(real64), intent(in) :: delta
    real(real64), intent(out) :: weights_new(size(weights))
    integer :: n, count_rem
    real(real64) :: tgt, lo, hi, need, cur, s
    logical, allocatable :: over(:), under(:), rem(:)

    n = size(weights)
    tgt = one / real(n, real64)
    if (delta <= zero) then
      weights_new = tgt
      return
    end if

    lo = tgt - delta
    hi = tgt + delta
    weights_new = weights

    allocate(over(n), under(n), rem(n))
    over = weights_new > hi
    under = weights_new < lo

    where (over)
      weights_new = hi
    end where
    where (under)
      weights_new = lo
    end where

    rem = .not.(over .or. under)
    need = one - sum(weights_new, mask = .not. rem)
    cur = sum(weights, mask = rem)
    count_rem = count(rem)

    if (count_rem > 0) then
      if (cur > zero) then
        where (rem)
          weights_new = weights * (need / cur)
        end where
      else
        where (rem)
          weights_new = need / real(count_rem, real64)
        end where
      end if
    end if

    weights_new = max(weights_new, zero)
    s = sum(weights_new)
    if (s <= zero) then
      weights_new = tgt
    else
      weights_new = weights_new / s
    end if

    deallocate(over, under, rem)
  end subroutine project_to_band

  subroutine cholesky_lower(a, L, info)
    !!* Compute the lower-triangular Cholesky factor of a symmetric matrix.
    !  Returns info = 0 on success; otherwise the failing pivot index.
    real(real64), intent(in) :: a(:,:)
    real(real64), intent(out) :: L(size(a,1), size(a,2))
    integer, intent(out) :: info
    integer :: n, i, j, k
    real(real64) :: sum_val

    n = size(a, 1)
    L = zero
    info = 0

    do j = 1, n
      sum_val = a(j, j)
      do k = 1, j - 1
        sum_val = sum_val - L(j, k) * L(j, k)
      end do
      if (sum_val <= zero) then
        info = j
        return
      end if
      L(j, j) = sqrt(sum_val)
      do i = j + 1, n
        sum_val = a(i, j)
        do k = 1, j - 1
          sum_val = sum_val - L(i, k) * L(j, k)
        end do
        L(i, j) = sum_val / L(j, j)
      end do
    end do
  end subroutine cholesky_lower

  subroutine standard_normals(z)
    !!* Generate independent standard normal variates via Box-Muller.
    real(real64), intent(out) :: z(:)
    integer :: n, i
    real(real64) :: u1, u2, radius, angle

    n = size(z)
    i = 1
    do while (i <= n)
      call random_number(u1)
      call random_number(u2)
      if (u1 <= tiny_val) cycle
      radius = sqrt(-2.0_real64 * log(u1))
      angle = 2.0_real64 * pi_val * u2
      z(i) = radius * cos(angle)
      if (i + 1 <= n) then
        z(i + 1) = radius * sin(angle)
      end if
      i = i + 2
    end do
  end subroutine standard_normals

  subroutine mv_normal_draw(mean_vec, chol, out)
    !!* Draw a multivariate normal sample with mean mean_vec and Cholesky factor chol.
    real(real64), intent(in) :: mean_vec(:)
    real(real64), intent(in) :: chol(:,:)
    real(real64), intent(out) :: out(:)
    integer :: n
    real(real64), allocatable :: z(:)

    n = size(mean_vec)
    allocate(z(n))
    call standard_normals(z)
    out = mean_vec + matmul(chol, z)
    deallocate(z)
  end subroutine mv_normal_draw

  real(real64) function mean_val(x)
    !!* Compute the arithmetic mean of a vector.
    real(real64), intent(in) :: x(:)
    if (size(x) <= 0) then
      mean_val = zero
    else
      mean_val = sum(x) / real(size(x), real64)
    end if
  end function mean_val

  real(real64) function sample_std(x)
    !!* Compute the sample standard deviation with Bessel's correction.
    real(real64), intent(in) :: x(:)
    real(real64) :: m, var_acc
    integer :: n

    n = size(x)
    if (n <= 1) then
      sample_std = zero
      return
    end if

    m = mean_val(x)
    var_acc = sum((x - m)**2)
    sample_std = sqrt(max(var_acc / real(n - 1, real64), zero))
  end function sample_std

  real(real64) function nanmean(x)
    !!* Compute the mean while ignoring IEEE NaN values.
    real(real64), intent(in) :: x(:)
    integer :: i, cnt
    real(real64) :: acc

    acc = zero
    cnt = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i))) then
        acc = acc + x(i)
        cnt = cnt + 1
      end if
    end do

    if (cnt > 0) then
      nanmean = acc / real(cnt, real64)
    else
      nanmean = ieee_value(one, ieee_quiet_nan)
    end if
  end function nanmean

  subroutine simulate_strategy(n, periods, sims, mu, sigma, rho, cost_rate, pre, post, trade, &
      band_width, seed, tol)
    !!* Simulate proportional-cost rebalancing with an optional no-trade band.
    !  Returns pre-/post-cost means, vols, and Sharpes plus turnover statistics.
    integer, intent(in) :: n, periods, sims
    real(real64), intent(in) :: mu, sigma, rho, cost_rate
    real(real64), intent(out) :: pre(stat_count), post(stat_count), trade(trade_stat_count)
    real(real64), intent(in), optional :: band_width
    integer, intent(in), optional :: seed
    real(real64), intent(in), optional :: tol

    real(real64) :: tol_val, delta
    real(real64), allocatable :: cov(:,:), chol(:,:), mean_vec(:)
    real(real64), allocatable :: pre_means(:), pre_vols(:), pre_srs(:)
    real(real64), allocatable :: post_means(:), post_vols(:), post_srs(:)
    real(real64), allocatable :: avg_turnovers(:), trade_freqs(:)
    real(real64), allocatable :: pre_rets(:), post_rets(:), per_turn(:), per_flag(:)
    real(real64), allocatable :: holdings(:), weights(:), weights_new(:)
    real(real64), allocatable :: returns_draw(:), holdings_post(:)
    real(real64) :: wealth, wealth_pre, wealth_post, traded_notional, per_turn_val, trade_flag
    real(real64) :: holdings_sum, pre_m, pre_s, post_m, post_s
    integer :: sim, t, info
    logical :: has_band

    if (present(tol)) then
      tol_val = tol
    else
      tol_val = 1.0e-12_real64
    end if

    has_band = present(band_width)
    if (has_band) then
      delta = max(band_width, zero)
    else
      delta = -one
    end if

    if (present(seed)) call seed_rng(seed)

    allocate(cov(n, n), chol(n, n), mean_vec(n))
    call make_cov(n, sigma, rho, cov)
    call cholesky_lower(cov, chol, info)
    if (info /= 0) then
      error stop "cholesky decomposition failed"
    end if
    mean_vec = mu

    allocate(pre_means(sims), pre_vols(sims), pre_srs(sims))
    allocate(post_means(sims), post_vols(sims), post_srs(sims))
    allocate(avg_turnovers(sims), trade_freqs(sims))

    allocate(pre_rets(periods), post_rets(periods), per_turn(periods), per_flag(periods))
    allocate(holdings(n), weights(n), weights_new(n), returns_draw(n), holdings_post(n))

    do sim = 1, sims
      holdings = one / real(n, real64)
      wealth = one

      do t = 1, periods
        call mv_normal_draw(mean_vec, chol, returns_draw)
        holdings = holdings * (one + returns_draw)
        wealth_pre = sum(holdings)
        if (wealth_pre <= tiny_val) wealth_pre = tiny_val
        weights = holdings / wealth_pre

        pre_rets(t) = (wealth_pre - wealth) / wealth

        if (has_band .and. delta >= zero) then
          call project_to_band(weights, delta, weights_new)
        else
          weights_new = weights
        end if

        holdings_post = wealth_pre * weights_new
        traded_notional = sum(abs(holdings_post - holdings))
        per_turn_val = traded_notional / max(wealth_pre, tiny_val)

        if (per_turn_val <= tol_val) then
          per_turn_val = zero
          traded_notional = zero
          trade_flag = zero
        else
          trade_flag = one
        end if
        per_turn(t) = per_turn_val
        per_flag(t) = trade_flag

        wealth_post = sum(holdings_post) - cost_rate * traded_notional
        post_rets(t) = (wealth_post - wealth) / wealth

        holdings_sum = sum(holdings_post)
        if (holdings_sum > zero) then
          holdings = holdings_post * (wealth_post / holdings_sum)
        else
          holdings = holdings_post
        end if
        wealth = wealth_post
      end do

      pre_m = mean_val(pre_rets)
      pre_s = sample_std(pre_rets)
      post_m = mean_val(post_rets)
      post_s = sample_std(post_rets)

      pre_means(sim) = pre_m
      pre_vols(sim) = pre_s
      if (pre_s > zero) then
        pre_srs(sim) = sqrt(days_year) * pre_m / pre_s
      else
        pre_srs(sim) = ieee_value(one, ieee_quiet_nan)
      end if

      post_means(sim) = post_m
      post_vols(sim) = post_s
      if (post_s > zero) then
        post_srs(sim) = sqrt(days_year) * post_m / post_s
      else
        post_srs(sim) = ieee_value(one, ieee_quiet_nan)
      end if

      avg_turnovers(sim) = mean_val(per_turn)
      trade_freqs(sim) = mean_val(per_flag)
    end do

    pre(1) = mean_val(pre_means)
    pre(2) = mean_val(pre_vols)
    pre(3) = nanmean(pre_srs)

    post(1) = mean_val(post_means)
    post(2) = mean_val(post_vols)
    post(3) = nanmean(post_srs)

    trade(1) = mean_val(avg_turnovers)
    trade(2) = mean_val(trade_freqs)

    deallocate(cov, chol, mean_vec, pre_means, pre_vols, pre_srs)
    deallocate(post_means, post_vols, post_srs, avg_turnovers, trade_freqs)
    deallocate(pre_rets, post_rets, per_turn, per_flag)
    deallocate(holdings, weights, weights_new, returns_draw, holdings_post)
  end subroutine simulate_strategy

  real(real64) function cube_root_band(sigma, rho, cost_rate, gamma)
    !!* Compute the cube-root approximation for the optimal band half-width.
    real(real64), intent(in) :: sigma, rho, cost_rate, gamma
    real(real64) :: var_pull

    var_pull = sigma * sigma * max(one - rho, 1.0e-12_real64)
    cube_root_band = (cost_rate / (gamma * var_pull)) ** (one / 3.0_real64)
  end function cube_root_band

  subroutine calibrate_k_for_sharpe(n, periods, sims, mu, sigma, rho, cost_rate, gamma, &
      best_sharpe, best_k, best_delta, base, grid_table, kmin, kmax, ksteps, seed)
    !!* Grid-search the scale factor k for the cube-root band using after-cost Sharpe.
    !  Returns optimal k, its delta, and a grid result table mirroring the Python output.
    integer, intent(in) :: n, periods, sims
    real(real64), intent(in) :: mu, sigma, rho, cost_rate, gamma
    real(real64), intent(out) :: best_sharpe, best_k, best_delta, base
    real(real64), allocatable, intent(out) :: grid_table(:,:)
    real(real64), intent(in), optional :: kmin, kmax
    integer, intent(in), optional :: ksteps
    integer, intent(in), optional :: seed

    integer :: steps, i
    real(real64) :: min_k, max_k, k_val, delta
    real(real64) :: pre_stats(stat_count), post_stats(stat_count)
    real(real64) :: trade_stats(trade_stat_count)

    if (present(kmin)) then
      min_k = kmin
    else
      min_k = 0.2_real64
    end if
    if (present(kmax)) then
      max_k = kmax
    else
      max_k = 2.0_real64
    end if
    if (present(ksteps)) then
      steps = ksteps
    else
      steps = 19
    end if
    if (steps < 1) steps = 1

    allocate(grid_table(steps, 7))
    base = cube_root_band(sigma, rho, cost_rate, gamma)

    best_sharpe = -1.0e30_real64
    best_k = min_k
    best_delta = base * best_k

    do i = 1, steps
      if (steps == 1) then
        k_val = min_k
      else
        k_val = min_k + real(i - 1, real64) * (max_k - min_k) / real(steps - 1, real64)
      end if
      delta = k_val * base
      if (present(seed)) then
        call simulate_strategy(n, periods, sims, mu, sigma, rho, cost_rate, pre_stats, post_stats, &
             trade_stats, band_width = delta, seed = seed)
      else
        call simulate_strategy(n, periods, sims, mu, sigma, rho, cost_rate, pre_stats, post_stats, &
             trade_stats, band_width = delta)
      end if
      grid_table(i, 1) = k_val
      grid_table(i, 2) = delta
      grid_table(i, 3) = post_stats(1)
      grid_table(i, 4) = post_stats(2)
      grid_table(i, 5) = post_stats(3)
      grid_table(i, 6) = trade_stats(1)
      grid_table(i, 7) = trade_stats(2)
      if (post_stats(3) > best_sharpe) then
        best_sharpe = post_stats(3)
        best_k = k_val
        best_delta = delta
      end if
    end do
  end subroutine calibrate_k_for_sharpe

  subroutine compare_strategies(n, periods, sims, mu, sigma, rho, cost_rate, delta_ref, &
      names, pre_tab, post_tab, trade_tab, seed)
    !!* Evaluate a collection of band widths and compile summary statistics.
    integer, intent(in) :: n, periods, sims
    real(real64), intent(in) :: mu, sigma, rho, cost_rate, delta_ref
    character(len=strategy_name_len), intent(out) :: names(:)
    real(real64), intent(out) :: pre_tab(:, :)
    real(real64), intent(out) :: post_tab(:, :)
    real(real64), intent(out) :: trade_tab(:, :)
    integer, intent(in), optional :: seed

    real(real64), dimension(nstrat) :: deltas
    logical, dimension(nstrat) :: use_band
    real(real64) :: pre_stats(stat_count), post_stats(stat_count)
    real(real64) :: trade_stats(trade_stat_count)
    integer :: i

    if (size(names) /= nstrat) error stop "names array wrong size"
    if (size(pre_tab, 1) /= stat_count .or. size(pre_tab, 2) /= nstrat) &
         error stop "pre_tab shape mismatch"
    if (size(post_tab, 1) /= stat_count .or. size(post_tab, 2) /= nstrat) &
         error stop "post_tab shape mismatch"
    if (size(trade_tab, 1) /= trade_stat_count .or. size(trade_tab, 2) /= nstrat) &
         error stop "trade_tab shape mismatch"

    names(1) = "100*delta"
    names(2) = "10*delta"
    names(3) = "3*delta"
    names(4) = "delta"
    names(5) = "delta/3"
    names(6) = "delta/10"
    names(7) = "delta/100"
    names(8) = "full_rebal (0)"
    names(9) = "buy_and_hold"

    deltas = [100.0_real64 * delta_ref, 10.0_real64 * delta_ref, 3.0_real64 * delta_ref, &
      delta_ref, delta_ref / 3.0_real64, 0.1_real64 * delta_ref, &
      0.01_real64 * delta_ref, zero, zero]

    use_band = [.true., .true., .true., .true., .true., .true., .true., .true., .false.]

    do i = 1, nstrat
      if (use_band(i)) then
        if (present(seed)) then
          call simulate_strategy(n, periods, sims, mu, sigma, rho, cost_rate, &
               pre_stats, post_stats, trade_stats, band_width = deltas(i), seed = seed)
        else
          call simulate_strategy(n, periods, sims, mu, sigma, rho, cost_rate, &
               pre_stats, post_stats, trade_stats, band_width = deltas(i))
        end if
      else
        if (present(seed)) then
          call simulate_strategy(n, periods, sims, mu, sigma, rho, cost_rate, &
               pre_stats, post_stats, trade_stats, seed = seed)
        else
          call simulate_strategy(n, periods, sims, mu, sigma, rho, cost_rate, &
               pre_stats, post_stats, trade_stats)
        end if
      end if
      pre_tab(:, i) = pre_stats
      post_tab(:, i) = post_stats
      trade_tab(:, i) = trade_stats
    end do
  end subroutine compare_strategies
  subroutine print_table(names, headers, data)
    !!* Render a metrics table similar to the Python pandas output.
    character(len=*), intent(in) :: names(:)
    character(len=*), intent(in) :: headers(:)
    real(real64), intent(in) :: data(:, :)
    integer :: nrows, ncols, name_width, i, j, header_len
    integer, parameter :: col_width = 12
    character(len=128) :: name_field
    character(len=col_width) :: header_field

    nrows = size(names)
    ncols = size(headers)
    if (size(data, 1) /= ncols .or. size(data, 2) /= nrows) error stop "print_table shape mismatch"

    name_width = len_trim("strategy")
    do i = 1, nrows
      name_width = max(name_width, len_trim(names(i)))
    end do

    name_field = adjustl("strategy")
    if (len_trim(name_field) < name_width) then
      name_field = trim(name_field) // repeat(" ", name_width - len_trim(name_field))
    end if
    write(*,'(A)', advance = 'no') name_field(1:name_width)
    do j = 1, ncols
      header_len = min(len_trim(headers(j)), col_width)
      header_field = repeat(" ", col_width)
      header_field(1:header_len) = headers(j)(1:header_len)
      write(*,'(A)', advance = 'no') "  " // header_field
    end do
    write(*,'(A)') ""

    do i = 1, nrows
      name_field = adjustl(names(i))
      if (len_trim(name_field) < name_width) then
        name_field = trim(name_field) // repeat(" ", name_width - len_trim(name_field))
      end if
      write(*,'(A)', advance = 'no') name_field(1:name_width)
      do j = 1, ncols
        write(*,'(A)', advance = 'no') "  "
        write(*,'(F12.6)', advance = 'no') data(j, i)
      end do
      write(*,'(A)') ""
    end do
  end subroutine print_table

  subroutine write_metrics_csv(filename, names, headers, data)
    !!* Persist tabular strategy metrics to CSV for optional reporting.
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: names(:)
    character(len=*), intent(in) :: headers(:)
    real(real64), intent(in) :: data(:, :)
    integer :: unit, i, j, nrows, ncols
    character(len=32) :: num_buf

    ncols = size(headers)
    nrows = size(names)
    if (size(data, 1) /= ncols .or. size(data, 2) /= nrows) error stop "write_metrics_csv shape mismatch"

    open(newunit = unit, file = filename, status = 'replace', action = 'write')
    write(unit,'(A)', advance = 'no') "strategy"
    do j = 1, ncols
      write(unit,'(A)', advance = 'no') "," // trim(headers(j))
    end do
    write(unit,'(A)') ""

    do i = 1, nrows
      write(unit,'(A)', advance = 'no') trim(names(i))
      do j = 1, ncols
        write(num_buf,'(F12.6)') data(j, i)
        write(unit,'(A)', advance = 'no') "," // trim(adjustl(num_buf))
      end do
      write(unit,'(A)') ""
    end do

    close(unit)
  end subroutine write_metrics_csv

end module trading_band_mod
