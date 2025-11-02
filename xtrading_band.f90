program xtrading_band
  !!* Driver for the proportional-cost trading band simulation.
  !  Mirrors the Python workflow: runs the core simulation, optional calibration, and reporting.
  use iso_fortran_env, only: real64
  use trading_band_mod, only: cube_root_band, calibrate_k_for_sharpe, compare_strategies
  use trading_band_mod, only: print_table, write_metrics_csv, days_year
  use trading_band_mod, only: stat_count, trade_stat_count, nstrat, strategy_name_len
  implicit none

  integer, parameter :: asset_count = 5
  integer, parameter :: simulation_periods = 200
  integer, parameter :: path_count = 100
  integer, parameter :: rng_seed = 7
  integer, parameter :: k_steps = 15
  real(real64), parameter :: mean_return = 0.001_real64
  real(real64), parameter :: volatility = 0.1_real64
  real(real64), parameter :: correlation = 0.30_real64
  real(real64), parameter :: transaction_cost = 0.005_real64
  real(real64), parameter :: risk_aversion = 3.0_real64
  real(real64), parameter :: k_min = 0.4_real64
  real(real64), parameter :: k_max = 1.8_real64
  character(len=*), parameter :: output_prefix = ""
  logical, parameter :: do_calibrate = .true.

  real(real64) :: delta_base, delta_ref
  real(real64) :: best_sharpe, best_k, best_delta
  real(real64), allocatable :: grid_table(:,:)
  real(real64) :: pre_tab(stat_count, nstrat)
  real(real64) :: post_tab(stat_count, nstrat)
  real(real64) :: trade_tab(trade_stat_count, nstrat)
  character(len=strategy_name_len) :: strategy_names(nstrat)
  character(len=16) :: pre_headers(3), post_headers(3)
  character(len=20) :: trade_headers(2)
  character(len=256) :: param_line
  character(len=*), parameter :: param_fmt = "('n=',I0,', T=',I0,', sims=',I0,', mu=',F9.6,', sigma=',F9.6,', rho=',F9.6,', c=',F9.6,', gamma=',F9.6,', seed=',I0)"
  logical :: have_prefix

  pre_headers = [character(len=16) :: "pre_mean", "pre_vol", "pre_sharpe"]
  post_headers = [character(len=16) :: "after_mean", "after_vol", "after_sharpe"]
  trade_headers = [character(len=20) :: "avg_turnover", "trade_freq"]

  delta_base = cube_root_band(volatility, correlation, transaction_cost, risk_aversion)

  write(*,"(A)") "parameters"
  write(*,*) "days_year", days_year
  write(param_line, param_fmt) asset_count, simulation_periods, path_count, &
       mean_return, volatility, correlation, transaction_cost, risk_aversion, rng_seed
  write(*,"(A)") trim(param_line)
  write(*,"(A,F12.8)") "delta_base (cube-root) = ", delta_base

  if (do_calibrate) then
    write(*,"(A)") ""
    write(*,"(A)") "sharpe-oriented calibration"
    call calibrate_k_for_sharpe(asset_count, simulation_periods, path_count, mean_return, &
         volatility, correlation, transaction_cost, risk_aversion, best_sharpe, best_k, &
         best_delta, delta_base, grid_table, kmin = k_min, kmax = k_max, ksteps = k_steps, &
         seed = rng_seed)
    write(*,"(A,F9.6)") "best k = ", best_k
    write(*,"(A,F12.8)") "best delta = ", best_delta
    write(*,"(A,F9.6)") "best after-cost sharpe = ", best_sharpe
    delta_ref = best_delta
  else
    delta_ref = delta_base
  end if

  call compare_strategies(asset_count, simulation_periods, path_count, mean_return, volatility, &
       correlation, transaction_cost, delta_ref, strategy_names, pre_tab, post_tab, trade_tab, &
       seed = rng_seed)

  write(*,"(A)") ""
  write(*,"(A)") "pre-cost results"
  call print_table(strategy_names, pre_headers, pre_tab)

  write(*,"(A)") ""
  write(*,"(A)") "after-cost results"
  call print_table(strategy_names, post_headers, post_tab)

  write(*,"(A)") ""
  write(*,"(A)") "trade statistics (turnover and trading frequency)"
  call print_table(strategy_names, trade_headers, trade_tab)

  have_prefix = len_trim(output_prefix) > 0
  if (have_prefix) then
    call write_metrics_csv(trim(output_prefix)//"_pre_cost.csv", strategy_names, pre_headers, &
         pre_tab)
    call write_metrics_csv(trim(output_prefix)//"_after_cost.csv", strategy_names, post_headers, &
         post_tab)
    call write_metrics_csv(trim(output_prefix)//"_trade_stats.csv", strategy_names, &
         trade_headers, trade_tab)
    if (do_calibrate) then
      call write_k_grid_csv(trim(output_prefix)//"_k_grid.csv", grid_table)
    end if
    if (do_calibrate) then
      write(*,"(A)") ""
      write(*,"(A)") "saved: "//trim(output_prefix)//"_pre_cost.csv, "//trim(output_prefix)//"_after_cost.csv, "//trim(output_prefix)//"_trade_stats.csv, "//trim(output_prefix)//"_k_grid.csv"
    else
      write(*,"(A)") ""
      write(*,"(A)") "saved: "//trim(output_prefix)//"_pre_cost.csv, "//trim(output_prefix)//"_after_cost.csv, "//trim(output_prefix)//"_trade_stats.csv"
    end if
  end if

  if (allocated(grid_table)) deallocate(grid_table)

contains

  subroutine write_k_grid_csv(filename, table)
    !!* Write the calibration grid (k, delta, metrics) to CSV if requested.
    character(len=*), intent(in) :: filename
    real(real64), intent(in) :: table(:,:)
    integer :: unit, i

    open(newunit = unit, file = filename, status = "replace", action = "write")
    write(unit,"(A)") "k,delta,after_mean,after_vol,after_sharpe,avg_turnover,trade_freq"
    do i = 1, size(table, 1)
      write(unit,"(F9.6,',',F12.8,*(:,',',F12.6))") table(i, :)
    end do
    close(unit)
  end subroutine write_k_grid_csv

end program xtrading_band
