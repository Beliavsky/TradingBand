@echo off
setlocal

rem Build with gfortran using the requested diagnostics.
gfortran -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function ^
  -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private ^
  trading_band_mod.f90 xtrading_band.f90 -o xtrading_band_gfortran.exe
if errorlevel 1 goto :end

rem Build with ifx, pointing explicitly at MSVC's linker.
ifx /nologo /traceback /check:bounds /warn:all /warn:unused ^
  /gen-interfaces /warn:interfaces /F512000000 ^
  /Qlocation,link,"C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Tools\MSVC\14.42.34433\bin\Hostx64\x64" ^
  trading_band_mod.f90 xtrading_band.f90 /exe:xtrading_band_ifx.exe

:end
endlocal
