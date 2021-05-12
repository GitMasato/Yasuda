@echo off
REM setlocal enabledelayedexpansion

REM "~dp0" path to this bat file
REM "pushed...popd" change the current directly and return to the previous one
pushd %~dp0..\..
  set DEM_EXE=%CD%\build\emgsim.exe
popd

pushd %~dp0
  set INI=%CD%\test_mag.inidem
  set PTCL=%CD%\input\test_ptcl.csv
  set MAG_1=%CD%\input\B_1.magdem
  set MAG_2=%CD%\input\B_2.magdem
  set MAG_3=%CD%\input\B_3.magdem
  set MAG_4=%CD%\input\B_4.magdem
  set OBJ_1=%CD%\input\obj_1.objdem
  set OBJ_2=%CD%\input\obj_2.objdem
  set STDOUT=%CD%\stdout.txt
  set STDERR=%CD%\stderr.txt
popd

call %DEM_EXE% %INI% %PTCL% %MAG_1% %MAG_2% %MAG_3% %MAG_4% %OBJ_1% %OBJ_2% > %STDOUT% 2> %STDERR%
REM call %DEM_EXE% %INI% %PTCL% %MAG_1% %MAG_2% %MAG_3% %MAG_4% %OBJ_1% %OBJ_2% > %STDOUT%
REM call %DEM_EXE% %INI% %PTCL% %MAG_1% %MAG_2% %MAG_3% %MAG_4% %OBJ_1% %OBJ_2% 2> %STDERR%
REM call %DEM_EXE% %INI% %PTCL% %MAG_1% %MAG_2% %MAG_3% %MAG_4% %OBJ_1% %OBJ_2%
