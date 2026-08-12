@echo off
setlocal
if not defined R_SCRIPT (
  where Rscript >nul 2>&1
  if errorlevel 1 (
    set "R_SCRIPT=C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
  ) else (
    set "R_SCRIPT=Rscript"
  )
)
set "RENV_CONFIG_AUTOLOADER_ENABLED=FALSE"
"%R_SCRIPT%" "%~dp0tests\run_example_test.R"
if errorlevel 1 exit /b %errorlevel%
type "%~dp0test_outputs\TEST_PASS.txt"
