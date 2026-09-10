@echo off
rem Wrapper batch script to execute microbiome covariate analysis in R
Rscript "%~dp0run_microbiome_covariate_models.R" > "%~dp0rout.txt" 2>&1
echo Exit code: %ERRORLEVEL% >> "%~dp0rout.txt"
