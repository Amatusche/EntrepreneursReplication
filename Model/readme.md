# Model replication

The model code is written in Julia. Install the packages imported by the scripts before running them. The median-return figure also requires `Data/out/medianreturnNWdec.txt` from the empirical analysis.

1. Set `leverage_mode = :none` in `Model/src/options.jl` to reproduce the baseline model. Leave `calibrate = false` to use the supplied calibrated parameters.
2. From the replication package's root directory, run:

   ```sh
   cd Model/src
   julia main.jl
   julia TABLES.jl
   julia PLOTS.jl
   ```

   Alternatively, run `julia main_parallel.jl` instead of `julia main.jl` to use parallel processes.

3. For the leverage extension, set `leverage_mode = :riskfree` in `options.jl` and repeat the three Julia commands. The current default is `:riskfree`.

Results are saved in `Model/out/stickyprices/no_leverage/` for the baseline and `Model/out/stickyprices/riskfree_leverage/` for the leverage extension:

- Model solutions and statistics: `.jld2` files in the result directory.
- Paper tables: `.tex` files in `Tables/`.
- Paper figures: `.tex` files in `Figures/`.
- Diagnostic plots: `Diagnostics/`.

Table and figure filenames follow the numbering in the JME submission, for example `Table_2_CalibrationExternal.tex`, `Figure_4_OutputDec_all3.tex`, and `Figure_C.2_OutputDec_all3.tex`. The two panels of Figure 8 are saved separately as `Figure_8a_BeqInvestmentResponse.tex` and `Figure_8b_BeqInvestmentResponse_direct.tex`. Both model versions export the full set in their respective result directories. Leverage outputs use the baseline filenames, except for results shown in the leverage appendix, which use their appendix numbers.

Keep the working directory at `Model/src` when running the scripts, and finish each model run before generating its tables and figures.
