# COVID19-SIR-model
## Disclaimer
This is a toy model that *may* be extended later into something more accurate.
## Description
A very simple SIR model [1] implementation and Levenberg-Marquardt [2,3] fit with the [Johns Hopkins University CSSE data for COVID-19](https://github.com/CSSEGISandData/COVID-19) [4].
Looking at the data for COVID-19 infections in the US reveals a sharp peak after the reported numbers have remained flat for the first 40 days.  

![CSSE data with linear fit](images/Infection_data.png)  

Attempting a linear fit of the data obviously fails.  
A well studied model for predicting infectious diseases is the so-called SIR model. 
It is a set of coupled ODEs for the Susceptible, Infected and Recovered population. A test run yields the following characteristic output:  

![SIR test](images/SIR_model_output.png)

Trying to fit the naive SIR model with constant coefficients to the CSSE data also proves difficult.  

Only a subset of the input data is used for the fit that is still less than ideal. It was not possible to capture the lack of infections in the first 40 days and sharp rise thereafter with a single set of constant coefficients.

![SIR model fit to CSSE data](images/SIR_model_fit.png)

## Requirements
The R script requires the following libraries:
* ggplot2
* deSolve
* lubridate
* FME
* reshape
* minpack.lm

## License
Please see the `LICENSE` file.

## References
[1] [Alemi, A. A., M. Bierbaum, C. R. Myers, and J. P. Sethna (2015): You Can Run, You Can Hide: The Epidemiology and Statistical Mechanics of Zombies.](https://arxiv.org/abs/1503.01104)

[2] Levenberg, K. (1944): A Method for the Solution of Certain Non-Linear Problems in Least Squares. Quarterly of Applied Mathematics. 2 (2): 164–168. doi:10.1090/qam/10666.

[3] Marquardt, D. (1963). An Algorithm for Least-Squares Estimation of Nonlinear Parameters. SIAM Journal on Applied Mathematics. 11 (2): 431–441. doi:10.1137/0111030

[4] https://github.com/CSSEGISandData/COVID-19
