#
# para_estim.r
#
# Description:
# ============
# Investigation of the Johns Hopkins University CSSE COVID-19 data for the US.
# A SIR model is defined and a Levenberg-Marquardt algorithm is used to fit the
# SIR parameters 'b' and 'nu' to the CSSE data.
#
#
# Copyright © 2020 Patrick Stegmann
#
# This file is part of COVID19-SIR-model.
#
# COVID19-SIR-model is free software:
# you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 

# Read CSV COVID-19 data into R:
InfData <- read.csv("./data/time_series_covid19_confirmed_global.csv",header = F,sep = ",")
RecData <- read.csv("./data/time_series_covid19_recovered_global.csv",header = F,sep = ",")
summary(InfData)
summary(RecData)

# Convert data into data frame
desc <- c('time','Infected')
datum <- c(as.Date(t(InfData[1,5:70]), "%m/%d/%Y"))
summary(datum)
infected <- as.double(t(InfData[227,5:70]))
recovered <- as.double(t(RecData[227,5:70]))
dfi <- data.frame(datum,infected,stringsAsFactors=FALSE)
dfr <- data.frame(datum,recovered,stringsAsFactors=FALSE)
names(dfi) = c("time","infected")

# Plot the infection data:
plot(dfi)
# Attempt a linear fit of the data:
fit <-lm(infected~datum,data=dfi)
# Realize that the linear fit fails:
abline(fit)
plot(dfr)

# Make a pretty plot of the data
library(ggplot2)
ggplot(dfi,aes(datum,infected)) + geom_point()
ggplot(dfr,aes(datum,recovered)) + geom_point()

# Set up the SIR ODE model:
library(deSolve)

# SIR model ODE
SIR <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -b*S*I
    dI <- b*S*I - nu*I
    dR <- nu*I
    list(c(dS, dI, dR))
  })
}

# Set up model parameters and initial conditions:
library(lubridate)

tiem <- ymd(datum)
tiem <- yday(tiem) - 1.0
parameters <- c(b = 0.00001, nu = 0.1)
state      <- c(S = 100000, I = 1., R = 0.0)
times      <- seq(tiem[1], tiem[66], by = 0.2)
# Do a first test calculation with the SIR model and check the results:
out <- ode(y = state, times = times, func = SIR, parms = parameters)
head(out)
plot(out)

# Define the residual for the Levenberg-Marquardt minimization:
library(FME)
library(reshape2)
ssq = function(parms){
  # initial population
  state      <- c(S = 100000, I = 1., R = 0.0)
  times      <- seq(1, 26, by = 1.0)
  b = parms[1]
  nu = parms[2]
  out <- ode(y = state, times = times, func = SIR, parms = list(b=b,nu=nu))
  
  # Filter data that contains time points where data are available
  outdf=data.frame(out)
  #outdf=outdf[outdf$time %in% seq(1, 66, by = 1.0),]
  
  # Evaluate predicted vs. experimental residual
  preddf=melt(outdf,id.var="time",value.name="I")
  expdfi=melt(dfi,id.var="time",value.name="infected")
  ssqres=sqrt((outdf$I - as.double(dfi$infected[40:66]) )**2 + ( outdf$R - as.double(dfr$recovered[40:66]) )**2 )
  # return predicted vs experimental residual
  return(ssqres)
}

# parameter fitting using Levenberg-Marquart algorithm
# initial guess for parameters
parms=c(b = 0.000006, nu = 0.05)
# actual fitting
library(minpack.lm) # library for least squares fit using Levenberg-Marquart algorithm
fitval=nls.lm(par=parms,fn=ssq)
# Check the output of the fit:
summary(fitval)
parest=as.list(coef(fitval))
parest

# degrees of freedom: # data points - # parameters
dof=3*nrow(df)-2
dof

# mean error
ms=sqrt(deviance(fitval)/dof)
ms

# variance Covariance Matrix
S=vcov(fitval)
S

# Compare the data and the ODE solution with optimized parameters:
times <- seq(1, 26, by = 1.0)
out <- ode(y = state, times = times, func = SIR, parms = parest)
outdf = data.frame(out)
gra = ggplot(data = dfi, aes(x=seq(1, 66, by = 1.0),y=dfi$infected)) + geom_point()
gra = gra + geom_line(data = outdf, aes(x=times+40,y=I))
gra = gra + labs(title="SIR model COVID-19 prediction") + xlab("Time [a]") + ylab("Number of COVID-19 infections [1 human]")
print(gra)

#End of Script.
