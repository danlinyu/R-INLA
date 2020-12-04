# Workspace:

setwd("/mnt/d/OneDrive - mail.montclair.edu/MyDocs/Collaboration/RUC/COVID19/model")

remove(list=ls())
my.dir <- paste(getwd(),"/",sep="")

# Load necessary packages:

require(INLA)

inla.setOption(pardiso.license = "/mnt/d/OneDrive - mail.montclair.edu/MyDocs/learningR/RINLA/licenses/pardiso.lic")

require(gstat)
require(geoR)
require(fields)
require(maptools)
require(lattice)
require(spdep)
require(rgdal)
library(plm)

# Load the data:

load('cncovidinladata.RData')

cn.adj <- paste(my.dir, "cnNew.graph", sep = "")

###########################################################################################################################
# Data Explanation:
# cncovid19.alldata includes all the data of China's COVID19 daily new cases and relevant variables:
# Only relevant variables are explained here:
# 1. Daily varying variables (social distancing measures that are hypothesized to have curbing effects on the daily new cases:
#
#		new.cases: new cases for that day.
#		travelcity: an index for inner city travel intensity (continuous)
#		school: status of school opening: 1 is closed, 0 is open.
#		work: status of workplaces: 1 is close, 0 is open.
#		tourists: a mined data from "yelp" like app in China, it is simply a count of clicks on popular tourists' sites
#		quarantine: status of whether the city is strictly quarantined: 1 yes, 0 no.
#		intercity: status of inter-city transportation: 1 is shutting down, 0 no.
#		addmig: daily migration index mined from China's search engine, Baidu (kind of the Chinese version of Google)
#		netmin: daily net migrants into the city, mined from Baidu.
#
# 2. Time invariant variables (background information):
#
#		tpop: total population, in 10,000 people
#		percapitalGDP: per capita GDP, measured in Yuan
#		secondaryInd: percentage of secondary industry in GDP
#		tertiaryind: percentage of tertiary industry in GDP
#		finincpc: local financial income per capita
#		finexppc: local financial expenses per capita
#		hospper10000: number of hospitals per 10,000 people in that city
#		hsbedper10000: number of hosptal beds per 10,000 people in that city
#		docper10000: number of doctors per 10,000 people in that city
#		jiaotongmi: road network density (total length of roads in that city divided by the city's land area)
#		TravelTime: network travel time from the center of the city to the nearest high speed railway station.
#
####################################################################################################################################
#
# There are two ending dates for building up the data sets:
# The first is from 1-25-2020 (that's 14 days after the first date with available data, to honor the 14 day quarantine policy) to 3-5-2020
# This is when the curve of COVID19 in China started to flatten
# The second is from 1-25-2020 to 05-04-2020

# The model with Inversed Normal Transformed variables (both dependent and independent variables are regressed on time invariant variables, and residuals taken as new
# dependent/independent variables. This model is again a "monster" model that takes forever to run.

mod.cncovid.intIV.full <- inla(fm.full,family="gaussian",
                       data=cn.inla.data.full, 
                       control.predictor=list(compute=TRUE),
                       control.fixed = list(prec.intercept = 1),
                       control.compute=list(return.marginals = TRUE, dic=TRUE, cpo=TRUE),
                       control.inla=list(cmin = 0, int.strategy = "eb"),
                       ## you'll need to change this one
                       num.threads="18:2", 
                       blas.num.threads = 1, 
                       verbose = FALSE)
					   
# If we are to model the relationships between social distancing measures and the daily new cases of COVID19, without going through the Inversed Normal Transformation,
# I guess we are going to model with the "possion" family, but how will we include the time-invariant variables to the model (the background information)? Please advise.
#