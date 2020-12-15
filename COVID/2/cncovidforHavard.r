require(INLA)

inla.setOption(scale.model.default=TRUE)
inla.setOption(num.threads="8:-8")
inla.setOption(inla.call = INLA:::inla.call.builtin())
inla.setOption(inla.call = "remote")

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

cn.adj <- "cnNew.graph"

# PC priors set up
cont.group <- list(model = "ar1", hyper = list(theta = list(prior = "pc.cor1", param = c(0.9, 0.9), initial = 3)))
hyper.prec = list(prec = list(prior="pc.prec", param=c(U.vcm,0.01)))

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


###############################################################################################
# Update 12 - 09 - 2020
###############################################################################################

# Based on the model run:

# The ones that group shall be removed: pref.id (the intercept, groupRho = 0.979), work.id
# (groupRho = 0.985),

# The ones "group" shall be replace with "replicate": addmig.id (groupRho = 0.101), school.id
# (groupRho = -0.222), tourists.id (groupRho = -0.056), quarantine.id (groupRho = -0.047)

d <- 0.0001

fm.full.update <- fm.full <- nc ~ -1 +
		f(pref.id, model = "besag", graph = cn.adj, group = intercept.day, constr = FALSE, hyper = hyper.prec, diagonal = d) +
		f(travelcity.id, travelcity, model = "besag", graph = cn.adj, group = travelcity.day, constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) +
		f(addmig.id, addmig, model = "besag", graph = cn.adj, replicate = addmig.day, constr = FALSE, hyper = hyper.prec, diagonal = d) +
		f(netmin.id, netmin, model = "besag", graph = cn.adj, group = netmin.day, constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) +
		f(school.id, school, model = "besag", graph = cn.adj, replicate = school.day, constr = FALSE, hyper = hyper.prec, diagonal = d) +
		f(work.id, work, model = "besag", graph = cn.adj, group = work.day, constr = FALSE, hyper = hyper.prec, diagonal = d) +
		f(tourists.id, tourists, model = "besag", graph = cn.adj, replicate = tourists.day, constr = FALSE, hyper = hyper.prec, diagonal = d) +
		f(quarantine.id, quarantine, model = "besag", graph = cn.adj, replicate = quarantine.day, constr = FALSE, hyper = hyper.prec, diagonal = d) +
		f(intercity.id, intercity, model = "besag", graph = cn.adj, group = intercity.day, constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) 

mod.cncovid.intIV.full.update <- inla(fm.full.update,
                                      family="gaussian",
                                      data=cn.inla.data.full, 
                                      control.mode = list(theta = "2.602 -0.293 7.554 0.310 4.278 1.863 0.562 3.760 -0.039 1.130 9.506 9.119 -1.213 0.757 3.531", restart = FALSE),
                                      ##control.mode = list(theta = "2.602 -0.293 7.554 0.310 4.278 1.863 0.562 3.760 -0.039 1.130 9.506 9.119 -1.213 0.757 3.531", fixed = TRUE),
                                      control.predictor=list(compute=TRUE),
                                      control.fixed = list(prec.intercept = 1),
                                      control.compute=list(return.marginals = TRUE, dic=TRUE, cpo=TRUE),
                                      control.inla=list(cmin = 0, int.strategy = "eb"),
                                      blas.num.threads = 1, 
                                      verbose = TRUE)
