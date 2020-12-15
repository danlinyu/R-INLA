# Workspace:

##setwd("/mnt/d/OneDrive - mail.montclair.edu/MyDocs/Collaboration/RUC/COVID19/model")
##remove(list=ls())
my.dir <- paste(getwd(),"/",sep="")

# Load necessary packages:

require(INLA)

##inla.setOption(pardiso.license = "/mnt/d/OneDrive - mail.montclair.edu/MyDocs/learningR/RINLA/licenses/pardiso.lic")

require(gstat)
require(geoR)
require(fields)
require(maptools)
require(lattice)
require(spdep)
require(rgdal)
library(plm)

load('cncovidinladata.RData')
cn.adj <- "cnNew.graph"

###########################################################################################################################
## Data Explanation:
## cncovid19.alldata includes all the data of China's COVID19 daily new cases and relevant variables:
## Only relevant variables are explained here:
##
## 1. Daily varying variables (social distancing measures that are hypothesized to have curbing effects on the daily new cases:
##
##		new.cases: new cases for that day.
##		travelcity: an index for inner city travel intensity (continuous)
##		school: status of school opening: 1 is closed, 0 is open.
##		work: status of workplaces: 1 is close, 0 is open.
##		tourists: a mined data from "yelp" like app in China, it is simply a count of clicks on popular tourists' sites
##		quarantine: status of whether the city is strictly quarantined: 1 yes, 0 no.
##		intercity: status of inter-city transportation: 1 is shutting down, 0 no.
##		addmig: daily migration index mined from China's search engine, Baidu (kind of the Chinese version of Google)
##		netmin: daily net migrants into the city, mined from Baidu.
##
## 2. Time invariant variables (background information):
##
##		tpop: total population, in 10,000 people
##		percapitalGDP: per capita GDP, measured in Yuan
##		secondaryInd: percentage of secondary industry in GDP
##		tertiaryind: percentage of tertiary industry in GDP
##		finincpc: local financial income per capita
##		finexppc: local financial expenses per capita
##		hospper10000: number of hospitals per 10,000 people in that city
##		hsbedper10000: number of hosptal beds per 10,000 people in that city
##		docper10000: number of doctors per 10,000 people in that city
##		jiaotongmi: road network density (total length of roads in that city divided by the city's land area)
##		TravelTime: network travel time from the center of the city to the nearest high speed railway station.
##
#####################################################################################################################################
##
## There are two ending dates for building up the data sets:
## The first is from 1-25-2020 (that's 14 days after the first date with available data, to honor the 14 day quarantine policy) to 3-5-2020
## This is when the curve of COVID19 in China started to flatten
## The second is from 1-25-2020 to 05-04-2020

## The model with Inversed Normal Transformed variables (both dependent and independent
## variables are regressed on time invariant variables, and residuals taken as new
## dependent/independent variables. This model is again a "monster" model that takes forever to
## run.

d <- 0.01
cont.group <- list(model = "ar1", hyper = list(theta = list(prior = "pc.cor1", param = c(0.9, 0.9), initial = 3)))
hyper.prec <- list(prec = list(prior = "pc.prec", param = c(1.6129, 0.01)))

resolution <- 1L
if (resolution > 1) {
    cn.inla.data.full$intercept.day <- cn.inla.data.full$intercept.day %/% resolution + 1
    cn.inla.data.full$travelcity.day <- cn.inla.data.full$travelcity.day %/% resolution + 1
    cn.inla.data.full$school.day <- cn.inla.data.full$school.day %/% resolution + 1
    cn.inla.data.full$netmin.day <- cn.inla.data.full$netmin.day %/% resolution + 1
    cn.inla.data.full$work.day <- cn.inla.data.full$work.day %/% resolution + 1
    cn.inla.data.full$tourists.day <- cn.inla.data.full$tourists.day %/% resolution + 1
    cn.inla.data.full$quarantine.day <- cn.inla.data.full$quarantine.day %/% resolution + 1
    cn.inla.data.full$intercity.day <- cn.inla.data.full$intercity.day %/% resolution + 1
    cn.inla.data.full$addmig.day <- cn.inla.data.full$addmig.day %/% resolution + 1
}

fm.full <- nc ~ -1 +
    f(pref.id, model = "besag", graph = cn.adj, group = intercept.day, 
      constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) +
    f(travelcity.id, travelcity, model = "besag", graph = cn.adj, group = travelcity.day,
      constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) +
    f(addmig.id, addmig, model = "besag", graph = cn.adj, group = addmig.day, constr = FALSE, 
      control.group = cont.group, hyper = hyper.prec, diagonal = d) + 
    f(netmin.id, netmin, model = "besag", graph = cn.adj, group = netmin.day, 
      constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) +
    f(school.id, school, model = "besag", graph = cn.adj, group = school.day, constr = FALSE,
      control.group = cont.group, hyper = hyper.prec, diagonal = d) +
    f(work.id, work, model = "besag", graph = cn.adj, group = work.day, constr = FALSE,
      control.group = cont.group, hyper = hyper.prec, diagonal = d) +
    f(tourists.id, tourists, model = "besag", graph = cn.adj, group = tourists.day, constr = FALSE, 
      control.group = cont.group, hyper = hyper.prec, diagonal = d) + 
    f(quarantine.id, quarantine, model = "besag", graph = cn.adj, group = quarantine.day,
      constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) +
    f(intercity.id, intercity, model = "besag", graph = cn.adj, group = intercity.day, constr = FALSE, 
      control.group = cont.group, hyper = hyper.prec, diagonal = d)

mod.cncovid.intIV.full <- inla(
    fm.full,family="gaussian",
    data=cn.inla.data.full, 
    control.predictor=list(compute=TRUE),
    control.fixed = list(prec.intercept = 1),
    control.compute=list(return.marginals = TRUE, dic=TRUE, cpo=TRUE),
    control.inla=list(cmin = 0, int.strategy = "eb"),
    control.mode = list(
        ##theta = "1.456 -0.602 4.613 2.189 3.726 1.328 2.134 2.736 2.799 3.810 2.739 3.056 3.171 9.498 7.202 1.175 0.897 3.053 2.465",
        theta = "2.364 -0.750 4.357 6.712 4.837 1.622 0.316 6.915 3.544 1.407 0.324 5.399 4.053 13.625 6.573 -1.190 -0.192 5.489 1.524",
        restart = TRUE), 
    ## you'll need to change this one
    num.threads="10:6", 
    blas.num.threads = 1, 
    inla.call = "remote", 
    verbose = TRUE)
					   
## If we are to model the relationships between social distancing measures and the daily new
## cases of COVID19, without going through the Inversed Normal Transformation, I guess we are
## going to model with the "possion" family, but how will we include the time-invariant
## variables to the model (the background information)? Please advise.
