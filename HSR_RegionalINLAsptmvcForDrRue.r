# INLA working directory:

#setwd("/mnt/d/sampledatafolder")

setwd("/mnt/d/OneDrive - mail.montclair.edu/MyDocs/Collaboration/RUC/HSR/AnalysisData")

# Load the necessary libraries:

remove(list=ls())
my.dir <- paste(getwd(),"/",sep="")


require(INLA)
inla.setOption(scale.model.default=FALSE)
inla.setOption(pardiso.license = "/mnt/d/OneDrive - mail.montclair.edu/MyDocs/learningR/RINLA/licenses/pardiso.lic")

require(gstat)
require(geoR)
require(fields)
require(maptools)
require(lattice)
require(spdep)
require(rgdal)
library(plm)
library(sf)

# Load the data:

load("justdata.RData")

# This is how the SOI graph nb is produced:

cnct.soi.nb<-graph2nb(soi.graph(tri2nb(coordinates(cnctsp)), coordinates(cnctsp)))

# And then the graph:

nb2INLA("cnct.graph", cnct.soi.nb)

# Read the graph to a graph object and examine if there are islands:

g <- inla.read.graph('cnct.graph')

# get the ids:
ids <- list()
if (length(g$cc$node) > 1) {
	for (i in 2:length(g$cc$node)) {
		ids[[i]] <- g$cc$node[[i]]
	}
}

#the k-nearest neighbor:
knn <- knearneigh(coordinates(cnctsp), k = 2)

#get the neighbor list for each of the pairs in the graph:
nblist <- knn$nn[unlist(ids),]

#assign the region id to the nblist:
row.names(nblist) <- unlist(ids)

# Alter the graph file based on the nblist's results and get the new graph updated.

cnct.adj <- paste(my.dir, "cnct.graph", sep = "")

d = 0.01
U.vcm <- 0.5/0.31

# Model 1:

fm.travletime.1 <- lnGDP ~ -1  +
			   f(county.id, model = "besag", graph = cnct.adj, group = id.year.int,
				 control.group = list(model = "ar1", hyper = "pc.cor1"), scale.model = TRUE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(urbanization.id, urbanization, model = "besag", graph = cnct.adj, group = id.year.urbanization,
				 control.group = list(model = "ar1", hyper = "pc.cor1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(laborPR.id, laborPR, model = "besag", graph = cnct.adj, group = id.year.laborPR,
			     control.group = list(model = "ar1", hyper = "pc.cor1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnFAI.id, lnFAI, model = "besag", graph = cnct.adj, group = id.year.lnFAI,
			     control.group = list(model = "ar1", hyper = "pc.cor1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnRevn.id, lnRevn, model = "besag", graph = cnct.adj, group = id.year.lnRevn,
			     control.group = list(model = "ar1", hyper = "pc.cor1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnExp.id, lnExp, model = "besag", graph = cnct.adj, group = id.year.lnExp,
			     control.group = list(model = "ar1", hyper = "pc.cor1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnDeposit.id, lnDeposit, model = "besag", graph = cnct.adj, group = id.year.lnDeposit,
			     control.group = list(model = "ar1", hyper = "pc.cor1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnStudSec.id, lnStudSec, model = "besag", graph = cnct.adj, group = id.year.lnStudSec,
			     control.group = list(model = "ar1", hyper = "pc.cor1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnPopDen.id, lnPopDen, model = "besag", graph = cnct.adj, group = id.year.lnPopDen,
			     control.group = list(model = "ar1", hyper = "pc.cor1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(Indstr.id, Indstr, model = "besag", graph = cnct.adj, group = id.year.Indstr,
			     control.group = list(model = "ar1", hyper = "pc.cor1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(TravelTime.id, TravelTime, model = "besag", graph = cnct.adj, group = id.year.TravelTime,
			     control.group = list(model = "ar1", hyper = "pc.cor1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d)


# Only pc.prior for the "besag" models:

fm.travletime.pc.besag <- lnGDP ~ -1  +
			   f(county.id, model = "besag", graph = cnct.adj, group = id.year.int,
				 control.group = list(model = "ar1"), scale.model = TRUE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(urbanization.id, urbanization, model = "besag", graph = cnct.adj, group = id.year.urbanization,
				 control.group = list(model = "ar1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(laborPR.id, laborPR, model = "besag", graph = cnct.adj, group = id.year.laborPR,
			     control.group = list(model = "ar1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnFAI.id, lnFAI, model = "besag", graph = cnct.adj, group = id.year.lnFAI,
			     control.group = list(model = "ar1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnRevn.id, lnRevn, model = "besag", graph = cnct.adj, group = id.year.lnRevn,
			     control.group = list(model = "ar1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnExp.id, lnExp, model = "besag", graph = cnct.adj, group = id.year.lnExp,
			     control.group = list(model = "ar1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnDeposit.id, lnDeposit, model = "besag", graph = cnct.adj, group = id.year.lnDeposit,
			     control.group = list(model = "ar1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnStudSec.id, lnStudSec, model = "besag", graph = cnct.adj, group = id.year.lnStudSec,
			     control.group = list(model = "ar1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(lnPopDen.id, lnPopDen, model = "besag", graph = cnct.adj, group = id.year.lnPopDen,
			     control.group = list(model = "ar1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(Indstr.id, Indstr, model = "besag", graph = cnct.adj, group = id.year.Indstr,
			     control.group = list(model = "ar1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d) +
			   f(TravelTime.id, TravelTime, model = "besag", graph = cnct.adj, group = id.year.TravelTime,
			     control.group = list(model = "ar1"), scale.model=TRUE, constr = FALSE,
				 hyper = list(prec=list(prior="pc.prec",
                              param=c(U.vcm,0.01),
                              initial=0)), diagonal = d)

mod.cnct.intIV.1 <- inla(fm.travletime.1,family="gaussian",data=cnct.inla.data, 
                  control.predictor=list(compute=TRUE),
				  control.fixed = list(prec.intercept = 1),
                  control.compute=list(dic=TRUE,cpo=TRUE),
				  control.inla=list(strategy="adaptive", int.strategy ="eb", cmin = 0), num.threads="22:1", verbose = FALSE)
				  
mod.cnct.intIV.pc.besag <- inla(fm.travletime.pc.besag,family="gaussian",data=cnct.inla.data, 
                  control.predictor=list(compute=TRUE),
				  control.fixed = list(prec.intercept = 1),
                  control.compute=list(dic=TRUE,cpo=TRUE),
				  control.inla=list(strategy="adaptive", int.strategy ="eb", cmin = 0), num.threads="8:1", verbose = FALSE)
			
