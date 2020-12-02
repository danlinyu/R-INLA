require(INLA)

require(gstat)
require(geoR)
require(fields)
require(maptools)
require(lattice)
require(spdep)
require(rgdal)
library(plm)
library(sf)

load("sampledataforDrRue.RData")
cnct.adj <- "cnctNEW.graph"

d = 0
U.vcm <- 0.1/0.31
cont.group <- list(model = "ar1", hyper = list(theta = list(prior = "pc.cor1", param = c(0.9, 0.9), initial = 3)))
hyper.prec = list(prec = list(prior="pc.prec", param=c(U.vcm,0.01)))

inla.setOption(scale.model.default=TRUE)

## no time-effect: county.id (the intercept), lnFAI.id, lnRevn.id, lnExp.id, lnDeposit.id,
## lnStudSec.id, lnPopDen.id

fm.travletime <- lnGDP ~ -1  +
    f(county.id, model = "besag", graph = cnct.adj, group = id.year.int,
      constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) +
    f(urbanization.id, urbanization, model = "besag", graph = cnct.adj, group = id.year.urbanization,
      constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) +
    f(laborPR.id, laborPR, model = "besag", graph = cnct.adj, group = id.year.laborPR,
      constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) +
    f(Indstr.id, Indstr, model = "besag", graph = cnct.adj, group = id.year.Indstr,
      constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) +
    f(TravelTime.id, TravelTime, model = "besag", graph = cnct.adj, group = id.year.TravelTime,
      constr = FALSE, control.group = cont.group, hyper = hyper.prec, diagonal = d) +
    f(lnFAI.id, lnFAI, model = "besag", graph = cnct.adj, constr = FALSE, hyper = hyper.prec, diagonal = d) +
    f(lnRevn.id, lnRevn, model = "besag", graph = cnct.adj, constr = FALSE, hyper = hyper.prec, diagonal = d) +
    f(lnExp.id, lnExp, model = "besag", graph = cnct.adj, constr = FALSE, hyper = hyper.prec, diagonal = d) +
    f(lnDeposit.id, lnDeposit, model = "besag", graph = cnct.adj, constr = FALSE, hyper = hyper.prec, diagonal = d) +
    f(lnStudSec.id, lnStudSec, model = "besag", graph = cnct.adj, constr = FALSE, hyper = hyper.prec, diagonal = d) +
    f(lnPopDen.id, lnPopDen, model = "besag", graph = cnct.adj, constr = FALSE, hyper = hyper.prec, diagonal = d) 

mod.cnct.intIV <- inla(fm.travletime,family="gaussian",
                       data=cnct.inla.data, 
                       control.predictor=list(compute=TRUE),
                       control.fixed = list(prec.intercept = 1),
                       control.compute=list(return.marginals = FALSE, dic=FALSE, cpo=FALSE),
                       control.inla=list(cmin = 0, int.strategy = "eb"),
                       ## you'll need to change this one
                       num.threads="18:2", 
                       blas.num.threads = 1, 
                       verbose = TRUE)

