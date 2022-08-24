########################################################################
#### Evolvability of 3-dimensional floral morphology in Dalechampia ####
########################################################################

rm(list=ls())

# Load packages
library(MCMCglmm)
library(evolvability)
library(reshape2)
library(plyr)

nms = c("GAD", "GSD", "ASD") # Focal traits

# 1. Descriptive plots ####
field = read.csv("data/field.csv")
field = na.omit(field)

# Compute species and population means
spmeans = ddply(field, .(species), summarize,
                n = length(species),
                npop = length(unique(pop)),
                mGSD = round(mean(GSD, na.rm=T), 2),
                GSDse = round(sd(GSD, na.rm=T)/sqrt(sum(!is.na(GSD))), 2),
                mGAD = round(mean(GAD, na.rm=T), 2),
                GADse = round(sd(GAD, na.rm=T)/sqrt(sum(!is.na(GAD))), 2),
                mASD = round(mean(ASD, na.rm=T), 2),
                ASDse = round(sd(ASD, na.rm=T)/sqrt(sum(!is.na(ASD))), 2))

spmeans

popmeans = ddply(field, .(species, pop), summarize,
                 n = length(pop),
                 GSD = mean(GSD, na.rm=T),
                 GAD = mean(GAD, na.rm=T),
                 ASD = mean(ASD, na.rm=T))

popmeans[order(popmeans$n),]
drop = which(popmeans$n<3) # Drop populations with very small sample size
popmeans = popmeans[-drop,]

table(popmeans$species) # Pops per species

scameans = popmeans[popmeans$species=="S", ]

# Figure 2
x11(height=5, width=9)

#cairo_pdf("Fig2.pdf", height=5, width=9, family = "Times")

par(mfrow=c(1,2))

plot(popmeans$GSD, popmeans$GAD, cex=sqrt(popmeans$n)/2, las=1,
     xlab="",
     ylab="",
     xlim=c(2,14),
     ylim=c(2,14))
mtext("Population-mean GSD (mm)", 1, line=2.5)
mtext("Population-mean GAD (mm)", 2, line=2.5)

lines(-10:20, -10:20, lty=2)

#points(spmeans$GSD, spmeans$GAD, pch=16)
points(scameans$GSD, scameans$GAD, cex=sqrt(scameans$n)/2, 
       xlab="", ylab="", las=1, col="blue")

plot(popmeans$ASD, popmeans$GSD-popmeans$GAD, cex=sqrt(popmeans$n)/2, 
     xlab="", ylab="", las=1)
abline(h=0, lty=2)
mtext("Population-mean ASD (herkogamy, mm)", 1, line=2.5)
mtext("Departure from optimum (GSD-GAD, mm)", 2, line=2.5)

points(scameans$ASD, scameans$GSD-scameans$GAD, cex=sqrt(scameans$n)/2, 
       xlab="", ylab="", las=1, col="blue")
#points(tapply(popmeans$ASD, popmeans$species, mean), tapply(popmeans$GSD-popmeans$GAD, popmeans$species, mean), 
#       cex=1, pch=16, 
#       xlab="", ylab="", las=1, col="black")

dev.off()

# 2. Estimating the G matrix for the S21 population ####

# Read data
dt = read.table("data/S21_diallel.txt", header = T)

dt$ind = dt$PlantID
dt$animal = dt$PlantID
dt$GSD = apply(subset(dt, select=c(GSDl, GSDc, GSDr)), 1, mean, na.rm=T)

# Select variables, mean-scale and multiply by 100 for better mixing
X = subset(dt, select=c(GAD, GSD, ASD, animal, ind, Mother, Date, stage))
X[,1:length(nms)] = apply(X[,1:length(nms)], 2, function(x) 100*x/mean(x, na.rm=T))

# Prepare pedigree
parentped = read.table("data/S21_parentped.txt", header=T)
ped = subset(dt, select = c("animal", "Mother", "Father"))
names(ped) = c("animal", "dam", "sire")

ped2 = MCMCglmm::prunePed(ped, X$animal)
ped2 = rbind(parentped, ped2)
ped2$sire = as.factor(ped2$sire)

invA = inverseA(ped2)$Ainv

# Set sampling parameters for MCMCglmm analysis
samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

# Set prior
n = length(nms)
alpha.mu = rep(0, n)
alpha.V = diag(n)*400
prior = list(R=list(V=diag(n), nu=n+0.002-1), 
             G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                    G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                    G3=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

# Run MCMC
a = Sys.time()
mod<-MCMCglmm(cbind(GAD, GSD, ASD) ~ -1+trait:factor(stage),
              random= ~us(trait):animal + us(trait):ind + us(trait):Date,
              rcov= ~us(trait):units,
              data=X, ginverse=list(animal = invA),
              family=rep("gaussian", n), prior=prior, nitt=nitt, burnin=burnin, thin=thin)
Sys.time() - a
# Thin 500 = 50 min

filename = paste0("GMatFit_", "thin_", thin, ".RData")

#save(mod, file=filename)

# Postprocessing
load(file="GMatFit_thin_500.RData")

n = length(nms)
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/10000 # Divide by 10000 because data were multiplied by 100
gmat = gmat*100 # Percent units
colnames(gmat) = rownames(gmat) = nms
round(gmat, 3) #G-matrix
round(cov2cor(gmat), 2) #Correlation matrix

# Credible intervals
q = apply(mod$VCV, 2, quantile, c(0.025, 0.975))
t(round(q/100, 3)[,1:9])

round(evolvabilityMeans(gmat[1:3, 1:3]), 3) # Means over the G-matrix

ci = evolvabilityMeansMCMC(mod$VCV[,1:(n*n)]/10000)
round(ci$post.medians*100, 3)

# Random selection gradients
set.seed(99)
betas = randomBeta(1000, 3)

evolvabilityBeta(gmat[1:3, 1:3], Beta=c(1,0,0))$a # GAD: a = 80.4
evolvabilityBeta(gmat[1:3, 1:3], Beta=c(0,1,0))$a #GSD: a = 77.9
evolvabilityBeta(gmat[1:3, 1:3], Beta=c(0,0,1))$a #ASD: a = 84.8

out = evolvabilityBeta(gmat[1:3, 1:3], betas)
str(out)

# Figure 3
x11(height=3.5, width=9)

#cairo_pdf("Fig3.pdf", height=3.5, width=9, fam="Times")
par(mfrow=c(1,3))

hist(out$e, las=1, main="", xlab="", col="lightblue", ylab="")
mtext("Evolvability (%)", 1, line=2.5)

hist(out$a*100, las=1, main="", xlab="", col="lightblue", ylab="")
mtext("Autonomy (%)", 1, line=2.5)

hist(out$c, las=1, main="", xlab="", col="lightblue", ylab="")
mtext("Conditional evolvability (%)", 1, line=2.5)

dev.off()

# 3. Estimating the population-level D matrix ####

# D and P matrix for Costa Rican populations
dat = read.table("data/costarica_greenhouse.txt", header=T)

dat$GSD = apply(subset(dat, select=c(GSDl, GSDc, GSDr)), 1, mean, na.rm=T)
dat$ind = dat$ID
dat$ind = factor(dat$ind)
dat$population = factor(dat$Pop)

unique(dat$population)
length(unique(dat$population))
table(dat$population)

# Mean P matrix
pdat = dat[dat$population %in% c("S1","S2","S6","S7","S8","S9","S11","S12","S20","S21","S22","S26"),]
pdat = na.omit(subset(pdat, select=c(GAD, GSD, ASD, population)))

pops = unique(pdat$population)
length(pops)

Plist = list()
for(i in 1:length(pops)){
  sub = pdat[pdat$population==pops[i],]
  pm = cov(sub[,1:3])
  Plist[[i]] = meanStdG(pm, colMeans(sub[,1:3]))
}

MeanP = apply(simplify2array(Plist), 1:2, mean)

# D matrix
dat = dat[dat$population %in% c("S1","S2","S3","S6","S7","S8","S9","S11","S12","S13","S16","S18","S19","S20","S21","S22","S26"),]
dat = dat[dat$ASD>0,]

df = ddply(dat, .(population), summarize,
           n = sum(!is.na(ASD)),
           GAD = mean(GAD, na.rm=T),
           GSD = mean(GSD, na.rm=T),
           ASD = mean(ASD, na.rm=T))
df = df[-17,]
rownames(df) = df[,1]           
df = df[,-1]
df

X = subset(dat, select=c(GAD, GSD, ASD, population, ind))
X[,1:length(nms)] = apply(X[,1:length(nms)], 2, function(x) 10*log(x)) #Ln-transform and multiply by 10 for mixing reasons
X = na.omit(X)
head(X)
str(X)

samples = 1000
thin = 100
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

n = length(nms)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*50
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

a = Sys.time()
mod<-MCMCglmm(cbind(GAD, GSD, ASD) ~ -1+trait,
              random= ~us(trait):population + us(trait):ind,
              rcov=~us(trait):units,
              data=X,
              family=rep("gaussian", n), prior=prior, nitt = nitt, burnin = burnin, thin = thin)
Sys.time() - a

#save(mod, file="Dmat_CostaRica_thin_100.RData")

# 4. Estimating the species-level D matrix ####
field = read.csv("data/field2.csv")
dat = field

dat = dat[-which(dat$species=="S2"),] # Exclude D. scandens taxon with near-zero ASD
dat$ASD[which(dat$ASD==0)] = 0.1 # Setting a single ASD=0 blossom to ASD=0.1

X = subset(dat, select=c(GAD, GSD, ASD, species, pop, patch))
X[,1:length(nms)] = apply(X[,1:length(nms)], 2, function(x) 10*log(x))
X = na.omit(X)
X$patch = as.factor(paste0(X$pop, X$patch))

head(X)
str(X)
table(dat$species)
length(unique(dat$pop))
length(unique(dat$species))

hist(X$GAD, breaks=20)
hist(X$GSD)
hist(X$ASD)

samples = 1000
thin = 100
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

n = length(nms)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*50
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G3=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

a = Sys.time()
mod<-MCMCglmm(cbind(GAD, GSD, ASD) ~ -1+trait,
              random= ~us(trait):species + us(trait):pop + us(trait):patch,
              rcov=~us(trait):units,
              data=X,
              family=rep("gaussian", n), prior=prior, nitt = nitt, burnin = burnin, thin = thin)
Sys.time() - a

#save(mod, file="Dmat_Species_thin_100.RData")

# 5. Divergence analysis ####

# G-matrix
load(file="GMatFit_thin_500.RData")
n = length(nms)
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = names(X)[1:length(nms)]
gmat

# D_P matrix
load(file="Dmat_CostaRica_thin_100.RData")
n = length(nms)
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = names(X)[1:length(nms)]
round(dmat, 3)
round(cov2cor(dmat), 2)

# Eigenvectors etc.
first_ev = eigen(gmat)$vectors[,1]
gmax = evolvabilityBeta(gmat, Beta = first_ev)$e
dmax = evolvabilityBeta(dmat, Beta = first_ev)$e

last_ev = eigen(gmat)$vectors[,nrow(gmat)]
gmin = evolvabilityBeta(gmat, Beta = last_ev)$e
dmin = evolvabilityBeta(dmat, Beta = last_ev)$e

# Evolvabilities, autonomies, and divergence along random selection gradients
set.seed(99)
nbeta = 1000
betas = randomBeta(nbeta, nrow(gmat))
ebeta = evolvabilityBeta(gmat, betas)$e
cbeta = evolvabilityBeta(gmat, betas)$c
abeta = evolvabilityBeta(gmat, betas)$a
dbeta = evolvabilityBeta(dmat, betas)$e

# Figure 4
x11(height = 9, width = 9)

#cairo_pdf("Fig4.pdf", height=8, width=8, fam="Times")
par(mfrow=c(2,2), mar=c(4,4,1,1))

xminvals = na.omit(c(log10(gmin), log10(diag(gmat)), log10(ebeta)))
xmaxvals = na.omit(c(log10(gmax), log10(diag(gmat)), log10(ebeta)))
yminvals = na.omit(c(log10(dmin), log10(diag(dmat)), log10(dbeta)))
ymaxvals = na.omit(c(log10(dmax), log10(diag(dmat)), log10(dbeta)))

plot(log10(ebeta), log10(dbeta), col="grey", las = 1,
     xlab="", xaxt="n", yaxt="n",
     ylab="",
     xlim=c(min(xminvals[xminvals>-Inf]), max(xmaxvals[xmaxvals>-Inf])),
     ylim=c(min(yminvals[yminvals>-Inf]), max(ymaxvals[ymaxvals>-Inf])))
#legend("topleft", "(a)", bty="n")
axis(1, at=seq(-0.4, 0.8, 0.2), signif(10^seq(-0.4, 0.8, 0.2), 2))

x3at = seq(-0.5, 1, .3)
x3 = (exp(sqrt(((10^x3at)/100)*(2/pi))))-1 # Y-axis given as proportional divergence
axis(2, at=x3at, signif(100*x3, 2), las=1)

mtext("Population divergence (%)", side=2, line=2.5)

points(log10(diag(gmat)), log10(diag(dmat)), pch=16)
points(log10(gmax), log10(dmax), col="red", pch=16)
points(log10(gmin), log10(dmin), col="blue", pch=16)

legend("bottomright", c(expression(paste(g[max])), expression(paste(g[min]))), pch=16, col=c("red", "blue"))

xminvals = na.omit(c(log10(gmin), log10(diag(gmat)), log10(cbeta)))
xmaxvals = na.omit(c(log10(gmax), log10(diag(gmat)), log10(cbeta)))
yminvals = na.omit(c(log10(dmin), log10(diag(dmat)), log10(dbeta)))
ymaxvals = na.omit(c(log10(dmax), log10(diag(dmat)), log10(dbeta)))

plot(log10(cbeta), log10(dbeta), col="grey", las = 1,
     xlab="", xaxt="n",
     ylab="", yaxt="n",
     xlim=c(min(xminvals[xminvals>-Inf]), max(xmaxvals[xmaxvals>-Inf])),
     ylim=c(min(yminvals[yminvals>-Inf]), max(ymaxvals[ymaxvals>-Inf])))
#legend("topleft", "(b)", bty="n")

axis(1, at=seq(-0.4, 0.8, 0.2), signif(10^seq(-0.4, 0.8, 0.2), 2))

x3at = seq(-0.5, 1, .3)
x3 = (exp(sqrt(((10^x3at)/100)*(2/pi))))-1
axis(2, at=x3at, signif(100*x3, 2), las=1)

cvals = NULL # Conditional evolvabilities per trait (cond. on all other traits)
for(i in 1:ncol(gmat)){
  b = rep(0, ncol(gmat))
  b[i] = 1
  cvals[i] = evolvabilityBeta(gmat, b)$c
}

points(log10(diag(gmat)), log10(diag(dmat)), pch=1, col="black")
points(log10(cvals), log10(diag(dmat)), pch=16, col="black")
arrows(log10(diag(gmat)), log10(diag(dmat)), log10(cvals), log10(diag(dmat)), code=2, length=0)

points(log10(gmax), log10(dmax), col="red", pch=16)
points(log10(gmin), log10(dmin), col="blue", pch=16)

# D_S matrix
load(file="Dmat_Species_thin_100.RData")

dsmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dsmat) = rownames(dsmat) = names(X)[1:length(nms)]
round(dsmat, 3)
round(cov2cor(dsmat), 2)

# Compute eigenvectors etc.
dsmax = evolvabilityBeta(dsmat, Beta = first_ev)$e
dsmin = evolvabilityBeta(dsmat, Beta = last_ev)$e
dsbeta = evolvabilityBeta(dsmat, betas)$e # Divergence along the random selection gradients

xminvals = na.omit(c(log10(gmin), log10(diag(gmat)), log10(ebeta)))
xmaxvals = na.omit(c(log10(gmax), log10(diag(gmat)), log10(ebeta)))
yminvals = na.omit(c(log10(dsmin), log10(diag(dsmat)), log10(dsbeta)))
ymaxvals = na.omit(c(log10(dsmax), log10(diag(dsmat)), log10(dsbeta)))

plot(log10(ebeta), log10(dsbeta), col="grey", las = 1,
     xlab="", xaxt="n",
     ylab="", yaxt="n",
     xlim=c(min(xminvals[xminvals>-Inf]), max(xmaxvals[xmaxvals>-Inf])),
     ylim=c(min(yminvals[yminvals>-Inf]), max(ymaxvals[ymaxvals>-Inf])))
#legend("topleft", "(c)", bty="n")

axis(1, at=seq(-0.4, 0.8, 0.2), signif(10^seq(-0.4, 0.8, 0.2), 2))
x3at = seq(-0.8, 2, 0.2)
x3 = (exp(sqrt(((10^x3at)/100)*(2/pi))))-1
axis(2, at=x3at, signif(100*x3, 2), las=1)

mtext("Evolvability (%)", side=1, line=2.5)
mtext("Species divergence (%)", side=2, line=2.5)

points(log10(diag(gmat)), log10(diag(dsmat)), pch=16)
points(log10(gmax), log10(dsmax), col="red", pch=16)
points(log10(gmin), log10(dsmin), col="blue", pch=16)

xminvals = na.omit(c(log10(gmin), log10(diag(gmat)), log10(cbeta)))
xmaxvals = na.omit(c(log10(gmax), log10(diag(gmat)), log10(cbeta)))
yminvals = na.omit(c(log10(dsmin), log10(diag(dsmat)), log10(dsbeta)))
ymaxvals = na.omit(c(log10(dsmax), log10(diag(dsmat)), log10(dsbeta)))

plot(log10(cbeta), log10(dsbeta), col="grey", las = 1,
     xlab="", xaxt="n",
     ylab="", yaxt="n",
     xlim=c(min(xminvals[xminvals>-Inf]), max(xmaxvals[xmaxvals>-Inf])),
     ylim=c(min(yminvals[yminvals>-Inf]), max(ymaxvals[ymaxvals>-Inf])))
#legend("topleft", "(d)", bty="n")

axis(1, at=seq(-0.4, 0.8, 0.2), signif(10^seq(-0.4, 0.8, 0.2), 2))
x3at = seq(-0.8, 2, 0.2)
x3 = (exp(sqrt(((10^x3at)/100)*(2/pi))))-1
axis(2, at=x3at, signif(100*x3, 2), las=1)

mtext("Conditional evolvability (%)", 1, line=2.5)

points(log10(diag(gmat)), log10(diag(dsmat)), pch=1, col="black")
points(log10(cvals), log10(diag(dsmat)), pch=16, col="black")
arrows(log10(diag(gmat)), log10(diag(dsmat)), log10(cvals), log10(diag(dsmat)), code=2, length=0)

points(log10(gmax), log10(dsmax), col="red", pch=16)
points(log10(gmin), log10(dsmin), col="blue", pch=16)

dev.off()
