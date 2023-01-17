#######################################################
########## Simulation Results for Manuscript ##########
#######################################################

######## Call on Functions and Load Packages ##########

### Function to name the matrix of all combinations of raters
name.pairs = function(n.pairs,mat){
  pairs.name = c(rep(NA,n.pairs))
  for(i in 1:n.pairs){
    pairs.name[i] = paste0("R",mat[1,i],"R",mat[2,i])
  }
  return(pairs.name)
}

source("SampAllSubjAllRatersWithCovar.R") 
source("ObsAgreeSum.R")

##### load packages #####
library(dplyr); library(tidyr); library(psych); library(ggplot2)
library(Matrix);library(limSolve); library(MASS); library(prodlim)
library(purrr); library(lme4); library(irr)

##### define parameters #####

# number of raters
num.rater = 6

# sd of rater random effects
sigma.rater = 0.5

# rater pairs
# matrix of raters to allow for pairwise comparisons
mat.raters = combn(x=num.rater, m=2)
# name the columns of the rater-pairs by rater numbers
colnames(mat.raters) = name.pairs(ncol(mat.raters),mat.raters)

BinOut.list = vector(mode="list",length=num.rater)
for(j in 1:num.rater){
  BinOut.list[[j]] = c(1,0)
}

BinOut.grid = expand.grid(BinOut.list)
MaxJointOutcomes = BinOut.grid[,order(ncol(BinOut.grid):1)]
colnames(MaxJointOutcomes) = paste0("Rater",seq(1,num.rater,1))

# number of subjects
num.subj = 200

# proportion of subjects who are cases
p.case = 0.5
n.case = round(num.subj*p.case)

# subject ID
id <- seq(1,num.subj,1)

# subject measures
group <- c(rep("Control", num.subj-n.case),rep("Case",n.case))
caseInd = group
caseInd = ifelse(caseInd=="Case", 1, 0)

age1 <- c(runif(num.subj-n.case,40,75),runif(n.case,55,90)) #(round(runif(num.subj,40,90)) - 25) / 50 #centered age
age <- (age1 - mean(age1)) / sd(age1)

# beta for subject measures
beta.case = 1.5
beta.age = 1 


subj.eta = beta.case*caseInd + beta.age*age 

# parameters that control the 
# probability of a positive evaluation 
# for the raters in the population
set.seed(4128796)
rater.blup = rnorm(num.rater,0,sigma.rater)

raterExperience = c(rep("Attending",3),rep("Fellow",3))
raterExpInd = raterExperience 
raterExpInd = ifelse(raterExpInd == "Fellow",1,0)

beta.fellow = -0.75

dat.raterExpInd = rep(raterExpInd,each=num.subj)

rater.eta = rater.blup + beta.fellow*raterExpInd

# get matrix of probabilities of a positive evaluation
# for each subject by all raters
ProbSubjRater = matrix(data = NA, nrow=num.subj, ncol=num.rater)
for(i in 1:ncol(ProbSubjRater)){
  ProbSubjRater[,i] = exp(rater.eta[i] + subj.eta) /
    (1 + exp(rater.eta[i] + subj.eta))
}
colnames(ProbSubjRater) = paste0("ProbRater",seq(1,num.rater,1))

### see what the maximum value of kappa can be for 
### each pair of raters and each subject
max.kappa.list = vector(mode = "list", length = ncol(mat.raters)); names(max.kappa.list) = paste0("Raters ", colnames(mat.raters))
max.kappa = rep(NA,length=ncol(mat.raters))
sum.kappa.list = vector(mode = "list", length = ncol(mat.raters))

for(i in 1:ncol(mat.raters)){
  prob1 = exp(rater.eta[mat.raters[1,i]] + subj.eta) /
    (1 + exp(rater.eta[mat.raters[1,i]] + subj.eta))
  prob2 = exp(rater.eta[mat.raters[2,i]] + subj.eta) /
    (1 + exp(rater.eta[mat.raters[2,i]] + subj.eta))
  
  exp.subj = prob1*prob2 + (1-prob1)*(1-prob2)
  
  Posprobs = cbind(prob1,prob2)
  Negprobs = cbind(1-prob1,1-prob2)
  
  p11.subj = apply(Posprobs, 1, min)
  p22.subj = apply(Negprobs, 1, min)
  
  ObsAgree.subj  = p11.subj+p22.subj
  
  m.kap.subj = (p11.subj + p22.subj - exp.subj) / (1 - exp.subj)
  
  sum.kappa.list[[i]] = summary(m.kap.subj)
  
  max.kappa.out = data.frame(prob1,prob2,p11.subj,p22.subj,ObsAgree.subj,exp.subj,m.kap.subj)
  max.kappa.list[[i]] = max.kappa.out
  
  max.kappa[i] = m.kap.subj[1]
}

# select the kappa values upon which to simulate the Bernoulli observations
kappaParm.pair = max.kappa - 0.05

# get theoretical probabilities of agreement based on the selected kappa values
true.probList = vector(mode="list", length=ncol(mat.raters)); names(true.probList) = paste0("Raters ", colnames(mat.raters))
TrueExpAgree = rep(NA,ncol(mat.raters))

for(i in 1:ncol(mat.raters)){
  probI = exp(rater.eta[mat.raters[1,i]] + subj.eta) /
    (1 + exp(rater.eta[mat.raters[1,i]] + subj.eta))
  probJ = exp(rater.eta[mat.raters[2,i]] + subj.eta) /
    (1 + exp(rater.eta[mat.raters[2,i]] + subj.eta))
  
  exp = probI*probJ + (1-probI)*(1-probJ)
  
  obs = kappaParm.pair[i]*(1-exp) + exp
  
  dis = 1-obs
  
  p11 = (probI + probJ - dis) / 2
  p22 = obs - p11
  p12 = probI - p11
  p21 = probJ - p11
  
  RaterI=paste0("Rater",mat.raters[1,i])
  RaterJ=paste0("Rater",mat.raters[2,i])
  
  kap = kappaParm.pair[i]
  
  p11a = p11
  p22a = p22
  
  out = data.frame(id,RaterI,RaterJ,kap,probI,probJ,
                   p11,p22,p12,p21,obs,exp,ProbSubjRater,p11a,p22a)
  
  names(out) = c("id","RaterI","RaterJ","kap","probI","probJ",
                 "p11","p22","p12","p21","obs","exp",
                 paste0("ProbRater",seq(1,num.rater,1)),
                 paste0("p11_R",mat.raters[1,i],"R",mat.raters[2,i]),
                 paste0("p22_R",mat.raters[1,i],"R",mat.raters[2,i]))
  
  true.probList[[i]] = out
  TrueExpAgree[i] = mean(exp)
}

# KappaM
kappa = sum(kappaParm.pair * (1-TrueExpAgree)) / (ncol(mat.raters) - sum(TrueExpAgree))

## Obtain A Matrix used in NNLS (Supplemental Appendix C)
# design matrix for prob of pos agreement for each pair of raters
A1 = matrix(NA, nrow=ncol(mat.raters), ncol=(2^num.rater), byrow = T)
for(k in 1:nrow(A1)){
  MaxJointOutcomes.sub = MaxJointOutcomes[,c(mat.raters[1,k],mat.raters[2,k])]
  MaxJointOutcomes.sum = apply(MaxJointOutcomes.sub,1,sum)
  
  A1.out = replace(MaxJointOutcomes.sum, MaxJointOutcomes.sum==1,0)
  A1.out = replace(A1.out, A1.out==2,1)
  A1[k,] = A1.out
}

# design matrix for prob of neg agreement for each pair of raters
A2 = t(apply(A1,1,rev)) 
# design matrix for the marginal probability of a + response for each rater
A3 = t(MaxJointOutcomes) 

A = as.matrix(rbind(A1,A2,A3))

A.list = vector(mode="list", length=num.rater-2)
# get the design matrices to get all necessary joint probabilities of agreement
for(i in 1:length(A.list)){
  RaterNumber = num.rater-i
  
  numrows = 2^RaterNumber
  numcols = 2^(RaterNumber+1)
  
  AAA = matrix(0,nrow=numrows,ncol=numcols)
  a=1;b=2
  for(j in 1:numrows){
    AAA[j,a:b]=1
    a=a+2
    b=b+2
  }
  
  A.list[[i]] = AAA
}

# get the joint probabilities of agreement
# have to do this by subject AND sequential combinations of raters (i.e. R1, R1+R2, R1+R2+R3,...)
JointProbs.list = vector(mode="list", length=num.rater)

b.list = vector(mode="list", length=length(true.probList))
for(i in 1:length(b.list)){
  b.list[[i]] = data.frame(true.probList[[i]][c("id",paste0("p11_R",mat.raters[1,i],"R",mat.raters[2,i]),
                                                paste0("p22_R",mat.raters[1,i],"R",mat.raters[2,i]),
                                                paste0("ProbRater",seq(1,num.rater,1)))])
}
b.all = reduce(b.list,full_join, by = c("id",paste0("ProbRater",seq(1,num.rater,1))))

b.subj = b.all %>%
  gather(ProbType,b,c(paste0("p11_",colnames(mat.raters)),
                      paste0("p22_",colnames(mat.raters)),
                      paste0("ProbRater",seq(1,num.rater,1)))) %>%
  arrange(id)

b.split = split(b.subj, f=b.subj$id)

MaxProbs.mat = sapply(b.split, function(x) nnls(A,x$b)$X)
colnames(MaxProbs.mat) = paste0("Subject",seq(1,num.subj,1))

JointProbs.list[[1]] = MaxProbs.mat

for(i in 1:length(A.list)){
  A.new = A.list[[i]]
  b.new = JointProbs.list[[i]]
  
  probNew.mat = apply(b.new,2, function(x) A.new %*% x)
  
  JointProbs.list[[i+1]] = probNew.mat
}

JointProbs.list[[num.rater]] = rbind(ProbSubjRater[,1],1-ProbSubjRater[,1])

### Get the conditional probabilities for each rater
### I.e. Prob(R2 | R1) or Prob(R4 | R1, R2,R3)

CondProbs.list = vector(mode="list", length=num.rater)
CondProbs.list[[1]] = JointProbs.list[[num.rater]]

j=2
for(i in num.rater:2){
  denom = JointProbs.list[[i]]
  
  # save only the prob where the next rater says POS
  numerator = apply(JointProbs.list[[i-1]],2,function(x) x[c(T,F)]) 
  
  prob.out = numerator / denom
  
  # some of the higher level probabilities are 0 
  # so the next level probabilities are also 0
  prob.out[is.na(prob.out)] = 0
  
  CondProbs.list[[j]] = prob.out
  
  j=j+1
}


### Now simulate the Bernoulli observations and use these
### observations to estimate kappa using 
### GLMM, Cohen's kappa and Fleiss's kappa

# number of simulations
S = 1000

# store output
## Cohen's kappa
KAPPA.PAIR.UNADJ = matrix(NA, nrow=S, ncol=ncol(mat.raters))
## Kappa under GLMM for each pair of 2 raters
KAPPA.PAIR.ADJ = matrix(NA, nrow=S, ncol=ncol(mat.raters))
## Kappa under GLMM across all raters
KAPPA.POP = matrix(NA, nrow=S, ncol=1)
## Fleiss's kappa
KAPPA.FLEISS = matrix(NA, nrow=S, ncol=1)

# set seed to reproduce results
set.seed(20211211)

for(s in 1:S){
  
  # empty matrix for bernoulli observations
  y.mat = matrix(NA,nrow=num.subj,ncol=num.rater+1)
  y.mat[,1] = seq(1, num.subj, 1)
  
  obs.mat1 = matrix(NA,nrow=num.subj, ncol=num.rater)
  y.mat[,2:ncol(y.mat)] = SampAllSubjAllRatersWithCovar(NumSubj=num.subj,
                                                        NumRater=num.rater, 
                                                        r=1, 
                                                        CondProbList = CondProbs.list,
                                                        obs.mat = obs.mat1,
                                                        BinOutList = BinOut.list)
  
  dat = data.frame(y.mat,caseInd,age)  
  
  colnames(dat) = c("Subj",paste0("Rater",seq(1,num.rater,1)),"Group","Age")
  
  dat.long = dat %>%
    gather(Rater,Y,paste0("Rater",seq(1,num.rater,1)))
  dat.long$Fellow = dat.raterExpInd
  
  # get observed probability of agreement
  cols.obsAgree = c("RaterI","RaterJ","ObsAgree","ObsDisagree")
  
  OBS.PROB.ALL = data.frame(matrix(nrow=0, ncol=length(cols.obsAgree)))
  
  for(r in 1:ncol(mat.raters)){
    prob.obs.all = ObsAgreeSum(dat[,paste0("Rater",mat.raters[1,r])], dat[,paste0("Rater",mat.raters[2,r])])
    
    out.ObsProb.all = cbind(paste0("Rater",mat.raters[1,r]),
                            paste0("Rater",mat.raters[2,r]),
                            prob.obs.all[1],
                            prob.obs.all[2])
    
    OBS.PROB.ALL = rbind(OBS.PROB.ALL, out.ObsProb.all)
  }
  colnames(OBS.PROB.ALL) = c(cols.obsAgree)
  ObsAgree = as.numeric(OBS.PROB.ALL[,"ObsAgree"])
  
  # get expected probability of agreement from GLMM
  mod.rater <- glmer(Y ~ Group + Age + Fellow + (1|Rater),
                    family = binomial("logit"),
                    data = dat.long,
                    control = glmerControl(optimizer = "bobyqa"))
  
  pred.dat = data.frame(dat.long$Rater,dat.long$Subj)
  colnames(pred.dat) = c("Rater","Subj")
  pred.dat$ExpPosProb = predict(mod.rater, type = "response")

  ExpAgree = matrix(NA,nrow=nrow(OBS.PROB.ALL),ncol=1)
  for(e in 1:nrow(OBS.PROB.ALL)){
    ExpPosProb1 = pred.dat[which(pred.dat$Rater == OBS.PROB.ALL[e,"RaterI"]),"ExpPosProb"]
    ExpPosProb2 = pred.dat[which(pred.dat$Rater == OBS.PROB.ALL[e,"RaterJ"]),"ExpPosProb"]
    
    expagree = ExpPosProb1*ExpPosProb2 + (1-ExpPosProb1)*(1-ExpPosProb2)
    ExpAgree[e,] = mean(expagree)
  }
  
  # calculate kappa for each pair of raters based on GLMM
  kappa.pair.adj = (ObsAgree - ExpAgree) / (1 - ExpAgree)
  kappa.est.adj = sum(kappa.pair.adj * (1-ExpAgree)) / (ncol(mat.raters) - sum(ExpAgree))
  
  # get kappa for each pair of raters - without adjustment (i.e., Cohen's kappa)
  kappa.pair.unadj = matrix(NA,nrow=ncol(mat.raters),ncol=1)
  
  for(j in 1:ncol(mat.raters)){
    x1 = dat[,paste0("Rater",mat.raters[1,j])]
    x2 = dat[,paste0("Rater",mat.raters[2,j])]
    
    kapest = cohen.kappa(cbind(x1,x2))
    kappa.pair.unadj[j,] = kapest$kappa
    
  }
  
  # calculate fleiss's kappa 
  fleiss.kappa = kappam.fleiss(dat[,paste0("Rater",seq(1,num.rater,1))])
  KAPPA.FLEISS[s,1] = fleiss.kappa$value
  
  KAPPA.PAIR.UNADJ[s,] = kappa.pair.unadj
  KAPPA.PAIR.ADJ[s,] = kappa.pair.adj
  
  KAPPA.POP[s,1] = kappa.est.adj
  
  if(s %% 10 == 0){print(paste0("S=",s))}
}








