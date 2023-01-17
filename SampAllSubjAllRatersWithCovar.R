### Function used for generating bernoulli observations
### for a number of raters >1
### under a kappa agreement structure,
### and this code will yield estimates of all pairwise kappas between the Raters

### Last Updated 2021.11.17 by Katelyn McKenzie

# NumRater = number of raters
# r = current rater we are sampling
# CondProbList = conditional probabilities of a positive eval given all previous raters
# obs.mat = matrix(NA,nrow=NumSubj, ncol=NumRater)- will contain generated observations
# BinOutList = BinOut.list -> used to get 
#               matrix of observations to obtain the correct Conditional Probability

# function will return "obs.mat"
# NumSubj = num.subj
# NumRater = num.rater
# r = 1
# CondProbList = CondProbs.list
# obs.mat = matrix(NA,nrow=num.subj, ncol=num.rater)
# BinOutList = BinOut.list

SampAllSubjAllRatersWithCovar <- function(NumSubj, NumRater, r, 
                                          CondProbList, obs.mat, BinOutList){
  if(r == (NumRater+1)){
    return(obs.mat)
  }else{
    if(r == 1){
      # get prob of a positive response
      pos.prob1 = CondProbList[[1]][1,]
      # draw observation
      yyy = rbinom(NumSubj,1,pos.prob1)
      # save observation
      obs.mat[,1] = yyy
      
      # update which rater we are drawing
      r=r+1 
      
      # move on to the next rater
      SampAllSubjAllRatersWithCovar(NumSubj, NumRater, r, CondProbList, obs.mat, BinOutList)
      
    }else{
      # all of the conditional probabilities for Rater r given
      # the previous rater's responses
      probs = CondProbList[[r]]
      
      bin.mat = expand.grid(BinOutList[1:(r-1)])
      if(ncol(bin.mat) == 1){
        bin.mat = bin.mat
        prev.obs = matrix(obs.mat[,1:(r-1)],nrow=NumSubj)
      }else{
        bin.mat = bin.mat[,order(ncol(bin.mat):1)]
        prev.obs = obs.mat[,1:(r-1)]
      }
      
      indices = matrix(c(apply(prev.obs, 1, function(x) row.match(x,bin.mat)),
                         seq(1,NumSubj,1)),
                       nrow=NumSubj)
      
      # get prob of a positive response
      pos.probVec = apply(indices, 1, function(x) CondProbList[[r]][x[1],x[2]])
      # draw observations
      yyy = rbinom(NumSubj,1,pos.probVec)
      # save observations
      obs.mat[,r] = yyy
      
      # update which rater we are drawing
      r=r+1 
      
      # move on to the next rater
      SampAllSubjAllRatersWithCovar(NumSubj, NumRater, r, CondProbList, obs.mat, BinOutList)
    }
  }
}