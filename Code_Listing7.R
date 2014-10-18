# Code Listing 7: Age and gender matching of single primary melanoma and control samples for CDKN2A sequencing.
# Copyright Sebastian Burgstaller-Muehlbacher, 2014. Distributed under GPLv3 or later, without warranty.
# Automatic matcher for age and gender, for a given set of samples


# getAgeDist returns a matrix containing the age distribution and the percentage,
# each Age in years contributed to the whole Dataset.
getAgeDistr <- function(inputSet){
  count_age<- matrix(0,ncol=2)
  total_rows <- nrow(inputSet)
  
  while(nrow(inputSet) > 0){
    count <- nrow(inputSet)
    age <- round(inputSet[count,6],0)
    count_age <- rbind(count_age,c(age, nrow(subset(inputSet,round(inputSet$Age_Dg,0) == age))/100))
    
    inputSet <- subset(inputSet, round(inputSet$Age_Dg,0) != age)
    
  }
  cat(nrow(count_age))
  return (count_age[2:nrow(count_age),])
}

addToMatches <- function(matches,addedStuff, total_req){
  random_sample <- sample(1:nrow(addedStuff), total_req)
    
  for(i in 1:length(random_sample)){
    if(total_req > 0)
      matches <- rbind(matches,addedStuff[random_sample[i],])
    
  }
  return(matches)
}

# this function does the matching
match_smpls <- function(age_Distr, fraction, match_pool, gender){
  
  match_pool <- subset(match_pool, match_pool$Geschlecht == gender)
  
  matches <- data.frame()
  age_group_chkr <- matrix(ncol=4, nrow=100)
  
  # sort the age distribution according to age and convert to data.frame
  age_Distr <- data.frame(age_Distr[order(age_Distr[,1]),])
  
  #fill age_group_chkr from 1 to 100 years
 
  for(i in 1:100){
    req_nr <- 0
    for(ii in 1:nrow(age_Distr)){
      if(age_Distr[ii,1] == i){
        req_nr <- round(age_Distr[ii,2]*fraction*100,0)
      } 
    }
    
    age_group_chkr[i,1] <- i
    age_group_chkr[i,2] <- nrow(subset(match_pool, match_pool$Alter == i))
    age_group_chkr[i,3] <- req_nr #required # of samples
    
    age_group_chkr[i,4] <- age_group_chkr[i,2]-round(req_nr,0) #inititally missing per age group
    
  }
  
  for(i in 1:100){
    if(age_group_chkr[i,4] >= 0 & age_group_chkr[i,2] > 0 & age_group_chkr[i,3] > 0) {
      all_from_AgeGroup <- subset(match_pool, match_pool$Alter == i)
      ran_sel <- sample(1:nrow(all_from_AgeGroup),age_group_chkr[i,3])
      for(ii in 1:length(ran_sel)){
        matches <- rbind(matches, all_from_AgeGroup[ran_sel[ii],])
        match_pool <- subset(match_pool, all_from_AgeGroup[ran_sel[ii],1] != match_pool$PatCode)
        age_group_chkr[i,2] <- age_group_chkr[i,2] - 1
        age_group_chkr[i,3] <- age_group_chkr[i,3] - 1
      }
    }
  }
  
  for(i in 1:100) {
    if(age_group_chkr[i,4] < 0 & age_group_chkr[i,2] > 0 & age_group_chkr[i,3] > 0){
      all_from_AgeGroup <- subset(match_pool, match_pool$Alter == i)
      matches <- rbind(matches, all_from_AgeGroup)
      age_group_chkr[i,2] <- 0 # set to zero, as all have been used up
      age_group_chkr[i,3] <- 0
      
      #remove from match pool
      for(i in 1:nrow(all_from_AgeGroup)){
        match_pool <- subset(match_pool, all_from_AgeGroup[i,1] != match_pool$PatCode)
      }
    }
  }
  
  # pool one up and one down and use them, do that as long as there are elements available
  for(i in 1:100){
    if(age_group_chkr[i,4] < 0 & age_group_chkr[i,3] > 0){
      req_count <- age_group_chkr[i,4] *(-1)
      
      pos <- i-1
      if(pos < 1) pos <- 1
      
      while(req_count > 0 & pos > 0){
        cat(pos,"\n")
        if(age_group_chkr[pos,2]>=req_count){
          all_from_AgeGroup <- subset(match_pool, match_pool$Alter == pos)
          ran_sel <- sample(1:nrow(all_from_AgeGroup),req_count)
          for(ii in 1:length(ran_sel)){
            matches <- rbind(matches, all_from_AgeGroup[ran_sel[ii],])
            match_pool <- subset(match_pool, match_pool$PatCode != all_from_AgeGroup[ran_sel[ii],1])
            age_group_chkr[pos,2] <- age_group_chkr[pos,2] - 1
            req_count <- 0
          }
        }
        else if(age_group_chkr[pos,2]<req_count){
          
          all_from_AgeGroup_d <- subset(match_pool, match_pool$Alter == pos)
          all_from_AgeGroup_u <- subset(match_pool, match_pool$Alter == (i-pos+i))
          
          cat("binhier",i,pos,(i-pos+i),nrow(all_from_AgeGroup_d),nrow(all_from_AgeGroup_u),"\n")
          
          totalCount <- nrow(all_from_AgeGroup_d) + nrow(all_from_AgeGroup_u)
          
          while(req_count > 0 & totalCount > 0){
            if(nrow(all_from_AgeGroup_d) > 0){
              for(ii in 1:nrow(all_from_AgeGroup_d)){
                matches <- rbind(matches, all_from_AgeGroup_d[ii,])
                match_pool <- subset(match_pool, match_pool$PatCode != all_from_AgeGroup_d[ii,1])
                req_count <- req_count - 1
                totalCount <- totalCount - 1
                age_group_chkr[pos,2] <- age_group_chkr[pos,2] - 1
              }
            }
            else if(nrow(all_from_AgeGroup_u) > 0){
              for(ii in 1:nrow(all_from_AgeGroup_u)){
                matches <- rbind(matches, all_from_AgeGroup_u[ii,])
                match_pool <- subset(match_pool, match_pool$PatCode != all_from_AgeGroup_u[ii,1])
                req_count <- req_count - 1
                totalCount <- totalCount - 1
                age_group_chkr[pos,2] <- age_group_chkr[pos,2] - 1
              }
            }
          }
          pos <- pos - 1 #go down one further
            
        }
      }
    }
  }
  return(matches)
}


#getReducedSet returns a matrix/table which contains a defined number of rows, determined
#by the varialbe finalCount. The idea is to randomly remove gender matched individuals
#with the same age, only leaving one individual with matched age and gender in the matrix
getReducedSet <- function(inputSet,finalCount){
  totalNumber <- nrow(inputSet)
  while(nrow(inputSet) >= finalCount){
    rnum <- sample(1:nrow(inputSet))
    sampleNr <- inputSet[rnum,inputSet$PatCode]
    age <- round(inputSet[rnum,inputSet$Age_Dg])
  }
}

  #load data from raw data file
  matching_raw <- read.table("matching_raw.txt", header=T, sep=",", stringsAsFactors=F)
  
  # split in cases and controls
  cases <- subset(matching_raw, matching_raw$ID =="C")
  #controls <- subset(matching_raw, matching_raw$ID =="C")

  #load samples to match for
  match_to <- read.csv("fam+mpm_list.csv", header=T, stringsAsFactors=F)
  
  #split in male and female and get length of each matrix
  match_to_m <- match_to[match_to$Geschlecht == "m",]
  count_m <- nrow(match_to_m)
  cat(count_m,"\n")

  match_to_f <- match_to[match_to$Geschlecht == "f",]
  count_f <- nrow(match_to_f)
cat(count_f,"\n")
  
  ageDistr_F <- getAgeDistr(inputSet = match_to_f)

  ageDistr_M <- getAgeDistr(inputSet = match_to_m)
  
  
  # reduce cases for those who should be matched for
  for(i in 1:nrow(match_to)){
    cases <- subset(cases, cases$PatCode != match_to[i,1])
  }

# match for cases and controls
cases_matched_m<- match_smpls(ageDistr_M,fraction=(count_m/nrow(match_to)*1.1),cases,gender="m")
cases_matched_f<- match_smpls(ageDistr_F,fraction=(count_f/nrow(match_to)*2),cases,gender="f")

library("sm")

tmp1<- cbind(cases_matched_m$Alter, cases_matched_m$Alter)
tmp1[,1] <- 1

tmp2 <- cbind(match_to_m$Age_Dg, match_to_m$Age_Dg)
tmp2[,1] <- 2

age_dens_m <- rbind(tmp1,tmp2)
sm.density.compare(age_dens_m[,2],age_dens_m[,1] , xlab="Age")
title(main="Age Distr. FM/MPM vs Ctrl male")

grp.f <- factor(age_dens_m[,1], levels=c(1,2), labels=c("FM/MPM","Ctrl"))

colfill<-c(2:(2+length(levels(grp.f))))
legend(locator(1), levels(grp.f), fill=colfill) 

#for females
  
tmp1<- cbind(cases_matched_f$Alter, cases_matched_f$Alter)
tmp1[,1] <- 1

tmp2 <- cbind(match_to_f$Age_Dg, match_to_f$Age_Dg)
tmp2[,1] <- 2

age_dens_m <- rbind(tmp1,tmp2)
sm.density.compare(age_dens_m[,2],age_dens_m[,1] , xlab="Age")
title(main="Age Distr. FM/MPM vs Ctrl female")

grp.f <- factor(age_dens_m[,1], levels=c(1,2), labels=c("FM/MPM","Ctrl"))

colfill<-c(2:(2+length(levels(grp.f))))
legend(locator(1), levels(grp.f), fill=colfill) 

