##Assignment 4 – RDD.
##Causal Inference and Research Design.
##Stefany Peña Montenegro.
##Sunday 14th, 2020.
-----------------------------------------------
  
  ##··PART II:Replication··##
  
  ##··My packages··##
  
  load.lib <- c('PerformanceAnalytics','xtable','stargazer', 'ggpubr','xts',"readxl","data.table","lubridate",'tidyr','rddensity','estimatr','rdrobust',
                'dplyr','ggplot2','forecast','RColorBrewer','quadprog','NMOF','fBonds','cvar','mFilter',"x13binary","seasonal","lmtest","RCurl" , "tidyverse","doBy","gdata","ggforce","ggpubr","haven","Hmisc","lubridate","rdd","readxl",
                "sandwich","stargazer","dagitty","rdrobust","lmtest")
  install.lib <- load.lib[!load.lib %in% installed.packages()]
  for(lib in install.lib) install.packages(lib)
  sapply(load.lib,require,character=TRUE)
  
  
  ##··My Data··##
  
  HansenDt <- read.csv("hansen_dwi.csv")
  HansenDt <- HansenDt %>% mutate(Date=as.Date(Date, format = "%d%b%Y"))
  str(HansenDt)
  
  ##··Point 3··##
     ##TreatmentDummy=DUI(DrivingUndertheInfluence)
  
  TreatmentD=ifelse(HansenDt$bac1>=0.08,1,0)
  HansenDt <- HansenDt %>% mutate(TreatmentD)
  stargazer(HansenDt,
            title = "Treatment variable descriptive statistics", 
            summary = T, 
            out = "TreatmentStats.txt",
            keep = "TreatmentD")
  
  
  ##··Point 4··## 
  
  McCrManipulationTest<- rdd::DCdensity(runvar = HansenDt$bac1,cutpoint = 0.08,
                                        htest = T, 
  )
  
  Results1McTest <- McCrManipulationTest[1:7]
  Results1McTest
  
  ggplot(HansenDt,aes(x=`bac1`))+ geom_histogram(aes(y= ..density.. ),color="grey", fill="NA", binwidth=0.001)+
    labs(title = "McCrary density for BAC", x = "Blood alcohol content", y = "Density")+
    geom_vline(xintercept = 0.08, linetype = "dashed",color="blue") +theme(axis.title.x =element_text(size=15),
                                                                           panel.background = element_blank())
  ##··Point 5··##
  
  DUI <- TreatmentD
  
  depvrs <- c("white","male","aged","acc")
  results <- list()
  for (i in 1:4) 
    {
    mdl <- lm(formula = HansenDt[,depvrs[i]] ~  HansenDt[,"DUI"] +  HansenDt[,"bac1_cent"] +  HansenDt[,"bac1_cent"]* HansenDt[,"DUI"], weights =  HansenDt[,"kernel_r"])
    mdl_r <- coeftest(mdl, vcov = vcovHC(mdl,"HC1"))
    results(depvrs) <- mdl_r
  }
  
  
  #Output Table
  stargazer(results[1],results[2],results[3],results[4],
            title = "Regression Discontinuity Estimates for the effect of exceeding BAC threshold o predetermined characteristics.",
            out = "RDD_out.tex",
            out.header = T,
            header = F,
            column.labels = c("White","Male","Age","Accident"),
            covariate.labels = c("DUI","BAC","BAC$*$DUI"),
            keep = "DUI", nobs = T, mean.sd = T, label = "RDDout"
  )
  
