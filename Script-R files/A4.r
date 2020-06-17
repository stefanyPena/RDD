##Assignment 4 – Regression Discontinuity Design.
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
 
  Dt <- HansenDt %>% mutate(bac1_cent=HansenDt$bac1-0.08) %>% mutate(kernel_r=ifelse(abs(bac1_cent)<0.05,1/(2*0.05),0)) %>% mutate(kernel_t=ifelse(abs(bac1_cent)<0.05,1-(abs(bac1_cent)/0.05),0))
  bac1_cent <- HansenDt$dbac1*HansenDt$bac1
  reg_1 <- lm_robust(male ~ bac1 + 'HansenDt$dbac1*HansenDt$bac1' + bac1_cent*HansenDt$bac1, data = HansenDt[(HansenDt$bac1 >= 0.05 & HansenDt$bac1 <= 0.13),])
  summary(reg_1)
  dbac1*bac1
  reg_2 <- lm_robust(age ~ HansenDt$bac1 + bac1_cent + bac1_cent*bac1, data = HansenDt[(HansenDt$bac1 >= 0.05 & HansenDt$bac1 <= 0.13),])
  summary(reg_2)
  
  reg_3 <- lm_robust(white ~ HansenDt$bac1 + bac1_cent + bac1_cent*HansenDt$bac1, data = HansenDt[(HansenDt$bac1 >= 0.05 & HansenDt$bac1 <= 0.13),])
  summary(reg_3)
  
  reg_4 <- lm_robust(acc ~ HansenDt$bac1 + bac1_cent + bac1_cent*bac1, data = HansenDt[(HansenDt$bac1 >= 0.05 & HansenDt$bac1 <= 0.13),])
  summary(reg_4)
  
  #Output Table
  stargazer( reg_1, reg_2, reg_3, reg_4,
            title = "Regression Discontinuity Estimates for the effect of exceeding BAC threshold o predetermined characteristics.",
            out = "RDD_out.tex",
            out.header = T,
            header = F,
            column.labels = c("White","Male","Age","Accident"),
            covariate.labels = c("DUI","BAC","BAC$*$DUI"),
            keep = "DUI", nobs = T, mean.sd = T, label = "RDDout"
  )
  
  ##··Point 6··##
  Datos=as.data.table(read.table("hansen_dwi.csv",header = T,sep=","))
  
  bin <- cut(HansenDt$bac1, seq(min(HansenDt$bac1),max(HansenDt$bac1),0.002))
  categories = Datos[bac1<=0.16]
   data_bin <- aggregate(HansenDt,list(bin),function(x) { return(c(mean(x),length(x)))})
 
    #MALE
   ggacc <- ggplot() + 
    geom_point(Datos = data_bin, aes(x=bac1[,1],y=acc[,1], alpha=.5)) + 
    stat_smooth(Datos = dt, aes(x = bac1, y = acc, color = factor(DUI), group = factor(DUI)), size = 0.5, method = lm)+
    geom_vline(xintercept = 0.08, linetype = "longdash")+
    scale_x_continuous(name = "BAC", limits = c(0,0.2))+
    scale_y_continuous(name = "")+
    coord_cartesian(ylim=c(0.05,0.25)) +
    scale_alpha_continuous(name = "BAC", breaks = "0.5", labels = "")+
    scale_color_manual(values = c("royalblue4","dodgerblue"), name = "", breaks = c("0","1"), labels = c("Control Fit","Treatment Fit"))+
    labs(title = "Panel A. Accident at scene")+
    theme_classic()+
    theme(legend.position = "bottom")
    
     malemeans <- split(categories$male, cut(categories$bac1, 100)) %>% 
    lapply(mean) %>% 
    unlist()
     agg_male_data <- data.frame(male = malemeans, bac1 = seq(0.0016,0.16, by = 0.0016))
     male_l=ggplot(categories, aes(bac1, male)) +
       geom_point(aes(x = bac1, y = male), data= agg_male_data) +
       stat_smooth(aes(bac1, male), method = "lm") +
       geom_vline(xintercept = 0.08)
     male_l
     male_q=ggplot(categories, aes(bac1, male)) +
       geom_point(aes(x = bac1, y = male), data = agg_male_data) +
       stat_smooth(aes(bac1, male,), formula = y ~ x + I(x^2),method = "lm") +
       geom_vline(xintercept = 0.08)
     
     ##WHITE
     whitemeans <- split(categories$white, cut(categories$bac1, 100)) %>% 
       lapply(mean) %>% 
       unlist()
     agg_white_data <- data.frame(white = whitemeans, bac1 = seq(0.0016,0.16, by = 0.0016))
     white_l=ggplot(categories, aes(bac1, white)) +
       geom_point(aes(x = bac1, y = white), data= agg_white_data) +
       stat_smooth(aes(bac1, white), method = "lm") +
       geom_vline(xintercept = 0.08)
     white_q=ggplot(categories, aes(bac1, white)) +
       geom_point(aes(x = bac1, y = white), data = agg_white_data) +
       stat_smooth(aes(bac1, white), formula = y ~ x + I(x^2),method = "lm") +
       geom_vline(xintercept = 0.08)
     
     ##AGED
     agedmeans <- split(categories$aged, cut(categories$bac1, 100)) %>% 
       lapply(mean) %>% 
       unlist()
     agg_aged_data <- data.frame(aged = agedmeans, bac1 = seq(0.0016,0.16, by = 0.0016))
     aged_l=ggplot(categories, aes(bac1, aged)) +
       geom_point(aes(x = bac1, y = aged), data = agg_aged_data) +
       stat_smooth(aes(bac1, aged), method = "lm") +
       geom_vline(xintercept = 0.08)
     aged_l
     aged_q=ggplot(categories, aes(bac1, aged)) +
       geom_point(aes(x = bac1, y = aged), data = agg_aged_data) +
       stat_smooth(aes(bac1, aged), formula = y ~ x + I(x^2), method = "lm") +
       geom_vline(xintercept = 0.08)
     
     
     ##ACC
     accmeans <- split(categories$acc, cut(categories$bac1, 100)) %>% 
       lapply(mean) %>% 
       unlist()
     agg_acc_data <- data.frame(acc = accmeans, bac1 = seq(0.0016,0.16, by = 0.0016))
     acc_l=ggplot(categories, aes(bac1, acc)) +
       geom_point(aes(x = bac1, y = acc), data = agg_acc_data) +
       stat_smooth(aes(bac1, acc), method = "lm") +
       geom_vline(xintercept = 0.08)
     acc_q=ggplot(categories, aes(bac1, acc)) +
       geom_point(aes(x = bac1, y = acc), data = agg_acc_data) +
       stat_smooth(aes(bac1, acc), formula = y ~ x + I(x^2),method = "lm") +
       geom_vline(xintercept = 0.08)
     
     ggarrange(male_l,white_l,aged_l,acc_l,ncol=2,nrow=2)
     ggarrange(male_q,white_q,aged_q,acc_q,ncol=2,nrow=2)
