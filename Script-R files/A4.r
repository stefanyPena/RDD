load.lib <- c('PerformanceAnalytics','xtable','stargazer', 'ggpubr','xts',"readxl","data.table","lubridate",'tidyr','rddensity','estimatr','rdrobust',
              'dplyr','ggplot2','forecast','RColorBrewer','quadprog','NMOF','fBonds','cvar','mFilter',"x13binary","seasonal","lmtest","RCurl") 
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib)
sapply(load.lib,require,character=TRUE)

url_robust <- "https://raw.githubusercontent.com/IsidoreBeautrelet/economictheoryblog/master/robust_summary.R"
eval(parse(text = getURL(url_robust, ssl.verifypeer = FALSE)),
     envir=.GlobalEnv)
