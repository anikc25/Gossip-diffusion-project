rm(list = ls())

library('data.table')
library('dplyr')
library('clubSandwich')
library('stargazer')
library('Hmisc')
library('readstata13')
library('car')
library('miceadds')
library('multiwayvcov')

globdir <- "C:/Users/asankar/Dropbox (MIT)/HIM DATASETS/Analysis/ArunEstherPaper/"
datadir<-paste0(globdir,"data/")
tablesdir <- paste0(globdir,"tables/")


# working at the Village-Month Level Restricted to Seeds Risk Set
villagexmonth_level <- fread(paste0(datadir,"haryana_rct_Figure3.csv"),header = TRUE, sep = ",", data.table = FALSE)

villagexmonth_level <-villagexmonth_level[villagexmonth_level$created_year>2016 & villagexmonth_level$created_year <2019,]

#Month vaccine given
villagexmonth_level$created_month=as.character(villagexmonth_level$created_month)
villagexmonth_level$created_month[nchar(villagexmonth_level$created_month)==1]=paste("0",villagexmonth_level$created_month[nchar(villagexmonth_level$created_month)==1],sep="")
villagexmonth_level$created_month=paste(villagexmonth_level$created_month,paste(villagexmonth_level$created_year),sep="-")
villagexmonth_level$created_month=paste("01-",villagexmonth_level$created_month,sep="")
villagexmonth_level=villagexmonth_level[villagexmonth_level$created_month!="01-07-2022",]
villagexmonth_level=villagexmonth_level[villagexmonth_level$created_month!="01-06-2020",]
villagexmonth_level$created_month=as.Date(villagexmonth_level$created_month,format="%d-%m-%Y")
villagexmonth_level$fes=paste(villagexmonth_level$district,villagexmonth_level$created_month,paste="_")
villagexmonth_level=villagexmonth_level[order(villagexmonth_level$created_month),]

month=villagexmonth_level$created_month[!duplicated(villagexmonth_level$created_month)]
month=as.character(month)
month=month[-1]



num_months = length(month)
res=data.frame(c(rep("nchildren",num_months)))

colnames(res)="vaccine"

n=1
for( i in month){
  model<-lm.cluster(villagexmonth_level[villagexmonth_level$created_month==i,"nchildren"]~gossip+trusted+slope+flat+trustgossip+district,"id_village_grp",data = villagexmonth_level[villagexmonth_level$created_month==i,]) %>% invisible()
  res$coef[n] <- summary(model)[2,1] %>% invisible() #gossipcoef
  res$sd[n] <- summary(model)[2,2] %>% invisible() #gossipstdev
  res$month[n] <- i
  n=n+1
}

res$month=as.Date(res$month,format="%Y-%m-%d")

#Technical parameters for graph
pd <- position_dodge(0.1) # move them .05 to the left and right

j=6  
p6=ggplot(data = res[res$vaccine=="nchildren",],aes(month,coef))+ geom_point(colour = "red")+geom_errorbar(aes(ymin=coef-sd, ymax=coef+sd), width=.1, position=pd,colour="red")+geom_hline(aes(yintercept=0),size=0.2)
p6+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(lineheight=.8, face="bold",hjust=0.5),legend.position = "none")+ scale_y_continuous(breaks=seq(-6,10,2))+ylab("Coefficient of Gossip Treatment")+xlab("Month")+scale_x_date(date_labels = "%b %y",date_breaks = "1 month",limits=c(as.Date("2017-02-01"),as.Date("2018-03-01")))
ggsave(paste0(tablesdir,"Figure3.png"))
dev.off()
