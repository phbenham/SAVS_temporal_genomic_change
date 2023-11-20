library(ggplot2)
library(gridExtra)

HabLoss<-read.csv("~/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/GeneticDiversityAnalyses/Final_GeneticDiversityAnalyses_26Aug2021/TidalMarsh_HabitatLoss_stats.csv", header=TRUE)

#z-transform all pop GenStats
#plot delta pi, theta, tajd on single figure 
#make three figures percent hab lost, modern ha, historic ha.


#1 more figure how does historic levels of diversity predict change in diversity?

pi<-subset(HabLoss, stat=="pi")
theta<-subset(HabLoss, stat=="theta")
tajd<-subset(HabLoss, stat=="tajd")

#pi.z<-scale(pi[,c(7:9)],center=TRUE,scale=TRUE)
#theta.z<-scale(theta[,c(7:9)],center=TRUE,scale=TRUE)
#tajd.z<-scale(tajd[,c(7:9)],center=TRUE,scale=TRUE)

#pi<-as.data.frame(cbind(pi[,c(1:6)],pi.z))
#theta<-as.data.frame(cbind(theta[,c(1:6)],theta.z))
#tajd<-as.data.frame(cbind(tajd[,c(1:6)],tajd.z))

pi[,c(2,3)]<-log(pi[,c(2,3)])
theta[,c(2,3)]<-log(theta[,c(2,3)])
tajd[,c(2,3)]<-log(tajd[,c(2,3)])



#pi
lm.pi1<-lm(delta~percent_loss, data=pi)
lm.pi2<-lm(delta~historic_ha, data=pi)
lm.pi3<-lm(delta~current_ha, data=pi)
lm.pi4<-lm(historic~historic_ha, data=pi)
lm.pi5<-lm(modern~current_ha,data=pi)
lm.pi6<-lm(modern~ historic_ha,data=pi)
print(summary(lm.pi1))
print(summary(lm.pi2))
print(summary(lm.pi3))
print(summary(lm.pi4))
print(summary(lm.pi5))
print(summary(lm.pi6))

#theta
lm.theta1<-lm(delta~percent_loss, data=theta)
lm.theta2<-lm(delta~historic_ha, data=theta)
lm.theta3<-lm(delta~current_ha, data=theta)
lm.theta4<-lm(historic~historic_ha, data=theta)
lm.theta5<-lm(modern~current_ha,data=theta)
lm.theta6<-lm(modern~historic_ha,data=theta)
print(summary(lm.theta1))
print(summary(lm.theta2))
print(summary(lm.theta3))
print(summary(lm.theta4))
print(summary(lm.theta5))
print(summary(lm.theta6))

#tajd
lm.tajd1<-lm(delta~percent_loss, data=tajd)
lm.tajd2<-lm(delta~historic_ha, data=tajd)
lm.tajd3<-lm(delta~current_ha, data=tajd)
lm.tajd4<-lm(historic~historic_ha, data=tajd)
lm.tajd5<-lm(modern~current_ha,data=tajd)
print(summary(lm.tajd1))
print(summary(lm.tajd2))
print(summary(lm.tajd3))
print(summary(lm.tajd4))
print(summary(lm.tajd5))

HabLoss.Z<-as.data.frame(rbind(pi,theta))
print(HabLoss.Z)

PerHab_plot<-ggplot(HabLoss.Z, aes(x=percent_loss, y=delta, color=stat, fill=stat)) + theme_bw() + geom_point(aes(shape=stat), size=4) + geom_hline(yintercept=0,color="gray") + geom_smooth(method="lm", se=TRUE, linetype='dashed', alpha=0.3, size=0.5) + labs(x="Percent tidal marsh habitat lost", y="Temporal change in diversity (modern-historic)") + theme(legend.position="none")

HistHab_plot<-ggplot(HabLoss.Z, aes(x=historic_ha, y=historic, color=stat, fill=stat)) + theme_bw() + geom_point(aes(shape=stat), size=4) + geom_hline(yintercept=0,color="gray") + geom_smooth(method="lm", se=TRUE, linetype='dashed', alpha=0.3, size=0.5) + labs(x="Historic tidal marsh extent (log-transformed)", y="Historic genetic diversity") + theme(legend.position="none")

ModHab_plot<-ggplot(HabLoss.Z, aes(x=current_ha, y=modern, color=stat, fill=stat)) + theme_bw() + geom_point(aes(shape=stat), size=4) + geom_hline(yintercept=0,color="gray") + geom_smooth(method="lm", se=TRUE, linetype='dashed', alpha=0.3, size=0.5)  + labs(x="Modern tidal marsh extent (log-transformed)", y="Modern genetic diversity") + theme(legend.position="none")

ModGen_HistHab<-ggplot(HabLoss.Z, aes(x= historic_ha, y=modern, color=stat, fill=stat)) + theme_bw() + geom_point(aes(shape=stat), size=4)  + geom_smooth(method="lm", se=TRUE, alpha=0.3, size=0.5)  + labs(x="Historic tidal marsh extent (log-transformed)", y="Modern genetic diversity") + theme(legend.position="none")

HabPlots<-grid.arrange(PerHab_plot, HistHab_plot, ModHab_plot, ncol=1,nrow=3)
HabPlots2<-grid.arrange(PerHab_plot, ModGen_HistHab, ncol=1,nrow=2)
print(HabPlots2)

Admix<-ggplot(HabLoss.Z, aes(x=admix, y=delta, color=stat, fill=stat))  + theme_bw() + geom_point(aes(shape=stat), size=4) + geom_hline(yintercept=0,color="gray") + geom_smooth(method="lm", se=TRUE, linetype='dashed', alpha=0.15, size=0.5)

#print(Admix)

lm.admix.pi<-lm(delta~admix, data=pi)
lm.admix.tajd<-lm(delta~admix, data=tajd)

print(summary(lm.admix.pi))
print(summary(lm.admix.tajd))




