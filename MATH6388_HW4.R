library(ggplot2)
library(tidyverse)
library(dplyr)
library(forcats)
library(ggrepel)

#Number 1
data<-read.table(file="stability.opti_mcc.groups.rarefaction", header=T)
data1 <- na.omit(data)
#plot one OTU
plot(x=data$numsampled, y=data$X0.03.F3D0, xlab="Sequences Sampled",ylab="OTUs", type="l", col="black", font.lab=3)
points(x=data$numsampled, y=data$lci.F3D0, type="l", col="blue")
points(x=data$numsampled, y=data$hci.F3D0, type="l", col="red")
legend(x=10000, y=150, c("0.03.F3D0", "lci", "hci"), c("black", "blue", "red"))

#plot multiple OTUs
plot(x=data$numsampled, y=data$X0.03.F3D0, xlab="Sequences Sampled",ylab="OTUs", type="l", col="black", font.lab=3)
points(x=data$numsampled, y=data$X0.03.F3D1, type="l", col="blue")
points(x=data$numsampled, y=data$X0.03.F3D141, type="l", col="red")
points(x=data$numsampled, y=data$X0.03.F3D142, type="l", col="green")
legend(x=10000, y=150, c("0.03.F3D0", "0.03.F3D1", "03.F3D141", "03.F3D142"), c("black", "blue", "red", "green"))

ggplot(data = data, aes(x=numsampled)) +
  geom_line(aes(y = X0.03.F3D0), color = "black") + 
  geom_line(aes(y = lci.F3D0), color="red", linetype="twodash") +
  geom_line(aes(y = hci.F3D0), color="green", linetype="twodash") +
  theme_minimal() + labs(y="OTUs", x = "Sequences Sampled") 

ggplot(data = data, aes(x=numsampled)) +
  geom_line(aes(y = X0.03.F3D0, color = "X0.03.F3D0")) + 
  geom_line(aes(y = lci.F3D0, color="lci.F3D0")) +
  geom_line(aes(y = hci.F3D0, color="hci.F3D0")) +
  theme_minimal() + labs(y="OTUs", x = "Sequences Sampled") +
  scale_colour_manual("", 
                      breaks = c("X0.03.F3D0", "lci.F3D0", "hci.F3D0"),
                      values = c("black", "red", "green"))

ggplot(data = data1, aes(x=numsampled)) +
  geom_line(aes(y = X0.03.F3D0), color = "black") + 
  geom_ribbon(aes(ymin=lci.F3D0,ymax=hci.F3D0), alpha=0.3) +
  theme_minimal() + labs(y="OTUs", x = "Sequences Sampled") 


#Number 2
summ <- read.table(file="stability.opti_mcc.groups.ave-std.summary", header=T)

#Rearrange group order
sum_avg <- summ[c(1,2,13:19,3:12),]
ggplot(data = sum_avg, aes(x=fct_inorder(group), y =invsimpson)) + geom_point(size = 2) +
  geom_errorbar(aes(ymax = invsimpson_hci, ymin = invsimpson_lci)) + theme_minimal()+
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x="Time Point", y="Inverse Simpson Diversity Index")

#Number 3
pcoa_a <- read.table(file="stability.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes", header=T)
pcoa_l <- read.table(file="stability.opti_mcc.thetayc.0.03.lt.ave.pcoa.loadings", header=T)
plot(pcoa_a$axis1, pcoa_a$axis2)
ggplot(data = pcoa_a, aes(x=axis1, y=axis2, label = group)) + geom_point() + theme_minimal() +
  geom_text_repel(size=2)
  

#Number 4
nmds <- read.table(file="stability.opti_mcc.thetayc.0.03.lt.ave.nmds.axes", header=T)
ggplot(data = nmds, aes(x=axis1, y=axis2, label = group)) + geom_point() + theme_minimal() +
  geom_text_repel(size=2)

#Number 6
spear <- read.table(file="stability.opti_mcc.0.03.subsample.spearman.corr.axes", header=T)
ggplot(data = spear, aes(x=axis1, y=axis2)) + geom_point() + 
  geom_segment(aes(x = 0, y = 0, xend = axis1, yend = axis2), color = "coral1", alpha=.25, linetype=1, arrow = arrow(length = unit(0.03, "npc")))+ theme_minimal()

#Number 7
phyl <- read.table(file="stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.1.phylodiv.rarefaction", header=T)
ggplot(data = phyl, aes(x=numSampled)) +
  geom_line(aes(y = F3D0), color = "black") + 
  geom_line(aes(y = F3D1), color="cornflowerblue") +
  geom_line(aes(y = F3D2), color="coral") +
  geom_line(aes(y = F3D3), color = "chartreuse1") + 
  geom_line(aes(y = F3D5), color="plum4") +
  geom_line(aes(y = F3D6), color="yellow3") +
  geom_line(aes(y = F3D7), color = "firebrick2") + 
  geom_line(aes(y = F3D8), color="dodgerblue2") +
  geom_line(aes(y = F3D9), color="darkorange2") +
  geom_line(aes(y = F3D141), color = "hotpink1") + 
  geom_line(aes(y = F3D142), color="mediumorchid1") +
  geom_line(aes(y = F3D143), color="brown2") +
  geom_line(aes(y = F3D144), color="darkturquoise") +
  geom_line(aes(y = F3D145), color = "darkolivegreen3") + 
  geom_line(aes(y = F3D146), color="gold") +
  geom_line(aes(y = F3D147), color="springgreen2") +
  geom_line(aes(y = F3D148), color = "tomato2") + 
  geom_line(aes(y = F3D149), color="royalblue2") +
  geom_line(aes(y = F3D150), color="violetred") +
  theme_minimal() + labs(y="OTUs", x = "Sequences Sampled") 

ggplot(data = data1, aes(x=numsampled)) +
  geom_line(aes(y = X0.03.F3D0), color = "black") + 
  geom_line(data = na.omit(phyl), aes(x=numSampled, y=F3D0, color = "red")) +theme_minimal()+
  labs(x="Sequences Sampled", y='F3D0 OTUs') + theme(legend.position = "none") 

