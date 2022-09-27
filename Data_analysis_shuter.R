library(dplyr)
library(ggplot2)
library(ggtree)
library(ggpubr)

## LOADING DATA ---------------------------------------------------------------

setwd("~/dimbo/")
growth_data <- read.csv("dimbo_vs_shuter.csv")

## PLOTTING --------------------------------------------------------------------

p1 <- ggplot(growth_data, 
       aes(x = Genome.Length.NCBI, 
           y = Genome.Size.Shuter)) +
  geom_point(size=5) + 
  scale_y_log10() + 
  scale_x_log10() +
  geom_smooth(method="lm",color="gray") +
  theme_pubclean() +
  xlab("Genome Length (bp)") +
  ylab("Genome Size (pg; Shuter et al.)")


p2 <- ggplot(growth_data, 
       aes(x = Doubling.Time.Shuter, 
           y = Doubling.Time.Dimbo)) +
  geom_point(size=5) + 
  scale_y_log10() + 
  scale_x_log10() +
  geom_smooth(method="lm",color="gray") +
  theme_pubclean() +
  xlab("Minimum Doubling Time (hours; Shuter et al.)") +
  ylab("Minimum Doubling Time (hours)") +
  geom_abline(slope=1,intercept=0,lty=2)


cor.test(growth_data$Genome.Size.Shuter%>%log10(),growth_data$Genome.Length.NCBI%>%log10())
cor.test(growth_data$Doubling.Time.Shuter%>%log10(),growth_data$Doubling.Time.Dimbo%>%log10())

png("Dimbo_vs_Shuter.png",width=1000,height=300)
ggarrange(p2,p1,labels=c("(a)","(b)"))
dev.off()


pdf("Dimbo_vs_Shuter.pdf",width=10,height=3)
ggarrange(p2,p1,labels=c("(a)","(b)"))
dev.off()
