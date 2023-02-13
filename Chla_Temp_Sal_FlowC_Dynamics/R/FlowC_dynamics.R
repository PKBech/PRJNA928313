################
### Packages ###
################

library(tidyr)
library(ggplot2)
library(dplyr)
library(scales)
library(ggbeeswarm)

#################
### Load data ###
#################

#Temperature(HOBO+Logger), Salinity(Logger), chla, abundance.
Metadata_w.FlowC <- read.csv("Chla_Temp_Sal_FlowC_Dynamics/Data/Metadata_w.FlowC.csv")
str(Metadata_w.FlowC)
Metadata_w.FlowC$Phase = as.factor(Metadata_w.FlowC$Phase)

# Keep only bioelements without negatives
metadata.bioelements = Metadata_w.FlowC %>%
  filter(element.type == "bioelement") %>%
  filter(sample.type == "sample")

# get T.abundance as cells per bioelement
metadata.bioelements$T.abundance.per.element = log10(metadata.bioelements$T.abundance*2)
str(metadata.bioelements)

#factor day
#metadata.bioelements$day = (metadata.bioelements$day)

#summarise
metadata.bioelements.stat = metadata.bioelements %>%
  dplyr::group_by(day, Phase) %>%
  dplyr::summarise(mean.t.abund = mean(T.abundance.per.element), 
            sd = sd(T.abundance.per.element))
str(metadata.bioelements.stat)

###################################
### Total abundance (flow data) ###
###################################
color_phase = c("#240785", "#e0b62b", "#f21395")
metadata.bioelements.stat$Phase = factor(metadata.bioelements.stat$Phase, levels= c("Early", "Peak", "Late"))
metadata.bioelements$Phase = factor(metadata.bioelements$Phase, levels= c("Early", "Peak", "Late"))


Abundance = ggplot(metadata.bioelements.stat, aes(x = day, y = mean.t.abund, col = Phase)) + 
  geom_quasirandom(data = metadata.bioelements, aes(y=T.abundance.per.element),
                   dodge.width=.1, cex=1, alpha = 0.5)+
  geom_point(aes(x = day, y = mean.t.abund), size = 2.5, alpha = 1) +
  geom_errorbar(aes(ymax = mean.t.abund-sd, ymin = mean.t.abund+sd),width = 0, 
                size = 0.7, alpha = 0.5)+
  xlab("")+ ylab("Log10(cells/element)\n")+
  scale_x_continuous(breaks = c(3,7,10,15,23,29,44,57,71,85,99,113))+
  theme_bw(base_size = 10)+
  scale_color_manual(values = color_phase)+
  scale_fill_manual(values = color_phase)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(face = "bold", colour = "black"),
        legend.position = "top",
        legend.background = element_rect(linetype = 1, size = 0.5, colour = "lightgrey"))

Abundance

color_phase = c("#240785", "#e0b62b", "#f21395")

Flow.C <- metadata.bioelements.stat %>%
  ggplot(aes(x = day, y = mean.t.abund, col = Phase)) + 
  geom_quasirandom(data = metadata.bioelements, aes(y=T.abundance.per.element, col = Phase), 
                  ,size = 2, dodge.width=0, alpha = 0.6, stroke = NA) + 
  geom_errorbar(aes(ymax = mean.t.abund-sd, ymin = mean.t.abund+sd, col = Phase), width = 0, linewidth = 0.7, alpha = 0.2)+
  geom_line(size = 1, alpha = 0.6, col = "black") +
  labs(x= "\nTime (days)", y = "Log10(cells/element)\n")+
  scale_x_continuous(breaks = c(3,7,10,15,23,29,44,57,71,85,99,113))+
  theme_bw(base_size = 10) +
  #facet_wrap(.~Measure, scales = "free_y") +
  scale_color_manual(values = color_phase)+
  scale_fill_manual(values = color_phase)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.position = "none",
        #legend.background = element_rect(linetype = 1, size = 0.5, colour = "lightgrey"),
        axis.text.x = element_text(angle = 45, vjust =1, hjust = 0.5))

ggsave(Flow.C, file = "Chla_Temp_Sal_FlowC_Dynamics/Figures/FlowC.svg", width=3.5, height=3)



