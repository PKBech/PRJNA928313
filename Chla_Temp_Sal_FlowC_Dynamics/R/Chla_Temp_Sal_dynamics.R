################
### Packages ###
################

library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)

#################
### Load data ###
#################

#Temperature(HOBO+Logger), Salinity(Logger) and Light (HOBO)
Metadata <- readxl::read_excel("Chla_Temp_Sal_FlowC_Dynamics/Data/Metadata.xlsx")

#Chla measurements
Chla_longForma <- read.csv("Chla_Temp_Sal_FlowC_Dynamics/Data/Chla_longFormat.csv")


#Chl a #####
Chla_longForma <- as.data.frame(Chla_longForma)[,c(1,2,3,4,5,7,8)]
colnames(Chla_longForma) <- c("Timepoint", "Date", "Day", "Chl-a (ug)", "Experiment", "Temperature", "Bryozoan appearence")
Chla_longForma$Timepoint <- as.factor(Chla_longForma$Timepoint)

#Summarize
chla.sum = Chla_longForma %>%
  dplyr::group_by(Day, Experiment) %>%
  dplyr::summarise(mean.chla = mean(`Chl-a (ug)`),
            sd.chla = sd(`Chl-a (ug)`))

# #Only succession study
# chla.sum.suc = chla.sum %>%
#   filter(Experiment == "Succession")

#Change day to match HOBO data
# chla.sum.suc$Day = chla.sum.suc$Day+0.5
chla.sum$Day = chla.sum$Day+0.5

#Temp data ####

#From wide to long format
metadata.longf = gather(Metadata, 
                        measurement, value, Temperature:Salinity, factor_key = TRUE)

###############
### Combine ###
###############


# combine datasets and remove all observations with Light = 0 and merge only matching observation point in chl a measurement
metadata.combined = merge(filter(Metadata, Light != 0), 
                          chla.sum, by = "Day", all.y = TRUE)


# Figure combined
MD = ggplot(metadata.combined, aes(x = Day)) + 
  geom_line(size = 1.5, aes(y = Temperature, color = "Temperature (C°)")) +
  geom_line(data= filter(metadata.combined, Experiment %in% c("Succession")),
            size = 1.5, aes(y = mean.chla*10, color = "Chl-a")) +
  geom_ribbon(data= filter(metadata.combined, Experiment %in% c("Succession")),
              aes(ymax = ((mean.chla + sd.chla)*10), ymin = ((mean.chla - sd.chla)*10)),
              alpha = 0.2, color = NA, fill = "#5674a8")+
  scale_y_continuous(name = "Temperature (C°)\n", breaks = c(0,5,10,15,20,25,30), 
                     sec.axis = sec_axis(~ . /10, name = "Chl-a (ug/element) \n"))+
  scale_x_continuous(breaks = c(3,7,10,15,23,29,44,57,71,85,99, 113))+
  scale_color_manual(values = c("Temperature (C°)" = "#d1bb8a",
                                "Chl-a" = "#5674a8"))+
  scale_fill_manual(values = c("Temperature (C°)" = "#d1bb8a",
                               "Chl-a" = "#5674a8"))+
  xlab("\nTime (days)")+
    theme_bw(base_size = 10)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", angle = 45, vjust =1, hjust = 0.5),
        axis.title.x = element_text(face = "bold", colour = "black"),
        axis.title.y = element_text(face = "bold", colour = "black"),
        legend.text = element_text(face ="bold", colour ="black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        #legend.background = element_rect(linetype = 1, size = 0.5, colour = "lightgrey")
        )

MD   
ggsave(MD, file = "Chla_Temp_Sal_FlowC_Dynamics/Figures/Chla_Temp.svg", width=3.45, height=2.6)
