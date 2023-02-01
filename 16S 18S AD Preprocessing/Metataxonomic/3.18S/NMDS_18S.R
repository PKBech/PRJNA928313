# NMDS 18S all element.type

NMDS_ord_bray = metaMDS(P18S_dist_bray, k=3)
#NMDS_ord_bray_spec = metaMDS(sqrt_genus_norm_dat_t, distance = "bray", autotransform = FALSE, k=3)
str(NMDS_ord_bray)

stressplot(NMDS_ord_bray)

#build a data frame with NMDS coordinates and metadata
NMDS1 = NMDS_ord_bray$points[,1]
NMDS2 = NMDS_ord_bray$points[,2]
NMDS3 = NMDS_ord_bray$points[,3]


# str(NMDS_ord_bray)
# NMDS_ord_bray$species


NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, NMDS2 = NMDS3, 
                  Chla_suc = P18S_filter_norm_meta_dat$chla, 
                  Day = P18S_filter_norm_meta_dat$day,
                  Subject = P18S_filter_norm_meta_dat$element.type)

#Set color theme
#fef0d9
#fdcc8a
#fc8d59
#d7301f
# cols_time <- c("0" = "#f7f7f7", "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("0", "1", "4","10"))

# cols_time <- c( "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("1", "4","10"))




p_time_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_suc_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Chla_suc)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_ass_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Subject)) +
  #stat_ellipse() +
  #scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)
# scale_colour_manual(
#   values = cols_time,
#   aesthetics = c("colour", "fill") 
# )  + xlim(-0.58,0.58) + ylim(-0.4,0.4) + theme_bw() +
# theme(#axis.line = element_line(color='black'),
#   plot.background = element_blank(),
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(), legend.position = "none") 


p_time_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)


p_Chla_suc_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Chla_suc)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_ass_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Chla_ass)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)


grid.arrange(p_time_NMDS12, p_time_NMDS23, p_Chla_suc_NMDS12, p_Chla_suc_NMDS23, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,4)),
             widths = c(2.7, 2.7), heights = c(2.7, 2.7))


# NDMS 18S only bioelements
NMDS_ord_bray = metaMDS(P18S_succession_dist_bray, k=3)
#NMDS_ord_bray_spec = metaMDS(sqrt_genus_norm_dat_t, distance = "bray", autotransform = FALSE, k=3)
str(NMDS_ord_bray)

stressplot(NMDS_ord_bray)

#build a data frame with NMDS coordinates and metadata
NMDS1 = NMDS_ord_bray$points[,1]
NMDS2 = NMDS_ord_bray$points[,2]
NMDS3 = NMDS_ord_bray$points[,3]


# str(NMDS_ord_bray)
# NMDS_ord_bray$species


NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, NMDS2 = NMDS3, 
                  Chla_suc = P18S_succession_filter_norm_meta_dat$chla, 
                  Day = P18S_succession_filter_norm_meta_dat$day,
                  Subject = P18S_succession_filter_norm_meta_dat$element.type)

#Set color theme
#fef0d9
#fdcc8a
#fc8d59
#d7301f
# cols_time <- c("0" = "#f7f7f7", "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("0", "1", "4","10"))

# cols_time <- c( "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("1", "4","10"))




p_time_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_suc_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Chla_suc)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_ass_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Subject)) +
  #stat_ellipse() +
  #scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)
# scale_colour_manual(
#   values = cols_time,
#   aesthetics = c("colour", "fill") 
# )  + xlim(-0.58,0.58) + ylim(-0.4,0.4) + theme_bw() +
# theme(#axis.line = element_line(color='black'),
#   plot.background = element_blank(),
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(), legend.position = "none") 


p_time_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)


p_Chla_suc_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Chla_suc)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_ass_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Chla_ass)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)


grid.arrange(p_time_NMDS12, p_time_NMDS23, p_Chla_suc_NMDS12, p_Chla_suc_NMDS23, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,4)),
             widths = c(2.7, 2.7), heights = c(2.7, 2.7))
