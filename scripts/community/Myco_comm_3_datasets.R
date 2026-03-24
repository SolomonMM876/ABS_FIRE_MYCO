#########################################
#source('Soil_ITS_Scripts/Jeff_proc/Soil_myco_community.R')


source('Bag_data_seq/Bag_myco_communties.R')

source('ABS_Second_Rnd/Bag_myco_comm_2nd_rnd.R')



library(patchwork)
################################
#Plot RA genera

genera<-(Interaction_Genus_myco_bag/Interval_Genus_myco_bag2)


genera



ggsave(filename = "plots/Community_genera_all_3.png", plot = genera, dpi=300,
       device = "png", width = 30, height = 50, units = "cm")



#plot top OTUs

OTU<-(TOP_OTU_bag/TOP_OTU_bag_2)+  plot_layout(heights = c(1, 0.5))  # the third plot gets half the vertical space


OTU

ggsave(filename = "plots/Community_OTU_TOP_all_3.png", plot = OTU, dpi=300,
       device = "png", width = 30, height = 70, units = "cm")



#plot most frequent OTUs
freq<-(top_freq_bag/top_freq_bag_2) +  plot_layout(heights = c(1, 0.5))  # the third plot gets half the vertical space

freq

ggsave(filename = "plots/Community_OTU_frew_3.png", plot = freq, dpi=300,
       device = "png", width = 30, height = 70, units = "cm")

