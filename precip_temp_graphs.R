#load libraries
library(ggplot2)
library(dplyr)

#source('Temp_Precip_import.R')

#precip per day#########

#facet all sites
p<-ggplot(Sites_Precip,aes(x = date, y = chirps)) +
         geom_line() +
  facet_wrap(~Site, scales = "free_y")

ggsave(p, file=paste0("Precip_Facet", i,".png"),path='precip_graph/', width = 14, height = 10, units = "cm")

#loop for each site individually
plots_list <- list()
for (i in unique(Sites_Precip$Site)) {
  plot_data <- subset(Sites_Precip, i == Site)
  
 
  
  current_plot<-ggplot(plot_data, aes(date, chirps)) +
    geom_line()+
    scale_x_date(date_breaks = "1 year", date_labels = "%b %Y") +
    ylab('Precipitation (mm)') +
    ggtitle(i)
    #facet_wrap(~Burn_Treatment) not needed
  
  plots_list[[i]] <- current_plot
}

#print graphs
for (i in unique(Sites_Precip$Site)) {
  ggsave(plots_list[[i]], file=paste0("plot_", i,".png"),path='precip_graph/', width = 14, height = 10, units = "cm")
}

####for monthly averages####


#calc monthly precip
Sites_Precip <-
  Sites_Precip %>%
  mutate(MonthYear = format(date, "%Y-%m")) %>%
  group_by(MonthYear) %>%
  mutate( AvgPrecipitation = mean(chirps, na.rm = TRUE))

p<-ggplot(Sites_Precip,aes(x = date, y = AvgPrecipitation)) +
  geom_col() +
  facet_wrap(~Site, scales = "free_y")
ggsave(p, file=paste0("plot_year_avg_",".png"),path='precip_graph/', width = 14, height = 10, units = "cm")



m_plots_list <- list()
for (i in unique(Sites_Precip$Site)) {
  plot_data <- subset(Sites_Precip, i == Site)
  
  
  current_plot<-ggplot(plot_data, aes(date, AvgPrecipitation)) +
    geom_col()+
    scale_x_date(date_breaks = "1 year", date_labels = "%b %Y") +
    ylab('Precipitation (mm)') +
    ggtitle(i)
  #facet_wrap(~Burn_Treatment) not needed
  
  m_plots_list[[i]] <- current_plot
}
for (i in unique(Sites_Precip$Site)) {
  ggsave(plots_list[[i]], file=paste0("plot_month_avg_", i,".png"),path='precip_graph/', width = 14, height = 10, units = "cm")
}





#daily temp max
filter(Tmax, between(date, as.Date("2018-01-01"), as.Date("2023-11-24"))) %>%
ggplot(aes(x = date, y = Maximum.temperature..Degree.C.)) +
  geom_line() +
  ggtitle('Tmax daily')->p

ggsave(p, file=paste0("plot_", p,".png"),path='precip_graph/', width = 14, height = 10, units = "cm")




###monthly temp avg

Tmax <-Tmax %>%
  mutate(MonthYear = format(date, "%Y-%m")) %>%
  group_by(MonthYear) %>%
  mutate( AvgTmax= mean(Maximum.temperature..Degree.C., na.rm = TRUE))

filter(Tmax, between(date, as.Date("2018-01-01"), as.Date("2023-11-24"))) %>%
ggplot(aes(x = date, y = AvgTmax)) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  geom_col() +  ggtitle('Tmax monthly')->p

ggsave(p, file=paste0("plot_Tmax", p,".png"),path='precip_graph/', width = 14, height = 10, units = "cm")

