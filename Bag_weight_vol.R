library(readxl)
library(tidyr)
library(dplyr)

library(ggplot2)


Fresh_Bag_weights <- read_excel("C:/Users/90957135/OneDrive - Western Sydney University/Bag_data/Fresh_Bag_weights.xlsx")



Fresh_Bag_weights<-Fresh_Bag_weights%>%
  group_by(Condition)%>%
  mutate(mean_weight= mean(weight_g, na.rm=T),
         std_err_weight= sd(weight_g, na.rm = TRUE) / sqrt(n()),
         mean_vol = mean(Volume,na.rm = T),
         std_err_vol= sd(Volume, na.rm = TRUE) / sqrt(n()),
         percentage_error_weight = (std_err_weight / mean_weight) * 100,
         percentage_error_vol = (std_err_vol / mean_vol) * 100)



# Box plots for the two bag conditions
 ggplot(Fresh_Bag_weights, aes(x = Condition, y = weight_g)) +
  geom_boxplot() +
  labs(x = "Condition", y = "Weight (g)", title = "Box Plot of Bag Weights")+
   theme_minimal()

# Bar graph with standard error
ggplot(Fresh_Bag_weights, aes(x = Condition, y=weight_g)) +
  geom_bar(stat="summary")+
  geom_errorbar(aes(ymin = mean_weight - std_err_weight, ymax = mean_weight + std_err_weight),
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Condition", y = "Mean Weight (g)", title = "Bar Graph of Bag Conditions with Standard Error")+
  theme_minimal()


Fresh_Bag_weights%>%
  ggplot(aes(x=Condition,y=weight_g))+
  geom_bar(stat="summary")+
  geom_errorbar(aes(ymin = mean_weight - std_err_weight, ymax = mean_weight + std_err_weight),
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Condition", y = "Mean Volumes cm^3", title = "Bar Graph of Bag Conditions with Standard Error")+
  theme_minimal()