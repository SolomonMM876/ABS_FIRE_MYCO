
#All meta data from 12 sites with bags collected
Bag_data<-read_excel('Processed_data/All_Bag_data.xlsx')
Seq_data<-read_excel('Processed_data/cleaned_seq_dat.xlsx')%>%
  mutate(Site=as.factor(Site), 
         Transect=as.factor(Transect),
         Location=as.factor(Location))%>%
  filter(guild=='Ectomycorrhizal')


#the difference between these two dfs is useful in adapting it to Jeffs' script
Bag_Seq_long<-left_join(Bag_data,Seq_data)
Bag_Seq_wide <- left_join(Bag_data, long_seq)

#taxa table useful in future analyses 
tax<-Seq_data%>%
  select(OTU,kingdom:species, guild)%>%
  distinct()


long_seq<-Seq_data  %>%
  pivot_wider(
    names_from = OTU,         # Column containing OTU names that will become new columns
    values_from = count,      # Column containing values that will fill the new columns
    id_cols = Tube_ID,  # Column(s) to keep as identifier
    values_fill = 0
  )

