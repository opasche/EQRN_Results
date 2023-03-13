library(tidyverse)
## This is all data, a list of length 68, each list entre has Date and Val time series
load("Data/data_all.Rdata")
## This is all summer data as a tibble object
load("Data/river_dat.Rdata")
## This is the meta-info for all stations
load("Data/info_all.Rdata")
head(info_all)
## This is the matrix of flow-connections (not equal to adjacency matrix of the graph)
load("Data/flow_con_all.Rdata")
## This is the weekly maxima data including NAs
load("Data/river_week.Rdata")
head(river_week)
## This is the daily maxima data of all common dates (in matrix form)
load("Data/river_dat_common.Rdata")
## This is the weekly maxima data of all common dates (in matrix form)
load("Data/river_week_common.Rdata")
## The coordinates of the flow connection edges for plotting
load("Data/coord_edges.Rdata")

Disch <- as_tibble(data_all[[1]])
colnames(Disch)[2] <- "station_1"
for (i in 2:length(data_all)){
  tb <- as_tibble(data_all[[i]])
  colnames(tb)[2] <- paste("station_",i,sep="")
  Disch <- full_join(Disch, tb, by="Date")
}
Disch <- Disch %>% arrange(Date)

write.csv(info_all, file="data_wrangled/info_discharges_initial.csv")
write_csv(Disch, path="data_wrangled/Discharges_all.csv")
write_csv(river_dat, path="data_wrangled/Discharges_summer.csv")



