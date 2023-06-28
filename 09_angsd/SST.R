library(rerddap)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(lubridate) # to use ymd_hms
library(viridis) 

# Function to get SST 
getSST <- function(dates, lat, lon){
  SST <- griddap(info("jplMURSST41"), 
                    time = dates,
                    latitude = lat,
                    longitude = lon,
                    fields = "analysed_sst",
                    fmt = "csv")
  SST$time <- as.Date(ymd_hms(SST$time))
  return(SST)
}


# SST <- griddap(info("jplMURSST41"), 
#                time = dates,
#                latitude = nap.lat,
#                longitude = nap.lon,
#                fields = "analysed_sst",
#                fmt = "csv")
# SST <- SST[,4 ]

#  get SST for all 4 sites 
dates= c('2018-01-01','2020-01-01')

nap.lat <-  c(34.42197,34.42197)
nap.lon <-  c(-119.95228,-119.95228)

mon.lat <-  c(36.6182,36.6182) # rounded from 36.61817
mon.lon <-  c(-121.897,-121.897)

dic.lat <- c(35.22445,35.22445)
dic.lon <- c(-120.87748,-120.87748)

pol.lat <- c(32.66533,32.66533)
pol.lon <- c(-117.26152,-117.26152)


nap.sst <- getSST(dates, nap.lat, nap.lon)
mon.sst <- getSST(dates, mon.lat, mon.lon)
dic.sst <- getSST(dates, dic.lat, dic.lon)
pol.sst <- getSST(dates, pol.lat, pol.lon)


SST=NULL
SST$nap.sst <- nap.sst$analysed_sst
SST$mon.sst <- mon.sst$analysed_sst
SST$dic.sst <- dic.sst$analysed_sst
SST$pol.sst <- pol.sst$analysed_sst
SST$dates <- nap.sst$time  
SST <- data.frame(SST)
# Make ggplot 

saveRDS(SST, "SST.rds")

ggplot(data=SST, aes(x=as.Date(dates)))+
  geom_line(aes(y=mon.sst, color ="Monterey")) +
  # geom_line(aes(y=dic.sst, color ="Diablo Canyon")) + 
  geom_line(aes(y=nap.sst, color ="Naples")) +
  geom_line(aes(y=pol.sst, color ="Point Loma")) +
  # scale_color_manual(name="Collection Sites", values=c("Monterey"="gold", "Diablo Canyon"="goldenrod","Naples"="firebrick3","Point Loma"="firebrick")) +
  # scale_color_bamako(reverse = TRUE, discrete = TRUE) +
  scale_color_viridis(option="magma", discrete=TRUE, begin=0, end=0.85, direction = 1) +
  labs(x="Year", y="Sea Surface Temperature (\u00B0C)") + 
  scale_x_date(labels=date_format ("%b %Y"), breaks=date_breaks("6 months"), limits = c(as.Date("2018-01-01"), as.Date("2020-01-01"))) + 
  theme_classic(base_size=18) +
  theme(axis.text.x = element_text(angle = 45, hjust =1)) 
  
# Get color code 

p <- ggplot(data=SST, aes(x=as.Date(dates)))+
  geom_line(aes(y=mon.sst, color ="Monterey")) +
  # geom_line(aes(y=dic.sst, color ="Diablo Canyon")) + 
  geom_line(aes(y=nap.sst, color ="Naples")) +
  geom_line(aes(y=pol.sst, color ="Point Loma")) +
  # scale_color_bamako(reverse = TRUE, discrete = TRUE) +
  # scale_color_viridis(option="magma", discrete=TRUE, begin=0, end=0.85, direction = 1) +
  scale_color_manual(values=c("#000004FF", "#D1426FFF","#FEB77EFF"), name= "Sites") +
  labs(x="Year", y="Sea Surface Temperature (\u00B0C)") + 
  scale_x_date(labels=date_format ("%b %Y"), breaks=date_breaks("6 months"), limits = c(as.Date("2018-01-01"), as.Date("2020-01-01"))) + 
  theme_classic(base_size=18) +
  theme(axis.text.x = element_text(angle = 45, hjust =1)) 
p


pdf("~/KW/figures/figure 1/sst.pdf", width =8.5, height = 5)
  p
dev.off()

tiff("~/KW/figures/figure 1/sst.tiff", width=12, height=5, units = "in", res=300)
p
dev.off()



build <- ggplot_build(p)
build$data[[3]]

"#972C80FF"
"#FEB77EFF"
"#000004FF"
"#D1426FFF"
