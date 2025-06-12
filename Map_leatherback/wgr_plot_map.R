############################
### WGR Sample Locations ###
############################
library("rnaturalearth")
library("rnaturalearthdata")
library("ggplot2")
library("ggpubr")
library("ggspatial")

setwd('/Users/osipova/Documents/LabDocs/Leatherback_popgen/')
sites = read.csv(file="Map_leatherback//WGR_collection_sites.csv")

world <- ne_countries(scale = "medium", returnclass = "sf")
pal=get_palette("Set3", 12)
sites=sites[order(sites$Location),]
#sites<-sites[sites$Location != "Indonesia",]

wgr_map=ggplot(data = world) + theme_classic() +
  geom_sf(color=NA, fill = "grey60") +
  geom_point(data = sites, aes(x = Long, y = Lat), size = 4, 
             shape = 21, fill = pal) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude")

wgr_map

ggsave("Map_leatherback/WGR_map.pdf", wgr_map, height=6, width=8)
