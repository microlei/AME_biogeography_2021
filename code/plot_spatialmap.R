# code for figure 1A (Doesn't work. Used photoshop instead)
library(ggplot2)
library(maps)
library(mapproj)

stt <- data.frame(name = c("Brewer's Bay", "Black Point", "Flat Cay"), # change these names!
                  lat = c(18.343027, 18.344616, 18.316767),
                  lon = c(-64.980451, -64.986030, -64.987961))
fla <- data.frame(name = c("FLA 015", "FLA 049", "FLA 073"), # change these names!
                  lat = c(24.722388, 25.110120, 25.386340),
                  lon = c(-82.828323, -80.303820, -80.162940))
stt.bg <- map_data("world", region = "Virgin Islands", subregion = "US")
fla.bg <- map_data("state", region = "florida")

stt.map <- ggplot() +
  geom_polygon(data=stt.bg, aes(x=long, y=lat, group=group), fill="grey", alpha=0.3) +
  geom_point(data=stt, aes(x=lon, y=lat, color=name)) +
  coord_map() +
  theme_light()

#fla.map
ggplot() +
  geom_polygon(data=fla.bg, aes(x=long, y=lat, group=group), fill="grey", alpha=0.3) +
  geom_point(data=fla, aes(x=lon, y=lat, color=name)) +
  coord_map() +
  theme_light()

ggplot() +
  geom_polygon(data = map_data("world"), aes(x=long, y=lat, group=group), fill="green", alpha=1) +
  coord_map(xlim=c(-84,-50), ylim=c(15,30)) +
  geom_point(data = bind_rows(stt, fla), aes(x=lon, y=lat, color=name)) +
  theme_light()
