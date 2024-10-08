#0. Set WD and load libraries ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(broom)
library(gapminder)
library(geojsonio)
library(patchwork)
library(plotly)
library(rgeos)
library(sf)
library(sp)
library(tidyverse)


#1. Define functions & themes ----
plottable <- function(x) {
  x <- tidy(x)
  x <- mutate(x,lat*-1)
  x <- x[,-c(2)]
  colnames(x)[7] = "latitude"
  colnames(x)[1] = "longitude"
  return(x)
}
originate <- function(y) {
  y <- select(y,longitude,latitude)
  y_centre <- as.data.frame(cbind(mean(y$longitude),mean(y$latitude)))
  colnames(y_centre) <- c("longitude","latitude")
  transform_y <- 0 - y_centre
  y$longitude <- y$longitude + transform_y$longitude
  y$latitude <- y$latitude + transform_y$latitude
  return(y)
}
originate_points <- function(a){
  hull_indices <- chull(a$longitude, a$latitude)
  hull_points <- a[hull_indices, ]
  centroid <- gCentroid(SpatialPolygons(list(Polygons(list(Polygon(hull_points)), ID = "1"))))
  centroid_coords <- coordinates(centroid)
  colnames(centroid_coords) <- c("longitude","latitude")
  transform_a <- 0 - centroid_coords
  transform_a <- as.data.frame(transform_a)
  vector_name <- paste("vector", deparse(substitute(a)), sep = "_")
  assign(vector_name, transform_a, envir = .GlobalEnv)
  a$longitude <- a$longitude + transform_a$longitude
  a$latitude <- a$latitude + transform_a$latitude
  return(a)
}
overlap <- function(block,ihc){
  block <- as.data.frame(block)
  ihc <- as.data.frame(ihc)
  
  block <- rbind(block,block[1, ])
  ihc <- rbind(ihc,ihc[1, ])
  
  block <- st_polygon(list(cbind(block$longitude, block$latitude))) %>%
    st_sfc() %>%
    st_sf()
  
  ihc <- st_polygon(list(cbind(ihc$longitude, ihc$latitude))) %>%
    st_sfc() %>%
    st_sf()
  intersection <- st_area(st_intersection(block, ihc))
  ihc_area <- st_area(ihc)
  return(intersection / ihc_area)
}
rotation <- function(b,c){
  rotation_angle_radians <- c*(3.141/180)
  
  for (j in 1:nrow(b)){
    old_long <- b[j,1]
    old_lat <- b[j,2]
    new_long <- old_long * cos(rotation_angle_radians) + old_lat * sin(rotation_angle_radians)
    new_lat <- -1 * old_long * sin(rotation_angle_radians) + old_lat * cos(rotation_angle_radians)
    b[j,1] <- new_long
    b[j,2] <- new_lat
  }
  return(b)
}
internal_select <- function(annotation,cells){
  coords_polygon <- cbind(annotation$longitude, annotation$latitude)
  polygon <- Polygon(coords_polygon)
  polygons <- Polygons(list(polygon), ID = "1")
  sp_polygon <- SpatialPolygons(list(polygons))
  coords_points <- cbind(cells$longitude, cells$latitude)
  sp_points <- SpatialPoints(coords_points)
  
  points_inside <- !is.na(over(sp_points, sp_polygon))
  cells <- cells[points_inside, ]
  return(cells)
  print(nrow(cells))
}
theme_spatial <- theme(panel.background = element_rect(fill = "black"),
                       plot.background = element_rect(fill = "black"),
                       panel.grid.major = element_line(color = "gray20"),
                       panel.grid.minor = element_line(color = "gray20"),
                       axis.title.x = element_text(color = "white"), 
                       axis.title.y = element_text(color = "white"),   
                       axis.text.x = element_text(color = "white"),  
                       axis.text.y = element_text(color = "white"),
                       plot.title = element_text(color = "white", face = "bold"),
                       plot.subtitle = element_text(color = "white", face = "italic"),
                       legend.background = element_rect(fill = "black", color = NA),
                       legend.text = element_text(color = "white"),
                       legend.title = element_text(color = "white"))

#2. Import block annotation ----

block <- read_delim("07199_block.txt")
colnames(block)[3] <- "longitude"
colnames(block)[4] <- "latitude"
block$latitude <- block$latitude * -1
block <- select(block, longitude,latitude)

p1 <- ggplot(block,aes(x=longitude,y=latitude)) + 
  geom_polygon(fill="violet",alpha=0.5) + 
  theme_spatial + 
  labs(x="Longitude",
       y="Latitude",
       title = "Block annotation")
ggplotly(p1) # imported block annotation

#3. Import IHC annotations ----

CD3_annotation <- read_delim("07199_CD3_annotation.txt")
colnames(CD3_annotation)[3] <- "longitude"
colnames(CD3_annotation)[4] <- "latitude"
CD3_annotation$latitude <- CD3_annotation$latitude * -1
CD3_annotation <- select(CD3_annotation, longitude,latitude)


p2 <- ggplot() + 
  geom_polygon(data=block,aes(x=longitude,y=latitude),fill="violet",alpha=0.5) + 
  geom_polygon(data=CD3_annotation,aes(x=longitude,y=latitude),fill="cyan",alpha=0.5) + 
  theme_spatial + 
  labs(x="Longitude",
       y="Latitude",
       title = "Block + IHC annotation (untransformed)")

ggplotly(p2) # imported IHC annotation

#4. Originate each annotation ----

block <- originate_points(block)
CD3_annotation <- originate_points(CD3_annotation)


#5. Compare centered outlines ----

p5 <- ggplot() + 
  geom_polygon(data = block,aes(x=longitude,y=latitude),fill="violet",alpha=0.5) + 
  geom_polygon(data = CD3_annotation,aes(x=longitude,y=latitude),fill="cyan",alpha=0.5) + 
  theme_spatial + 
  labs(x="Longitude",
       y="Latitude",
       title = "Block + IHC annotations (originated)")
ggplotly(p5)

#6. Calculate a scale factor to size up IHC annotation----

print(paste0("Proportion of block area: ",
             round(100*overlap(CD3_annotation, block),2),
             "%"
             ))

print(paste0("Proportion of CD3 area: ",
             round(100*overlap(block, CD3_annotation),2),
             "%"
             ))

CD3_scale_factor <- sqrt(1 / overlap(CD3_annotation,block)) # calculates CD3 scale factor
print(paste0("CD3_Scale Factor: ",round(CD3_scale_factor,2),"x"))



CD3_annotation <- CD3_annotation * CD3_scale_factor



p6 <- ggplot() + 
  geom_polygon(data = block,aes(x=longitude,y=latitude),fill="violet",alpha=0.5) + 
  geom_polygon(data = (CD3_annotation),aes(x=longitude,y=latitude),fill="cyan",alpha=0.5) + 
  theme_spatial + 
  labs(x="Longitude",
       y="Latitude",
       title = paste0("Block + IHC annotations (scaled) ",round(100*overlap(CD3_annotation, block),2),
                      "% overlap"))
ggplotly(p6)

print(paste0("CD3 Annotation overlap: ",
             round(100*overlap(CD3_annotation, block),2),
             "%"
))


#7. Calculate a CD3 rotation angle ----



# create blank df for results and set first angle range (1:360)
overlap_trial <- data.frame()
overlap_list <- list()
angle_range <- c(1, seq(10, 360, by = 10))



for (j in angle_range){
  backup <- CD3_annotation
  CD3_annotation <- rotation(CD3_annotation,j)
  overlap_list[[j]] <- data.frame(
    angle = j, 
    overlap = overlap(block,CD3_annotation))
  overlap_trial <- do.call(rbind, overlap_list)
  progress <- round(((j/360)*100),2)
  print(paste0(j,
               " - ",
               progress,
               "% complete"))
  CD3_annotation <- backup
}
means <- rowMeans(embed(overlap_trial$overlap, 2))
max_mean_index <- which.max(means)
new_angle_range <- overlap_trial$angle[c(max_mean_index, max_mean_index + 1)]
angle_range <- seq(new_angle_range[1],new_angle_range[2],length.out=10)
range_anno_1 <- annotate("rect", xmin=c(new_angle_range[1]), 
                   xmax=c(new_angle_range[2]), 
                   ymin=c(min(overlap_trial$overlap)*100) , 
                   ymax=c(100), alpha=0.5, fill="cyan")
ggplotly(
  ggplot(overlap_trial,aes(x=angle,y=overlap*100)) + 
           theme_spatial + 
           geom_point(color="white") + 
           geom_line(color="white") + 
           labs(x="Angle of rotation (degrees)",
                y="Overlap between sections (%)",
                title="Rotation optimisation") + 
    range_anno_1
)
# calculate angle pair to continue with and run 2nd pass
for (j in angle_range){
  backup <- CD3_annotation
  CD3_annotation <- rotation(CD3_annotation,j)
  overlap_list[[j]] <- data.frame(
    angle = j, 
    overlap = overlap(block,CD3_annotation))
  overlap_trial <- do.call(rbind, overlap_list)
  progress <- round(((j/new_angle_range[2])*100),2)
  print(paste0(j,
               " - ",
               progress,
               "% complete"))
  CD3_annotation <- backup
}
means <- rowMeans(embed(overlap_trial$overlap, 2))
max_mean_index <- which.max(means)
new_angle_range <- overlap_trial$angle[c(max_mean_index, max_mean_index + 1)]
angle_range <- seq(new_angle_range[1],new_angle_range[2],length.out=10)
max_overlap <- which.max(overlap_trial$overlap)
CD3_max_angle <- overlap_trial$angle[max_overlap]
ggplotly(
  ggplot(overlap_trial,aes(x=angle,y=overlap*100)) + 
    theme_spatial + 
    labs(x="Angle of rotation (degrees)",
         y="Overlap between sections (%)",
         title="Rotation optimisation") + 
    range_anno_1 + 
    geom_vline(xintercept=CD3_max_angle, 
               color="cyan", 
               linewidth=0.75) + 
    geom_point(color="white") + 
    geom_line(color="white")
  )




# Print the result
print(CD3_max_angle)
print(overlap(CD3_annotation,block))
CD3_annotation <- rotation(CD3_annotation,CD3_max_angle)
print(overlap(CD3_annotation,block))

#9. Plot all scaled, and rotated annotations ----

p6 <- ggplot() + 
  geom_polygon(data = block,aes(x=longitude,y=latitude),fill="violet",alpha=0.5) + 
  geom_polygon(data = (CD3_annotation),aes(x=longitude,y=latitude),fill="cyan",alpha=0.5) + 
  theme_spatial + 
  labs(x="Longitude",
       y="Latitude",
       title = "Block + IHC annotations (scaled + rotated)")
ggplotly(p6)
print(overlap(block,CD3_annotation)*100)
#10. Import cell locations from IHC data ----

CD3_cells <- read.csv("CD3_positive.csv",check.names = F)
colnames(CD3_cells)[7] <- "longitude"
colnames(CD3_cells)[8] <- "latitude"
CD3_cells <- select(CD3_cells,longitude,latitude)
CD3_cells$latitude <- CD3_cells$latitude * -1

#11. Transform CD3 cell coordinates using IHC transformation vectors ----
CD3_cells$longitude <- CD3_cells$longitude + vector_CD3_annotation$longitude
CD3_cells$latitude <- CD3_cells$latitude + vector_CD3_annotation$latitude
CD3_cells <- CD3_cells * CD3_scale_factor
CD3_cells <- rotation(CD3_cells,CD3_max_angle)
CD3_cells <- internal_select(CD3_annotation,CD3_cells)



p7 <- ggplot() + 
  labs(x="Longitude",
       y="Latitude",
       title = "Block + IHC annotation + CD3+ cells") + 
  geom_polygon(data = block,aes(x=longitude,y=latitude),fill="violet",alpha=0.5) + 
  geom_polygon(data = (CD3_annotation),aes(x=longitude,y=latitude),fill="cyan",alpha=0.5) + 
  geom_point(data = CD3_cells,aes(x=longitude,y=latitude),color="white",alpha=0.1) + 
  theme_spatial
ggplotly(p7)





#12. ----
CD8_cells <- read.csv("CD8_positive.csv",check.names = F)
colnames(CD8_cells)[7] <- "longitude"
colnames(CD8_cells)[8] <- "latitude"
CD8_cells <- select(CD8_cells,longitude,latitude)
CD8_cells$latitude <- CD8_cells$latitude * -1


CD8_cells$longitude <- CD8_cells$longitude + vector_CD8_annotation$longitude
CD8_cells$latitude <- CD8_cells$latitude + vector_CD8_annotation$latitude
CD8_cells <- CD8_cells * CD8_scale_factor
CD8_cells <- rotation(CD8_cells,CD8_max_angle)



CD8_cells <- internal_select(CD8_annotation,CD8_cells)

p8 <- ggplot() + 
  labs(x="Longitude",
       y="Latitude",
       title = "Block + IHC annotation + T cells") + 
  geom_polygon(data = block,aes(x=longitude,y=latitude),fill="violet",alpha=0.5) + 
  geom_polygon(data = (CD3_annotation),aes(x=longitude,y=latitude),fill="cyan",alpha=0.5) + 
  geom_point(data = CD3_cells,aes(x=longitude,y=latitude),color="white",alpha=0.1) +
  theme_spatial
p8




final_plot <- (p1 + p2) / (p5 + p6) / (p7 + p8)

#final_plot



pdf_file <- paste0("Spatial Plot.pdf")
pdf(pdf_file)
print(final_plot)
dev.off()

if (file.exists(pdf_file)) {
  cat("PDF successfully created at", pdf_file)
} else {
  cat("Failed to create PDF.")
}


# step 1 originate
# step 2 scale (4.51)

#13. Import the sampling grid annotation ----

grid <- read_delim("07199_grid.txt")
colnames(grid)[3] <- "longitude"
colnames(grid)[4] <- "latitude"
grid$latitude <- grid$latitude * -1
grid <- select(grid,longitude,latitude)
grid <- rbind(grid,grid[1,])

p8 <- p1 + geom_polygon(data = grid,aes(x=longitude,y=latitude),color="gray",fill = "gray",alpha=0.7)
ggplotly(p8)
#14. Transform the grid location using the block vectors ----

grid$longitude <- grid$longitude + vector_block$longitude
grid$latitude <- grid$latitude + vector_block$latitude

p7 <- ggplot() + 
  labs(x="Longitude",
       y="Latitude",
       title = "CD3+ cells and grid sample") + 
  geom_point(data = CD3_cells,aes(x=longitude,y=latitude),color="white",alpha=0.1) +
  theme_spatial + 
  geom_polygon(data = grid, aes(x=longitude,y=latitude),color="gray",fill = "gray",alpha=0.75) + 
  theme_spatial
ggplotly(p7)

target_cells <- CD3_cells


cells_in_grid <- point.in.polygon(
  target_cells$longitude, target_cells$latitude,
  grid$longitude, grid$latitude
)

print(paste0("cells present in grid location: ",
             sum(cells_in_grid > 0)))
print(overlap(block,grid))





