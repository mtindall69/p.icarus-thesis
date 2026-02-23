#Blueness Quantification: Color Segmentation with Recolorize

install.packages("recolorize") # CRAN release

library(recolorize)

## TESTS

# # load an image that comes with the package:
# img <- system.file("extdata/corbetti.png", package = "recolorize")
# rc <- recolorize2(img, cutoff = 10)
# 
# img <- readImage("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/1.12d_1.png")
# 
# #bins >= 4
# 
# # get list of all PNGs that come with the package:
# images <- dir(system.file("extdata", package = "recolorize"),
#               pattern = ".png", full.names = TRUE)
# images <- dir("C:\Users\maddi\demo-grounded-sam\data_raw\segmentation_masks-butterflies", 
#               pattern = ".png", full.names = TRUE)
# 
# # for every image...
# for (i in 1:length(images)) {
#   
#   # get an initial fit with generic clustering
#   init_fit <- recolorize2(images[i], method = "hist", bins = 6, cutoff = 30, plotting = FALSE)
#   
#   # drop small patches
#   refined_fit <- thresholdRecolor(init_fit, pct = 0.01, plotting = FALSE)
#   
#   # store in an output variable
#   if (i == 1) {
#     colormap_list <- list(refined_fit)
#   } else {
#     colormap_list[[i]] <- refined_fit
#   }
# }
# 
# # compare original to recolored images:
# layout(matrix(1:10, nrow = 2, byrow = TRUE))
# par(mar = rep(0, 4))
# o <- lapply(colormap_list, function(i) plot(i$original_img))
# r <- lapply(colormap_list, function(i) plotImageArray(recoloredImage(i)))


# BATCH PROCESSING
library(dplyr)
# img_dir <- "C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/segmentation_masks-butterflies/test"
# images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
# rc_list <- vector("list", length = length(images))
# for (i in 1:length(images)) {
#   rc <- suppressMessages(recolorize2(images[i], bins = 6,
#                                      cutoff = 30, plotting = TRUE))
#   
#   
#   rc_list[[i]] <- data.frame(
#     image_id = basename(images[i]),
#     size = rc$sizes,
#     R = round(rc$centers[,1] * 255),
#     G = round(rc$centers[,2] * 255),
#     B = round(rc$centers[,3] * 255)
#   )
# }
# combined_df <- bind_rows(rc_list)
# write.csv(combined_df, "C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/segmentation_masks-butterflies/test/test2.csv", 
#           row.names = FALSE)
# 
# # CONVERT RGB to HSV
# install.packages("colorspace")
# library(colorspace)
# 
# rgb6bins <- read.csv("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/segmentation_masks-butterflies/test/test.csv")
# 
# # Create sRGB object (0–255 input)
# srgb <- sRGB(R = rgb6bins$R/255, G = rgb6bins$G/255, B = rgb6bins$B/255)
# 
# # Convert to HSV
# hsv_obj <- as(srgb, "HSV")
# hsv_coords <- coords(hsv_obj)     # columns: H, S, V
# 
# rgb6bins$H <- hsv_coords[, "H"]   # in degrees
# rgb6bins$S <- hsv_coords[, "S"]   # 0–1
# rgb6bins$V <- hsv_coords[, "V"]   # 0–1
# 
# write.csv(rgb6bins, "C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/segmentation_masks-butterflies/test/6bin_hsv.csv", 
#           row.names = FALSE)


# Whole dataset parallel batching
library(colorspace)
library(parallel)
library(recolorize)
rm(list = ls())

img_dir <- "C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/segmentation_masks-butterflies"
images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
n_cores <-  15  
cl <- makeCluster(n_cores)
#apply every package to the core
clusterEvalQ(cl, {
  library(recolorize)
  library(colorspace)
})

result <- parLapply(cl, images, function(img) {
  
  rc <- suppressMessages(
    recolorize2(img, bins = 6, cutoff = 30, plotting = FALSE)
  )
  
  hsv_colors <- as(RGB(rc$centers), "HSV")@coords
  size <- rc$sizes
  H <- hsv_colors[,1]
  S <- hsv_colors[,2]
  V <- hsv_colors[,3]
  
  total_size <- sum(size, na.rm = TRUE)
  
  # filter blue by using a true/false vector
  blue_index <- (H >= 170 & H <= 270) &
    (S >= 0.2) &
    (V >= 0.35)
  
  blue_size <- sum(size[blue_index])
  
  data.frame(
    image_id = basename(img),
    blue_size = blue_size,
    total_size = total_size
  )
})
stopCluster(cl)
# combine all
final_result <- do.call(rbind, result)
write.csv(final_result, "C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/segmentation_masks-butterflies/blueness.csv", 
          row.names = FALSE)


