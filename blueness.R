#Blueness Quantification: Color Segmentation with Recolorize

install.packages("recolorize") # CRAN release

library(recolorize)

# load an image that comes with the package:
img <- system.file("extdata/corbetti.png", package = "recolorize")
rc <- recolorize2(img, cutoff = 10)

img <- readImage("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/1.12d_1.png")

#bins >= 4

# get list of all PNGs that come with the package:
images <- dir(system.file("extdata", package = "recolorize"),
              pattern = ".png", full.names = TRUE)
images <- dir("C:\Users\maddi\demo-grounded-sam\data_raw\segmentation_masks-butterflies", 
              pattern = ".png", full.names = TRUE)

# for every image...
for (i in 1:length(images)) {
  
  # get an initial fit with generic clustering
  init_fit <- recolorize2(images[i], method = "hist", bins = 6, cutoff = 30, plotting = FALSE)
  
  # drop small patches
  refined_fit <- thresholdRecolor(init_fit, pct = 0.01, plotting = FALSE)
  
  # store in an output variable
  if (i == 1) {
    colormap_list <- list(refined_fit)
  } else {
    colormap_list[[i]] <- refined_fit
  }
}

# compare original to recolored images:
layout(matrix(1:10, nrow = 2, byrow = TRUE))
par(mar = rep(0, 4))
o <- lapply(colormap_list, function(i) plot(i$original_img))
r <- lapply(colormap_list, function(i) plotImageArray(recoloredImage(i)))


#Chen code
library(dplyr)
img_dir <- "C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/segmentation_masks-butterflies/test"
images <- list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
rc_list <- vector("list", length = length(images))
for (i in 1:length(images)) {
  rc <- suppressMessages(recolorize2(images[i], bins = 6,
                                     cutoff = 30, plotting = FALSE))
  
  
  rc_list[[i]] <- data.frame(
    image_id = basename(images[i]),
    size = rc$sizes,
    R = round(rc$centers[,1] * 255),
    G = round(rc$centers[,2] * 255),
    B = round(rc$centers[,3] * 255)
  )
}
combined_df <- bind_rows(rc_list)
write.csv(combined_df, "C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/segmentation_masks-butterflies/test/test.csv", 
          row.names = FALSE)
