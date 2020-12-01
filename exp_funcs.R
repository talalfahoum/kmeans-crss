library(imager) 
library(tidyverse) 
library(tidymodels) 
library(sp) 
library(scales)
library(cowplot)
devtools::install_github("sharlagelfand/dmc") 
library(dmc)

# to ensure reproducibility
set.seed(2087)


# HELPER 1 - helper for process_image
tidy_color <- function(cluster) {
  centres <- tidy(cluster)
  centres <- centres %>% mutate(col = rgb(R,G,B))
  dmc_near <- map(centres$col, ~dmc(.x))
  centres <- centres %>% mutate(dmc_n = dmc_near)
  return(centres)
}

# FUNCTION 1 - process_image
process_image <- function(image_file_name, k_list) {
  ## process_image(image_file_name, k_list) takes in an image of some resolution
  ## and and list of potential cluster centre values; ultimately it returns a tibble
  ## containing various pieces of information that includes a  complete a k-means clustering 
  ## for each cluster centre value in k_list. 
  ##
  ## Input:
  ## - image_file_name: a downloaded PNG or JPEG image 
  ## - k_list: a list with each value comprising the cluster size in the clustering
  ## we will complete within the function
  ##
  ## Output:
  ## - cluster_info: a tibble of tibbles that includes the following columns:
  ## kclust - the clustering done for each centre, augmented - the augmented
  ## cluster for each kclust entry, glanced - our kclust passed through the glance function,
  ## which is important for scree plots, and tidied - the tidied clusters that also 
  ## contain  their associated RGB values and their nearest DMC thread colour information
  ##
  ## Example:
  ##   library(imager)
  ##   library(tidyverse) 
  ##   library(tidymodels)
  ##   devtools::install_github("sharlagelfand/dmc") 
  ##   library(dmc)
  ##   library(cowplot)
  ##   NOTE: for all the functions below that use functions above them, you will need
  ##   the aforementioned libraries
  ##
  ##   file_name <- "a1_image.png" # make sure image is in working directory
  ##   k_l = c(2,6,4)
  ##   cluster_info <- process_image(file_name, k_l)
  ##
  ##   this should return the aforementuioned tibble
  
  
  im <- imager::load.image(image_file_name)
  tidy_dat <- as.data.frame(im, wide = "c") %>% rename(R = c.1, G = c.2, B = c.3)
  dat <- select(tidy_dat,c(-x,-y))
  
  kclusts = tibble(k_list) %>%
    mutate(
      kclust = map(k_list, ~kmeans(x = dat, centers = .x, nstart = 4)),
      augmented = map(kclust, augment, tidy_dat),
      glanced = map(kclust, glance),
      tidied = map(kclust, ~tidy_color(.x)) # need to account for RGB and DMC
    )
  
  
  return(kclusts)
  
}




# FUNCTION 2 - scree_plot
scree_plot <- function(cluster_info) {
  ## scree_plot(cluster_info) takes in a cluster_info tibble that contains
  ## columns that are a result of a k-means clustering; these columns are 
  ## detailed above
  ##
  ## Input:
  ##  - cluster_info: a tibble of tibbles that contains information/columns
  ## about k-means clusterings for one or more cluster centres.As it relates to our
  ## scree_plot function,cluster_info harbors a valuable 'glanced' column that we
  ## will use to generate the scree plot.
  ## Output:
  ##  a ggplot object that comprises a scree plot that tries to estimate the best 
  ## cluster centre number-from 1 to 4 - for conducting a k_means clustering for our image.
  ## Example:
  ##  cluster_info <- process_image('my_image.png', c(1:4))
  ##  my_scree <- scree_plot(cluster_info)
  ##  my_scree
  
  clusterings <- cluster_info %>%
    unnest(cols = c(glanced))
  
  ggplot(clusterings, aes(k_list, tot.withinss)) + geom_line() +
    geom_point()
}

# HELPER 2 - helper for color_strips
square <- function(x, label_size) { 
  ggplot()  + 
    coord_fixed(xlim=c(0,1), ylim = c(0,1)) + theme_void() + 
    theme(plot.background = element_rect(fill = x)) + 
    geom_text(aes(0.5,0.5),label = x , size = label_size)
}


# FUNCTION 3 - color_strips
color_strips <- function(cluster_info) {
  
  ## color_strips(cluster_info) takes in a tibble that is cluster_info that harbors
  ## information about various k-means clusterings, and it generates one or more color
  ## color strips that convey hex-specific information about our cluster centres - per
  ## clustering completed.
  ## Input:
  ##  - cluster_info:  a tibble of tibbles that contains information/columns
  ## about k-means clusterings for one or more cluster centres. As it pertains to our
  ## color_strips function, cluster_info contains valuable information about the hex code
  ## designation of our cluster centre color.
  ##
  ## Output:
  ##  - this function returns one or more color strips that depict the hex 
  ## codes for the nearest DMC color to our cluster centre designated hex
  ## code for that particular cluster centre. In other words, a clustering of
  ## 2 centres would return a strip with two blocks, one for each color.
  ##  
  ## Example:
  ##  cluster_info <- process_image('my_image.png', c(1:4))
  ##  my_strip <- color_strips(cluster_info)
  ##  my_strip
  
  clusterings_ <- cluster_info %>%
    unnest(cols = c(tidied))
  
  dmc_table <- clusterings_ %>%
    unnest(cols = c(dmc_n))
  
  k_unique = unique(dmc_table$k_list)
  for (k_val in k_unique) {
    dmc_k_table = subset(dmc_table, k_list == k_val)
    dmc_col = dmc_k_table$hex
    n_col = length(dmc_col)
    
    rect_dat <- tibble(x1 = c(0:(n_col-1)), x2 = c(1:n_col), y1 = rep(0,n_col),
                       y2 =rep(1,n_col), colour = dmc_col)
    curr_plot = rect_dat %>% ggplot() + coord_fixed() + 
      geom_rect(aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=colour), color="black") +
      geom_text(aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=colour), size=24/n_col) + 
      scale_fill_manual(values = rect_dat$colour)+ theme_void() + theme(legend.position = "none")
    
    print(curr_plot)
  }
  
}

# HELPER 3 - helper for make_pattern
change_resolution <- function(image_df, x_size)
{
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ## a lower resolution image. Any non-coordinate columns in the data
  ## frame are summarized with their most common value in the larger
  ## grid cell.
  ##
  ## Input:
  ## - image_df: A data frame in wide format. The x-coordinate column MUST
  ##             be named 'x' and the y-coordinate column MUST be named 'y'.
  ##             Further columns have no naming restrictions.
  ## - x_size:   The number of cells in the x-direction. The number of cells
  ##             in the vertical direction will be computed to maintain the 
  ##             perspective. There is no guarantee that the exact number
  ##             of cells in the x-direction is x_size
  ##
  ## Output:
  ## - A data frame with the same column names as image_df, but with fewer 
  ##   entries that corresponds to the reduced resolution image.
  ##
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager') 
  ##   im <- load.image(fpath)
  ##   im_dat<- as.data.frame(im,wide = "c") %>% rename(R = c.1, G = c.2, B = c.3) %>%
  ##            select(x,y,R,G,B)
  ##   agg_image <- change_resolution(im_dat, 50)
  
  if(!require(sp)) {
    stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
  }
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }
  
  sp_dat <- image_df 
  gridded(sp_dat) = ~x+y
  
  persp = (gridparameters(sp_dat)$cells.dim[2]/gridparameters(sp_dat)$cells.dim[1])
  y_size = floor(x_size*persp)
  orig_x_size = gridparameters(sp_dat)$cells.dim[1]
  orig_y_size = gridparameters(sp_dat)$cells.dim[2]
  
  x_res = ceiling(orig_x_size/x_size)
  y_res = ceiling(orig_y_size/y_size)
  
  gt = GridTopology(c(0.5,0.5), c(x_res, y_res),
                    c(floor(orig_x_size/x_res), floor(orig_y_size/y_res)))
  SG = SpatialGrid(gt)
  agg = aggregate(sp_dat, SG, function(x) names(which.max(table(x)))[1] )
  agg@grid@cellsize <- c(1,1)
  df <- agg %>% as.data.frame %>% rename(x = s1, y = s2)  %>% select(colnames(image_df))
  
  return(df)
  
}

# FUNCTION 4 - make_pattern

make_pattern <- function(cluster_info, k, x_size, black_white=FALSE, background_color=NULL){
  
  ## make_pattern(cluster_info, k, x_size, black_white=FALSE, background_color=NULL)
  ## uses information about a certain image obtained from 'cluster_info' and generates
  ## a cross-stitch pattern whose appearance depends on our other parameters explained
  ## below
  ## 
  ## Input:
  ## - cluster_info:  A tibble containing valuable information on different k-means clusterings
  ## - k:             The chosen cluster size, 
  ## - x_size:        The (approximate) total number of possible stitches in the horizontal 
  ##                   direction of our eventual pattern, 
  ## - black_white: Whether or not we want toPrint the pattern in black and white
  ##                (TRUE) or colour (FALSE, default)
  ## - background_color: The colour of the background, which should not be stitched in the
  ##                     pattern. (Default is to not have a colour)
  ##
  ## Output:
  ## - a cross-stitch pattern of our image at hand that can be followed, complete 
  ##   with a legend that has thread colour, and a guide grid.
  ## Example:
  ##   cluster_info <- process_image('my_image.png', c(1:6))
  ##   my_pattern <- make_pattern(cluster_info, 6, black_white=TRUE)
  ##   my_pattern
  ##
  
  
  
  if (is.null(background_color)){
    
    clusterings <- cluster_info %>%
      unnest(cols = c(tidied))
    
    dmc_table <- clusterings %>%
      unnest(cols = c(dmc_n))
    
    dmc_sub <- subset(dmc_table, k_list == k)
    dmc_col = dmc_sub$hex
    dmc_name = dmc_sub$name
    
    aug_initial <- cluster_info %>%
      unnest(cols = c(augmented))
    
    augment_final <- subset(aug_initial, k_list == k)
    
    res_input <- select(augment_final,c(-tidied,-glanced, -kclust))
    
    
    im_df <- change_resolution(res_input, x_size)
    
    
    if(black_white == FALSE) { # handles the absence of a background
      ggplot(im_df, aes(x=x, y = y, fill = .cluster)) + 
        geom_point(aes(col = .cluster, shape = .cluster)) + 
        scale_colour_manual(values = dmc_col,
                            label =  dmc_name) +
        scale_y_reverse() + theme_void() +background_grid(color.major = "black",
                                                          color.minor = "black")
    } else if(black_white == TRUE){
      bw_vec = rep(c("#000000"), k)
      ggplot(im_df, aes(x=x, y = y, fill = .cluster)) + 
        geom_point(aes(col = .cluster, shape = .cluster)) + 
        scale_colour_manual(values = bw_vec,
                            label =  dmc_name) +
        scale_y_reverse() + theme_void() +background_grid(color.major = "black",
                                                          color.minor = "black")
    } 
    
    
  } else if(not(is.null(background_color))) { # handles the presence of a background
    
    clusterings <- cluster_info %>%
      unnest(cols = c(tidied))
    
    dmc_table2 = clusterings %>%
      unnest(cols = c(dmc_n))
    dmc_sub2 = subset(dmc_table2, k_list == k)
    dmc_new = subset(dmc_sub2, col!= background_color)
    dmc_col2 = dmc_new$hex
    dmc_name2 = dmc_new$name
    
    row_background = which(grepl(background_color, dmc_sub2$col))
    clust_remove = dmc_sub2$cluster[row_background]
    
    aug_initial <- cluster_info %>%
      unnest(cols = c(augmented))
    
    augment_next <- subset(aug_initial, k_list == k)
    
    augment_background <- subset(augment_next, augment_next$.cluster!= clust_remove)
    
    res_input <- select(augment_background,c(-tidied,-glanced, -kclust))
    
    
    im_df <- change_resolution(res_input, x_size)
    
    
    if(black_white == FALSE) {
      ggplot(im_df, aes(x=x, y = y, fill = .cluster)) + 
        geom_point(aes(col = .cluster, shape = .cluster)) + 
        scale_colour_manual(values = dmc_col2,
                            label =  dmc_name2) +
        scale_y_reverse() + theme_void() +background_grid(color.major = "black",
                                                          color.minor = "black")
    } else if(black_white == TRUE){
      bw_vec = rep(c("#000000"), k)
      ggplot(im_df, aes(x=x, y = y, fill = .cluster)) + 
        geom_point(aes(col = .cluster, shape = .cluster)) + 
        scale_colour_manual(values = bw_vec,
                            label =  dmc_name2) +
        scale_y_reverse() + theme_void() +background_grid(color.major = "black",
                                                          color.minor = "black")
      
      
    }
    
  }
  
  
}
