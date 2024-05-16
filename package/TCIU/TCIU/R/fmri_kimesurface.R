#' @title interactive graph object of 3D kime-series
#' @description Use \code{plotly} to display in 3D the kime-series as 2D manifolds (kimesurface) over the cartesian domain.
#'
#' @param fmridata a 4d array which contains the spatial and temporal record of fMRI result or a single real valued vector.
#' @param voxel_location a 3d array indicating the spatial location of the brain.
#' @param is.4d The default is true. If change to false, need to input a vector instead of array.
#'
#' @details The function \code{fmri_kimesurface} is display in 3D the kime-series as 2D manifolds (kimesurface) over the Cartesian domain. It helps transform the fMRI time-series data at a fixed voxel location into a kimesurface (kime-series). User can choose to provide the 4D array of the fMRI spacetime image and the voxel_location or a single time-series vector, then a 3D visualization will be shown.
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return an interactive plot in 3D kimesurface
#' @export
#'
#' @import plotly DT scales
#' @importFrom extraDistr rlaplace
#' @importFrom spatstat.core blur
#' @importFrom spatstat.geom as.im
#' 
#' @examples
#' # sample fMRI time-series vector of a single voxel
#' sample_voxel = sample[[5]]
#' fmri_kimesurface(sample_voxel, is.4d = FALSE)[[1]]
#' fmri_kimesurface(sample_voxel, is.4d = FALSE)[[2]]
#' fmri_kimesurface(sample_voxel, is.4d = FALSE)[[3]]
#' fmri_kimesurface(sample_voxel, is.4d = FALSE)[[4]]
#' 

fmri_kimesurface <- function(fmridata, 
                             voxel_location = NULL,
                             is.4d = TRUE) {
  # randomly generate 8 phi kime-phases for each of the 10 time
  phi_8_vec <- matrix(NA, ncol=10, nrow = 8)
  #if(rand_opt=="laplace"){
    for (t in 1:10) { 
      # for a given t, generate 8 new phases
      set.seed(t);
      phi_8_vec[ ,t] <-
        extraDistr::rlaplace(8,mu=0,sigma=0.5)
      # rank-order the phases for consistency
      # within the same foliation leaf
      phi_8_vec[ ,t] <- sort(phi_8_vec[ ,t])
      # force phases in [-pi: pi)
      for (i in 1:8) {
        if (phi_8_vec[i,t] < -pi) 
          phi_8_vec[i,t] <- -pi
        if (phi_8_vec[i,t] >= pi) 
          phi_8_vec[i,t] <- pi
      }
    }
  #}
  if (is.4d == TRUE & is.null(voxel_location) == FALSE ){
    Voxel = fmridata[voxel_location[1],
                     voxel_location[2],
                     voxel_location[3], ]
  }else{
    Voxel = fmridata
  }

  fMRI_ON<-Voxel[c(rep(TRUE,10),rep(FALSE,10))]
  fMRI_OFF<-Voxel[c(rep(FALSE,10),rep(TRUE,10))]
  
  # construct the 160 (time) by 3 (fesatures) DF
  df3D_ON <- data.frame(time=1:10, phi=c(phi_8_vec[1,],phi_8_vec[1,],phi_8_vec[2,],phi_8_vec[2,],phi_8_vec[3,],phi_8_vec[3,],
                                         phi_8_vec[4,],phi_8_vec[4,],phi_8_vec[5,],phi_8_vec[5,],phi_8_vec[6,],phi_8_vec[6,],
                                         phi_8_vec[7,],phi_8_vec[7,],phi_8_vec[8,],phi_8_vec[8,]), switch=c(rep(TRUE,10),rep(FALSE,10)), fMRI=Voxel)
  # dim(df3D_ON); head(df3D_ON, 15)
  
  # Convert the long DF representing fMRI_ON and fMRI_OFF from polar coordinates to Cartesian coordinates
  matrix_ON <- matrix(0, nrow = 21, ncol = 21) 
  matrix_OFF <- matrix(0, nrow = 21, ncol = 21) 
  for (t in 1:10) {
    for (p in 1:8) {
      x = 11+t*cos(phi_8_vec[p,t])
      y = 11+t*sin(phi_8_vec[p,t])
      matrix_ON[x,y]  <- fMRI_ON[(p-1)*10 +t]
      matrix_OFF[x,y] <- fMRI_OFF[(p-1)*10 +t]
    }
  }
  
  # fix the plot_ly Text Lables
  x <- vector()
  y <- vector()
  i <- 1
  for (t in 1:10) {
    for (p in 1:8) {
      x[i] = 11+t*cos(phi_8_vec[p,t])
      y[i] = 11+t*sin(phi_8_vec[p,t])
      i <- i+1
    }
  }
  
  matrix_ON_smooth <- (1/10000)*as.matrix(spatstat.explore::blur(spatstat.geom::as.im(matrix_ON), sigma=0.5))
  matrix_OFF_smooth <- (1/10000)*as.matrix(spatstat.explore::blur(spatstat.geom::as.im(matrix_OFF), sigma=0.5))
  

  xx2 <- 11 + seq(0,10,1/2) %o% cos(seq(-pi, pi, 2*pi/20))
  yy2 <- 11 + seq(0,10,1/2) %o% sin(seq(-pi, pi, 2*pi/20))
 
  #plot 2D into 3D and make the text of the diameter (time), height (r), and phase (phi)
  f <- list(family = "Courier New, monospace", size = 18, color = "black")
  x <- list(title = "k1", titlefont = f)
  y <- list(title = "k2", titlefont = f)
  z <- list(title = "fMRI Kime-series", titlefont = f)
  zd <- list(title = "fMRI Kime-ON/OFF difference", titlefont = f)
  
  # Added: Reinterpolation
  ON_transformed <- matrix(0, nrow = 21, ncol = 21)
  OFF_transformed <- matrix(0, nrow = 21, ncol = 21)
  ON_OFF_transformed <- matrix(0, nrow = 21, ncol = 21)
  
  cart_x <- as.vector(rep(seq(1,21),21))
  cart_y <- as.vector(sort(rep(seq(1,21),21)))
  cart_z_on <- as.vector(matrix_ON_smooth)
  cart_z_diff <- as.vector(matrix_ON_smooth-matrix_OFF_smooth)
  cart_z_off <- as.vector(matrix_OFF_smooth)
  
  for(i in 1:21){
    for(j in 1:21){
      #interpolations
      int_res_on <- akima::interp(cart_x,cart_y,cart_z_on,xx2[i,j],yy2[i,j])
      int_res_off <- akima::interp(cart_x,cart_y,cart_z_off,xx2[i,j],yy2[i,j])
      int_res_diff <- akima::interp(cart_x,cart_y,cart_z_diff,xx2[i,j],yy2[i,j])
      #insert data
      ON_transformed[i,j] = as.numeric(int_res_on$z)
      OFF_transformed[i,j] = as.numeric(int_res_off$z)
      ON_OFF_transformed[i,j] = as.numeric(int_res_diff$z)
      # if None set to 0
      if(is.na(ON_transformed[i,j])){
        ON_transformed[i,j] = 0
      }
      if(is.na(ON_OFF_transformed[i,j])){
        ON_OFF_transformed[i,j] = 0
      }
      if(is.na(OFF_transformed[i,j])){
        OFF_transformed[i,j] = 0
      }
    }
  }
  #Custom hovering texts
  custom_txt <- matrix(NA, nrow=21, ncol=21)
  custom_txtOFF <- matrix(NA, nrow=21, ncol=21)
  custom_txt_DIFF <- matrix(NA, nrow=21, ncol=21)
  custom_txt_ON_OFF <- matrix(NA, nrow=21, ncol=21)
  
  for (xdir in 1:21) {
     for (ydir in 1:21) {
       custom_txt[xdir,ydir] <- paste(' fMRI: ', round(ON_transformed[xdir,ydir], 3),
                      '\n time: ', round((xdir-1)/2, 0),
                      '\n phi: ', round(-pi+pi/10*(ydir-1), 2))
       custom_txtOFF[xdir,ydir] <- paste(' fMRI: ', round(OFF_transformed[xdir,ydir], 3),
                      '\n time: ', round((xdir-1)/2, 0),
                      '\n phi: ', round(-pi+pi/10*(ydir-1), 2))
       custom_txt_DIFF[xdir,ydir] <- paste(' fMRI: ', round(ON_transformed[xdir,ydir]-OFF_transformed[xdir,ydir], 3),
                      '\n time: ', round((xdir-1)/2, 0),
                      '\n phi: ', round(-pi+pi/10*(ydir-1), 2))
       # custom_txt_ON_OFF[xdir,ydir] <- paste(' fMRI: ', round(ON_OFF_transformed[xdir,ydir], 3),
       #                '\n time: ', round((xdir-1)/2, 0),
       #                '\n phi: ', round(-pi+pi/10*(ydir-1), 2))
     }
  }
  # Start Plotting
  plot1<-plot_ly(x = ~xx2, y = ~yy2, z = ~ON_transformed, type = "surface", colors=c("#FFFFFF","#0000FF"), # scatterpolar
          text = custom_txt, hoverinfo = "text", showlegend = FALSE) %>% 
    # trace the main Z-axis
    add_trace(x=11, y=11, z=0:0.15, type="scatter3d", mode="lines", 
              line = list(width = 10, color="red"), name="Space(x)", 
              hoverinfo="none", showlegend = FALSE) %>%
    layout(dragmode = "turntable", title = "ON Kime-Surface/Kime-Series at a fixed voxel location",
           scene = list(xaxis = x, yaxis = y, zaxis = z), showlegend = FALSE)
  
  # Plot OFF kime-surface
  plot2<-plot_ly(x = ~xx2, y = ~yy2, z = ~OFF_transformed, type = "surface",colors=c("#FFFFFF","#0000FF"),   # scatterpolar
          text = custom_txtOFF, hoverinfo = "text", showlegend = FALSE) %>% 
    # trace the main Z-axis
    add_trace(x=11, y=11, z=0:0.15, type="scatter3d", mode="lines", 
              line = list(width = 10, color="red"), name="Space(x)", 
              hoverinfo="none", showlegend = FALSE) %>%
    layout(dragmode = "turntable", title = "OFF Kime-Surface/Kime-Series at a fixed voxel location",
           scene = list(xaxis = x, yaxis = y, zaxis = z), showlegend = FALSE)
  
  plot3<-plot_ly(x = ~xx2, y = ~yy2, z = ~ON_transformed-OFF_transformed, type = "surface",colors = fmri_split_ab_bl(ON_transformed-OFF_transformed), #colors=c("#FFFF00","#FFFFFF","#0000FF"),  # scatterpolar
          text = custom_txt_DIFF, hoverinfo = "text", showlegend = FALSE) %>% 
    # trace the main Z-axis
    add_trace(x=11, y=11, z=-0.15:0.15, type="scatter3d", mode="lines", 
              line = list(width = 10, color="red"), name="Space(x)", 
              hoverinfo="none", showlegend = FALSE) %>%
    layout(dragmode = "turntable", title = "Difference for ON & OFF Kime-Surface/Kime-Series at a fixed voxel location",
           scene = list(xaxis = x, yaxis = y, zaxis = zd), showlegend = FALSE)
  
  return(list(df3D_ON,plot1,plot2,plot3))
}
