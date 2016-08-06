#' @title Plots of value-edged networks.
#' @description Generates a visualization of a value-edged network.
#'
#' @param sociomatrix A square numeric matrix (socimatrix) with real valued edges
#'(no NA's).
#' @param threshold The threshold for removing edges from the network in order to
#' calculate the positions for the nodes using the futcherman reingold algorithm.
#' The value is multiplied against max(abs(sociomatrix)) to determine the
#' threshold. Defaults to 0.5.
#' @param save_pdf Logical indicating whether the plot should be saved to a PDF.
#' @param pdf_name The name we would like to give to the output file. Be sure to
#' include a ".pdf" extension.
#' @param output_directory The directory where the user would like to output the
#' PDF if save_pdf == TRUE.
#' @param comparison_network An optional argument providing a second square
#' numeric matrix (socimatrix) with real valued edges '(no NA's) to be visually
#' compared to sociomatrix. The second network will be Procrustes transformed so
#' that it appears most similar without chaning hte relativel positions of
#' nodes. Defualts to NULL.
#' @param comparison_names An optional string vector of length two providing
#' titles for each of the two networks to be compared. Defaults to NULL.
#' @param seed Optional argument to set the seed for the network layout
#' algorithm so that plots look the same across multiple runs. Defaults to NULL
#' but can be a positive integer (eg. 12345).
#' @param white_background Defaults to FALSE. If TRUE, then network is plotted
#' on a white background with black lettering.
#' @param show_legend Logical indicating whether a legend with extremal edge
#' values should be shown. Defualts to TRUE.
#' @param title The title we wish to give our plot.
#' @param identical_node_positions Logical indicating whether node positions
#' should be fixed to be the same when comparing networks. Defualts to FALSE.
#' @examples
#' set.seed(12345)
#' sociomatrix <- matrix(rnorm(400,0,20),20,20)
#' colnames(sociomatrix) <- rownames(sociomatrix) <- letters[1:20]
#' plot_network(sociomatrix)
#' @export
plot_network <- function(sociomatrix,
                         threshold = 0.5,
                         save_pdf = FALSE,
                         pdf_name = "Test.pdf",
                         output_directory = "./",
                         comparison_network = NULL,
                         comparison_names = NULL,
                         seed = NULL,
                         white_background = FALSE,
                         show_legend = TRUE,
                         title = "",
                         identical_node_positions = FALSE
                         ){

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # check input
  if (class(sociomatrix) != "matrix" & class(sociomatrix) != "data.frame") {
    stop("You must provide the network as a numeric matrix.")
  }

  if (nrow(sociomatrix) != ncol(sociomatrix)) {
    stop("You must provide a square matrix.")
  }

  if (white_background) {
    # generate edge colors
    negcolors <- colorRampPalette(c('red','white'))
    poscolors <- colorRampPalette(c('white','blue'))
    negcolors <- negcolors(25)
    poscolors <- poscolors(25)
  } else {
    # generate edge colors
    negcolors <- colorRampPalette(c('red','black'))
    poscolors <- colorRampPalette(c('black','blue'))
    negcolors <- negcolors(25)
    poscolors <- poscolors(25)
  }

  # check to see if we provided a comparison network, and if so, deal with it.
  COMPARISON <- FALSE
  if (!is.null(comparison_network)) {

    if (is.null(comparison_names)) {
      comparison_names <- c("","")
    } else {
      if(length(comparison_names) != 2) {
        stop("You must provide a comparison_names object as a vector containing two strings.")
      }
    }

    # check optional input
    if (class(comparison_network) != "matrix" &
       class(comparison_network) != "data.frame") {
      stop("You must provide the network as a numeric matrix.")
    }

    if (nrow(comparison_network) != ncol(comparison_network)) {
      stop("You must provide a square matrix.")
    }

    if (nrow(comparison_network) != nrow(sociomatrix)) {
      stop("You must provide two matrices with the same dimensions.")
    }
    COMPARISON <- TRUE

    # get the network ready for plotting

    diag(comparison_network) <- 0

    # create temporary matrices that can be altered
    temp <- temp2 <- matrix(comparison_network[,],
                            nrow(comparison_network),
                            ncol(comparison_network))

    # determine the threshold for removing edges
    cutoff <- max(abs(temp))*threshold

    # remove edges
    temp[which(abs(temp ) < cutoff)] <- 0

    # create a network object using adjacency matrix with edges removed
    net3 <- igraph::graph.adjacency(temp ,mode="directed",
                                   weighted=TRUE,diag=FALSE)

    #create layout with Fuchterman Reingold
    layout_c <- igraph::layout_with_fr(net3, weights = igraph::E(net3)$weight)

    # create a second network object with the un-truncated network
    net4 <- igraph::graph.adjacency(temp2,mode="directed",
                                    weighted=TRUE,diag=FALSE)

    # get an edgelist
    edgelist_c <- igraph::get.edgelist(net4)
    # get the edge weights
    weights_c <- igraph::E(net4)$weight

    # order edgeweights from smallest absolute value to largest
    ordering <- order(abs(weights_c), decreasing = F)
    edgelist_c <- edgelist_c[ordering,]
    weights_c <- weights_c[ordering]

    # generate edge widths
    negbreaks_c <- seq(min(weights_c), 0, length.out = 26)
    posbreaks_c <- seq(0, max(weights_c), length.out = 26)
    widbreaks_c <- seq(0,max(abs(weights_c)),length.out = 50)
    widths_c <- seq(0,5,length.out = 50)
  }

  diag(sociomatrix) <- 0

  # create temporary matrices that can be altered
  temp <- temp2 <- matrix(sociomatrix[,],nrow(sociomatrix),ncol(sociomatrix))

  # determine the threshold for removing edges
  cutoff <- max(abs(temp))*threshold

  # remove edges
  temp[which(abs(temp ) < cutoff)] <- 0

  # create a network object using adjacency matrix with edges removed
  net <- igraph::graph.adjacency(temp ,mode="directed",
                                 weighted=TRUE,diag=FALSE)

  #create layout with Fuchterman Reingold
  layout <- igraph::layout_with_fr(net, weights = igraph::E(net)$weight)

  # create a second network object with the un-truncated network
  net2 <- igraph::graph.adjacency(temp2,mode="directed",
                                 weighted=TRUE,diag=FALSE)

  # get an edgelist
  edgelist <- igraph::get.edgelist(net2)
  # get the edge weights
  weights <- igraph::E(net2)$weight

  # order edgeweights from smallest absolute value to largest
  ordering <-order(abs(weights), decreasing = F)
  edgelist <- edgelist[ordering,]
  weights <- weights[ordering]

  # generate edge widths
  negbreaks <- seq(min(weights), 0, length.out = 26)
  posbreaks <- seq(0, max(weights), length.out = 26)
  widbreaks <- seq(0,max(abs(weights)),length.out = 50)
  widths <- seq(0,5,length.out = 50)

  if (COMPARISON) {
    if (identical_node_positions) {
      # if we are using the same positions, then just make the layouts identical
      layout_c <- layout
    } else {
      layout_c <- vegan::procrustes(layout, layout_c, scale = F)$Yrot
    }
  }

  ##### If we are saving a PDF
  if(save_pdf) {
    #get current working directory
    cur_directory <- getwd()
    setwd(output_directory)

    # if we are making a comparison, two plots next to eachother
    if (COMPARISON) {
      pdf(file = pdf_name, width = 24, height = 12)
      #start plot
      if (white_background) {
        par(bg = "white", mar = c(2,2,2,2), xpd=TRUE, mfrow = c(1,2))
        plot(layout,pch = 20, cex = 1, col = "white", axes = F,
             xlab = "", ylab = "", main = comparison_names[1],col.main = "black",
             xlim = c((min(layout[,1])-2), (max(layout[,1])+2)),
             ylim = c((min(layout[,2])-2), (max(layout[,2])+2)))
      } else {
        par(bg = "black", mar = c(2,2,2,2), xpd=TRUE, mfrow = c(1,2))
        plot(layout,pch = 20, cex = 1, col = "black", axes = F,
             xlab = "", ylab = "", main = comparison_names[1],col.main = "white",
             xlim = c((min(layout[,1])-2), (max(layout[,1])+2)),
             ylim = c((min(layout[,2])-2), (max(layout[,2])+2)))
      }

      # add in edges
      for(i in 1:length(weights)){
        cur1 <- layout[edgelist[i,1],]
        cur2 <- layout[edgelist[i,2],]
        curweight <- weights[i]

        # find edge color
        nf <- TRUE
        counter <- 1
        bin <- 1
        while(nf){
          if(curweight > 0){
            if(posbreaks[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }else{
            if(negbreaks[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }
          counter <- counter +1
        }

        # find edge width
        nf <- TRUE
        counter <- 1
        wid <- 1
        while(nf){
          if(widbreaks[counter] >= abs(curweight)){
            wid <- counter
            nf <- FALSE
          }
          counter <- counter +1
        }
        if(curweight > 0){
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = poscolors[bin], lwd = widths[wid])
        }else{
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = negcolors[bin], lwd = widths[wid])
        }
      }
      if (white_background) {
        text(layout,labels = rownames(sociomatrix), col = "black")
        if (show_legend) {
          legend("bottom", inset = 0, title = "Edge Values",title.col = "black",
                 legend = c(round(min(sociomatrix),2), round(max(sociomatrix),2)),
                 fill = c("red","blue"), horiz = T, bg = "white",text.col = "black",
                 box.col = "white")
        }
      } else {
        text(layout,labels = rownames(sociomatrix), col = "white")
        if (show_legend) {
          legend("bottom", inset = 0, title = "Edge Values",title.col = "white",
                 legend = c(round(min(sociomatrix),2), round(max(sociomatrix),2)),
                 fill = c("red","blue"), horiz = T, bg = "black",text.col = "white")
        }
      }


      # now for the comparison network

      if (white_background) {
        plot(layout_c ,pch = 20, cex = 1, col = "white", axes = F,
             xlab = "", ylab = "", main = comparison_names[2],col.main = "black",
             xlim = c((min(layout_c[,1]) - 2), (max(layout_c[,1]) + 2)),
             ylim = c((min(layout_c[,2]) - 2), (max(layout_c[,2]) + 2)))
      } else {
        plot(layout_c ,pch = 20, cex = 1, col = "black", axes = F,
             xlab = "", ylab = "", main = comparison_names[2],col.main = "white",
             xlim = c((min(layout_c[,1]) - 2), (max(layout_c[,1]) + 2)),
             ylim = c((min(layout_c[,2]) - 2), (max(layout_c[,2]) + 2)))
      }

      # add in edges
      for(i in 1:length(weights_c)){
        cur1 <- layout_c[edgelist_c[i,1],]
        cur2 <- layout_c[edgelist_c[i,2],]
        curweight <- weights_c[i]

        # find edge color
        nf <- TRUE
        counter <- 1
        bin <- 1
        while(nf){
          if(curweight > 0){
            if(posbreaks_c[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }else{
            if(negbreaks_c[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }
          counter <- counter +1
        }

        # find edge width
        nf <- TRUE
        counter <- 1
        wid <- 1
        while(nf){
          if(widbreaks_c[counter] >= abs(curweight)){
            wid <- counter
            nf <- FALSE
          }
          counter <- counter +1
        }
        if(curweight > 0){
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = poscolors[bin], lwd = widths_c[wid])
        }else{
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = negcolors[bin], lwd = widths_c[wid])
        }
      }
      if (white_background) {
        text(layout_c,labels = rownames(comparison_network), col = "black")
        if (show_legend) {
          legend("bottom", inset = 0, title = "Edge Values",title.col = "black",
                 legend = c(round(min(comparison_network),2),
                            round(max(comparison_network),2)),
                 fill = c("red","blue"), horiz = T, bg = "white",text.col = "black",
                 box.col = "white")
        }
      } else {
        text(layout_c,labels = rownames(comparison_network), col = "white")
        if (show_legend) {
          legend("bottom", inset = 0, title = "Edge Values",title.col = "white",
                 legend = c(round(min(comparison_network),2),
                            round(max(comparison_network),2)),
                 fill = c("red","blue"), horiz = T, bg = "black",text.col = "white")
        }
      }

      dev.off()
    } else {
      # for only a single plot
      pdf(file = pdf_name, width = 12, height = 12)
      #start plot
      if (white_background) {
        par(bg = "white", mar = c(2,2,2,2), xpd = TRUE)
        plot(layout,pch = 20, cex = 1, col = "white", axes = F,
             xlab = "", ylab = "", main = title,
             xlim = c((min(layout[,1]) - 2), (max(layout[,1]) + 2)),
             ylim = c((min(layout[,2]) - 2), (max(layout[,2]) + 2)))
      } else {
        par(bg = "black", mar = c(2,2,2,2), xpd = TRUE)
        plot(layout,pch = 20, cex = 1, col = "black", axes = F,
             xlab = "", ylab = "", main = title,
             xlim = c((min(layout[,1]) - 2), (max(layout[,1]) + 2)),
             ylim = c((min(layout[,2]) - 2), (max(layout[,2]) + 2)))
      }

      # add in edges
      for(i in 1:length(weights)){
        cur1 <- layout[edgelist[i,1],]
        cur2 <- layout[edgelist[i,2],]
        curweight <- weights[i]

        # find edge color
        nf <- TRUE
        counter <- 1
        bin <- 1
        while(nf){
          if(curweight > 0){
            if(posbreaks[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }else{
            if(negbreaks[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }
          counter <- counter +1
        }

        # find edge width
        nf <- TRUE
        counter <- 1
        wid <- 1
        while(nf){
          if(widbreaks[counter] >= abs(curweight)){
            wid <- counter
            nf <- FALSE
          }
          counter <- counter +1
        }
        if(curweight > 0){
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = poscolors[bin], lwd = widths[wid])
        }else{
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = negcolors[bin], lwd = widths[wid])
        }
      }
      if (white_background) {
        text(layout,labels = rownames(sociomatrix), col = "black")
        if (show_legend) {
          legend("bottom", inset = 0, title = "Edge Values",title.col = "black",
                 legend = c(round(min(sociomatrix),2), round(max(sociomatrix),2)),
                 fill = c("red","blue"), horiz = T, bg = "white",text.col = "black",
                 box.col = "white")
        }

      } else {
        text(layout,labels = rownames(sociomatrix), col = "white")
        if (show_legend) {
          legend("bottom", inset = 0, title = "Edge Values",title.col = "white",
                 legend = c(round(min(sociomatrix),2), round(max(sociomatrix),2)),
                 fill = c("red","blue"), horiz = T, bg = "black",text.col = "white")
        }
      }
      dev.off()

    } # end of comparison or not conditional else statement
    #reset working directory
    setwd(cur_directory)
  } else {


    # if we are making a comparison, two plots next to eachother
    if (COMPARISON) {
      #start plot
      if (white_background) {
        par(bg = "white", mar = c(2,2,2,2), xpd=TRUE, mfrow = c(1,2))
        plot(layout,pch = 20, cex = 1, col = "white", axes = F,
             xlab = "", ylab = "", main = comparison_names[1],col.main = "black",
             xlim = c((min(layout[,1])-2), (max(layout[,1])+2)),
             ylim = c((min(layout[,2])-2), (max(layout[,2])+2)))
      } else {
        par(bg = "black", mar = c(2,2,2,2), xpd=TRUE, mfrow = c(1,2))
        plot(layout,pch = 20, cex = 1, col = "black", axes = F,
             xlab = "", ylab = "", main = comparison_names[1],col.main = "white",
             xlim = c((min(layout[,1])-2), (max(layout[,1])+2)),
             ylim = c((min(layout[,2])-2), (max(layout[,2])+2)))
      }

      # add in edges
      for(i in 1:length(weights)){
        cur1 <- layout[edgelist[i,1],]
        cur2 <- layout[edgelist[i,2],]
        curweight <- weights[i]

        # find edge color
        nf <- TRUE
        counter <- 1
        bin <- 1
        while(nf){
          if(curweight > 0){
            if(posbreaks[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }else{
            if(negbreaks[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }
          counter <- counter +1
        }

        # find edge width
        nf <- TRUE
        counter <- 1
        wid <- 1
        while(nf){
          if(widbreaks[counter] >= abs(curweight)){
            wid <- counter
            nf <- FALSE
          }
          counter <- counter +1
        }
        if(curweight > 0){
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = poscolors[bin], lwd = widths[wid])
        }else{
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = negcolors[bin], lwd = widths[wid])
        }
      }
      if (white_background) {
        text(layout,labels = rownames(sociomatrix), col = "black")
        legend("bottom", inset = 0, title = "Edge Values",title.col = "black",
               legend = c(round(min(sociomatrix),2), round(max(sociomatrix),2)),
               fill = c("red","blue"), horiz = T, bg = "white",text.col = "black",
               box.col = "white")
      } else {
        text(layout,labels = rownames(sociomatrix), col = "white")
        legend("bottom", inset = 0, title = "Edge Values",title.col = "white",
               legend = c(round(min(sociomatrix),2), round(max(sociomatrix),2)),
               fill = c("red","blue"), horiz = T, bg = "black",text.col = "white")
      }

      # now for the comparison network

      if (white_background) {
        plot(layout_c ,pch = 20, cex = 1, col = "white", axes = F,
             xlab = "", ylab = "", main = comparison_names[2],col.main = "black",
             xlim = c((min(layout_c[,1]) - 2), (max(layout_c[,1]) + 2)),
             ylim = c((min(layout_c[,2]) - 2), (max(layout_c[,2]) + 2)))
      } else {
        plot(layout_c ,pch = 20, cex = 1, col = "black", axes = F,
             xlab = "", ylab = "", main = comparison_names[2],col.main = "white",
             xlim = c((min(layout_c[,1]) - 2), (max(layout_c[,1]) + 2)),
             ylim = c((min(layout_c[,2]) - 2), (max(layout_c[,2]) + 2)))
      }

      # add in edges
      for(i in 1:length(weights_c)){
        cur1 <- layout_c[edgelist_c[i,1],]
        cur2 <- layout_c[edgelist_c[i,2],]
        curweight <- weights_c[i]

        # find edge color
        nf <- TRUE
        counter <- 1
        bin <- 1
        while(nf){
          if(curweight > 0){
            if(posbreaks_c[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }else{
            if(negbreaks_c[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }
          counter <- counter +1
        }

        # find edge width
        nf <- TRUE
        counter <- 1
        wid <- 1
        while(nf){
          if(widbreaks_c[counter] >= abs(curweight)){
            wid <- counter
            nf <- FALSE
          }
          counter <- counter +1
        }
        if(curweight > 0){
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = poscolors[bin], lwd = widths_c[wid])
        }else{
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = negcolors[bin], lwd = widths_c[wid])
        }
      }

      if (white_background) {
        text(layout_c,labels = rownames(comparison_network), col = "black")
        if (show_legend) {
          legend("bottom", inset = 0, title = "Edge Values",title.col = "black",
                 legend = c(round(min(comparison_network),2),
                            round(max(comparison_network),2)),
                 fill = c("red","blue"), horiz = T, bg = "white",text.col = "black",
                 box.col = "white")
        }
      } else {
        text(layout_c,labels = rownames(comparison_network), col = "white")
        if (show_legend) {
          legend("bottom", inset = 0, title = "Edge Values",title.col = "white",
                 legend = c(round(min(comparison_network),2),
                            round(max(comparison_network),2)),
                 fill = c("red","blue"), horiz = T, bg = "black",text.col = "white")
        }
      }

    } else {
      #start plot
      if (white_background) {
        par(bg = "white", mar = c(2,2,2,2), xpd = TRUE)
        plot(layout,pch = 20, cex = 1, col = "white", axes = F,
             xlab = "", ylab = "", main = title,
             xlim = c((min(layout[,1]) - 2), (max(layout[,1]) + 2)),
             ylim = c((min(layout[,2]) - 2), (max(layout[,2]) + 2)))
      } else {
        par(bg = "black", mar = c(2,2,2,2), xpd = TRUE)
        plot(layout,pch = 20, cex = 1, col = "black", axes = F,
             xlab = "", ylab = "", main = title,
             xlim = c((min(layout[,1]) - 2), (max(layout[,1]) + 2)),
             ylim = c((min(layout[,2]) - 2), (max(layout[,2]) + 2)))
      }

      # add in edges
      for (i in 1:length(weights)) {
        cur1 <- layout[edgelist[i,1],]
        cur2 <- layout[edgelist[i,2],]
        curweight <- weights[i]

        # find edge color
        nf <- TRUE
        counter <- 1
        bin <- 1
        while(nf){
          if(curweight > 0){
            if(posbreaks[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }else{
            if(negbreaks[counter] >= curweight){
              bin <- counter
              nf <- FALSE
            }
          }
          counter <- counter +1
        }

        # find edge width
        nf <- TRUE
        counter <- 1
        wid <- 1
        while(nf){
          if(widbreaks[counter] >= abs(curweight)){
            wid <- counter
            nf <- FALSE
          }
          counter <- counter +1
        }
        if(curweight > 0){
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = poscolors[bin], lwd = widths[wid])
        }else{
          lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
                col = negcolors[bin], lwd = widths[wid])
        }
      }
      if (white_background) {
        text(layout,labels = rownames(sociomatrix), col = "black")
        if (show_legend) {
          legend("bottom", inset = 0, title = "Edge Values",title.col = "black",
                 legend = c(round(min(sociomatrix),2), round(max(sociomatrix),2)),
                 fill = c("red","blue"), horiz = T, bg = "white",text.col = "black",
                 box.col = "white")
        }
      } else {
        text(layout,labels = rownames(sociomatrix), col = "white")
        if (show_legend) {
          legend("bottom", inset = 0, title = "Edge Values",title.col = "white",
                 legend = c(round(min(sociomatrix),2), round(max(sociomatrix),2)),
                 fill = c("red","blue"), horiz = T, bg = "black",text.col = "white")
        }
      }
    }
  }

  par(bg = "white")
  # do not return anything.
}

