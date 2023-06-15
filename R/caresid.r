#' Correspondence Analysis plot with visualization of significant associations based on chi-square standardized residuals
#'
#' @description  Performs a Correspondence Analysis (CA) on a contingency table and creates a scatterplot
#' of the row and column points on the selected dimensions. Optionally, the function
#' can add segments to the plot to visualize significant associations between row and column
#' categories on the basis of positive (unadjusted) standardized residuals larger than a given threshold. The
#' segments can be optionally labelled with the corresponding residual value.\cr
#' Visit this \href{https://drive.google.com/file/d/1Z3jhiNgVk7jjhlwnNgH9aspnZxlLe_km/view?usp=share_link}{LINK} to access the package's vignette.\cr
#'
#' @param cross.tab A dataframe representing the input contingency table.
#' @param dim1 The first dimension to plot (default is 1).
#' @param dim2 The second dimension to plot (default is 2).
#' @param segments Logical. If TRUE, add segments to the plot to connect row to column
#'                 points (or viceversa) with positive standardized residuals larger than a given threshold
#'                 (default is FALSE).
#' @param category Character vector. If provided, only add segments from that/those row (or column) category(ies)
#'                 to the column (or row) categories where the corresponding standardised residuals are
#'                 positive and larger than a given threshold. If NULL (default) all the categories are considered.
#' @param mult.comp Logical. If TRUE, adjust the residuals' significance threshold for multiple comparisons using
#'                  Sidak's method (default is FALSE).
#' @param label.residuals Logical. If TRUE, the value of the positive standardised residual will be shown as a label
#'                  at the midpoint of every segment (default is FALSE).
#' @param residual.label.size Numeric. The size of the residuals' label (default is 2).
#' @param dot.size Numeric. The size of the scatterplot's points (default is 1).
#' @param dot.label.size Numeric. The size of the points' label (default is 2.5).
#' @param axis.label.size Numeric. The size of the axis labels (default is 9).
#' @param square Logical. If TRUE, set the ratio of y to x to 1 (default is FALSE).
#'
#' @return A list with two elements: \itemize{
#' \item{\code{stand.residuals} contains the unadjusted standardized residuals for all cells.}
#' \item{\code{resid.sign.thres} contains the threshold used to determine significant residuals.}
#'         }
#'
#' @details If the \code{segment} argument is \code{FALSE} (default), a regular symmetric CA biplot is rendered.
#'
#' If the \code{segment} argument is \code{TRUE}, the function adds segments to the plot to connect
#' row and column points with positive (unadjusted) standardized residuals larger than a given threshold, indicating
#' a significant association. The threshold is 1.96 if \code{mult.comp} is \code{FALSE}, and is
#' adjusted for multiple comparisons if \code{mult.comp} is \code{TRUE}.
#'
#' In the latter case, the threshold for significant residuals is calculated using the Sidak's method.
#' It is based on an adjusted 0.05 alpha level which is calculated as \code{1-(1 - 0.05)^(1/(nr*nc))},
#' where \code{nr} and \code{nc} are the number of rows and columns in the table respectively.
#' The adjusted alpha is then converted to a critical two-tailed z value (see Beasley-Schumacker 1995).\cr
#'
#' Please note, all the visualised associations (if any) are significant at least at alpha 0.05.\cr
#'
#' Optionally, the residual segments can be labelled with the corresponding residual value by setting
#' the \code{label.residuals} to \code{TRUE}.\cr
#'
#' The idea of connecting points in a CA plot based on the value of standardized residuals can serve
#' to visually highlight certain associations in your data. However, please note that while this function can help
#' visualize the associations in the contingency table, it does not replace other formal approaches for the
#' interpretation of the CA scatterplot and formal statistical tests for assessing the
#' significance and strength of the association.
#'
#' @references Beasley TM and Schumacker RE (1995), Multiple Regression Approach to Analyzing
#' Contingency Tables: Post Hoc and Planned Comparison Procedures,
#' The Journal of Experimental Education, 64(1): 86, 89.
#'
#' @import ggplot2
#' @import ggrepel
#' @import ca
#' @importFrom stats qnorm
#' @importFrom graphics par
#'
#' @export
#'
#' @examples
#'
#' # Create a toy dataset (famous Eye-color Hair-color dataset)
#'
#' mytable <- structure(list(BLACK_H = c(68, 20, 15, 5),
#' BROWN_H = c(119, 84, 54, 29),
#' RED_H = c(26, 17, 14, 14),
#' BLOND_H = c(7, 94, 10, 16)),
#' class = "data.frame",
#' row.names = c("Brown_E", "Blue_E", "Hazel_E", "Green_E"))
#'
#' # EXAMPLE 1
#' # Run the function:
#'
#' result <- caresid(mytable, segments=TRUE)
#'
#'
#' # EXAMPLE 2
#' # As above, but adjusting for multiple comparisons:
#'
#' result <- caresid(mytable, segments=TRUE, mult.comp=TRUE)
#'
#'
#' # EXAMPLE 3
#' # As in the first example, but selecting only 2 row categories;
#' # residual labels are shown:
#'
#' result <- caresid(mytable, segments=TRUE, category=c("Brown_E", "Green_E"), label.residuals=TRUE)
#'
#'
caresid <- function(cross.tab, dim1=1, dim2=2, segments = FALSE, category = NULL,
                    mult.comp = FALSE, label.residuals = FALSE, residual.label.size = 2,
                    dot.size = 1, dot.label.size = 2.5, axis.label.size = 9, square = FALSE) {

  Dim1=Dim2=ID=Type=x=y=xend=yend=label=angle=NULL

  # Perform the correspondence analysis
  ca_results <- ca(cross.tab)

  # Get the coordinates
  col_coords <- data.frame(Dim1 = ca_results$colcoord[,dim1],
                           Dim2 = ca_results$colcoord[,dim2],
                           ID = colnames(cross.tab),
                           Type = "Column")
  row_coords <- data.frame(Dim1 = ca_results$rowcoord[,dim1],
                           Dim2 = ca_results$rowcoord[,dim2],
                           ID = rownames(cross.tab),
                           Type = "Row")

  coor <- rbind(col_coords, row_coords)

  # Extract inertia and compute percentages
  summary_ca <- summary(ca_results)
  total_inertia <- sum(summary_ca$scree[,2])
  inertia1 <- summary_ca$scree[dim1,2]
  inertia2 <- summary_ca$scree[dim2,2]
  percentage1 <- round(inertia1 / total_inertia * 100, 2)
  percentage2 <- round(inertia2 / total_inertia * 100, 2)

  # Create dimension labels
  xlab <- paste0("Dimension ", dim1, " (Inertia: ", round(inertia1, 2), "; ", round(percentage1, 2), "%)")
  ylab <- paste0("Dimension ", dim2, " (Inertia: ", round(inertia2, 2), "; ", round(percentage2, 2), "%)")

  # Calculate expected frequencies and standardized residuals
  row.totals <- rowSums(cross.tab)
  col.totals <- colSums(cross.tab)
  grand.total <- sum(cross.tab)

  expected <- outer(row.totals, col.totals) / grand.total
  residuals <- (cross.tab - expected) / sqrt(expected)

  # Calculate the threshold for a residual significant at least at alpha 0.05
  if (mult.comp==TRUE) {
    alpha.adj = 1-(1 - 0.05)^(1/(nrow(cross.tab)*ncol(cross.tab)))
    resid.threshold <- qnorm(alpha.adj/2, mean = 0, sd = 1, lower.tail = F)
  } else {
    resid.threshold <- 1.96}

  segments_df <- NULL
  # Plot segments if required
  if (segments) {
    significant_pairs <- which(residuals > resid.threshold, arr.ind = TRUE)
    if(length(significant_pairs) > 0) {
      significant_pairs_df <- data.frame(row = rownames(cross.tab)[significant_pairs[,1]],
                                         col = colnames(cross.tab)[significant_pairs[,2]])

      segments_df <- data.frame(x = col_coords$Dim1[significant_pairs[,2]],
                                y = col_coords$Dim2[significant_pairs[,2]],
                                xend = row_coords$Dim1[significant_pairs[,1]],
                                yend = row_coords$Dim2[significant_pairs[,1]],
                                row = rownames(cross.tab)[significant_pairs[,1]],
                                col = colnames(cross.tab)[significant_pairs[,2]])

      if (!is.null(category)) {
        segments_df <- segments_df[(segments_df$row %in% category) | (segments_df$col %in% category),]
        significant_pairs_df <- significant_pairs_df[(significant_pairs_df$row %in% category) | (significant_pairs_df$col %in% category),]
      }
    }
  }

  plot <- ggplot() +
    geom_hline(yintercept=0, linetype="dashed", linewidth = 0.25, color = "gray") +
    geom_vline(xintercept=0, linetype="dashed", linewidth = 0.25, color = "gray") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.title = element_text(size = axis.label.size), plot.title = element_text(size = 10)) +
    labs(title="Correspondence Analysis Scatterplot", x=xlab, y=ylab) +
    scale_color_manual(values=c("blue", "red"))

  if (!is.null(segments_df)) {
    plot <- plot + geom_segment(data=segments_df, aes(x=x, y=y, xend=xend, yend=yend), linetype="dotted", linewidth = 0.25)

    if(label.residuals) {
      slopes <- with(segments_df, (yend - y) / (xend - x))
      angles <- atan(slopes) * (180 / pi)

      midpoint_x <- with(segments_df, (x + xend) / 2)
      midpoint_y <- with(segments_df, (y + yend) / 2)

      d <- 0.04
      offset_x <- d * sin(angles * pi / 180)
      offset_y <- d * cos(angles * pi / 180)

      midpoint_x <- midpoint_x + offset_x
      midpoint_y <- midpoint_y + offset_y

      residuals_df <- data.frame(x = midpoint_x,
                                 y = midpoint_y,
                                 label = round(residuals[as.matrix(significant_pairs_df)],2),
                                 angle = angles)

      plot <- plot + geom_text(data = residuals_df, aes(x = x, y = y, label = label, angle = angle),
                               size = residual.label.size, hjust = 1, color = "darkgrey")
    }
  }

  plot <- plot +
    geom_point(data=coor, aes(x=Dim1, y=Dim2, color=Type), size=dot.size) +
    geom_text_repel(data=coor, aes(x=Dim1, y=Dim2, label=ID, color=Type), size = dot.label.size)

  if (square) {
    plot <- plot + coord_fixed(ratio = 1)
  }

  print(plot)

  return(list(stand.residuals = residuals,
              resid.sign.thres = resid.threshold))
}
