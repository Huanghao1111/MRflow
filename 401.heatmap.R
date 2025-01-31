# Required libraries
library(readxl)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(data.table)

# Read and prepare data
hhh <- fread("data/Table1.csv") %>% as.data.frame()

# Set row names from the first column
rownames(hhh) = hhh[,1]

# Remove the first column as it's now used as row names
hhh = hhh[,-1]
hhh_transformed <- hhh
# Transform P-values to -log10 scale for better visualization
hhh_transformed[,1:4] <- -log10(as.matrix(hhh[,1:4]))
hhh = hhh_transformed

# Extract OR values for the forest plot
or_values = hhh[,5]

# Define track labels for different MR methods
track_labels <- c("IVW", "MR-Egger", "WMe", "WMo")

# Define color scale for P-values
# Colors transition from orange (#fc7a57) through white to teal (#83c7c1)
col_fun = colorRamp2(c(0, 1, 3),
                     c("#fc7a57", "white", "#83c7c1"))

# Create legend for the heatmap
lgd = Legend(
  col_fun = col_fun,
  title = "-log10(P value)",
  at = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
  labels = as.character(c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)),
  legend_height = unit(6, "cm"),
  title_gp = gpar(fontsize = 12, fontface = "bold"),
  labels_gp = gpar(fontsize = 10),
  direction = "vertical",
  legend_gp = gpar(fill = 0),
  border = TRUE,
  # Key parameter: set label position to right side
  grid_width = unit(3, "mm"),  # Adjust color block width
)

# Define color mapping for OR values
or_col_fun = colorRamp2(c(min(or_values), 1, max(or_values)),
                        c("blue", "white", "red"))

# Clear any existing plots and set up new device
if (dev.cur() > 1) dev.off()
dev.new(width = 10, height = 10)

# Initialize circular plot settings
circos.clear()
# Set gap between sectors
circos.par$gap.degree <- 40
# Set starting angle for the first sector
circos.par$start.degree <- 20
# Set track margins
circos.par$track.margin <- c(0, 0.001)

# Calculate degree for each sector
sector_degree <- (360 - 4 * circos.par$gap.degree) / 4

# Create circular heatmap
for (i in 1:4) {
  # Prepare data for current track
  data_tmp <- as.matrix(hhh[,i])
  if (i == 1) {
    rownames(data_tmp) <- rownames(hhh)
  }
  colnames(data_tmp) <- colnames(hhh)[i]

  # Draw heatmap track
  circos.heatmap(data_tmp,
                 col = col_fun,
                 rownames.side = "outside",
                 cluster = F,
                 cell.border = NA,
                 track.height = 0.04)

  # Get current track information
  current_sector <- get.cell.meta.data("sector.index")
  current_xlim <- get.cell.meta.data("xlim")
  current_ylim <- get.cell.meta.data("ylim")

  # Calculate sector positions
  sector_start <- circos.par$start.degree + (i-1) * (sector_degree + circos.par$gap.degree)
  sector_mid <- sector_start + sector_degree/2

  # Add method labels
  circos.text(
    x = current_xlim[2],
    y = mean(current_ylim),
    labels = track_labels[i],
    sector.index = current_sector,
    track.index = get.current.track.index(),
    facing = "inside",
    niceFacing = TRUE,
    adj = c(-0.5, 0),  # Adjust label position
    cex = 0.5,
    font = 1,
    family = "Times New Roman"
  )
}

# Initialize circular forest plot
n_genes <- nrow(hhh)
sector_degrees <- 360 - (circos.par$gap.degree * 4)  # Total available degrees

# Prepare data frame for forest plot
df <- data.frame(
  factors = rep("group", nrow(hhh)),
  x = 1:nrow(hhh),
  y = or_values
)

# Add OR value track (forest plot)
circos.track(
  factors = df$factors,
  x = df$x,
  y = df$y,
  panel.fun = function(x, y) {
    # Add reference line at OR = 1
    circos.lines(x, rep(1, length(x)), col = "gray", lty = 2)
    # Add points colored by OR value
    circos.points(x, y,
                  col = ifelse(y > 1, "red", "blue"),
                  pch = 16,
                  cex = abs(log2(y)) * 0.5)
  },
  track.height = 0.04,
  bg.border = "black",
)

# Draw legend
draw(lgd,
     x = unit(0.95, "npc"),
     y = unit(0.95, "npc"),
     just = c("right", "top"))

# Clear circular plot parameters
circos.clear()
