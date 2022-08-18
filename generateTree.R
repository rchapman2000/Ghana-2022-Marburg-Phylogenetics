library("ggtree")
library("treeio")
library("TDbook")
library("ggplot2")

# Sets the working directory to the location of the script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Parses the nexus tree produced by iqtree as a tree object.
data <- read.nexus("deduped.aln.timetree.nex")

# Parses in the annotation file as a csv file using the read.csv function()
df_annotationData <- read.csv("sample-annotations.tsv", sep = "\t", header = FALSE)
# Adds the column names to the annotation dataframe
colnames(df_annotationData) <- c("sample", "country", "name")

# Plots the tree using the ggtree package. 
# The %<+% operator ties the annotation data to the tree and allows for information
# like the shorter name and country to be displayed on the graph.
ggtree(data, mrsd = "2022-06-28") %<+% df_annotationData + geom_tiplab(aes(label=name), size=3) +
  geom_tippoint(aes(color=country), shape=16, size=2) +
  theme_tree2() + xlim(1700, 2070) + hexpand(.5, direction = 1) +
  theme(legend.position=c(0.32, 0.85), legend.title.align = 0.5, 
        legend.text = element_text(size=10),
        legend.box.background = element_rect(color="black", size=0.75),
        plot.margin = margin(1, 1, 1, 1, unit = "cm")) + 
  guides(color=guide_legend(title="Country", ncol=2))


# Saves the plot as a jpg image.
ggsave("Marburgvuris-phylogenetic-tree.pdf", height=22, width=26, dpi=300, unit="cm")
