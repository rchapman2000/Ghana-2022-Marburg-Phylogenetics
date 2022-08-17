library("ggtree")
library("treeio")
library("TDbook")

# Sets the working directory to the location of the script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Parses the nexus tree produced by iqtree as a tree object.
data <- read.nexus("sequences.aln.timetree.nex")

# Parses in the annotation file as a csv file using the read.csv function()
df_annotationData <- read.csv("sample-annotations.tsv", sep = "\t", header = FALSE)
# Adds the column names to the annotation dataframe
colnames(df_annotationData) <- c("sample", "country", "name")

# Plots the tree using the ggtree package. 
# The %<+% operator ties the annotation data to the tree and allows for information
# like the shorter name and country to be displayed on the graph.
ggtree(data, mrsd= "2022-06-25") %<+% df_annotationData + geom_tiplab(aes(label=name), size=3) +
  geom_tippoint(aes(color=country), shape=16, size=2) + theme_tree() + xlim(1700, 2050) +
  theme(legend.position="bottom", legend.box.background = element_rect(color="black", size=0.5)) + 
  guides(color=guide_legend(title="Country", nrow=2))

# Saves the plot as a jpg image.
ggsave("Marburgvuris-phylogenetic-tree.jpg", height=22, width=26, dpi=300, unit="cm")
