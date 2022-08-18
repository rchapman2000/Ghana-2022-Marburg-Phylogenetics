# Marburg Virus Phylogenetic Analysis - Ghana, August 2022
This repository houses intermediate data and workflow commands for the phylogenetic analysis of Marburg Virus Genomes sequenced in Ghana.

## Data Curation
Published Marburg Virus genomes were curated from the [Virus Pathogen Resource (Vipr) Database](https://www.viprbrc.org/brc/home.spg?decorator=vipr). Metadata for 103 total genomes was downloaded. The following parameters were used in the search:
  - Family: *Filoviridae*, Genus: *Marburgvirus*, Species: *Marburg marburgvirus*
  - Complete Genomes Only
  
From this metadata file, a ```metadata.tsv``` and ```sequences.fasta``` file (for input into Nextstrain), were produced using a custom python script (included in this repository). T This script filters out sequences missing a collection date and country of origin. As well, the script makes use of biopython's Entrez library to download the sequences and further metadata from genbank. To run the script, the following command was used:
```
python3 build_nextstrain.py -i Marburg-Genome-Metadata-Filtered.tsv --email EMAIL_ADDRESS
```

Additionally, genomes which were laboratory strains, Ravn virus genomes, and those that were tree outliers or fell onto the same branch as ravn virus sequences (based on a preliminary phylogenetic tree) were removed from the Initial metadata.

Finally, the sequences and metadata for the genomes produced in Ghana were added to the files.

## NextStrain Build
[Nextstrain Installation](https://docs.nextstrain.org/en/latest/install.html)

A next strain build for the sequences was constructed based on the **Zika Virus Tutorial** ([Github](https://github.com/nextstrain/zika-tutorial), [Documentation](https://docs.nextstrain.org/en/latest/tutorials/creating-a-workflow.html)). The Snakemake file and build information is present in this repository. Briefly the build workflow consits of the following steps: index, filter, align, generate tree, refine, ancestral, translate, traits, and export.

The GenBank Entry [NC_001608](https://www.ncbi.nlm.nih.gov/nuccore/NC_001608) was used as a reference genome for the nextstrain build.

The following file structure was created:
```
marburg-nextstrain/
├── config/
│   ├── auspice_config.json
│   ├── colors.tsv
│   ├── dropped_strains.txt
│   └── marburg-reference.gb
├── data/
│   ├── metadata.tsv
│   ├── sequences.fasta
│   └── countries.txt (this file is not necessary, just placed here for convenience)
└── Snakefile
```

Then the following command was run:
```
nextstrain build --cpu 8 .
```

Which produced the following final file structure:
```
marburg-nextstrain/
├── auspice/
│   └── marburg.json
├── config/
│   ├── auspice_config.json
│   ├── colors.tsv
│   ├── dropped_strains.txt
│   └── marburg-reference.gb
├── data/
│   ├── metadata.tsv
│   ├── sequences.fasta
│   └── countries.txt (this file is not necessary, just placed here for convenience)
├── results/
│   └── ...
└── Snakefile
```

The nextstrain was visualized using the following command:
```
nextstrain view auspice/
```

## Manual Phylogenetic Tree
A phylogenetic tree was also manually generated using command-line and R tools.
### Data Preparation
First, sequences deduplicated using ```dedupe.sh``` from the [BBTools Suite](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/)
```
dedupe.sh -in=InputData/sequences.fasta -out=deduped.fasta
```

Next, the deduplicated sequences were aligned using [Mafft](https://mafft.cbrc.jp/alignment/software/).
```
mafft --thread 8 --retree 2 --reorder deduped.fasta > deduped.aln
```

Then, the alignment was used to create a phylogenetic tree using [IQTree](http://www.iqtree.org/). To improve phylogenetic dating, IQTree can take collection dates for the samples in a ```dates.tsv`` file (a tab separated file containing sample name and collection date). To create this file, the following commands were run:
```
# Grabs the sample and date column from the metadata.tsv file using cut. 
# Uses tail to remove the header by grabbing only subsequent lines.
cut -f1,3 InputData/metadata.tsv | tail -n +2 > sample-dates.tsv

# Dates in the metadata file contain XX for ambiguous date information.
# IQTree does not take this format, thus it must be converted.
# The python3 shell was used to accomplish this:
python3

i = open("sample-dates.tsv", "r")
o = open("sample-dates-format.tsv", "w+")
for l in i:
    split = l.strip("\n").split("\t")
    date = split[1].replace("-XX", "")
    o.write("{0}\t{1}\n".format(split[0], date))

i.close()
o.close()
```
**NOTE:** A few dates were present in the format YYYY-DD-MM. These were manually changed to be in the format YYYY-MM-DD.

Next, IQTree was run using the General time reversible (GTR) substitution model and Ultrafast bootstrap approximation for assessing branch support. As well, the sample collection dates were supplied.
```
iqtree -s deduped.aln -m GTR -B 1000 --date sample-dates-format.tsv
```
### Visualization
The phylogenetic tree was imported into R using the [treeio package](10.18129/B9.bioc.treeio) and visualized using [ggtree](10.18129/B9.bioc.ggtree)

Initially, the sample name, country, and accession number were parsed from the ```metadata.tsv``` file:
```
cut -f1,4,6 InputData/metadata.tsv | tail -n +2 > sample-annotations.tsv
```

The country field was used to color tree tips, and the accession number was displayed as a shorter name for the sample (**Note:** Samples from Ghana did not have accession numbers, and were thus manually given names MV-22-### to be displayed on the tree)

Finally, the RScript ```generateTree.R``` (included in repository) was run to generate the tree.

## References:
1. Vipr Database: https://www.viprbrc.org/brc/home.spg?decorator=vipr
2. Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics, Volume 25, Issue 11, 1 June 2009, Pages 1422–1423, https://doi.org/10.1093/bioinformatics/btp163
3. James Hadfield, Colin Megill, Sidney M Bell, John Huddleston, Barney Potter, Charlton Callender, Pavel Sagulenko, Trevor Bedford, Richard A Neher, Nextstrain: real-time tracking of pathogen evolution, Bioinformatics, Volume 34, Issue 23, 01 December 2018, Pages 4121–4123, https://doi.org/10.1093/bioinformatics/bty407
4. Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution, 30(4), 772–780. https://doi.org/10.1093/molbev/mst010
5. Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular biology and evolution, 32(1), 268–274. https://doi.org/10.1093/molbev/msu300
6. Wang L, Lam TT, Xu S, Dai Z, Zhou L, Feng T, Guo P, Dunn CW, Jones BR, Bradley T, Zhu H, Guan Y, Jiang Y, Yu G (2020). “treeio: an R package for phylogenetic tree input and output with richly annotated and associated data.” Molecular Biology and Evolution, 37, 599-603. doi: 10.1093/molbev/msz240. 
7. Yu G, Smith D, Zhu H, Guan Y, Lam TT (2017). “ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data.” Methods in Ecology and Evolution, 8, 28-36. doi: 10.1111/2041-210X.12628, http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12628/abstract. 


