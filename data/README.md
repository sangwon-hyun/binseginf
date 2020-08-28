** Access to data
- The data for the Snijders analysis (primarily displayed in Section 4.1 and 4.2) is located within the GLAD package (see `data(snijders)`).

You can download the package using the instructors from https://www.bioconductor.org/packages/release/bioc/html/GLAD.html. Specifically, use the code below to access this data.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GLAD")
library("GLAD")
data(snijders)
```

The above code loads 15 data frames, from `gm00143` to `gm13330`. If you wish to cite these datasets for your work, please use the below reference (MLA format): Snijders, Antoine M., et al. "Assembly of microarrays for genome-wide measurement of DNA copy number." Nature genetics 29.3 (2001): 263-264.

- The data for the Botton analysis (primarily displayed in Section 4.3) is located in this GitHub repository. These files were generated in the following way. First, we describe the aCGH data. This file is directly found at https://github.com/etal/cnvkit-examples/blob/master/cell/CL_acgh.cnr. Second, we describe the sequencing data.

 -- First, download and install the CNVkit Python package (see https://github.com/etal/cnvkit). 

 -- Next, download the latest version of hg19 genome (see https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies#hg19).

 -- Then, run the following Makefiles in the command line. This will require fixing the paths in within these files so they can appropriately find your installed version of CNVkit as well as the downloaded hg19 genome.

 1. https://github.com/etal/cnvkit-examples/blob/master/Makefile
 2. https://github.com/etal/cnvkit-examples/blob/master/intervals/Makefile (You will need to grab the access-5k-mappable.hg19.bed file from the CNVkit package.)
 3. https://github.com/etal/cnvkit-examples/blob/master/compare/cnvkit-flat/Makefile

 This results in the CL_flat.cnr file, which is exactly what is shown in this GitHub repository folder. If you wish to cite these datasets for your work, please use the below reference (MLA format): Talevich, Eric, et al. "CNVkit: genome-wide copy number detection and visualization from targeted DNA sequencing." PLoS computational biology 12.4 (2016): e1004873.