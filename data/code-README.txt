** Access to data
- The data for the Snijders analysis (primarily displayed in Section 4.1 and 4.2) is located within the GLAD package (see `data(snijders)`).

- The data for the Botton analysis (primarily displayed in Section 4.3) is located within our binseginfhelper GitHub repository (see https://github.com/linnykos/binseginf_cnvkit/tree/master/data). These files were generated in the following way. First, we describe the aCGH data. This file is directly found at https://github.com/etal/cnvkit-examples/blob/master/cell/CL_acgh.cnr. Second, we describe the sequencing data.

 -- First, download and install the CNVkit Python package (see https://github.com/etal/cnvkit). 

 -- Next, download the latest version of hg19 genome (see https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies#hg19).

 -- Then, run the following Makefiles in the command line. This will require fixing the paths in within these files so they can appropriately find your installed version of CNVkit as well as the downloaded hg19 genome.

 1. https://github.com/etal/cnvkit-examples/blob/master/Makefile
 2. https://github.com/etal/cnvkit-examples/blob/master/intervals/Makefile (You will need to grab the access-5k-mappable.hg19.bed file from the CNVkit package.)
 3. https://github.com/etal/cnvkit-examples/blob/master/compare/cnvkit-flat/Makefile

 This results in the CL_flat.cnr file, which is exactly what is shown in this GitHub repository folder.