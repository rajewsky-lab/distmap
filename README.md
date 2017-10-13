# DistMap: single cell spatial distributed mapping 

## About
This `R` package is created and maintained by Nikos Karaiskos.
DistMap can be used to spatially map single cell RNA sequencing data
by using an existing reference database of in situs.

DistMap accompanies the following publication

*The Drosophila Embryo at Single Cell Transcriptome Resolution*, <br />
[*Science* 358, 194 (2017)](http://science.sciencemag.org/content/358/6360/194)

N. Karaiskos<sup>#</sup>, P. Wahle<sup>#</sup>, J. Alles, A. Boltengagen,
S. Ayoub, C. Kocks, N. Rajewsky<sup>&</sup> and R. Zinzen<sup>&</sup>

<sup>#</sup> Contributed equally <br />
<sup>&</sup> Corresponding authors: [N. Rajewsky](mailto:rajewsky@mdc-berlin.de), [R. Zinzen](mailto:robert.zinzen@mdc-berlin.de)

[Contact the author](mailto:nikolaos.karaiskos@mdc-berlin.de) in case you've found a bug. 

## Installation
The easiest way to install `DistMap` is through `devtools`

```
library(devtools)
install_github("rajewsky-lab/DistMap")
```

## Usage
The `DistMap` object is used to store the following structures:
* `raw.data`, the raw data (e.g. UMI counts) of the experiment, provided
by the user as a matrix with genes as rows and cells as columns.
* `data` is the normalized data, provided by the user as a matrix similar
to the raw data.
* `binarized.data` is the binarized version of the single cell data 
computed via the `binarizeSingleCellData` function.
* `insitu.matrix` is the matrix of the reference database, provided by
the user, with genes as columns and positions (bins) as rows. See the
included example used in the paper.
* `geometry`, a matrix containing the cartesian coordinates of each bin
in three dimensional space. Provided by the user, bins as rows and coordinates
as columns, see `geometry.txt.` provided as an example.

The first step is to initialize the `DistMap` object
```
dm = new("DistMap",
         raw.data=raw.data,
         data=normalized.data,
         insitu.matrix=insitu.matrix,
         geometry=geometry)
```
Then the binarized single cell data is computed and the cells are mapped
onto the reference atlas
```
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dm <- mapCells(dm)
```
Once the cells have been mapped, the `DistMap` functions can be used
to compute a vISH or a gradient of a gene and visualize the expression
pattern
```
computeVISH(dm, 'sna', threshold=0.75)
computeGeneGradient(dm, 'sna')
```


