HPGLseq: multi-comparisons Differential Expression Analysis for RNA-Seq
====================================================

The purpose of this package is to facilitate the automation of multi-comparison differential expression analysis .


## Installation

To begin, install [Bioconductor](http://www.bioconductor.org/) along with a
few dependencies that HPGLseq uses:

```r
source("http://bioconductor.org/biocLite.R")
biocLite(c('limma', 'preprocessCore', 'sva'))
```

Next, use [devtools](https://github.com/hadley/devtools) to install the latest
version of HPGLseq and HPGLseq from Github:
```r
require(devtools)
install_github("cbcbSEQ", user="kokrah")
install_github("HPGLseq", user="rtemanni")
```

If all went well you should now be able to load HPGLseq:
```r
require(HPGLseq)
vignette('HPGLseqDemo', package='HPGLseq')
```
