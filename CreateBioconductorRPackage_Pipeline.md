学到的重要一课，框架从开始就要搭好，别想着最后再改，真费时间
# Basically flowing link in Bioconductor website: https://bioconductor.org/developers/package-guidelines/
  http://bioconductor.org/developers/
  https://bioconductor.org/developers/how-to/buildingPackagesForBioc/
  http://bioconductor.org/developers/how-to/buildingPackagesForBioc/

# Details

## Create a new folder /Users/taofang/Documents/Bioconductor to install necessary tools
But package is stored for now in /Users/taofang/Documents/GeneModuleAnnotationPaper/code/GENEMABR

## Package developers better/should always use devel version of Bioconductor

### Using ‘bioc-devel’ during mid-October to mid-April not for now (25/04/2019)
1. install the devel version of R on macOS
   download R-3.5-branch-el-capitan-sa-x86_64.tar.gz   to Bioconductor folder and run :
   command : tar fvxz R-3.5-branch-el-capitan-sa-x86_64.tar.gz -C /
   Failed.

   Try to use pkg version instead: R-3.5-branch-el-capitan.pkg
   Succeed. Then restart to Rstudio. It will automatically use updated version R.

2. install correct verion of BiocManager by :
  if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
  BiocManager::install(version="devel")
  BiocManager::valid()

### Using ‘bioc-devel’ during mid-April to mid-October
1. use the release version of R
  use R-3.5.3.pkg

2. invoke the function install(version="devel") (from the BiocManager package):
  if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
  BiocManager::install(version = "devel")
  BiocManager::valid()              # checks for out of date packages
  failed. error : Error: Bioconductor version '3.9' requires R version '3.6'; see https://bioconductor.org/install
  reason: Package authors should develop against the version of R that will be available to users when the Bioconductor devel branch becomes the Bioconductor release branch.
  try release R verion 3.6. Failed as release version 3.6 not found

  Information from bioconductor website (https://www.bioconductor.org/developers/release-schedule/): This release(3.9) will use R-3.6.0 (“Planting of a Tree”). The official release date is schedule for Tuesday April 30th.

  The bioconductor verion now used is 3.8 instead of 3.9, not the lastest version.
  we will for now use this setting  
  information from http://bioconductor.org/install/
  The current release of Bioconductor is version 3.8; it works with R version 3.5.3. Users of older R and Bioconductor must update their installation to take advantage of new features and to access packages that have been added to Bioconductor since the last release.

  The development version of Bioconductor is version 3.9; it works with R version 3.6.0. More recent ‘devel’ versions of R (if available) will be supported during the next Bioconductor release cycle.

  Problem solved we don't find R 3.6 from UK mirrors. But its there in US mirrors :)

## DESCRIPTION
bioViews

## NAMESPACE

## NEWS

## CITATION
do this later !!!

## Including Data

Raw Data and the inst/extdata/ Directory
It is often desirable to show a workflow which involves parsing or loading of raw files. Bioconductor recommends finding existing raw data already provided in another package or the hubs, however if this is not applicable, raw data files should be included in the inst/extdata. Files of these type are often accessed utilizing system.file(). Bioconductor requires documentation on these files in an inst/script/ directory.

## automatically load R package
add command like this in vignettes
if(!require(glmnet)){
    install.packages("glmnet")
    library(glmnet)
}


## installation?
To install this package, start R (version "3.6") and enter:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scDD", version = "3.9")

also form github
https://kbroman.org/pkg_primer/pages/github.html https://cran.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html

then install by
library(devtools)
install_github("TaoDFang/GENEMABR")



## vignettes
learn some formats from: https://github.com/campbio/celda/blob/master/vignettes/DecontX-analysis.Rmd
Better to generate pdf format file, its easier to be distributed.
<!-- Need to isnatll Tex tool tinytex::install_tinytex() or http://www.tug.org/mactex/ -->

remember to add  VignetteEngine https://www.bioconductor.org/packages/devel/bioc/vignettes/BiocStyle/inst/doc/AuthoringRmdVignettes.html
latex_engine options: pdflatex, xelatex and lualatex,
at first lualatex only works for macos while pdflatex and xelatex failed after only install basic mactex.
But lualatex cause some errors in R CMD check process.
So maybe better to install full version MacTex and  use pdflatex engine 

## coding style
https://bioconductor.org/developers/how-to/coding-style/

## use external package
https://kbroman.org/pkg_primer/pages/depends.html  , http://r-pkgs.had.co.nz/namespace.html

things to be noticed here, to use other package, we need to add some thing in DESCRIPTION and NAMESPACE
file. While roxygen2 document() function will remove this automatically. So remember to add this thing after build
funtion. In this pacakge, we used add:

import(glmnet)
import(igraph)
import(Matrix)

## documentation
by roxygen2
http://r-pkgs.had.co.nz/data.html#documenting-data
see GeneModuleAnnotation.Rmd for detail


## .gitignore
hidden file(check them by ls -a )

## News

## CITATION
to be done

## register in support site
"Maintainer must register at the support site; visit https://support.bioconductor.org/accounts/signup/ ."

## unit test
http://bioconductor.org/developers/how-to/unitTesting-guidelines/
use testthat package http://r-pkgs.had.co.nz/tests.html

library(testthat)
library(GENEMABR)






## R CMD check environment
http://bioconductor.org/developers/how-to/buildingPackagesForBioc/

In linux/macOS commond line (go to folder /Users/taofang/Documents/GeneModuleAnnotationPaper/github_package):
R CMD build GENEMABR

R CMD check GENEMABR_0.99.0.tar.gz


In R:
library(BiocCheck)
pathToPkg <- file.path("/Users/taofang/Documents/GeneModuleAnnotationPaper/github_package/","GENEMABR")
BiocCheck(pathToPkg)

pass with information:
$error
character(0)

$warning
character(0)

$note
[1] "Potential intermediate files found:"                                                                                                                                             
[2] "Consider shorter lines; 140 lines (19%) are > 80 characters long."                                                                                                               
[3] "Consider multiples of 4 spaces for line indents, 90 lines(12%) are not."                                                                                                         
[4] "Cannot determine whether maintainer is subscribed to the bioc-devel mailing list (requires admin credentials).\nSubscribe here: https://stat.ethz.ch/mailman/listinfo/bioc-devel"

##  submission
https://bioconductor.org/developers/package-submission/

## keep tring and debuging until finished in the end
