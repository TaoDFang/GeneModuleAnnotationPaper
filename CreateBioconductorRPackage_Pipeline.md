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

  The bioconductor verion now used is 3.8 instead of 3.9, not the lastest version.
  we will for now use this setting  
  information from http://bioconductor.org/install/
  The current release of Bioconductor is version 3.8; it works with R version 3.5.3. Users of older R and Bioconductor must update their installation to take advantage of new features and to access packages that have been added to Bioconductor since the last release.

  The development version of Bioconductor is version 3.9; it works with R version 3.6.0. More recent ‘devel’ versions of R (if available) will be supported during the next Bioconductor release cycle.

## DESCRIPTION
bioViews

## NAMESPACE

## NEWS

## CITATION
do this later !!!

## Including Data
to be continue
