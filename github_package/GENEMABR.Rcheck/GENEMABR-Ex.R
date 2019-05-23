pkgname <- "GENEMABR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('GENEMABR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("find_root_ids")
### * find_root_ids

flush(stderr()); flush(stdout())

### Name: find_root_ids
### Title: find_root_ids
### Aliases: find_root_ids

### ** Examples

find_root_ids(selected_pathways=c("GO:0005834","R-HSA-111469"))



cleanEx()
nameEx("fisher_exact_test")
### * fisher_exact_test

flush(stderr()); flush(stdout())

### Name: fisher_exact_test
### Title: fisher_exact_test
### Aliases: fisher_exact_test

### ** Examples

fetRes <- fisher_exact_test(selected_pathways=c("GO:0007250","GO:0008625"),
  gene_input=c("TRPC4AP","CDC37","TNIP1","IKBKB","NKIRAS2","NFKBIA","TIMM50",
     "RELB","TNFAIP3","NFKBIB","HSPA1A","NFKBIE","SPAG9","NFKB2","ERLIN1",
     "REL","TNIP2","TUBB6","MAP3K8"),gene_pathway_matrix=NULL)



cleanEx()
nameEx("from_id2name")
### * from_id2name

flush(stderr()); flush(stdout())

### Name: from_id2name
### Title: from_id2name
### Aliases: from_id2name

### ** Examples

from_id2name((selected_pathways=list(c("GO:0032991#GO:0044425#GO:0044464"),"R-HSA-5357801")))



cleanEx()
nameEx("get_steps")
### * get_steps

flush(stderr()); flush(stdout())

### Name: get_steps
### Title: get_steps
### Aliases: get_steps

### ** Examples

get_steps(selected_pathways=c("GO:0005834","R-HSA-111469"))



cleanEx()
nameEx("regression_selected_pathways")
### * regression_selected_pathways

flush(stderr()); flush(stdout())

### Name: regression_selected_pathways
### Title: regression_selected_pathways
### Aliases: regression_selected_pathways

### ** Examples

rspResults <- regression_selected_pathways(gene_input=c("TRPC4AP","CDC37",
  "TNIP1","IKBKB","NKIRAS2", "NFKBIA","TIMM50","RELB","TNFAIP3","NFKBIB",
  "HSPA1A","NFKBIE","SPAG9","NFKB2","ERLIN1","REL","TNIP2",
  "TUBB6","MAP3K8"),
 gene_pathway_matrix=NULL,lambda=0.007956622,alpha=0.5)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
