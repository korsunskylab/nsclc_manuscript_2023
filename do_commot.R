#!/home/ik97/anaconda3/envs/r430/bin/Rscript
fnames = commandArgs(trailingOnly=TRUE)

library(furrr)
library(future)
library(dplyr)
library(glue)
library(data.table)

library(reticulate)
use_python('/home/ik97/anaconda3/envs/commot/bin/python')
commot = import('commot')
anndata = import('anndata')
nperm = 100L

## REVISIONS: update CD8 and CD4 names 
meta_data = fread('/n/data1/bwh/medicine/korsunsky/lab/ik97/NirMerfish/share/cells/meta_data.csv')[
    , .(cell, type_lvl2)
][
    grepl('TCF7', type_lvl2)
]
meta_data$cell = as.character(meta_data$cell)
## REVISIONS: update CD8 and CD4 names 


plan(multicore)
res = future_map(fnames, function(fname) {
    tryCatch({
        adata = anndata$read_h5ad(fname)        
        
        ## REVISIONS: update CD8 and CD4 names 
        adata$obs$fine_type = as.character(adata$obs$fine_type)
        .m = meta_data[cell %in% rownames(adata$obs)]
        i = match(.m$cell, rownames(adata$obs))
        adata$obs$fine_type[i]= .m$type_lvl2
        adata$obs$fine_type = factor(adata$obs$fine_type)
        ## REVISIONS: update CD8 and CD4 names         
        
        pathway = gsub('.*cache_(.*?)_(.*)_lr.h5ad.*', '\\1', fname)
        fname_out = gsub('cache', 'permrev', fname)
        message(glue('fname_out: {fname_out}'))
        
        if (!file.exists(fname_out)) {
            commot$tl$cluster_communication(adata, database_name='cellchat', pathway_name=pathway, clustering='fine_type', n_permutations=nperm)
            lr_pairs = adata$uns$get('commot-cellchat-info')$df_ligrec[, 1:2]    
            for (i in 1:nrow(lr_pairs)) {
                lr_pair = as.character(unlist(lr_pairs[i, ]))
                commot$tl$cluster_communication(adata, database_name='cellchat', lr_pair = lr_pair, clustering='fine_type', n_permutations=nperm)    
            }
            adata$write_h5ad(fname_out)    
            message('Successfully wrote out file')
        } else {
            message(glue('ALREADY PROCESSED FILE {fname_out}'))
            # message(glue('ALREADY PROCESSED FILE commot/perm_{pathway}_{lib}_lr.h5ad'))
        }
    }, error = function(e) {
        message(glue('FAILED ON {fname}'))
        print(e)
    })
})
