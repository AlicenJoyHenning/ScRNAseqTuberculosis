import anndata as ad
import pandas as pd
import scipy 
import os

os.getcwd()
os.chdir("/home/alicen/2024/PBMC/velocyto/")

def starsolo_velocity_anndata(input_dir, project_name):
    # Load Genes and Cells identifiers
    obs = pd.read_csv(os.path.join(input_dir,'barcodes.tsv'), header = None, index_col = 0)
    
    # Remove index column name to make it compliant with the anndata format
    obs.index.name = None

    var = pd.read_csv(os.path.join(input_dir,"features.tsv"), sep='\t',names = ('gene_ids', 'feature_types'), index_col = 1)
    var.index.name = None

    from scipy import io,sparse

    spliced=scipy.sparse.csr_matrix(scipy.io.mmread(os.path.join(input_dir,"spliced.mtx")).T)
    ambiguous=scipy.sparse.csr_matrix(scipy.io.mmread(os.path.join(input_dir,"ambiguous.mtx")).T)
    unspliced=scipy.sparse.csr_matrix(scipy.io.mmread(os.path.join(input_dir,"unspliced.mtx")).T)
    adata=ad.AnnData(X=spliced,obs=obs,var=var,layers={'spliced':spliced,"ambiguous":ambiguous,"unspliced":unspliced})
    adata.var_names_make_unique()

    # adata.write_h5ad(os.path.join(output_dir)) # add this back with output_dir specified if you want the hdf5 output
    
    # calculate nuclear fraction 
    
    exon_sum = adata.layers['spliced'].sum(axis=1)
    intron_sum = adata.layers['unspliced'].sum(axis=1)
    nuclear_fraction = intron_sum/(exon_sum + intron_sum)
    
    # save to adata object 
    adata.obs["nf"] = nuclear_fraction
    
    # extract with barcodes & nuclear fraction 
    df_nf = adata.obs
    df_nf.to_csv(project_name + ".csv", index = True )
    
    return adata


# NOTE: input directory should contain barcodes.tsv, features.tsv with 3 mtx from spliced, ambigious, unspliced
#	output directory must include full name to output for that sample "SRR1234.h5ad" 




# List of project names and their corresponding input directories
projects = [
    ("/home/alicen/2024/PBMC/post_alignment/PBMC_TB_1/velocyto", "PBMC_TB_1"),
    ("/home/alicen/2024/PBMC/post_alignment/PBMC_TB_2/velocyto", "PBMC_TB_2"),
    ("/home/alicen/2024/PBMC/post_alignment/PBMC_TB_3/velocyto", "PBMC_TB_3"),
    ("/home/alicen/2024/PBMC/post_alignment/PBMC_LTBI_1/velocyto", "PBMC_LTBI_1"),
    ("/home/alicen/2024/PBMC/post_alignment/PBMC_LTBI_2/velocyto", "PBMC_LTBI_2"),
    ("/home/alicen/2024/PBMC/post_alignment/PBMC_HC_1/velocyto", "PBMC_HC_1"),
   ("/home/alicen/2024/PBMC/post_alignment/PBMC_HC_2/velocyto", "PBMC_HC_2")
]

# Loop to run each sample 
for  input_dir, project_name in projects:
    starsolo_velocity_anndata(input_dir=input_dir, project_name=project_name)
  
