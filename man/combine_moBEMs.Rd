% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_moBEMs.R
\name{combine_moBEMs}
\alias{combine_moBEMs}
\title{Downloads data from provided data links and converts in binary event matrices. Returns one matrix combining all omics
Omics can include mutation, copy number, gene expression, proteomics and fusion events}
\usage{
combine_moBEMs(
  data_links,
  sanger_mut_event_id = "ensembl_gene_id",
  sanger_model_id = "model_id",
  broad_mut_event_id = "EnsemblGeneID",
  broad_model_mut_id = "ModelID",
  cn_event_id = "symbol",
  cn_group_gain = "Amplification",
  cn_group_loss = "Deletion",
  fusion_event_start = "gene_symbol_5prime",
  fusion_event_end = "gene_symbol_3prime",
  proteomics_thresh = 2.3,
  proteomics_event_id = "symbol",
  gene_expr_thresh = 2.3,
  gene_expr_event_id = "ensembl_gene_id",
  events_matrix = NULL
)
}
\arguments{
\item{data_links}{list where the names correspond to omic types. The values are links to the datasets to be downloaded from CMP. For utilising the Broad Mutation calls requires link to data, and variant conversion files etc.}

\item{sanger_mut_event_id}{Column name in sanger mutation data to use for event IDs. Default ensembl_gene_id}

\item{sanger_model_id}{Column name in sanger CMP download data to use for model IDs. Default model_id}

\item{broad_mut_event_id}{Column name in broad mutation data to use for event IDs. Default EnsemblGeneID}

\item{broad_model_mut_id}{Column name in sanger mutation data to use for model IDs. Default ModelID}

\item{cn_event_id}{Column name in copy number data to use for event IDs. Default symbol.}

\item{cn_group_gain}{Copy number category to use for generating amplification events, one of Gain or Amplification. Default Amplification}

\item{cn_group_loss}{Copy number category to use for generating loss events, one of Loss or Deletion. Default Deletion}

\item{fusion_event_start}{Column name to use for the first identifier of the fusion event. Default gene_symbol_5prime}

\item{fusion_event_end}{Column name to use for the second identifier of the fusion event. Default gene_symbol_3prime}

\item{proteomics_thresh}{Value to use for the z-score threshold. Default 2.3}

\item{proteomics_event_id}{Column name to use for the proteomics event identifier. Default symbol.}

\item{gene_expr_thresh}{Value to use for the z-score threshold. Default 2.3}

\item{gene_expr_event_id}{Column name to use for the gene expression event identifier. Default ensembl_gene_id.}

\item{events_matrix}{Optional. Matrix of events to include in final binary event matrices e.g. cancer driver genes. Each row should be an event
Columns should match event ids e.g. sanger_mut_event_id. Optionally include a column to specify events per omic. E.g. adding column "expr" with 1 for include 0 exclude will use only those events =1 for gene expression.}
}
\value{
Combined Binary Event Mutation matrix
}
\description{
Downloads data from provided data links and converts in binary event matrices. Returns one matrix combining all omics
Omics can include mutation, copy number, gene expression, proteomics and fusion events
}
\examples{
# example code
data_links<-list(mut="https://cog.sanger.ac.uk/cmp/download/mutations_summary_latest.csv.gz",
prot="https://cog.sanger.ac.uk/cmp/download/proteomics_latest.csv.gz",
fusion="https://cog.sanger.ac.uk/cmp/download/fusions_latest.csv.gz"
)
library(data.table)
omic_moBEM<-combine_moBEMs(data_links)
#example all omics available - not run due to time constraints:
#data_links<-list(mut="https://cog.sanger.ac.uk/cmp/download/mutations_summary_latest.csv.gz",
#cn="https://cog.sanger.ac.uk/cmp/download/cnv_summary_latest.csv.gz",
#prot="https://cog.sanger.ac.uk/cmp/download/proteomics_latest.csv.gz",
#fusion="https://cog.sanger.ac.uk/cmp/download/fusions_latest.csv.gz",
# expr="https://cog.sanger.ac.uk/cmp/download/rnaseq_latest.csv.gz",
# cmpmodel="https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz",
#drivermutfile="https://cog.sanger.ac.uk/cmp/download/driver_mutations_latest.csv",
# cancerdriverfile="https://cog.sanger.ac.uk/cmp/download/driver_genes_latest.csv",
# cpgmutfile="https://cog.sanger.ac.uk/cmp/download/cancer_predisposition_variants_latest.csv"
# )
# all_sanger_omic_moBEM<-combine_moBEMs(data_links)

}
