#'Convert CMP expression data into matrix format for discretisation.

#' @param expr_data copy number data from CMP download
#' @param gene_id_col The name of the column in copy number data to use as the gene identifier. Default ensembl_gene_id.
#' @param model_id_col The name of the column in copy number data to use as the model identifier. Default model_id.
#' @param expr_val Value to use to populate matrix. Default rsem_tpm.
#' @param source Name of the data source to use. One of Sanger or Broad.
#' @param log_transform Should data be log transformed. Default. TRUE.
#' @param filename Name of file to save matrix as in RDS format. Only used if saveRDS =TRUE. Default ExpressionMatrix.
#' @param save_RDS Save expression matrix to RDS file. Default TRUE.
#' @return  Gene expression matrix
#' @export
#' @examples
#' all_expr<-"https://cog.sanger.ac.uk/cmp/download/rnaseq_latest.csv.gz"
#' library(data.table)
#' #not run
#' #expr_data<-fread(all_expr)
#' #expr_tpm_mat<-get_gene_expression_input(expr_data,save_RDS=FALSE,expr_val="rsem_tpm",log_transform = FALSE)
#'
#'

get_gene_expression_input<-function(expr_data,filename="ExpressionMatrix",gene_id_col="ensembl_gene_id",model_id_col="model_id",expr_val="rsem_tpm",source="Sanger",log_transform=TRUE,save_RDS=TRUE){
  if(length(expr_data$data_source==source)>1){
    subset_data<-expr_data[expr_data$data_source==source,]
    ExprMatrix<-reshape2::dcast(subset_data,as.formula(paste0(gene_id_col,"~",model_id_col)),value.var=expr_val,fun.aggregate = mean,fill=0)
    rownames(ExprMatrix)<-ExprMatrix[,1]
    ExprMatrix<-ExprMatrix[,2:ncol(ExprMatrix)]

    if(log_transform){
      ExprMatrix<-log(ExprMatrix+1,2)
    }
    if(save_RDS){
      saveRDS(ExprMatrix,file=paste0(filename,".Rds"))
    }else{
      return(ExprMatrix)
    }
  }else{
    stop("No data available for that data source. Check inputs.")
  }
}

get_zscore_gene_expr<-function(expr_data,threshold=2,events=NULL){

  #Z-scale the expression data and then discretise:
  EXPscale<-t(scale(t(expr_data)))
  EXPlow<-EXPscale<(-threshold)+0
  if(!is.null(events)){
    EXPlow<-EXPlow[intersect(rownames(EXPlow),events),]
  }
  rownames(EXPlow)<-sapply(rownames(EXPlow),function(x) paste0(x,"_GeneExpr_Low"))
  EXPhigh<-EXPscale>(threshold)+0
  if(!is.null(events)){
    EXPhigh<-EXPhigh[intersect(rownames(EXPhigh),events),]
  }

  rownames(EXPhigh)<-sapply(rownames(EXPhigh),function(x) paste0(x,"_GeneExpr_High"))
  EXPz<-rbind(EXPlow,EXPhigh)
  return(EXPz)
}

