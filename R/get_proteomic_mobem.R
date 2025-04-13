#'Takes proteomics data download from CMP and discretises into a binary event matrix

#' @param proteomics_data from CMP download
#' @param z_thresh Threshold to use for discretising z-score, default 2.
#' @param direction Which direction to consider for discretising. One of higher or lower, z-score > z_thresh or z-score < z_thresh respectively
#' @param gene_id_col The name of the column in proteomics_data to use as the gene identifier. Default symbol.
#' @param model_id_col The name of the column in proteomics_data to use as the model identifier. Default model_id.
#' @param annot_rownames Optional suffix to add to all rownames e.g. High for high expressed proteins. Default NULL.
#' @param events Optional. Set of events to select for final discrete matrix. Default NULL.
#' @param filename Optional. Name of file to save matrix as in tsv format. Default NULL.
#' @return Proteomic Binary Event Matrix
#' @export
#' @examples
#' library(data.table)
#' all_prot<-"https://cog.sanger.ac.uk/cmp/download/proteomics_latest.csv.gz"
#' proteomics_data<-fread(all_prot)
#' #Proteomic binary matrices
#' prot_high<-get_proteomic_mobem(proteomics_data,z_thresh=2.3,annot_rownames = "High")
#' prot_low<-get_proteomic_mobem(proteomics_data,z_thresh=c(-2.3),direction="lower",annot_rownames = "Low")
#'
get_proteomic_mobem<-function(proteomics_data,z_thresh=2,direction=c("higher","lower"),gene_id_col="symbol",
                              model_id_col="model_id",annot_rownames=NULL,events=NULL,filename=NULL){
  direction<-match.arg(direction)
  proteomics_data$binary_category<-0
  if(direction=="higher"){
    proteomics_data[proteomics_data$zscore>z_thresh,"binary_category"]<-1
  }else{
    proteomics_data[proteomics_data$zscore<z_thresh,"binary_category"]<-1
  }
  BEM<-reshape2::dcast(proteomics_data,as.formula(paste0(gene_id_col,"~",model_id_col)),value.var="binary_category",fun.aggregate = mean,fill=0)
  rownames(BEM)<-BEM[,1]
  BEM<-BEM[,2:ncol(BEM)]

  BEM<-(BEM!=0)+0
  if(!is.null(events)){
    sel_events<-intersect(events,rownames(BEM))
    if(length(sel_events)==0){
      stop("No selected events in the proteomics matrix. Check inputs")
    }
    if(length(sel_events)<10){
      warning("Fewer than 10 events in proteomics matrix.")
    }
    BEM<-BEM[sel_events,]
  }
  if(!is.null(annot_rownames)){
    rnames<-rownames(BEM)
    rnames<-paste(rnames,annot_rownames,sep="_")
    rownames(BEM)<-rnames
  }
  if(!is.null(filename)){
    write.table(BEM,file=paste0(filename,".tsv"),quote=F,sep="\t")
  }
  return(BEM)
}



