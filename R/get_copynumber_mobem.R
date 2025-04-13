#'Takes copy number data download from CMP and discretises into a binary event matrix

#' @param cn_data copy number data from CMP download
#' @param cn_threshold Category to use for discretising, one of Gain, Amplification, Loss or Deletion.
#' @param gene_id_col The name of the column in copy number data to use as the gene identifier. Default symbol.
#' @param model_id_col The name of the column in copy number data to use as the model identifier. Default model_id.
#' @param data_types Name of the data type to use e.g. WES.
#' @param data_source Name of the data source to use. One of Sanger or Broad.
#' @param annot_rownames Optional suffix to add to all rownames e.g. Gain for copy number gain events. Default NULL.
#' @param events Optional. Set of events to select for final discrete matrix. Default NULL.
#' @param filename Optional. Name of file to save matrix as in tsv format. Default NULL.
#' @return Copy Number Binary Event Mutation matrix
#' @examples
#' driver_cn<-"https://cog.sanger.ac.uk/cmp/download/cnv_summary_latest.csv.gz"
#' library(data.table)
#' cn_data<-fread(driver_cn)
#' sanger_wes_Amp<-get_copynumber_mobem(cn_data,cn_threshold="Amplification",annot_rownames="Amplification")

#'
#' @export
#'

get_copynumber_mobem<-function(cn_data,cn_threshold=c("Gain","Amplification","Loss","Deletion"),
                               gene_id_col="symbol",model_id_col="model_id",data_types="WES",
                               data_source="Sanger",annot_rownames=NULL,events=NULL,filename=NULL){
  cn_threshold<-match.arg(cn_threshold)
  cn_data$binary_category<-0
  cn_data[cn_data$cn_category==cn_threshold,"binary_category"]<-1

  subset_data<-cn_data[cn_data$data_type==data_types&cn_data$source==data_source,]
  if(sum(cn_data$data_type==data_types&cn_data$source==data_source)>0){
    mutBEM<-reshape2::dcast(subset_data,as.formula(paste0(gene_id_col,"~",model_id_col)),value.var="binary_category",fun.aggregate = mean,fill=0)
    rownames(mutBEM)<-mutBEM[,1]
    mutBEM<-mutBEM[,2:ncol(mutBEM)]

    mutBEM<-(mutBEM!=0)+0
    if(!is.null(events)){
      sel_events<-intersect(events,rownames(mutBEM))
      if(length(sel_events)==0){
        stop("No selected events in the proteomics matrix. Check inputs")
      }
      if(length(sel_events)<10){
        warning("Fewer than 10 events in proteomics matrix.")
      }
      mutBEM<-mutBEM[sel_events,]
    }
    if(!is.null(annot_rownames)){
      rnames<-rownames(mutBEM)
      rnames<-paste(rnames,annot_rownames,sep="_")
      rownames(mutBEM)<-rnames
    }
    if(!is.null(filename)){
      write.table(mutBEM,file=paste0(filename,".tsv"),quote=F,sep="\t")
    }
    return(mutBEM)
  }else{

    warning("No copy number data for that combination of data type and source. Check inputs")
    return(NULL)
  }

}



