#'Takes mutation data download from CMP and discretises into a binary event matrix

#' @param mutation_data mutation data from CMP download
#' @param vaf_thresh Threshold to include mutation. Default 0.15.
#' @param gene_id_col The name of the column in mutation data to use as the gene identifier. Default ensembl_gene_id.
#' @param model_id_col The name of the column in mutation data to use as the model identifier. Default model_id.
#' @param data_types Name of the data type to use e.g. WES.
#' @param data_source Name of the data source to use. One of Sanger or Broad.
#' @param annot_rownames Optional suffix to add to all rownames e.g. mutation. Default NULL.
#' @param events Optional set of events (mutations) to include in final matrix.
#' @param filename Optional. Name of file to save matrix as in tsv format. Default NULL.
#' @return mutation Binary Event Mutation matrix
#' @export
#' @examples
#' library(data.table)
#' mutation_file_path<-"https://cog.sanger.ac.uk/cmp/download/mutations_summary_latest.csv.gz"
#' mutation_data<-fread(mutation_file_path)
#' sanger_wes<-get_mutation_mobem(mutation_data,annot_rownames = "mut")
#' broad_wgs<-get_mutation_mobem(mutation_data,data_source="Broad",data_type="WGS",annot_rownames = "mut")
#'
#'


get_mutation_mobem<-function(mutation_data,gene_id_col="ensembl_gene_id",model_id_col="model_id",vaf_thresh=0.15,
                             data_types="WES",data_source="Sanger",annot_rownames=NULL,events=NULL,filename=NULL){


  if(sum(mutation_data$data_type==data_types&mutation_data$source==data_source)>0){
    subset_data<-mutation_data[mutation_data$data_type==data_types&mutation_data$source==data_source,]
    subset_data<-subset_data[subset_data$vaf>vaf_thresh,]
    mutBEM<-reshape2::dcast(subset_data,as.formula(paste0(gene_id_col,"~",model_id_col)),value.var="vaf",fun.aggregate = mean,fill=0)
    rownames(mutBEM)<-mutBEM[,1]
    mutBEM<-mutBEM[,2:ncol(mutBEM)]

    mutBEM<-(mutBEM>=vaf_thresh)+0
    if(!is.null(events)){
      sel_events<-intersect(events,rownames(mutBEM))
      if(length(sel_events)==0){
        stop("Selected events not present in mutation matrix. Check inputs")
      }
      if(length(sel_events)<10){
        warning("Fewer than 10 selected events in mutation matrix.")
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
    warning("No mutation data for that combination of data type and source. Check inputs")
    return(NULL)
  }

}


