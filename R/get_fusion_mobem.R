#'Takes fusion data download from CMP and discretises into a binary event matrix

#' @param fusion_data from CMP download
#' @param start_symbol The name of the column in fusion_data to use as the start of fusion identifier. Default gene_symbol_5prime.
#' @param end_symbol The name of the column in fusion_data to use as the end of fusion identifier. Default gene_symbol_3prime.
#' @param model_id_col The name of the column in fusion_data to use as the model identifier. Default model_id.
#' @param annot_rownames Optional suffix to add to all row names e.g. Fusion to denote fusion events. Default NULL.
#' @param events Optional. Set of events to select for final discrete matrix. Default NULL.
#' @param filename Optional. Name of file to save matrix as in tsv format. Default NULL.
#' @return Fusion Binary Event Matrix
#' @export
#' @examples
#' all_fusion<-"https://cog.sanger.ac.uk/cmp/download/fusions_latest.csv.gz"
#' library(data.table)
#' fusion_data<-fread(all_fusion)
#' fusion<-get_fusion_mobem(fusion_data,annot_rownames = "fusion")
#'
#'
get_fusion_mobem<-function(fusion_data,start_symbol="gene_symbol_5prime",end_symbol="gene_symbol_3prime",
                           model_id_col="model_id",annot_rownames=NULL,events=NULL,filename=NULL){
  fusion_data$fusion_id<-paste(fusion_data[[start_symbol]],fusion_data[[end_symbol]],sep ="_")
  fusion_data$is_fusion<-1
  BEM<-reshape2::dcast(fusion_data,as.formula(paste0("fusion_id~",model_id_col)),value.var="is_fusion",fun.aggregate = mean,fill=0)
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




