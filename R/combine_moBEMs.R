#'Downloads data from provided data links and converts in binary event matrices. Returns one matrix combining all omics
#'Omics can include mutation, copy number, gene expression, proteomics and fusion events

#' @param data_links list where the names correspond to omic types. The values are links to the datasets to be downloaded from CMP. For utilising the Broad Mutation calls requires link to data, and variant conversion files etc.
#' @param sanger_mut_event_id Column name in sanger mutation data to use for event IDs. Default ensembl_gene_id
#' @param sanger_model_id Column name in sanger CMP download data to use for model IDs. Default model_id
#' @param broad_mut_event_id Column name in broad mutation data to use for event IDs. Default EnsemblGeneID
#' @param broad_model_mut_id Column name in sanger mutation data to use for model IDs. Default ModelID
#' @param cn_event_id Column name in copy number data to use for event IDs. Default symbol.
#' @param cn_group_gain Copy number category to use for generating amplification events, one of Gain or Amplification. Default Amplification
#' @param cn_group_loss Copy number category to use for generating loss events, one of Loss or Deletion. Default Deletion
#' @param fusion_event_start Column name to use for the first identifier of the fusion event. Default gene_symbol_5prime
#' @param fusion_event_end Column name to use for the second identifier of the fusion event. Default gene_symbol_3prime
#' @param proteomics_thresh Value to use for the z-score threshold. Default 2.3
#' @param proteomics_event_id Column name to use for the proteomics event identifier. Default symbol.
#' @param gene_expr_thresh Value to use for the z-score threshold. Default 2.3
#' @param gene_expr_event_id Column name to use for the gene expression event identifier. Default ensembl_gene_id.
#' @param events_matrix Optional. Matrix of events to include in final binary event matrices e.g. cancer driver genes. Each row should be an event
#' Columns should match event ids e.g. sanger_mut_event_id. Optionally include a column to specify events per omic. E.g. adding column "expr" with 1 for include 0 exclude will use only those events =1 for gene expression.
#' @return Combined Binary Event Mutation matrix
#' @examples
#' # example code
#' data_links<-list(mut="https://cog.sanger.ac.uk/cmp/download/mutations_summary_latest.csv.gz",
#' prot="https://cog.sanger.ac.uk/cmp/download/proteomics_latest.csv.gz",
#' fusion="https://cog.sanger.ac.uk/cmp/download/fusions_latest.csv.gz"
#' )
#' library(data.table)

#' omic_moBEM<-combine_moBEMs(data_links)

#' #example all omics available - not run due to time constraints:
#' #data_links<-list(mut="https://cog.sanger.ac.uk/cmp/download/mutations_summary_latest.csv.gz",
#' #cn="https://cog.sanger.ac.uk/cmp/download/cnv_summary_latest.csv.gz",
#' #prot="https://cog.sanger.ac.uk/cmp/download/proteomics_latest.csv.gz",
#' #fusion="https://cog.sanger.ac.uk/cmp/download/fusions_latest.csv.gz",
#'# expr="https://cog.sanger.ac.uk/cmp/download/rnaseq_latest.csv.gz",
#'# cmpmodel="https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz",
#' #drivermutfile="https://cog.sanger.ac.uk/cmp/download/driver_mutations_latest.csv",
#'# cancerdriverfile="https://cog.sanger.ac.uk/cmp/download/driver_genes_latest.csv",
#'# cpgmutfile="https://cog.sanger.ac.uk/cmp/download/cancer_predisposition_variants_latest.csv"
#'# )
#' # all_sanger_omic_moBEM<-combine_moBEMs(data_links)
#'
#' @export
#'



combine_moBEMs<-function(data_links,
                         sanger_mut_event_id="ensembl_gene_id",sanger_model_id="model_id",
                         broad_mut_event_id="EnsemblGeneID",broad_model_mut_id="ModelID",
                         cn_event_id="symbol", cn_group_gain="Amplification",cn_group_loss="Deletion",
                         fusion_event_start="gene_symbol_5prime",fusion_event_end="gene_symbol_3prime",
                         proteomics_thresh=2.3,proteomics_event_id="symbol",
                         gene_expr_thresh=2.3,gene_expr_event_id="ensembl_gene_id",
                         events_matrix=NULL
                         ){
  omics<-names(data_links)
  omic_mat<-NULL
  for(i in omics){

    if(i == "mut"){

      #sort the subset_meta according to a set of rules.
      order_list<-rbind(c("WGS","Sanger"),
                        c("WES","Sanger"),
                        c("TGS","Sanger"),
                        c("WGS","Broad"),
                        c("WES","Broad"),
                        c("TGS","Broad"))

      mutation_data<-fread(data_links[["mut"]])
      #set the metadata to specify the order that data should be included. Only include new data if it adds models.
      mut_mat<-NULL

      for(j in 1:nrow(order_list)){
        if(!is.null(events_matrix)){
          if(sanger_mut_event_id%in%colnames(events_matrix)&"mut"%in%colnames(events_matrix)){
            events=events_matrix[events_matrix[,"mut"]==1,sanger_mut_event_id]
          }else{
            if(sanger_mut_event_id%in%colnames(events_matrix)){
              events=events_matrix[,sanger_mut_event_id]
            }else{
              events=NULL
            }
          }
        }else{
          events=NULL
        }
        temp_mat<-get_mutation_mobem(mutation_data,data_types=order_list[j,1],data_source=order_list[j,2],annot_rownames="mut",
                                     gene_id_col=sanger_mut_event_id,model_id_col=sanger_model_id,events=events)
        print(sum(is.na(temp_mat)))
        if(!is.null(temp_mat)){
          #have data for selected combination. Check if new models can be added.
          if(is.null(mut_mat)){
            mut_mat<-temp_mat
          }else{
            new_models<-setdiff(colnames(temp_mat),colnames(mut_mat))
            print(paste("New mutation models",new_models))
            if(length(new_models)>1){
              #add new models to data set
              mut_mat<-combine_matrix(mut_mat,temp_mat[,new_models])
            }
            if(length(new_models)==1){
              #add new models to data set
              mut_mat<-combine_matrix(mut_mat,matrix(temp_mat[,new_models],ncol=1),model_name=new_models)
            }
          }
        }
      }

      if("broadmut"%in%names(data_links)){
        if(!is.null(events_matrix)){
          if(broad_mut_event_id%in%colnames(events_matrix)&"mut"%in%colnames(events_matrix)){
            events=events_matrix[events_matrix[,"mut"]==1,broad_mut_event_id]
          }else{
            if(broad_mut_event_id%in%colnames(events_matrix)){
              events=events_matrix[,broad_mut_event_id]
            }else{
              events=NULL
            }
          }
        }else{
          events=NULL
        }
        #do this one slightly differently and look for broadmut in the data links list and do within the i==mut one
        BroadBinary<-get_mutation_mobem_fromBroad(data_links[["broadmut"]],cancerdriverfile = data_links[["cancerdriverfile"]],
                                                  drivermutfile=data_links[["drivermutfile"]],cpgmutfile=data_links[["cpgmutfile"]],
                                                  MutCodeFile=data_links[["MutCodeFile"]],
                                                  gene_id_col=broad_mut_event_id,model_id_col=broad_model_mut_id,events=events,annot_rownames="mut")
        print(sum(is.na(BroadBinary)))
        if(!is.null(mut_mat)){
          #combine the Sanger and Broad mutation data so that we use Broad data where Sanger missing mutation data:
          cmplist<-fread(data_links[["cmpmodel"]])

          m1<-as.data.frame(mut_mat)

          model_update<-unlist(cmplist[base::match(colnames(BroadBinary),cmplist$BROAD_ID),"model_id"])
          selcols<-!is.na(model_update)
          BroadBinary<-BroadBinary[,selcols]

          colnames(BroadBinary)<-model_update[!is.na(model_update)]
          BroadBinary<-BroadBinary[,colnames(BroadBinary)%in%cmplist$model_id]
          m2<-as.data.frame(BroadBinary)

          AddBroad<-setdiff(colnames(BroadBinary),colnames(m1))
          mut_mat<-merge(m1,m2[,AddBroad],by="row.names", all.x=TRUE)
          rownames(mut_mat)<-mut_mat$Row.names
          mut_mat <- mut_mat[, -1]
          mut_mat[is.na(mut_mat)]<-0
          mut_mat<-as.matrix(mut_mat)
          if(!is.null(mut_mat)){
            omic_mat<-mut_mat
          }
        }else{

          omic_mat<-BroadBinary
        }
      }

    }
    if(i =="cn"){

      order_list<-rbind(c("WGS","Sanger"),
                        c("WES","Sanger"),
                        c("WGS","Broad"),
                        c("WES","Broad"))

      cn_data<-fread(data_links[["cn"]])
      cn_mat<-NULL
      for(j in 1:nrow(order_list)){
        if(!is.null(events_matrix)){
          if(cn_event_id%in%colnames(events_matrix)&"cn"%in%colnames(events_matrix)){
            events=events_matrix[events_matrix[,"cn"]==1,cn_event_id]
          }else{
            if(cn_event_id%in%colnames(events_matrix)){
              events=events_matrix[,cn_event_id]
            }else{
              events=NULL
            }
          }
        }else{
          events=NULL
        }
        temp_matA<-get_copynumber_mobem(cn_data,data_types=order_list[j,1],data_source=order_list[j,2],
                                        cn_threshold=cn_group_gain,annot_rownames=cn_group_gain,
                                        gene_id_col=cn_event_id,model_id_col=sanger_model_id,events=events)
        temp_matL<-get_copynumber_mobem(cn_data,data_types=order_list[j,1],data_source=order_list[j,2],
                                        cn_threshold=cn_group_loss,annot_rownames=cn_group_loss,
                                        gene_id_col=cn_event_id,model_id_col=sanger_model_id,events=events)
        temp_mat<-combine_matrix(temp_matA,temp_matL)
        if(!is.null(temp_mat)){
          #have data for selected combination. Check if new models can be added.
          if(is.null(cn_mat)){
            cn_mat<-temp_mat
          }else{
            new_models<-setdiff(colnames(temp_mat),colnames(cn_mat))
            if(length(new_models)>1){
              #add new models to data set
              cn_mat<-combine_matrix(cn_mat,temp_mat[,new_models])
            }
            if(length(new_models)==1){
              #add new models to data set
              cn_mat<-combine_matrix(cn_mat,matrix(temp_mat[,new_models],ncol=1),model_name=new_models)
            }
          }
        }
      }

      if(!is.null(cn_mat)){
        if(!is.null(omic_mat)){
          omic_mat<-combine_matrix(omic_mat,cn_mat)}else{
            omic_mat<-cn_mat
          }
      }
    }
    if(i =="prot"){
      proteomics_data<-fread(data_links[["prot"]])
      if(!is.null(events_matrix)){
        if(proteomics_event_id%in%colnames(events_matrix)&"prot"%in%colnames(events_matrix)){
          events=events_matrix[events_matrix[,"prot"]==1,proteomics_event_id]
        }else{
          if(proteomics_event_id%in%colnames(events_matrix)){
            events=events_matrix[,proteomics_event_id]
          }else{
            events=NULL
          }
        }
      }else{
        events=NULL
      }

      prot_high<-get_proteomic_mobem(proteomics_data,z_thresh=proteomics_thresh,annot_rownames = "Prot_High",gene_id_col = proteomics_event_id,model_id_col = sanger_model_id,events=events)

      prot_low<-get_proteomic_mobem(proteomics_data,z_thresh=proteomics_thresh,direction="lower",annot_rownames = "Prot_Low",gene_id_col = proteomics_event_id,model_id_col = sanger_model_id,events=events)

      prot_mat<-combine_matrix(prot_high,prot_low)
      if(!is.null(omic_mat)){
        omic_mat<-combine_matrix(omic_mat,prot_mat)}else{
          omic_mat<-prot_mat
        }
    }
    if(i == "fusion"){
      fusion_data<-fread(data_links[["fusion"]])

      fusion<-get_fusion_mobem(fusion_data,annot_rownames = "fusion",start_symbol=fusion_event_start,end_symbol=fusion_event_end,model_id_col=sanger_model_id)
      if(!is.null(omic_mat)){
        omic_mat<-combine_matrix(omic_mat,fusion)}else{
          omic_mat<-fusion
        }
    }
    if(i =="expr"){
      expr_data<-fread(data_links[["expr"]])
      if(!is.null(events_matrix)){
        if(gene_expr_event_id%in%colnames(events_matrix)&"expr"%in%colnames(events_matrix)){
          events=events_matrix[events_matrix[,"expr"]==1,gene_expr_event_id]
        }else{
          if(proteomics_event_id%in%colnames(events_matrix)){
            events=events_matrix[,gene_expr_event_id]
          }else{
            events=NULL
          }
        }
      }else{
        events=NULL
      }

      expr_tpm_mat<-get_gene_expression_input(expr_data,save_RDS=FALSE,expr_val="rsem_tpm",log_transform = FALSE,gene_id_col = gene_expr_event_id,model_id_col = sanger_model_id)
      expr_binary<-get_zscore_gene_expr(expr_tpm_mat,threshold = gene_expr_thresh,events)
      if(!is.null(omic_mat)){
        omic_mat<-combine_matrix(omic_mat,expr_binary)}else{
          omic_mat<-expr_binary
        }
    }
  }
  return(omic_mat)

}

combine_matrix<-function(matrix1,matrix2,model_name=NULL){
  rnames<-union(rownames(matrix1),rownames(matrix2))
  if(is.null(model_name)){
    cnames<-union(colnames(matrix1),colnames(matrix2))}else{
      cnames<-c(colnames(matrix1),model_name)
      colnames(matrix2)<-model_name
    }
  out_mat<-matrix(NA,nrow=length(rnames),ncol=length(cnames))
  rownames(out_mat)<-rnames
  colnames(out_mat)<-cnames
  out_mat[rownames(matrix1),colnames(matrix1)]<-matrix1
  out_mat[rownames(matrix2),colnames(matrix2)]<-matrix2
  return(out_mat)
}

