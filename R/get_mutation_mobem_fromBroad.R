#'Takes mutation data download from Broads depmap.org and discretises into a binary event matrix

#' @param broad_mutation_data mutation data from CMP download
#' @param vaf_thresh Threshold to include mutation. Default 0.15.
#' @param vafCol Name of the column in the mutation data containing the VAF values. Default AF.
#' @param gene_id_col The name of the column in mutation data to use as the gene identifier. Default EnsemblGeneID.
#' @param model_id_col The name of the column in mutation data to use as the model identifier. Default ModelID.
#' @param variantCol Name of the column that contains the variant annotation e.g missense, frameshift, intronic etc
#' @param GeneCol Name of the column that corresponds to identifiers used in the cancer driver data e.g. symbol
#' @param annot_rownames Optional suffix to add to all rownames e.g. mutation. Default NULL.
#' @param drivermutfile Name of CMP link to cancer driver variants file
#' @param cancerdriverfile Name of CMP link to cancer drivers annotation file, e.g names of cancer drivers and type Oncogene, LoF, Ambiguous.
#' @param cpgmutfile Name of CMP link to cancer predisposition mutation file
#' @param events Optional. Set of events to select for final discrete matrix. Default NULL.
#' @param MutCodeFile Name of csv file containing mappings between Broad variant annotation classes and the Sanger CMP variants annotation.
#' @return mutation Binary Event Mutation matrix
#' @examples
#' library(rtracklayer)
#' drivermutfile="https://cog.sanger.ac.uk/cmp/download/driver_mutations_latest.csv"
#' cancerdriverfile="https://cog.sanger.ac.uk/cmp/download/driver_genes_latest.csv"
#' cpgmutfile="https://cog.sanger.ac.uk/cmp/download/cancer_predisposition_variants_latest.csv"
#' geneannotfile="https://cog.sanger.ac.uk/cmp/download/gene_identifiers_latest.csv.gz"
#' MutCodeFile <- system.file("data", "MutationCodes.RData", package = "GenerateMoBEMs")
#' #subset of the data download from depmap.org used for example purposes:
#' omicsfile <- system.file("data", "random_sampleOSM.RData", package = "GenerateMoBEMs")
#' BroadBinary<-get_mutation_mobem_fromBroad(omicsfile,cancerdriverfile = cancerdriverfile,drivermutfile=drivermutfile,cpgmutfile=cpgmutfile,MutCodeFile=MutCodeFile)


#' @export
#'



get_mutation_mobem_fromBroad<-function(broad_mutation_data,gene_id_col="EnsemblGeneID",model_id_col="ModelID",variantCol="VariantInfo",GeneCol="HugoSymbol",vaf_thresh=0.15,vafCol="AF",cancerdriverfile,drivermutfile,cpgmutfile,MutCodeFile,events=NULL,annot_rownames=NULL){
  cancerDrivers<-read.csv(cancerdriverfile)
  # Get the file extension
  file_extension <- tools::file_ext(broad_mutation_data)

  # Check the extension and load accordingly
  if (file_extension == "csv") {
    mut_broad <- read.csv(file = broad_mutation_data)
  } else if (file_extension == "RData") {
    load(broad_mutation_data)  # Assumes the RData file contains a variable called 'mut_broad'
  } else {
    stop("Unsupported file type. Please provide a .csv or .RData file.")
  }

  #filter for those mutations passing vaf filter:
  mut_broad<-mut_broad[mut_broad[,vafCol]>=vaf_thresh,]
  # Get the file extension
  file_extension <- tools::file_ext(MutCodeFile)

  # Check the extension and load accordingly
  if (file_extension == "csv") {
    Mutcodes <- read.csv(file = MutCodeFile)
  } else if (file_extension == "RData") {
    load(MutCodeFile)  # Assumes the RData file contains a variable called 'mut_broad'
  } else {
    stop("Unsupported file type. Please provide a .csv or .RData file.")
  }

  #for LoF or Ambiguous we use the classification rather than the specific variants
  LoF_AmbF<-mut_broad[mut_broad[,variantCol]%in%Mutcodes[Mutcodes$SangerName!="None","BroadName"],]
  cGenes<-cancerDrivers[which(cancerDrivers[,2]%in%c("LoF","ambiguous")),"symbol"]
  LoF_AmbF<-LoF_AmbF[LoF_AmbF[,GeneCol]%in%cGenes,]

  Mutmatrix2<-reshape2::dcast(LoF_AmbF,as.formula(paste0(gene_id_col,"~",model_id_col)),value.var = variantCol,fun.aggregate = length)
  rownames(Mutmatrix2)<-Mutmatrix2[,1]
  Mutmatrix2<-Mutmatrix2[,2:ncol(Mutmatrix2)]
  Mutmatrix2<-(Mutmatrix2!=0)+0
  #use exact matches for GoF mutations
  cclehg38gr<-makeGRangesFromDataFrame(mut_broad,keep.extra.columns = TRUE,seqnames.field = "Chrom",start.field = "Pos",end.field = "Pos",strand.field = "Strand",
                                       starts.in.df.are.0based = FALSE)

  drivermutations<-read.csv(drivermutfile,header=T,stringsAsFactors = FALSE)
  df<-data.frame(chrom=paste0("chr",drivermutations$chr_name),start=drivermutations$chr_start,end=drivermutations$chr_end,
                 drvrinfo=drivermutations$effect,gene=drivermutations$symbol,ref_aa=drivermutations$ref_aa,pep_coord=drivermutations$pep_coord,alt_aa=drivermutations$alt_aa_list)

  cpgenes<-read.csv(cpgmutfile,header=T,stringsAsFactors = FALSE)

  dfc<-data.frame(chrom=paste0("chr",cpgenes$chr_name),start=cpgenes$chr_start,end=cpgenes$chr_end,drvrinfo=cpgenes$effect,
                  gene=cpgenes$symbol,ref_aa=cpgenes$ref_aa,pep_coord=cpgenes$pep_coord,alt_aa=cpgenes$alt_aa_list)
  df<-rbind(df,dfc)


  ##This needs sorting ###
  df<-df[df$gene%in%cancerDrivers[,"symbol"],]

  granges<-makeGRangesFromDataFrame(df,starts.in.df.are.0based = FALSE)
  seqlevelsStyle(granges) = "UCSC"  # necessary

  Overlaps<-as.data.frame(findOverlaps(granges,cclehg38gr))

  ccleSel<-as.data.frame(cclehg38gr[Overlaps$subjectHits,],row.names=NULL)
  sangerSel<-df[Overlaps$queryHits,]


  ccle_Mcode<-sapply(ccleSel[,variantCol],function(x) Mutcodes[match(unlist(strsplit(x,"&",fixed=T)),Mutcodes$BroadName),"BroadCode"])
  sanger_Mcode<-Mutcodes[match(sangerSel$drvrinfo,Mutcodes$SangerName),"SangerCode"]
  ccleSel$CheckVar<-as.logical(unlist(sapply(1:length(ccle_Mcode),function(x) sum(ccle_Mcode[[x]]%in%sanger_Mcode[x])>0)))


  ccleUse<-ccleSel[ccleSel$CheckVar,]
  Mutmatrix<-reshape2::dcast(ccleUse,as.formula(paste0(gene_id_col,"~",model_id_col)),value.var = "width",fun.aggregate = length)
  rownames(Mutmatrix)<-Mutmatrix[,1]
  Mutmatrix<-Mutmatrix[,2:ncol(Mutmatrix)]
  Mutmatrix<-Mutmatrix!=0+0

  #Combine mutations from LoF/Amb and GoF/point mutation matrices.
  usecl<-union(colnames(Mutmatrix),colnames(Mutmatrix2))
  usegenes<-union(rownames(Mutmatrix),rownames(Mutmatrix2))
  m1<-matrix(0,nrow=length(usegenes),ncol=length(usecl),dimnames=list(usegenes,usecl))
  m2<-matrix(0,nrow=length(usegenes),ncol=length(usecl),dimnames=list(usegenes,usecl))
  m1[rownames(Mutmatrix),colnames(Mutmatrix)]<-Mutmatrix
  m2[rownames(Mutmatrix2),colnames(Mutmatrix2)]<-Mutmatrix2
  ComCCLES<-((m1+m2)!=0)+0
  if(!is.null(events)){
    sel_events<-intersect(events,rownames(ComCCLES))
    if(length(sel_events)==0){
      stop("Selected events not present in mutation matrix. Check inputs")
    }
    if(length(sel_events)<10){
      warning("Fewer than 10 selected events in mutation matrix.")
    }

    ComCCLES<-ComCCLES[sel_events,]
  }
  if(!is.null(annot_rownames)){
    rownames(ComCCLES)<-paste(rownames(ComCCLES),annot_rownames,sep="_")
  }

  return(ComCCLES)
}



