
getSeurat_all <- function(CosMx.in, data.in, default.cosmx = "RNA", default.data.in = "RNA"){
  DefaultAssay(CosMx.in) <- default.cosmx
  DefaultAssay(data.in) <- default.data.in
  
  
  if(default.data.in == "SCT"){
    data.anchors <- FindTransferAnchors(reference = CosMx.in, query = data.in, dims = 1:30, #k.anchor = 20,
                                        reference.reduction = "pca")
  }else{
    data.anchors <- FindTransferAnchors(reference = CosMx.in, query = data.in, dims = 1:30, 
                                        reference.reduction = "pca")
  }
 
  

  predictions <- TransferData(anchorset = data.anchors, refdata = CosMx.in$cell_type, dims = 1:30)
  predictions_subtype <- TransferData(anchorset = data.anchors, refdata = CosMx.in$cell_type_subcluster, dims = 1:30)
  
  
ggplot(predictions, aes(x=prediction.score.max))+
    geom_histogram()+
    ggplot(predictions_subtype, aes(x=prediction.score.max))+
    geom_histogram()
  
  predictions$predicted.id <- ifelse(predictions$prediction.score.max < 0.4,"Unassigned",predictions$predicted.id)
  predictions_subtype[which(predictions$predicted.id == "Unassigned"),"predicted.id"] <- "Unassigned"
  #table(rownames(predictions) == rownames(predictions_subtype))
  
  
  out <- data.frame(data.in@meta.data, mainClass = predictions$predicted.id, subClass = predictions_subtype$predicted.id)
  out$mainClass_fromSub <- reshape2::colsplit(out$subClass,"_", LETTERS[1:2])[,1]
  px <- ggplot(out, aes(x=mainClass_fromSub, fill=mainClass))+
    geom_bar(position = "fill")+
    scale_fill_manual(values = MetBrewer::met.brewer("Renoir", length(unique(out$mainClass))))
  print(px)
  idx <- which(out$mainClass != out$mainClass_fromSub)
  
  print(sprintf("%s cells have been set to unassigned because of incorrect parent child status", length(idx)))
  out[out$mainClass_fromSub != out$mainClass,"subClass"] <- "Unassigned"
  out
}
