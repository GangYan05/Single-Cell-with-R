# existing reference
ref <- BlueprintEncodeData()
ref
# cell type annotation 
pred <- SingleR(test=sce, ref=ref, labels=ref$label.main)
table(pred$labels)
# plot the confidence score heatmap
plotScoreHeatmap(pred)



