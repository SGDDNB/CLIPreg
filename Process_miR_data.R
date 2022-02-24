# Processing of miRNA file

library(data.table)

miR=fread("C:/Users/e0545037/Downloads/Predicted_Targets_Context_Scores.default_predictions.txt/Predicted_Targets_Context_Scores.default_predictions.txt")
miR$`Predicted relative KD`[miR$`Predicted relative KD`=="NULL"]=0
miR$`Predicted relative KD`=as.numeric(miR$`Predicted relative KD`)
miR=miR[abs(miR$`Predicted relative KD`)>1]

miR_names=unique(miR$miRNA)

miR_data=list()

for (i in miR_names) {
  miR_i=miR[miR$miRNA==i,]
  miR_data[[i]]=unique(miR_i$`Gene ID`)
}

save(miR_data,file = "data/miR_data.RData")
