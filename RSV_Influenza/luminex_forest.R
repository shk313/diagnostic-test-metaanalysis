###################### Calc raw Se and Sp #########################
library(dplyr)
library(mada)
library(forestplot)
library(grid)
library(ggplot2)

file.path <- "C:/Users/SuzanneKeddie/OneDrive - London School of Hygiene and Tropical Medicine/diagnostic_test_meta_analysis/examples/RSV_Influenza/"
data <- read.csv(paste0(file.path, 'rsv_ready.csv'))

## Influenza Chan AH1 three zeros in 2x2 table

## try ordering by pathogen
data_order <- data[order(data$Pathogen),]

data_sub <- subset(data, Genotype == "A + B")
### Table of TP, FN, FP, TN
data_sub2 <- data_sub[,toupper(colnames(data_sub)) %in% c("TP","FP","TN","FN")]

sums <- rowSums(data_sub2)

test_data <- madad(data_sub2)
se_mean <- unlist(test_data$sens[1])
se_lower <- test_data$sens$sens.ci[,1]
se_upper <- test_data$sens$sens.ci[,2]
sp_mean <- unlist(test_data$spec[1])
sp_lower <- test_data$spec$spec.ci[,1]
sp_upper <- test_data$spec$spec.ci[,2]

## Save out this info?


########################

##########################################
## Forest plot

## Two plots on the same page, separate plots for Se and Sp
tabletext <- cbind(c("Study",data_sub$study_name), c("TP", data_sub$TP), c("FP", data_sub$FP), c("FN", data_sub$FN), c("TN", data_sub$TN), c("Pathogen", data_sub$Genotype))

means <- c(NA, se_mean)
lowers <- c(NA,se_lower)
uppers <- c(NA,se_upper)
means <- round(means,2)
lowers <- round(lowers,2)
uppers <- round(uppers,2)
tabledata <- structure(list(
  mean = means, lower = lowers, upper=uppers),
  .Names = c("mean","lower",'upper'),
  row.names = c(NA, -12), ## 23 if the whole dataset
  class = "data.frame")
tabledata$upper[tabledata$upper == 1] <- 1.00001 ## can't plot as upper is the same as mean

vals <- c(paste0(means,"(",lowers,",",uppers,")"))
vals[1] <- "Sensitivity(CI)"
tabletext <- cbind(tabletext, vals)
grid.newpage()
#grid.show.layout(grid.layout(1,2, widths = c(2.3,1)))
pushViewport(viewport(layout=grid.layout(1,2, widths = c(2.5,1))))  ## 2.85 for INfluenza
pushViewport(viewport(layout.pos.col = 1))
forestplot(tabletext,
           tabledata,
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=.8)),
           #hrzl_lines = gpar(lwd=1),
           zero=0.5,
           clip=c(0.5,1.0001),
           boxsize = 0.25,
           graph.pos = 7,
           xticks = c(0.5,0.75,1.0),
           #xlab = "Sensitivity",
           vertices = TRUE,
           new_page = FALSE)
popViewport()
pushViewport(viewport(layout.pos.col = 2))
means <- c(NA, sp_mean)
lowers <- c(NA,sp_lower)
uppers <- c(NA,sp_upper)
means <- round(means,2)
lowers <- round(lowers,2)
uppers <- round(uppers,2)
tabledata <- structure(list(
  mean = means, lower = lowers, upper=uppers),
  .Names = c("mean","lower",'upper'),
  row.names = c(NA, -12),
  class = "data.frame")
tabledata$upper[tabledata$upper == 1] <- 1.00001
vals <- c(paste0(means,"(",lowers,",",uppers,")"))
vals[1] <- "Specificity(CI)"

forestplot(vals,#tabletext,
           tabledata,
          txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=.8)),
            boxsize = 0.25,
          graph.pos = 1,
          zero=0.5,
            clip = c(0.5,1.0001),
            #xlab = "Specificity",
            xticks = c(0.5, 0.75, 1),
          vertices = TRUE,
            new_page = FALSE)





####### ROC scatter plot

plot(1-c1,s1,xlim=c(0,1),ylim=c(0,1))
lines(1-c_mean,s_mean,type="p",col="red")
plt<-data.frame(c=sp_mean,s=se_mean,N=sums,path = data_order$Pathogen)

plot <- ggplot(plt)+
  geom_point(aes(x=c,y=s,col=path,size=N), shape = 22)+
  scale_x_reverse("Specificity", limits = c(1,0.9),expand=c(0,0))+
  scale_y_continuous("Sensitivity", limits = c(0.7,1), expand=c(0,0))+
  scale_color_discrete("Pathogen")+
  scale_size_continuous("Sample Size")+
  theme_bw()
  #xlim(0,1)+
  #ylim(0,1))
#+scale_size(range=c(0.4,3))





