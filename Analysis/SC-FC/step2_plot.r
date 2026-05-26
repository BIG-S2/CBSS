cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/fMRI_coupling/PopulationAnalysis
module load r/4.4.0;R
library(tidyverse)
library(pals)         # for jet colormap
library(scales)       # for label_number
library(ggplot2)
METHOD='CBSS_Line_mask_' #''
fileT=Sys.glob(paste0('G360_',METHOD,'out/*_metrics.txt'))
temp1=gsub(paste0('G360_',METHOD,'out/'),'',gsub('_metrics.txt','',fileT))
file1=unlist(lapply(strsplit(temp1,'_'),'[',1))
file2=unlist(lapply(strsplit(temp1,'_'),'[',2))
net0=unique(file1);net0=net0[c(1,2,4,3,5:9,11:12,10)]
net00=c("AUD","CON","DAN","DMN","FPN","LAN","OAN","PMN","SMN","Vis1","Vis2","VMN")
M=matrix(NA,length(net00),length(net00));rownames(M)=colnames(M)=net0
for(ii in 1:length(fileT)){
filetemp=read.table(fileT[ii])
M[file1[ii],file2[ii]]=M[file2[ii],file1[ii]]=round(filetemp$avgcorr_demean,2)
}
rownames(M)=colnames(M)=net00
heatM<-function(M,out='SC2FC_acc.pdf',lmin=0,lmax=1){
	dfM <- as.data.frame(M)
	dfM$a <- rownames(dfM)
	df_long <- dfM %>%
	  pivot_longer(-a, names_to = "b", values_to = "FiberCount")

	nets <- rownames(M)

	gg <- ggplot(df_long, aes(x = a, y = b, fill = FiberCount)) +
	  geom_tile(color = "white", linewidth = 0.3) +
	  geom_text(aes(label = round(FiberCount, 2)), size = 7, color = "black") +
	  scale_fill_gradientn(
		colors = pals::jet(256)[-c(1:30)],
		limits = c(lmin,lmax),                 # ← fixed range
		oob = scales::squish,               # clamp values outside [0, 0.3]
		na.value = "black",
		name = "Value",
		labels = scales::label_number()
	  ) +
	  coord_fixed() +
	  labs(x = "", y = "", title = "Fiber counts between network pairs") +
	  theme_minimal(base_size = 12) +
	  theme(
		axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 18),
		axis.text.y = element_text(size = 18),
		panel.grid = element_blank(),
		plot.title = element_text(face = "bold", size = 20),
		legend.position = "right",
		legend.text  = element_text(size = 16),
		legend.title = element_text(size = 18)
	  )
	ggsave(out, gg, width = 10, height = 9, dpi = 300)
}
heatM(M,out=paste0(METHOD,'SC2FC_acc.pdf'),lmin=0,lmax=0.25)



#
fileT=Sys.glob(paste0('G360_',METHOD,'out/*_metrics.txt'))
temp1=gsub(paste0('G360_',METHOD,'out/'),'',gsub('_metrics.txt','',fileT))
file1=unlist(lapply(strsplit(temp1,'_'),'[',1))
file2=unlist(lapply(strsplit(temp1,'_'),'[',2))
net0=unique(file1);net0=net0[c(1,2,4,3,5:9,11:12,10)]
net00=c("AUD","CON","DAN","DMN","FPN","LAN","OAN","PMN","SMN","Vis1","Vis2","VMN")
M=matrix(NA,length(net00),length(net00));rownames(M)=colnames(M)=net0
for(ii in 1:length(fileT)){
filetemp=read.table(fileT[ii])
M[file1[ii],file2[ii]]=M[file2[ii],file1[ii]]=round(filetemp$avgrank,2)
}
rownames(M)=colnames(M)=net00
heatM(M,out=paste0(METHOD,'SC2FC_rank.pdf'),lmin=0.5,lmax=0.9)

####################################################################################
##Calculate weight
####################################################################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/fMRI_coupling/PopulationAnalysis
module load r/4.4.0;R
library(data.table)
B_test=fread('../../../PCA_ELReg/FA_UKB_Phase3_proj/Mean.txt',data.table=F)
rownames(B_test)=B_test[,1];B_test=B_test[,-1]
fileT=Sys.glob(paste0('G360_out/*_pred.txt'))
filetemp=fread(fileT[1],data.table=F)
UID=intersect(filetemp[,1],rownames(B_test))
B_test=B_test[UID,]
B_test <- scale(B_test)
temp1=gsub('G360_out/','',gsub('_pred.txt','',fileT))
file1=unlist(lapply(strsplit(temp1,'_'),'[',1))
file2=unlist(lapply(strsplit(temp1,'_'),'[',2))
net0=unique(file1);net0=net0[c(1,2,4,3,5:9,11:12,10)]
net00=c("AUD","CON","DAN","DMN","FPN","LAN","OAN","PMN","SMN","Vis1","Vis2","VMN")
M=matrix(NA,dim(B_test)[2],length(file1));rownames(M)=1:dim(B_test)[2];colnames(M)=paste0(file1,'_',file2)
M1=M
for(ii in 1:length(fileT)){
print(length(fileT)-ii)
filetemp=fread(fileT[ii],data.table=F)
rownames(filetemp)=filetemp[,1];filetemp=filetemp[,-1];filetemp=filetemp[UID,]
filetemp <-scale(filetemp)
cortemp=abs(t(B_test)%*%filetemp/dim(B_test)[1])
M[,paste0(file1[ii],'_',file2[ii])]=apply(cortemp,1,mean)
M1[,paste0(file1[ii],'_',file2[ii])]=apply(cortemp,1,max)
}
AA=read.csv('/overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/template/Fiber_atlas_100points_FullName_qc_LenCurv.csv')
AAM=as.data.frame(cbind(AA,M))
colnames(AAM)=c(colnames(AA),colnames(M))
AAM1=as.data.frame(cbind(AA,M1))
colnames(AAM1)=c(colnames(AA),colnames(M1))
write.csv(AAM,file='SC2FC_MeanWeight.csv',quote=F,row.names=F)
write.csv(AAM1,file='SC2FC_MaxWeight.csv',quote=F,row.names=F)
########################################
#Plot
########################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/fMRI_coupling/PopulationAnalysis
module load r/4.4.0;R
library(data.table)
library(tidyverse)
library(circlize)
library(RColorBrewer)
df <- read_csv("SC2FC_MaxWeight.csv", show_col_types = FALSE)
#df <- read_csv("SC2FC_MeanWeight.csv", show_col_types = FALSE)
df=df[,c(8,9,11,18,24,32,37,45,52,56,60,70,72,83,95)]
colnames(df)=c('From','To','Network_Pair','AUD','CON','DMN','DAN','FPN','LAN','OAN','PMN','SMN','VMN','Vis1','Vis2')
df_long <- df %>%
  select(Network_Pair, where(is.numeric)) %>%    # keep only numeric FC columns
  pivot_longer(-Network_Pair, 
               names_to = "FC_Network", 
               values_to = "Weight") %>%
  group_by(Network_Pair, FC_Network) %>%
  summarise(AvgWeight = mean(Weight, na.rm = TRUE), .groups = "drop")
df_long=as.data.frame(df_long)
df_long[,1]=gsub('Visual','Vis',gsub('Somatomotor','SMN',gsub('Ventral-Multimodal','VMN',df_long[,1])))
df_long[,1]=gsub('Auditory','AUD',gsub('Language','LAN',gsub('Posterior-Multimodal','PMN',df_long[,1])))
df_long[,1]=gsub('Dorsal-Attention','DAN',gsub('Cingulo-Opercular','CON',gsub('Orbito-Affective','OAN',df_long[,1])))
df_long[,1]=gsub('Default','DMN',gsub('Frontoparietal','FPN',df_long[,1]))
df=df_long
sc_pairs=df_long[,1]
nets <- c('AUD','CON','DAN','DMN','FPN','LAN','OAN','PMN','SMN','Sub','Vis1','Vis2','VMN')
get_links_for_FC <- function(fc_name, df) {
  df %>%
    filter(FC_Network == fc_name) %>%
    separate(Network_Pair, into = c("SC_Net1", "SC_Net2"), sep = "-") %>%
    mutate(SC_Pair = paste(SC_Net1, SC_Net2, sep = "_")) %>%
    group_by(SC_Pair, SC_Net1, SC_Net2) %>%
    summarise(weight = mean(AvgWeight, na.rm = TRUE), .groups = "drop") %>%
    mutate(FC = fc_name)
}
fc_list <- nets[-10]
n_fc <- length(fc_list)
all_nets <- nets
distinct_cols <- c(brewer.pal(12, "Set3"), "#000000")[1:13]  # add black if needed
names(distinct_cols) <- nets  # same order as your "nets" vector
distinct_cols[1:4]=c('red','pink','blue','purple')
distinct_cols[7:8]=c('green','yellow');distinct_cols[11:12]=c('lightgreen','brown')
angles <- seq(0, 2*pi, length.out = n_fc + 1)[- (n_fc + 1)]
radius <- 0.8
positions <- data.frame(
  x = radius * cos(angles),
  y = radius * sin(angles),
  angle = angles,
  FC = fc_list
)
pdf('SCFC_ChordDiagram_top20perc.pdf')
plot.new()
par(usr = c(-2, 2, -2, 2), xpd = TRUE)
for (i in seq_len(n_fc)) {
  x0 <- positions$x[i]
  y0 <- positions$y[i]
  r <- 0.33   # make each chord diagram larger
  fig_x <- (x0 - r + 1.6) / 3.2
  fig_y <- (y0 - r + 1.6) / 3.2
  fig_xend <- (x0 + r + 1.6) / 3.2
  fig_yend <- (y0 + r + 1.6) / 3.2
  fig_region <- pmin(pmax(c(fig_x, fig_xend, fig_y, fig_yend), 0), 1)
  par(fig = fig_region, new = TRUE)
  df_fc <- get_links_for_FC(positions$FC[i], df)
  if (nrow(df_fc) == 0) next
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = 10,
    track.margin = c(0, 0),
    cell.padding = c(0, 0, 0, 0),
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1, 1)
  )
  df_fc1=df_fc[df_fc$weight>quantile(df_fc$weight,0.8),]
	chordDiagram(
	  df_fc1[, c("SC_Net1", "SC_Net2", "weight")],
	  grid.col = distinct_cols,
	  transparency = 0.4,
	  directional = 0,
	  annotationTrack = "grid",
	  preAllocateTracks = list(track.height = 0.05)
	)
  title(main = positions$FC[i], cex.main = 0.9, line = -1.5)
}
dev.off()
pdf('SCFC_ChordDiagram_top10perc.pdf')
plot.new()
par(usr = c(-2, 2, -2, 2), xpd = TRUE)
for (i in seq_len(n_fc)) {
  x0 <- positions$x[i]
  y0 <- positions$y[i]
  r <- 0.33   # make each chord diagram larger
  fig_x <- (x0 - r + 1.6) / 3.2
  fig_y <- (y0 - r + 1.6) / 3.2
  fig_xend <- (x0 + r + 1.6) / 3.2
  fig_yend <- (y0 + r + 1.6) / 3.2
  fig_region <- pmin(pmax(c(fig_x, fig_xend, fig_y, fig_yend), 0), 1)
  par(fig = fig_region, new = TRUE)
  df_fc <- get_links_for_FC(positions$FC[i], df)
  if (nrow(df_fc) == 0) next
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = 10,
    track.margin = c(0, 0),
    cell.padding = c(0, 0, 0, 0),
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1, 1)
  )
  df_fc1=df_fc[df_fc$weight>quantile(df_fc$weight,0.9),]
	chordDiagram(
	  df_fc1[, c("SC_Net1", "SC_Net2", "weight")],
	  grid.col = distinct_cols,
	  transparency = 0.4,
	  directional = 0,
	  annotationTrack = "grid",
	  preAllocateTracks = list(track.height = 0.05)
	)
  title(main = positions$FC[i], cex.main = 0.9, line = -1.5)
}
dev.off()


# install.packages("plotrix") # if not installed
library(plotrix)

# -----------------------------
# Define network names and colors
# -----------------------------
nets <- c("AUD", "CON", "DAN", "DMN", "FPN", "LAN",
          "OAN", "PMN", "SMN", "Sub", "Vis1", "Vis2", "VMN")

cols <- c("red", "pink", "blue", "purple", "#80B1D3", "#FDB462",
          "green", "yellow", "#D9D9D9", "#BC80BD", "lightgreen", "brown", "#000000")

n <- length(nets)
angles <- seq(0, 2*pi, length.out = n + 1)[- (n + 1)]
r <- 1  # radius for circular layout
x <- r * cos(angles)
y <- r * sin(angles)
pdf("innercirc.pdf", width = 4.5, height = 4.5)
plot(NA, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5),
     asp = 1, axes = FALSE, xlab = "", ylab = "",
     main = "Network Legend (Labels Inside Circles)")
points(x, y, pch = 21, bg = cols, col = "black", cex = 5)
text(x, y, nets, col = c(rep("white",4),rep("black",7),rep("white",2)), font = 2, cex = 1)
draw.circle(0, 0, r, border = "gray80", lty = 2)
dev.off()

####################################################################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/fMRI_coupling/PopulationAnalysis
module load r/4.4.0;R
library(tidyverse)
library(pals)         # for jet colormap
library(scales)       # for label_number
library(ggplot2)
METHOD='CBSS_Line_mask_' #''
fileT=Sys.glob(paste0('G360_',METHOD,'out/*_metrics.txt'))
temp1=gsub(paste0('G360_',METHOD,'out/'),'',gsub('_metrics.txt','',fileT))
file1=unlist(lapply(strsplit(temp1,'_'),'[',1))
file2=unlist(lapply(strsplit(temp1,'_'),'[',2))
net0=unique(file1);net0=net0[c(1,2,4,3,5:9,11:12,10)]
net00=c("AUD","CON","DAN","DMN","FPN","LAN","OAN","PMN","SMN","Vis1","Vis2","VMN")
M=matrix(NA,length(net00),length(net00));rownames(M)=colnames(M)=net0
for(ii in 1:length(fileT)){
filetemp=read.table(fileT[ii])
M[file1[ii],file2[ii]]=M[file2[ii],file1[ii]]=round(filetemp$avgcorr_demean,2)
}
rownames(M)=colnames(M)=net00
M1=M
METHOD='CBSS_Line_nonmask_' #''
fileT=Sys.glob(paste0('G360_',METHOD,'out/*_metrics.txt'))
temp1=gsub(paste0('G360_',METHOD,'out/'),'',gsub('_metrics.txt','',fileT))
file1=unlist(lapply(strsplit(temp1,'_'),'[',1))
file2=unlist(lapply(strsplit(temp1,'_'),'[',2))
net0=unique(file1);net0=net0[c(1,2,4,3,5:9,11:12,10)]
net00=c("AUD","CON","DAN","DMN","FPN","LAN","OAN","PMN","SMN","Vis1","Vis2","VMN")
M=matrix(NA,length(net00),length(net00));rownames(M)=colnames(M)=net0
for(ii in 1:length(fileT)){
filetemp=read.table(fileT[ii])
M[file1[ii],file2[ii]]=M[file2[ii],file1[ii]]=round(filetemp$avgcorr_demean,2)
}
rownames(M)=colnames(M)=net00
Contrast=(M1/M)^2
indd=sort(-apply(Contrast,2,median),ind=T)$ix
#boxplot(Contrast[indd,indd])

library(ggplot2)
indd <- sort(-apply(Contrast, 2, median), index.return = TRUE)$ix
C <- as.matrix(Contrast[indd, indd])
df <- as.data.frame(as.table(C))
colnames(df) <- c("SC_net", "FC_net", "contrast")
df <- df[!is.na(df$contrast), ]
df$FC_net <- factor(df$FC_net, levels = colnames(C))
df$SC_net <- factor(df$SC_net, levels = rownames(C))
library(RColorBrewer)
pp <- ggplot(df, aes(x = 1, y = contrast)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.4,
    fill  = "grey90",
    color = "black"
  ) +
  geom_jitter(
    aes(color = SC_net),
    width  = 0.12,
    height = 0,
    size   = 2.2,
    alpha  = 0.9
  ) +
  facet_wrap(~ FC_net, nrow = 3) +
  scale_x_continuous(breaks = NULL, labels = NULL) +
  scale_color_brewer(palette = "Paired") +
  labs(
    x     = NULL,
    y     = "Coupled / non-coupled contrast index",
    color = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    # keep axis lines but remove all grid lines
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    # pull legend closer to the panels
	legend.key.width   = unit(0.85, "lines"),
    legend.spacing.x   = unit(0.8, "lines"),
    legend.text        = element_text(margin = margin(l = 0.2, r = 0)),
    legend.position    = "right",
    legend.box.margin  = margin(0, -30, 4, -5.5),
    legend.margin      = margin(0, 40, 4, -5.5)
  )

ggsave("CBSS_Line_couple_contrast.pdf", plot = pp, width = 5, height = 3.5)


####################################################################################

########################################
#Plot
########################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/fMRI_coupling/PopulationAnalysis
module load r/4.4.0;R
library(data.table)
library(tidyverse)
library(circlize)
library(RColorBrewer)
df <- read_csv("SC2FC_MaxWeight.csv", show_col_types = FALSE)
#df <- read_csv("SC2FC_MeanWeight.csv", show_col_types = FALSE)
df=df[,c(8,9,11,18,24,32,37,45,52,56,60,70,72,83,95)]
colnames(df)=c('From','To','Network_Pair','AUD','CON','DMN','DAN','FPN','LAN','OAN','PMN','SMN','VMN','Vis1','Vis2')
df_long <- df %>%
  select(Network_Pair, where(is.numeric)) %>%    # keep only numeric FC columns
  pivot_longer(-Network_Pair, 
               names_to = "FC_Network", 
               values_to = "Weight") %>%
  group_by(Network_Pair, FC_Network) %>%
  summarise(AvgWeight = mean(Weight, na.rm = TRUE), .groups = "drop")
df_long=as.data.frame(df_long)
df_long[,1]=gsub('Visual','Vis',gsub('Somatomotor','SMN',gsub('Ventral-Multimodal','VMN',df_long[,1])))
df_long[,1]=gsub('Auditory','AUD',gsub('Language','LAN',gsub('Posterior-Multimodal','PMN',df_long[,1])))
df_long[,1]=gsub('Dorsal-Attention','DAN',gsub('Cingulo-Opercular','CON',gsub('Orbito-Affective','OAN',df_long[,1])))
df_long[,1]=gsub('Default','DMN',gsub('Frontoparietal','FPN',df_long[,1]))
df=df_long
sc_pairs=df_long[,1]
nets <- c('AUD','CON','DAN','DMN','FPN','LAN','OAN','PMN','SMN','Sub','Vis1','Vis2','VMN')
get_links_for_FC <- function(fc_name, df) {
  df %>%
    filter(FC_Network == fc_name) %>%
    separate(Network_Pair, into = c("SC_Net1", "SC_Net2"), sep = "-") %>%
    mutate(SC_Pair = paste(SC_Net1, SC_Net2, sep = "_")) %>%
    group_by(SC_Pair, SC_Net1, SC_Net2) %>%
    summarise(weight = mean(AvgWeight, na.rm = TRUE), .groups = "drop") %>%
    mutate(FC = fc_name)
}
#fc_list <- nets[-10]
# "AUD"  "CON"  "DAN"  "DMN"  "FPN"  "LAN"  "OAN"  "PMN"  "SMN"  "Sub" "Vis1" "Vis2" "VMN"
#fc_list <- c("AUD", "DAN", "LAN", "SMN", "Vis1", "Vis2")
fc_list <- c("CON", "DMN", "FPN", "OAN", "PMN", "VMN")
n_fc <- length(fc_list)

angles <- seq(0, 2*pi, length.out = n_fc + 1)[-(n_fc + 1)]

radius <- 0.68   # smaller radius = centers closer together
positions <- data.frame(
  x = radius * cos(angles),
  y = radius * sin(angles),
  angle = angles,
  FC = fc_list
)

pdf('SCFC_ChordDiagram_top10perc_2.pdf', width = 9, height = 9)
#pdf('SCFC_ChordDiagram_top20perc_2.pdf', width = 9, height = 9)
plot.new()
par(usr = c(-2, 2, -2, 2), xpd = TRUE)

for (i in seq_len(n_fc)) {
  x0 <- positions$x[i]
  y0 <- positions$y[i]

  r <- 0.46   # larger subplot size

  fig_x <- (x0 - r + 1.6) / 3.2
  fig_y <- (y0 - r + 1.6) / 3.2
  fig_xend <- (x0 + r + 1.6) / 3.2
  fig_yend <- (y0 + r + 1.6) / 3.2
  fig_region <- pmin(pmax(c(fig_x, fig_xend, fig_y, fig_yend), 0), 1)

  par(fig = fig_region, new = TRUE)

  df_fc <- get_links_for_FC(positions$FC[i], df)
  if (nrow(df_fc) == 0) next

  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = 6,
    track.margin = c(0, 0),
    cell.padding = c(0, 0, 0, 0),
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1, 1)
  )

  #cutoff <- quantile(df_fc$weight, 0.8, na.rm = TRUE)
  cutoff <- quantile(df_fc$weight, 0.9, na.rm = TRUE)
  df_fc1 <- df_fc[df_fc$weight >= cutoff, ]

  chordDiagram(
    df_fc1[, c("SC_Net1", "SC_Net2", "weight")],
    grid.col = distinct_cols,
    transparency = 0.4,
    directional = 0,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.07)
  )

  title(main = positions$FC[i], cex.main = 3.0, line = -4.0)
}

dev.off()