#cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/Aging;module load r/4.4.0;R
#Voxel level BAD
library(data.table)
A=fread('variables/demo_20230308.csv',data.table=F)
B_traintemp=fread(paste0('../../PCA_ELReg/FA_UKB_Phase12/test_PCs_1.txt'),data.table=F)
B_testtemp=fread(paste0('../../PCA_ELReg/FA_UKB_Phase3_proj/test_PCs_1.txt'),data.table=F)
UID_train1=intersect(A[which(!is.na(A[,'Age_img'])),1],B_traintemp[,1])
UID_test1=intersect(A[which(!is.na(A[,'Age_img'])),1],B_testtemp[,1])
rownames(A)=A[,1]
A_train1=A[as.character(UID_train1),]
A_test1=A[as.character(UID_test1),]
mu0=mean(A_train1[,'Age_img']);sd0=sd(A_train1[,'Age_img'])
pred_age=fread('out/Net_cor_pred_60900_withMean/Pred_demo_20230308_Age_img.csv',data.table=F)
pred_age=t(t(pred_age)*sd0+mu0)
BAD=pred_age-A_test1[,'Age_img']
BAD2=BAD;BAD2[T]=NA;colnames(BAD2)=c('Vis1','Vis2','SMN','CON','DAN','LAN','FPN','AUD','PMN','DMN','OAN','VMN','Sub')
colnames(BAD2)=paste0('BAD_',colnames(BAD2))
set.seed(2025);validind=sample(1:dim(A_test1)[1],floor(dim(A_test1)[1])/2)
testind=setdiff(1:dim(A_test1)[1],validind)
for(jj in 1:dim(BAD)[2]){
A_test1$BAD=BAD[,jj]
BAD2[validind,jj] <- BAD[validind,jj]-predict(lm(BAD ~ Age_img, data = A_test1[testind,]), newdata = data.frame(Age_img = A_test1[validind, 'Age_img']))
BAD2[testind,jj] <- BAD[testind,jj]-predict(lm(BAD ~ Age_img, data = A_test1[validind,]), newdata = data.frame(Age_img = A_test1[testind, 'Age_img']))
}
A_test1=cbind(A_test1,BAD2)
boxplot(BAD~sex,data=A_test1)
C0=fread('variables/demo_20230308.csv',data.table=F)
C1=fread('variables/ukb45941_cognitive_coded.csv',data.table=F)
C2=fread('variables/mental_health_100060.csv',data.table=F)
C3=fread('variables/lifestyle_and_environment.csv',data.table=F)
rownames(C0)=C0[,1];C0=C0[as.character(UID_test1),-1]
rownames(C1)=C1[,1];C1=C1[as.character(UID_test1),-1]
rownames(C2)=C2[,1];C2=C2[as.character(UID_test1),-1]
rownames(C3)=C3[,1];C3=C3[as.character(UID_test1),-1]
C=cbind(C0[,c(1,3,6,8,9,10)],C1,C2,C3)
D=fread('variables/confound/UKB_Age_Site_Motion.txt',data.table=F)
rownames(D)=D[,1];D=D[,-1];D=D[as.character(UID_test1),]
data0=cbind(D,C);data0=as.data.frame(cbind(data0,scale(BAD2)))
#######
data0$QC_Imaging_Site=factor(data0$QC_Imaging_Site)
CorT=pT=sdT=matrix(NA,116-14+1,13);for(jj in 1:13){
print(13-jj)
for(ii in 14:116){
data0$y=scale(data0[,ii])
data0$BAD=scale(data0[,jj+116])
if(!colnames(data0)[ii]%in%c('sex','race','edu'))
mod0=lm(y~BAD+QC_Imaging_Site+Sex+brain_position1+brain_position2+brain_position3+brain_position4+volumetric_scaling+
Sex_Age+Age2+Sex_Age2+race+edu,data=data0)
if(colnames(data0)[ii]=='sex')mod0=lm(y~BAD+QC_Imaging_Site+brain_position1+brain_position2+brain_position3+brain_position4+volumetric_scaling+Age2+race+edu,data=data0)
if(colnames(data0)[ii]=='race')mod0=lm(y~BAD+QC_Imaging_Site+Sex+brain_position1+brain_position2+brain_position3+brain_position4+volumetric_scaling+
Sex_Age+Age2+Sex_Age2+edu,data=data0)
if(colnames(data0)[ii]=='edu')mod0=lm(y~BAD+QC_Imaging_Site+Sex+brain_position1+brain_position2+brain_position3+brain_position4+volumetric_scaling+
Sex_Age+Age2+Sex_Age2+race,data=data0)
CorT[ii-13,jj]=as.numeric(mod0$coef['BAD']);pT[ii-13,jj]=summary(mod0)$coef['BAD','Pr(>|t|)'];sdT[ii-13,jj]=summary(mod0)$coef['BAD','Std. Error']; 
}
}
rownames(CorT)=rownames(pT)=rownames(sdT)=colnames(data0)[14:116]
colnames(CorT)=colnames(pT)=colnames(sdT)=colnames(BAD2)
write.csv(CorT,file='out/Net_cor_pred_60900_withMean/BAD_rvalue.csv',quote=F)
write.csv(pT,file='out/Net_cor_pred_60900_withMean/BAD_pvalue.csv',quote=F)
write.csv(sdT,file='out/Net_cor_pred_60900_withMean/BAD_sd.csv',quote=F)
INDD=c('Smoking_status_20116','Average_weekly_red_wine_intake_1568','sex')
CorT=CorT[INDD,];pT=pT[INDD,];sdT=sdT[INDD,]

#top brain aging factors
                                Smoking_status_20116
                                          0.07408399
                       Alcohol_intake_frequency_1558
                                         -0.05953187
                 Average_weekly_red_wine_intake_1568
                                          0.06486884
                                 Digit_matches_23324
                                         -0.05756316
                                   Ever_smoked_20160
                                          0.04840754
          Average_weekly_beer_plus_cider_intake_1588
                                          0.04802515
                                   Numeric_path_6348
                                          0.04253103
                  Average_weekly_spirits_intake_1598
                                          0.04596022
                                  Cereal_intake_1458
                                         -0.04132755
                        Current_tobacco_smoking_1239
                                          0.04034741
                              Alphanumeric_path_6350
                                          0.03710494
                                     Tea_intake_1488
                                         -0.03556333
                                  Coffee_intake_1498
                                          0.03447731
Average_weekly_champagne_plus_white_wine_intake_1578
                                          0.03648454
                                 Puzzles_solved_6373
                                         -0.02906227


#write.csv(data0,file='out/cor_pred_60900_withMean/BAD_Net.csv',quote=F)
data0=data0[,c(indd+13,which(colnames(data0)=='BAD'))]
colnames(data0)=c('SmokingStatus','Digit_matches','Alcohol_freq_less','RedWine','Ever_smoke','Numeric_path','Weekly_spirits','Cereal',
'CurrentTobacco','Champagne_wine','Time_alpha_numeric','Tea','Beer_cider','Coffee','Sex','Puzzles_solved','BAD')
names(CorT)=colnames(data0)[-dim(data0)[2]]
names(sdT)=colnames(data0)[-dim(data0)[2]]
#write.csv(data0,file='out/cor_pred_60900_withMean/BAD.csv',quote=F,row.names=F)

######################################################################
#plot
######################################################################
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

# --- tidy inputs (same as before) ---
df_eff <- as.data.frame(CorT) %>%
  rownames_to_column("Variable") %>%
  pivot_longer(-Variable, names_to = "Network", values_to = "Effect")

df_se <- as.data.frame(sdT) %>%
  rownames_to_column("Variable") %>%
  pivot_longer(-Variable, names_to = "Network", values_to = "SE")

df <- df_eff %>%
  left_join(df_se, by = c("Variable", "Network")) %>%
  mutate(
    Network = str_replace(Network, "^BAD_", ""),
    CI_lo   = Effect - 1.96 * SE,
    CI_hi   = Effect + 1.96 * SE
  )

# order networks by median |effect| (keep or change as you like)
net_order <- df %>%
  group_by(Network) %>% summarize(med_abs = median(abs(Effect)), .groups="drop") %>%
  arrange(desc(med_abs)) %>% pull(Network)

df <- df %>%
  mutate(
    Network  = factor(Network, levels = net_order),
    Variable = str_replace_all(Variable, "_20116|_1568", "")
  )

# --- rotated forest plot (vertical CIs), darker colors, large text, tight gap ---
dodge_w <- 0.35  # smaller dodge -> less horizontal gap between the 3 predictors
df[df[,1]=='Average_weekly_red_wine_intake',1]='RedWine'
df[df[,1]=='Smoking_status',1]='Smoking'
p_forest_v <- ggplot(df, aes(x = Network, y = Effect, colour = Variable, group = Variable)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi),
                width = 0, size = 1.1,
                position = position_dodge(width = dodge_w)) +
  geom_point(size = 3.8, position = position_dodge(width = dodge_w)) +
  # darker palette
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(
    x = NULL,
    y = "",
    colour = NULL,
    title = ""
  ) + scale_y_continuous(breaks = c(0.025, 0.075,0.125), labels = c("0.025", "0.075","0.125")) +
  theme_minimal(base_size = 22) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x = element_text(size = 40, angle = 35, hjust = 1),  # bigger x labels
    axis.text.y = element_text(size = 40),
     legend.position      = c(0.4, 1.10),
    legend.justification = c(0.5, 1),
    legend.direction     = "horizontal",                   # laid out left→right
    legend.box.background = element_rect(fill = scales::alpha("white", 0.5), colour = NA),
    legend.text          = element_text(size = 40),
    legend.spacing.x     = unit(6, "pt"),
    legend.key.width     = unit(18, "pt"),         
    plot.title = element_text(face = "bold", hjust = 0.5, size = 36),
    plot.margin = margin(8, 14, 6, 14)
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))+
  coord_cartesian(ylim = c(0.025, 0.125), clip = "off")

ggsave("aging_effects_networks_forest_rotated.pdf",
       p_forest_v, width = 11, height = 5, units = "in")

##############################################################################################################################################################################
#Fiber level BAD
#cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/Aging;module load r/4.4.0;R
library(data.table)
A=fread('variables/demo_20230308.csv',data.table=F)
B_traintemp=fread(paste0('../../PCA_ELReg/FA_UKB_Phase12/test_PCs_1.txt'),data.table=F)
B_testtemp=fread(paste0('../../PCA_ELReg/FA_UKB_Phase3_proj/test_PCs_1.txt'),data.table=F)
UID_train1=intersect(A[which(!is.na(A[,'Age_img'])),1],B_traintemp[,1])
UID_test1=intersect(A[which(!is.na(A[,'Age_img'])),1],B_testtemp[,1])
rownames(A)=A[,1]
A_train1=A[as.character(UID_train1),]
A_test1=A[as.character(UID_test1),]
mu0=mean(A_train1[,'Age_img']);sd0=sd(A_train1[,'Age_img'])
pred_age=fread('out/Net_cor_pred/Pred_demo_20230308_Age_img.csv',data.table=F)
pred_age=t(t(pred_age)*sd0+mu0)
BAD=pred_age-A_test1[,'Age_img']
BAD2=BAD;BAD2[T]=NA;colnames(BAD2)=c('Vis1','Vis2','SMN','CON','DAN','LAN','FPN','AUD','PMN','DMN','OAN','VMN','Sub')
colnames(BAD2)=paste0('BAD_',colnames(BAD2))
set.seed(2025);validind=sample(1:dim(A_test1)[1],floor(dim(A_test1)[1])/2)
testind=setdiff(1:dim(A_test1)[1],validind)
for(jj in 1:dim(BAD)[2]){
A_test1$BAD=BAD[,jj]
BAD2[validind,jj] <- BAD[validind,jj]-predict(lm(BAD ~ Age_img, data = A_test1[testind,]), newdata = data.frame(Age_img = A_test1[validind, 'Age_img']))
BAD2[testind,jj] <- BAD[testind,jj]-predict(lm(BAD ~ Age_img, data = A_test1[validind,]), newdata = data.frame(Age_img = A_test1[testind, 'Age_img']))
}
A_test1=cbind(A_test1,BAD2)
boxplot(BAD~sex,data=A_test1)
C0=fread('variables/demo_20230308.csv',data.table=F)
C1=fread('variables/ukb45941_cognitive_coded.csv',data.table=F)
C2=fread('variables/mental_health_100060.csv',data.table=F)
C3=fread('variables/lifestyle_and_environment.csv',data.table=F)
rownames(C0)=C0[,1];C0=C0[as.character(UID_test1),-1]
rownames(C1)=C1[,1];C1=C1[as.character(UID_test1),-1]
rownames(C2)=C2[,1];C2=C2[as.character(UID_test1),-1]
rownames(C3)=C3[,1];C3=C3[as.character(UID_test1),-1]
C=cbind(C0[,c(1,3,6,8,9,10)],C1,C2,C3)
D=fread('variables/confound/UKB_Age_Site_Motion.txt',data.table=F)
rownames(D)=D[,1];D=D[,-1];D=D[as.character(UID_test1),]
data0=cbind(D,C);data0=as.data.frame(cbind(data0,scale(BAD2)))
#######
data0$QC_Imaging_Site=factor(data0$QC_Imaging_Site)
CorT=pT=sdT=matrix(NA,116-14+1,13);for(jj in 1:13){
print(13-jj)
for(ii in 14:116){
data0$y=scale(data0[,ii])
data0$BAD=scale(data0[,jj+116])
if(!colnames(data0)[ii]%in%c('sex','race','edu'))
mod0=lm(BAD~y+QC_Imaging_Site+Sex+brain_position1+brain_position2+brain_position3+brain_position4+volumetric_scaling+
Sex_Age+Age2+Sex_Age2+race+edu,data=data0)
if(colnames(data0)[ii]=='sex')mod0=lm(BAD~y+QC_Imaging_Site+brain_position1+brain_position2+brain_position3+brain_position4+volumetric_scaling+Age2+race+edu,data=data0)
if(colnames(data0)[ii]=='race')mod0=lm(BAD~y+QC_Imaging_Site+Sex+brain_position1+brain_position2+brain_position3+brain_position4+volumetric_scaling+
Sex_Age+Age2+Sex_Age2+edu,data=data0)
if(colnames(data0)[ii]=='edu')mod0=lm(BAD~y+QC_Imaging_Site+Sex+brain_position1+brain_position2+brain_position3+brain_position4+volumetric_scaling+
Sex_Age+Age2+Sex_Age2+race,data=data0)
CorT[ii-13,jj]=as.numeric(mod0$coef['y']);pT[ii-13,jj]=summary(mod0)$coef['y','Pr(>|t|)'];sdT[ii-13,jj]=summary(mod0)$coef['y','Std. Error']; 
}
}
rownames(CorT)=rownames(pT)=rownames(sdT)=colnames(data0)[14:116]
colnames(CorT)=colnames(pT)=colnames(sdT)=colnames(BAD2)
write.csv(CorT,file='out/Net_cor_pred_60900_withMean/BAD_rvalue.csv',quote=F)
write.csv(pT,file='out/Net_cor_pred_60900_withMean/BAD_pvalue.csv',quote=F)
write.csv(sdT,file='out/Net_cor_pred_60900_withMean/BAD_sd.csv',quote=F)
INDD=c('Smoking_status_20116','Average_weekly_red_wine_intake_1568','sex')
CorT=CorT[INDD,];pT=pT[INDD,];sdT=sdT[INDD,]

#top brain aging factors
                                Smoking_status_20116
                                          0.07408399
                       Alcohol_intake_frequency_1558
                                         -0.05953187
                 Average_weekly_red_wine_intake_1568
                                          0.06486884
                                 Digit_matches_23324
                                         -0.05756316
                                   Ever_smoked_20160
                                          0.04840754
          Average_weekly_beer_plus_cider_intake_1588
                                          0.04802515
                                   Numeric_path_6348
                                          0.04253103
                  Average_weekly_spirits_intake_1598
                                          0.04596022
                                  Cereal_intake_1458
                                         -0.04132755
                        Current_tobacco_smoking_1239
                                          0.04034741
                              Alphanumeric_path_6350
                                          0.03710494
                                     Tea_intake_1488
                                         -0.03556333
                                  Coffee_intake_1498
                                          0.03447731
Average_weekly_champagne_plus_white_wine_intake_1578
                                          0.03648454
                                 Puzzles_solved_6373
                                         -0.02906227



data0=data0[,c(indd+13,which(colnames(data0)=='BAD'))]
colnames(data0)=c('SmokingStatus','Digit_matches','Alcohol_freq_less','RedWine','Ever_smoke','Numeric_path','Weekly_spirits','Cereal',
'CurrentTobacco','Champagne_wine','Time_alpha_numeric','Tea','Beer_cider','Coffee','Sex','Puzzles_solved','BAD')
names(CorT)=colnames(data0)[-dim(data0)[2]]
names(sdT)=colnames(data0)[-dim(data0)[2]]
#write.csv(data0,file='out/cor_pred_60900_withMean/BAD.csv',quote=F,row.names=F)

######################################################################
#plot
######################################################################
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

# --- tidy inputs (same as before) ---
df1_eff <- as.data.frame(CorT) %>%
  rownames_to_column("Variable") %>%
  pivot_longer(-Variable, names_to = "Network", values_to = "Effect")

df1_se <- as.data.frame(sdT) %>%
  rownames_to_column("Variable") %>%
  pivot_longer(-Variable, names_to = "Network", values_to = "SE")

df1 <- df1_eff %>%
  left_join(df1_se, by = c("Variable", "Network")) %>%
  mutate(
    Network = str_replace(Network, "^BAD_", ""),
    CI_lo   = Effect - 1.96 * SE,
    CI_hi   = Effect + 1.96 * SE
  )

# order networks by median |effect| (keep or change as you like)
net_order <- df1 %>%
  group_by(Network) %>% summarize(med_abs = median(abs(Effect)), .groups="drop") %>%
  arrange(desc(med_abs)) %>% pull(Network)

df1 <- df1 %>%
  mutate(
    Network  = factor(Network, levels = net_order),
    Variable = str_replace_all(Variable, "_20116|_1568", "")
  )

# --- rotated forest plot (vertical CIs), darker colors, large text, tight gap ---
dodge_w <- 0.35  # smaller dodge -> less horizontal gap between the 3 predictors
df1[df1[,1]=='Average_weekly_red_wine_intake',1]='RedWine'
df1[df1[,1]=='Smoking_status',1]='Smoking'
df$Method  <- "Streamline"
df1$Method <- "Voxel-based"
dd <- bind_rows(df, df1)

# order networks as they appear (or reuse your previous net_order if desired)
dd <- dd %>% mutate(Network = factor(Network, levels = unique(Network)))

# mapping + 2-level dodge: colour by Variable, shape/linetype by Method
dodge_w <- 0.45
pos_dodge <- position_dodge(width = dodge_w)

p_forest_v <- ggplot(
  dd,
  aes(x = Network, y = Effect,
      colour = Variable,
      shape  = Method,
      linetype = Method,
      group  = interaction(Variable, Method))
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi),
                width = 0, size = 1.1, position = pos_dodge) +
  geom_point(size = 3.8, position = pos_dodge) +
  # darker colours for variables
  scale_color_brewer(palette = "Dark2") +
  # shapes/linetypes to distinguish methods
  scale_shape_manual(values = c("Streamline" = 16, "Voxel-based" = 17)) +
  scale_linetype_manual(values = c("Streamline" = "solid", "Voxel-based" = "dashed")) +
  # y-axis: only two ticks and zoom to [0.025, 0.125]
  scale_y_continuous(breaks = c(0.025,0.075, 0.125), labels = c("0.025","0.075","0.125")) +
  coord_cartesian(ylim = c(0.025, 0.15), clip = "off") +
  labs(
    x = NULL, y = "",
    colour = NULL, shape = NULL, linetype = NULL,
    title = ""
  ) +
  theme_minimal(base_size = 22) 
  p_forest_v=p_forest_v+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x = element_text(size = 34, angle = 35, hjust = 1),
    axis.text.y = element_text(size = 34),
    # legend: one row, slightly above panel, with background
    legend.position      = c(0.40, 1.05),
    legend.justification = c(0.5, 1),
    legend.direction     = "horizontal",
    legend.box.background = element_rect(fill = scales::alpha("white", 0.5), colour = NA),
    legend.text          = element_text(size = 38),
    legend.spacing.x     = unit(2, "pt"),
    legend.key.width     = unit(18, "pt"),
	legend.spacing.y     = unit(-4, "pt"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 36),
    plot.margin = margin(8, 14, 6, 14)
  ) 
ggsave("aging_effects_networks_forest_both_methods.pdf",
       p_forest_v, width = 11, height = 7, units = "in")