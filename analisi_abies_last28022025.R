################################################################################
# Set up working directory

# setwd("")

setwd("C:\\aaa_lavori\\lav_waed_abies")

##########################################################################################################
# load R packages  and set the random seeds

source("load_libraries_abies.R")

set.seed(123) # set seeds 

#######################################################################################################
# clean working area

cat("\014") 

################################################################################
# Loading data of GCMS data ( Areas) and seeds

dati_sel=readRDS("dati_sel_abies.rds") # data relative to seeds
dati_seeds=read.xlsx("dati_abies_seeds.xlsx",1) # data relative to seeds


cat("\014") 

################################################################################
# list variabile's name

names(dati_sel)  # vedere i nomi delle matrici 

# [1] "ID"                    "Species"               "date"                  "a-pinene"              "camphene"             
# [6] "b-pinene"              "myrcene"               "Limonene"              "β-Phellandrene"        "g-terpinene"          
# [11] "p-cymene"              "terpinolene"           "Limonene_oxide_cis"    "Pinocarvone"           "Unk_1"                
# [16] "bornylacetate"         "4-ol-terpinen"         "trans-Verbenol"        "verbenone"             "b-elemene"            
# [21] "sesquiterpeni.1"       "sesquiterpene.2"       "a-humulene"            "GermacreneD"           "sesquiterpene.3"      
# [26] "B.caryophillene-oxide" "1.Germacrene.D-4-ol"   "unk_.2"                "α-epi-Cadinol"         "Selina-6-en-4-ol"     
# [31] "α-Cubebene"            "a-copaene"             "b-copaene"             "TOT_Mono"              "TOT_mono&sesqui"    



# data numerosity

nrow(dati_sel) #  170

#######################################################################################################
# Data normalization

dati_sel_rel=100*dati_sel[,4:31]/dati_sel$`TOT_mono&sesqui`

mat_final_rel=data.frame(species=dati_sel$Species,dati_sel_rel)

write.xlsx(mat_final_rel,"dati_relativi_totali.xlsx",overwrite = T) 

################################################################################
#  monoTh AND  sesquiTh separated

dati_monosesqui=cbind(100*dati_sel[,4:18]/dati_sel$TOT_Mono,100*dati_sel[,19:31]/(dati_sel$`TOT_mono&sesqui`-dati_sel$TOT_Mono))

# only  monoTh on monothtot and sesquiTh on the sum to underweight sesqui


dati_monosesquimiche=cbind(100*dati_sel[,4:18]/dati_sel$TOT_Mono,100*dati_sel[,19:31]/(dati_sel$`TOT_mono&sesqui`))

##########################################################################
# Creazione delle matrici di analisi per PCA
##########################################################################

dati_sel_rel=dati_monosesquimiche # PCA on full data normalized matrix

X=dati_sel_rel
Y=dati_sel$Species

# X=asin(sqrt(X/100)) # arsin trasformation 



#########################################################################################
# PCA explore variable and outlier detection

data_pca=data.frame(X,Species=Y)


abies_pca=ordinate(data_pca , cols = 1:25, model = ~ prcomp(., scale. = TRUE)) 

abies_pca %>%
  augment_ord() %>%
  ggbiplot(axis.type = "predictive",aes(alpha = .3)) +
  theme_bw() +
  geom_rows_point(aes(color=Species,shape=Species))+
  stat_rows_center(alpha = .8, fun.center = "mean",size = 5,aes(color = Species))+
  geom_origin() +
  ggtitle("Principal Component Analisys Abies terpene data") +
  labs(color = "Species")

ggsave("PCA_biplot.png")



#########################################################################################
# data selection by non parametric ANOVA o Kruskal Wallis ( optional) 

a=col_kruskalwallis(X,Y) # by using uso il pacchetto MatrixTests

row.names(a)[which(a$pvalue<0.05)] # all coumpound discriminate!
row.names(a)[which(a$pvalue>0.05)] 

write.xlsx(data.frame(a,compound_terpene=row.names(a)),"test_kruskal.wallis.xlsx")

##########################################################################################################
# Faccio la Posthoc analysis

res_posthoc=list() # list postHoc with KS test
res_posthoc_dunn=list() # list postHoc with Dunn test PMCMRplus
res_GG=list()

#################

j=1  

for ( i in 2:29) { # run on the columns starting fron numeric values
  
formula_t=paste(names(mat_final_rel)[i], "~ species + 0")

MM=glm(formula_t,  data = mat_final_rel)

GG <- posthoc(MM)


dunn_df=as.data.frame(gsub("\\s+", " ", paste0(names(mat_final_rel[i])," ",capture.output(summary(kwAllPairsDunnTest(mat_final_rel[,i], factor(mat_final_rel$species),p.adjust.method = "bonferroni")))[2:4])))

names(dunn_df)=c("Compound Pair Z p_values sign")
res_posthoc_dunn[[j]]=dunn_df
res_posthoc[[j]]=data.frame(terpene_compound=names(mat_final_rel[i]),
                            GG$CI,groups=GG$Grouping,
                            test_sp=gsub("species","",row.names(GG$PvaluesMatrix)),
                            GG$PvaluesMatrix,
                            ks_pval=KruskalWallisAllPvalues(mat_final_rel[,i], factor(mat_final_rel$species))
                            )
res_GG[[j]]=GG

j=j+1

}

res_df=do.call("rbind",res_posthoc)
res_df_dunn=do.call("rbind",res_posthoc_dunn)

#################################################################
# write results in two separated sheets defined by list

write.xlsx(list(dunn_test=res_df_dunn,
                t_test=res_df),
                "posthoc_kruskal.wallis.xlsx")

##############################################################################################
# summary tables

my_desc=psych::describeBy(mat_final_rel[,-1], mat_final_rel$species)

write_excel_csv2(rbindlist(my_desc, idcol = 'species'), file = 'summary_terpene_by_species.xlsx')


######################################################################################
# boxplot main discriminant compounds


names(dati_sel_rel)[c(1,3,5,8,17,23,25)]

dati_sel_box=dati_sel_rel
dati_sel_box$Species=factor(dati_sel$Species)
dati_sel_box=janitor::clean_names(dati_sel_box) # clean names of compounds

dir.create("boxplot_sel")

setwd("boxplot_sel")

# [1] "a-pinene"            
ggboxplot(dati_sel_box,"species",names(dati_sel_box)[1],fill="red") +ylim(0,50)+ylab(paste0(names(dati_sel_box)[1]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box)[1],".png"))
# "b-pinene"            
ggboxplot(dati_sel_box,"species",names(dati_sel_box)[3],fill="red") +ylim(0,50)+ylab(paste0(names(dati_sel_box)[3]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box)[3],".png"))

# "Limonene"
ggboxplot(dati_sel_box,"species",names(dati_sel_box)[5],fill="red") +ylim(20,90)+ylab(paste0(names(dati_sel_box)[5]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box)[5],".png"))

# [4] "p-cymene"
ggboxplot(dati_sel_box,"species",names(dati_sel_box)[8],fill="red") +ylim(0,3)+ylab(paste0(names(dati_sel_box)[8]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box)[8],".png"))

# "sesquiterpeni.1"  
ggboxplot(dati_sel_box,"species",names(dati_sel_box)[17],fill="red") +ylim(0,5)+ylab(paste0(names(dati_sel_box)[17]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box)[17],".png"))

# "1.Germacrene.D-4-ol"
ggboxplot(dati_sel_box,"species",names(dati_sel_box)[23],fill="red") +ylim(0,50)+ylab(paste0(names(dati_sel_box)[23]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box)[23],".png"))

# [7] "Selina-6-en-4-o

ggboxplot(dati_sel_box,"species",names(dati_sel_box)[25],fill="red") +ylim(0,50)+ylab(paste0(names(dati_sel_box)[25]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box)[25],".png"))

setwd("..")


dir.create("boxplot_all")

setwd("boxplot_all")


for ( i in 1:25) {

ggboxplot(dati_sel_box,"species",names(dati_sel_box)[i],fill="red") +ylim(0,50)
  
ggsave(paste0("boxplot_",names(dati_sel_box)[i],".png"))

}

setwd("..")



####################################################################################
# Ordination Supervised methods : LDA & & KNN

# LDA
# Prior probabilities of groups: the proportion of training observations in each group. 
# For example, there are 31% of the training observations in the setosa group
# Group means: group center of gravity. Shows the mean of each variable in each group.
# Coefficients of linear discriminants: Shows the linear combination of predictor variables that are used to form the LDA

##############################################################################################################################
#############################################################################################################################
# checking collinearity by Condition Number for obtain data selected for next LDA procedure
# why?
# Multicollinearity means that your predictors are correlated. Why is this bad?
#   
# Because LDA, like regression techniques involves computing a matrix inversion, which is inaccurate 
# if the determinant is close to 0 (i.e. two or more variables are almost a linear combination of each other).
# More importantly, it makes the estimated coefficients impossible to interpret. If an increase in X1
# ,say, is associated with an decrease in X2 and they both increase variable Y, every change in X1 will 
# be compensated by a change in X2 and you will underestimate the effect of X1 on Y.
# In LDA, you would underestimate the effect of X1 on the classification.
########################################################################################################




dati_sel_rel=dati_monosesquimiche

CN(dati_sel_rel) # 
multiCol(dati_sel_rel, graf = TRUE)
names(dati_sel_rel)[c(20,16,23)]

dati_sel_rel_sel=dati_sel_rel[,-c(20,16,23)]
CN(dati_sel_rel_sel) 
multiCol(dati_sel_rel_sel, graf = TRUE)

###################################################################################
X=dati_sel_rel_sel
Y=dati_sel$Species

training.samples <- Y %>% createDataPartition(p = 0.8, list = FALSE)

dati_train=X[training.samples, ]
dati_test=X[-training.samples, ]
dati_test$Y=Y[-training.samples]
dati_train$Y=Y[training.samples]

model <- lda(Y~., data =dati_train)

predictions_lda <- model %>% predict(dati_test)
model_accuracy_lda=mean(predictions_lda$class==dati_test$Y)
res_lda_abiesGC=confusionMatrix(predictions_lda$class,factor(dati_test$Y))

print(xtable(res_lda_abiesGC$table), type = "html",file="confusion_matrix.docx")
print(xtable(res_lda_abiesGC$byClass), type = "html",file="index_classification.docx")
res_lda_abiesGC$overall
res_lda_abiesGC


###########################################################################################################################################
# plots LDA

abies_lda <- lda_ord(X, Y, axes.scale = "standardized")
abies_lda %>%
  as_tbl_ord() %>%
  augment_ord() %>%
  mutate_rows(
    species = grouping,
    data = ifelse(.element == "active", "centroid", "case")
  ) %>%
  ggbiplot() +
  theme_bw() +
  geom_rows_point(aes(
    color = grouping,
    size = data, 
    alpha = data
  )) +
  ggtitle("Standardized coefficient biplot of mediterranean Abies spp. ") +
  expand_limits(y = c(-3, 5))+
  labs(color = "Species")

ggsave("LDA_biplot.png")


################################################################################
# Now working on data seeds


Xseed=dati_seeds[,4:8]
Yseed=dati_seeds$species

write.xlsx(list(data.frame(Y,X),data.frame(Yseed,Xseed)),"dati_PCA_seeds.xlsx")

cat("\014") 

########################################################################################################
# Summary descriptive data 

my_desc_seed=psych::describeBy(Xseed,Yseed)
table_seeds=rbindlist(my_desc_seed, idcol = 'species')
table_seeds$vars=rownames(my_desc_seed$`Abies alba`)[table_seeds$vars]
write_excel_csv2(table_seeds, file = 'summary_seeds_by_species.xlsx')


########################################################################################################
# One way ANOVA

col_oneway_welch(Xseed, Yseed) # utilizzo il pacchetto MatrixTests per fare l'ANOVA a 1 via 

########################################################################################################

X=Xseed
Y=Yseed

################################################################################
# PCA

ord <- PCA(Xseed, , scale.unit = T,graph = T)
ggord(ord, Yseed,arrow=NULL,txt=NULL)
pc_seeds=ord$ind$coord # new dimensions

#################################################################################################
# 

data_pca=data.frame(X,Species=Y)

seeds_abies_pca=ordinate(data_pca , cols = 1:5, model = ~ prcomp(., scale. = TRUE)) 

seeds_abies_pca %>%
  augment_ord() %>%
  ggbiplot(axis.type = "predictive",aes(alpha = .3)) +
  theme_bw() +
  geom_rows_point(aes(color=Species,shape=Species))+
  stat_rows_center(alpha = .8, fun.center = "mean",size = 5,aes(color = Species))+
  geom_origin() +
  ggtitle("PCA of Seeds Abies data matrix") +
  labs(color = "Species")

ggsave("PCA_biplot_seeds.png")

#####################################################################
# References
# https://yangxiaozhou.github.io/data/2019/10/02/linear-discriminant-analysis.html
# https://cmdlinetips.com/2020/12/canonical-correlation-analysis-in-r/
# https://www.r-bloggers.com/2021/05/linear-discriminant-analysis-in-r/
# https://vitalflux.com/pca-vs-lda-differences-plots-examples/
# https://towardsai.net/p/data-science/lda-vs-pca
# https://stats.stackexchange.com/questions/23353/pca-lda-cca-and-pls
# https://www.geeksforgeeks.org/classifying-data-using-support-vector-machinessvms-in-r/
# https://mdatools.com/docs/pca.html
# https://rdrr.io/cran/mixOmics/man/plsda.html
# http://mixomics.org/methods/spls-da/