################################################################################
# Set up working directory

setwd("")

##########################################################################################################
# load R packages  and set the random seeds

source("load_libraries_abies.R")

set.seed(123) # set seeds 

################################################################################
# Loading data of GCMS data ( Areas) and seeds

################################################################################
# To achieve  reproducibility setup work dir where data are stored.

 setwd("")


set.seed(123)

##########################################################################################################
# Load R packages

source("load_libraries_abies.R")

cat("\014") 

################################################################################
# Loading data 

dati_sel=read.xlsx("dati_abies_GC.xlsx",1) # dati relativi alle analisi GC MS

dati_seeds=read.xlsx("dati_abies_seeds.xlsx",1) # dati relativi ai paraemtri morfologici dei semi


cat("\014") 

################################################################################
# eliminazione dei dati outlier e salvataggio su un file di archiviazione R formato rds 

saveRDS(dati_sel,"dati_sel_abies.rds")

################################################################################
# Visualizzazione della matrice di lavoro GC MS e dei  nomi delle variabili


# View(dati_sel)

names(dati_sel)

# [1] "ID"                    "Species"               "date"                  "a-pinene"              "camphene"             
# [6] "b-pinene"              "myrcene"               "Limonene"              "β-Phellandrene"        "g-terpinene"          
# [11] "p-cymene"              "terpinolene"           "Limonene_oxide_cis"    "Pinocarvone"           "Unk_1"                
# [16] "bornylacetate"         "4-ol-terpinen"         "trans-Verbenol"        "verbenone"             "b-elemene"            
# [21] "sesquiterpeni.1"       "sesquiterpene.2"       "a-humulene"            "GermacreneD"           "sesquiterpene.3"      
# [26] "B.caryophillene-oxide" "1.Germacrene.D-4-ol"   "unk_.2"                "α-epi-Cadinol"         "Selina-6-en-4-ol"     
# [31] "α-Cubebene"            "a-copaene"             "b-copaene"             "TOT_Mono"              "TOT_mono&sesqui"    


nrow(dati_sel) #  170

#######################################################################################################
# Normalizzazione dei dati sul totale dei mono&sesqui perchè parto dal dato delle aree

dati_sel_rel=100*dati_sel[,4:31]/dati_sel$`TOT_mono&sesqui`

mat_final_rel=data.frame(species=dati_sel$Species,dati_sel_rel)

 
# solo dei mono e dei sesqui separati

dati_monosesqui=cbind(100*dati_sel[,4:18]/dati_sel$TOT_Mono,
                      100*dati_sel[,19:31]/(dati_sel$`TOT_mono&sesqui`-dati_sel$TOT_Mono))

# solo dei mono e dei sesqui sommati per fare pesare meno i sesqui

dati_monosesquimiche=cbind(100*dati_sel[,4:18]/dati_sel$TOT_Mono,
                           100*dati_sel[,19:31]/(dati_sel$`TOT_mono&sesqui`))

# la scelta della normalizzazione dipende dalla situazione e dagli obbiettivi degli analisi

write.xlsx(mat_final_rel,"dati_relativi_totali.xlsx",overwrite = T) 


###############################################################################################################################
# controllo di collinearità opzionale in questo caso non utilizzato

dati_sel_rel=dati_monosesquimiche

CN(dati_sel_rel) # 
multiCol(dati_sel_rel, graf = TRUE)
a=multiCol(dati_sel_rel, graf = TRUE)

id_VIF=which(names(dati_sel_rel) %in% c("b-elemene","1.Germacrene.D-4-ol","GermacreneD") ==T)

dati_sel_rel_LDA=dati_sel_rel[, -id_VIF] # selected data for LDA

CN(dati_sel_rel_LDA) 
multiCol(dati_sel_rel_LDA, graf = TRUE)

# Collinearity well explained

# https://stats.stackexchange.com/questions/70899/what-correlation-makes-a-matrix-singular-and-what-are-implications-of-singularit
# https://math.stackexchange.com/questions/2120542/48-reasons-why-a-matrix-is-singular
# https://stackoverflow.com/questions/62573178/pca-lda-analysis-r
##########################################################################

dati_sel_rel=dati_monosesquimiche 


##########################################################################
# Creazione delle matrici di analisi per PCA


X=dati_sel_rel
Y=dati_sel$Species

# X=asin(sqrt(X/100)) # eventuale trasfomazione dei dati normalizzati 


Xseed=dati_seeds[,4:8]
Yseed=dati_seeds$species

write.xlsx(list(data.frame(Y,X),data.frame(Yseed,Xseed)),"dati_PCA.xlsx")

cat("\014") 

#########################################################################################
# PCA explore variable and outlier detection

data_pca=data.frame(X,Species=Y)

abies_pca=ordinate(data_pca , cols = 1:28, model = ~ prcomp(., scale. = TRUE)) 

abies_pca %>%
  augment_ord() %>%
  ggbiplot(axis.type = "predictive") +
  theme_bw() +
  geom_rows_point(alpha = .3,aes(color=Species,shape=Species))+
  stat_rows_center(alpha = .8, fun.center = "mean",size = 5,aes(color = Species))+
  geom_origin() +
  ggtitle("Biplot PCA - Abies spp.") +
  labs(color = "Species")

ggsave("PCA_biplot_ord.png")

######################################################################################

color_cluster_dendro <- c("#FF0000", #cluster 1
                          "#669900", #cluster 2
                          "#00CCCC", #cluster 3
                          "#9933FF") #cluster 4


type=c(rep("monoterpenes",15),rep("sesquiterpenes",13))


pca_prcomp_abies <- prcomp(X,scale = T) #scaled PCA 

PCA_final_abies <- factoextra::fviz_pca_biplot(pca_prcomp_abies, 
                                                   # fill individuals by groups
                                                   geom.ind = c("point"),
                                                   col.ind = "black",
                                                   fill.ind = Y,
                                                   pointshape = 21, #i numeri definiscono  uno stile di punti!
                                                   palette = "color_cluster_dendro", 
                                                   addEllipses = T, 
                                                   ellipse.level = 0.10,
                                                   ellipse.type = "convex",
                                                   geom.var = c("arrow", "text"), 
                                                   arrowsize = 0.3,
                                                   labelsize = 2,
                                                   col.var = type,
                                                   legend.title = list(fill = "Species", color = "Compounds"),
                                                   title = "Biplot PCA - Abies spp. data",
                                                   repel = T  
) +
  ggpubr::color_palette("color_cluster_dendro") +
  theme(legend.direction = "horizontal") +
  theme(legend.position = "bottom") +
  theme(legend.box = "vertical")

PCA_final_abies

ggsave(filename="PCA_abies_novel.png")


#########################################################################################
# data selection by non parametric ANOVA o Kruskal Wallis ( optional) 

a=col_kruskalwallis(X,Y) # uso il pacchetto MatrixTests

row.names(a)[which(a$pvalue<0.05)] # all coumpound discriminate!
row.names(a)[which(a$pvalue>0.05)] 

write.xlsx(data.frame(a,compound_terpene=row.names(a)),"test_kruskal.wallis.xlsx")

##########################################################################################################
# Doing Posthoc analysis

res_posthoc=list()
res_posthoc_dunn=list()
res_GG=list()

##########

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

write.xlsx(list(dunn_test=res_df_dunn,t_test=res_df),"posthoc_kruskal.wallis.xlsx")

##############################################################################################
# summary tables

my_desc=psych::describeBy(mat_final_rel[,-1], mat_final_rel$species)

readr::write_excel_csv2(rbindlist(my_desc, idcol = 'species'), file = 'summary_terpene_by_species.xlsx')

######################################################################################
# boxplot main discriminant compounds

names(dati_sel_rel)[c(1,3,5,8,17,23,25)] # today another 2 terpinolene and b.phellanderene

dati_sel_box=dati_sel_rel
dati_sel_box$Species=factor(dati_sel$Species)
dati_sel_box=janitor::clean_names(dati_sel_box) # clean names of compounds check if they are ok

dir.create("boxplot_sel")

setwd("boxplot_sel")


dati_sel_box_short=dati_sel_box
dati_sel_box_short$species
levels(dati_sel_box_short$species)<-c("A. Alba","A. nebrodensis","A. pinsapo") # change names for legend

names(dati_sel_box_short)[6]="b.phellandrene"
ggboxplot(dati_sel_box_short,"species",names(dati_sel_box_short)[6],fill="red") +ylim(0,1)+ylab(paste0(names(dati_sel_box_short)[6]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box_short)[6],"_a.png"),width = 4,height = 3.5)

names(dati_sel_box_short)[9]="terpinolene"
ggboxplot(dati_sel_box_short,"species",names(dati_sel_box_short)[9],fill="red") +ylim(0,0.5)+ylab(paste0(names(dati_sel_box_short)[9]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box_short)[9],"_a.png"),width = 4,height = 3.5)


names(dati_sel_box_short)[1]="a.Pinene"
ggboxplot(dati_sel_box_short,"species",names(dati_sel_box_short)[1],fill="red") +ylim(0,50)+ylab(paste0(names(dati_sel_box_short)[1]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box_short)[1],".png"),width = 4,height = 3.5)


names(dati_sel_box_short)[3]="b.Pinene"            
ggboxplot(dati_sel_box_short,"species",names(dati_sel_box_short)[3],fill="red") +ylim(0,50)+ylab(paste0(names(dati_sel_box_short)[3]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box_short)[3],".png"),width = 4,height = 3.5)

names(dati_sel_box_short)[5]="Limonene"
ggboxplot(dati_sel_box_short,"species",names(dati_sel_box_short)[5],fill="red") +ylim(20,90)+ylab(paste0(names(dati_sel_box_short)[5]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box_short)[5],".png"),width = 4,height = 3.5)

# [4] "p-Cymene"
ggboxplot(dati_sel_box_short,"species",names(dati_sel_box_short)[8],fill="red") +ylim(0,3)+ylab(paste0(names(dati_sel_box_short)[8]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box)[8],".png"),width = 3.5,height = 3.5)

names(dati_sel_box_short)[17]="Sesquiterpeni.1"  
ggboxplot(dati_sel_box_short,"species",names(dati_sel_box_short)[17],fill="red") +ylim(0,5)+ylab(paste0(names(dati_sel_box_short)[17]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box_short)[17],".png"),width = 4,height = 3.5)

# "1.Germacrene.D-4-ol"
names(dati_sel_box_short)[23]="Germacrene.D.4.ol"  
ggboxplot(dati_sel_box_short,"species",names(dati_sel_box_short)[23],fill="red") +ylim(0,50)+ylab(paste0(names(dati_sel_box_short)[23]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box_short)[23],".png"),width = 4,height = 3.5)

# [7] "Selina-6-en-4-o
names(dati_sel_box_short)[25]="Selina.6.en.4.ol"
ggboxplot(dati_sel_box_short,"species",names(dati_sel_box_short)[25],fill="red") +ylim(0,50)+ylab(paste0(names(dati_sel_box_short)[25]," (%)"))
ggsave(paste0("boxplot_",names(dati_sel_box_short)[25],".png"),width = 4,height = 3.5)

setwd("..")



dir.create("boxplot")


#############################################################################################################

setwd("boxplot")

dati_sel_box=dati_sel_rel
dati_sel_box$Species=factor(dati_sel$Species)
dati_sel_box=janitor::clean_names(dati_sel_box) # pulisco i nomi

for ( i in 1:25) {

ggboxplot(dati_sel_box,"species",names(dati_sel_box)[i],fill="red") +ylim(0,max(dati_sel_box[i])+2)
  
ggsave(paste0("boxplot_",names(dati_sel_box)[i],".png"))

}

setwd("..")


####################################################################################
# Ordination Supervised methods : LDA 

# LDA
# Prior probabilities of groups: the proportion of training observations in each group. For example, there are 31% of the training observations in the setosa group
# Group means: group center of gravity. Shows the mean of each variable in each group.
# Coefficients of linear discriminants: Shows the linear combination of predictor variables that are used to form the LDA d

#########################################################################
# retrive selected data

X=dati_sel_rel_LDA
Y=dati_sel$Species

training.samples <- Y %>% createDataPartition(p = 0.8, list = FALSE) # create training and test data with caret R package

dati_train=X[training.samples, ]
dati_test=X[-training.samples, ]
dati_test$Y=Y[-training.samples]
dati_train$Y=Y[training.samples]

model <- lda(Y~., data =dati_train)

variance_explained=model$svd^2 / sum(model$svd^2)

# Wilks lamda testing

summary(manova(as.matrix(X[training.samples, ])~factor(dati_train$Y)),tol=0,test="Wilks")

summary(manova(as.matrix(X[training.samples, ])~factor(dati_train$Y)),tol=0)

################################################################
# LD plot on test data

predictions_lda <- model %>% predict(dati_test)

dataset = data.frame(Y = as.factor(dati_test$Y),
                     lda = predictions_lda$x)


dataset$Y=as.factor(gsub("Abies","A.",as.character(dataset$Y)))

centroids <- dataset |> group_by(Y) |> dplyr::summarise(centroid1 = mean(lda.LD1),
                                                        centroid2 = mean(lda.LD2))

p=ggplot(dataset) + 
  geom_point(aes(lda.LD1, lda.LD2, colour = Y, shape = Y), size = 2.5) + 
  geom_point(aes(centroid1, centroid2, colour = Y, shape = Y), size = 8, data = centroids) +
  geom_text(aes(centroid1, centroid2, label = Y), size = 5, data = centroids)+
  labs(title="LDA Biplot",subtitle="Test data",x="LD1",y="LD2")+guides(color = guide_legend("Species"),
                                                      shape = guide_legend("Species"))+
  xlim(-4, 5)

p # to visualise and save with RStudio

# ggord(model, dati_train$Y,arrow=F,repel =T)+labs(title="LDA Biplot",subtitle="Test data",x="LD1",y="LD2")+xlim(-40, 15)+ylim(-10, 15)+guides(color = guide_legend("Species"))

#########################################################################################################################################################################################
# exporting accuracy indexs and confusion matrix

model_accuracy_lda=mean(predictions_lda$class==dati_test$Y)

res_lda_abiesGC=confusionMatrix(predictions_lda$class,factor(dati_test$Y))

kable(as.table(res_lda_abiesGC$table), "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  cat(., file = "confusion_matrix_pops_abies.docx")

kable(as.table(round(res_lda_abiesGC$overall,2)), "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  cat(., file = "confusion_accuracy_pops_abies.docx")

kable(as.table(format(res_lda_abiesGC$byClass,digits=1)), "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  cat(., file = "confusion_acc_byClass_pops_abies.docx")


###########################################################################################################################################
# plots LDA with ordr packages

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
  ggtitle("Standardized LDA Biplot of Abies spp.",subtitle="Training data") +
  expand_limits(y = c(-3, 5))+
  labs(color = "Species")

ggsave("LDA_biplot.png")


################################################################################
# Ora sto lavorando sui dati dei semi

# kruskal wallis

col_oneway_welch(Xseed, Yseed) # utilizzo il pacchetto MatrixTests per fare l'ANOVA a 1 via 

summarySE(Xseed,"length","Yseed")
summarySE(Xseed,"mc","Yseed") 
summarySE(Xseed,"weight100","Yseed") 
summarySE(Xseed,"length","Yseed") 
summarySE(Xseed,"width","Yseed")

X=Xseed
Y=Yseed

################################################################################
# PCA seeds

ord <- PCA(Xseed, graph = FALSE)
ggord(ord, Yseed,arrow=NULL,txt=NULL)



PCA_final_abies_seeds <- factoextra::fviz_pca_biplot(ord, 
                                               # fill individuals by groups
                                               geom.ind = c("point"),
                                               col.ind = "black",
                                               fill.ind = Yseed,
                                               pointshape = 21, #i numeri definiscono  uno stile di punti!
                                               palette = "color_cluster_dendro", 
                                               addEllipses = T, 
                                               ellipse.level = 0.10,
                                               ellipse.type = "convex",
                                               geom.var = c("arrow", "text"), 
                                               arrowsize = 0.3,
                                               labelsize = 2,
                                               title = "Biplot PCA - Abies seeds spp. data",
                                               repel = T  
) +
  ggpubr::color_palette("color_cluster_dendro") +
  theme(legend.direction = "horizontal") +
  theme(legend.position = "bottom") +
  theme(legend.box = "vertical")

ggsave(filename="PCA_abies_seeds novel.png")

################################################################################
# LDA seeds

training.samples <- Y %>% createDataPartition(p = 0.8, list = FALSE)
dati_train=X[training.samples, ]
dati_test=X[-training.samples, ]
dati_test$Y=Y[-training.samples]
dati_train$Y=Y[training.samples]

model <- lda(Y~., data =dati_train)
predictions_lda <- model %>% predict(dati_test)
model_accuracy_lda=confusionMatrix(predictions_lda$class,factor(dati_test$Y))

###########################################################################################################################################
# plots DA biplot seeds

ggord(model, Y[training.samples],
      arrow=NULL,
      txt=NULL,
      xlims=c(-8,8),
      ylims =c(-5,7),
      obslab =F)

#####################################################################
# References

# https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab11_LDA_Model-assessment.html
# https://dicook.github.io/mulgar_book/14-lda.html
# https://stats.stackexchange.com/questions/560375/how-to-interpret-the-output-of-lda-discriminant-analysis-in-r
# https://cmdlinetips.com/2020/12/canonical-correlation-analysis-in-r/
# https://www.r-bloggers.com/2021/05/linear-discriminant-analysis-in-r/
# https://vitalflux.com/pca-vs-lda-differences-plots-examples/
# https://towardsai.net/p/data-science/lda-vs-pca
# https://stats.stackexchange.com/questions/23353/pca-lda-cca-and-pls
# https://www.geeksforgeeks.org/classifying-data-using-support-vector-machinessvms-in-r/
# https://mdatools.com/docs/pca.html
# https://rdrr.io/cran/mixOmics/man/plsda.html
# http://mixomics.org/methods/spls-da/
# https://stackoverflow.com/questions/46977743/pca-analysis-remove-centroid
