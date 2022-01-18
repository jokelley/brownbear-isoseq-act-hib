---
title: "Bear Tissue - PacBio IsoSeq vs Illumina RNAseq"
author: "Joanna Kelley"
date: "10 Jun 2021"
output:
  html_notebook:
    fig_caption: yes
    toc: yes
    toc_depth: 3
    toc_float: yes
  pdf_document:
    toc: yes
    toc_depth: '3'
  html_document:
    keep_md: TRUE
---

```{r}
library(dplyr)
library(ggplot2)
library(reshape)
library(corrplot)
library(edgeR)
```


```{r}
# --------------------
# Read in files
# --------------------
isoseq = read.csv("hq.no5merge.collapsed.filtered.mapped_fl_count-namesFixed.csv", header = T)
illumina = read.csv("AllBears-Subset_combined_est_counts_shortreads.csv", header = T)

# --------------------
# Merge the two datasets
# --------------------
combo = merge(isoseq, illumina, by.x = "id", by.y = "id")

# --------------------
# Split the id column to pull out the gene number
# --------------------
tmp <- unlist(strsplit(combo$id, '\\.'))
combo$gene <- tmp[seq(2,length(tmp),3)]


# --------------------
# log cpm of isoform matrix 
# --------------------
logcpm <- cpm(as.matrix(combo[,2:36]), log=TRUE, prior.count = 0.1)

# --------------------
# Pearson correlations of log cpm of isoform matrix 
# --------------------
logcpm.cor = cor(logcpm, method = "pearson")
corrplot(logcpm.cor, method="circle", type = "upper", diag = FALSE, cl.lim = c(-0.1,1), title = "Isoform level correlation", tl.cex = 0.5, tl.col = c(rep("tomato3",18), rep("darkgrey",17))) #, order = "hclust")

# --------------------
# log cpm of gene count matrix 
# --------------------
combo.gene = aggregate(. ~ gene, combo[,2:37], sum)
logcpm.gene <- cpm(as.matrix(combo.gene[,2:36]), log=TRUE, prior.count = 0.1)


# --------------------
# Pearson correlations of log cpm of gene count matrix 
# --------------------
logcpm.gene.cor = cor(logcpm.gene, method = "pearson")
corrplot(logcpm.gene.cor, method="circle", type = "upper", diag = FALSE, cl.lim = c(-0.01,1), title = "Gene level correlation", tl.cex = 0.5) #, order = "hclust", hclust.method = "median")

```

```{r}
# --------------------
# Correlations within categories
# --------------------
platform = c(rep("PB",18), rep("IL",17))
tissue = c(rep("adipose",6), rep("liver",6), rep("muscle",6),rep("adipose",6), rep("liver",5), rep("muscle",6))
season = c(rep("hib",2), rep("act",3), rep("hib",1), rep("hib",3), rep("act",3), rep("hib",3), rep("act",3), rep("hib",3), rep("act",3), rep("hib",3), rep("act",2), rep("hib",3), rep("act",3))
bear = c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,3,1,2,3,1,2,3)
info = as.data.frame(cbind(colnames(logcpm.gene.cor),platform, tissue, season,bear))

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cr  =(cormat)[ut]
    )
}
flat = flattenCorrMatrix(logcpm.gene.cor)

#flat = flattenCorrMatrix(logcpm.cor) # isoform matrix

by.cat1 = merge(flat,info,by.x = "row",by.y = "V1")
by.cat = merge(by.cat1,info,by.x = "column", by.y = "V1")

by.cat %>%
   group_by(platform.x, platform.y) %>% 
   summarise_at(vars("cr"), mean)

by.cat %>%
   group_by(platform.x,platform.y,tissue.x,tissue.y) %>% 
   summarise_at(vars("cr"), mean) %>% print(n = Inf)

by.cat %>%
   group_by(bear.x,bear.y) %>% 
   summarise_at(vars("cr"), mean)

by.cat %>%
   group_by(tissue.x,tissue.y) %>% 
   summarise_at(vars("cr"), mean)

by.cat %>%
    group_by(platform.x, platform.y,tissue.x,tissue.y,season.x,season.y) %>% 
    summarise_at(vars("cr"), mean) %>% print(n = Inf)

```

```{r}
liver = read.table("all-shortread/liver-tappAS_DIUGene_Proteins.tsv", sep="\t",header=T, comment.char = "!")

fat = read.table("all-shortread/fat-tappAS_DIUGene_Proteins.tsv", sep="\t",header=T, comment.char = "!")

muscle = read.table("all-shortread/muscle-tappAS_DIUGene_Proteins.tsv", sep="\t",header=T, comment.char = "!")

aa = merge(fat, liver, by.x = "X.Gene", by.y = "X.Gene")
bb = merge(aa, muscle, by.x = "X.Gene", by.y = "X.Gene")

diu.all = bb[bb$DIU.Result.x == "DIU" & bb$DIU.Result.y == "DIU" & bb$DIU.Result == "DIU",]
dim(diu.all) 
diu.all.isoform = diu.all[diu.all$Major.Isoform.Switching.x == "YES" & diu.all$Major.Isoform.Switching.y == "YES" & diu.all$Major.Isoform.Switching == "YES",]

rngtt = as.data.frame(read.csv("all-shortread/RNGTT-counts.csv", header = T, row.names = 1))
t.rngtt = as.data.frame(t(rngtt))
t.rngtt$names <- rownames(t.rngtt)

library(stringr)
t.rngtt$names <- str_split_fixed(t.rngtt$names, "_", 3) 
t.rngtt$season <- t.rngtt$names[,3]
t.rngtt$tissue <- t.rngtt$names[,2]

t.rngtt$combo <- paste(t.rngtt$season,t.rngtt$tissue)

yy = aggregate(cbind(PB.17293.1,PB.17293.2,PB.17293.6,PB.17293.13,PB.17293.12,PB.17293.15,PB.17293.16,PB.17293.26, PB.17293.29) ~ combo, data = t.rngtt, FUN = mean, NA.rm = T)

zz = aggregate(cbind(PB.17293.1,PB.17293.2,PB.17293.6,PB.17293.13,PB.17293.12,PB.17293.15,PB.17293.16,PB.17293.26, PB.17293.29) ~ season + tissue, data = t.rngtt, FUN = mean, NA.rm = T)

```

```{r}
dim(liver[liver$DIU.Result == "DIU" & liver$Major.Isoform.Switching == "YES",])
dim(liver[liver$DIU.Result == "DIU" & liver$Major.Isoform.Switching == "NO",])
dim(liver[liver$DIU.Result != "DIU" & liver$Major.Isoform.Switching == "YES",])
dim(liver[liver$DIU.Result != "DIU" & liver$Major.Isoform.Switching == "NO",])
dim(liver)

dim(fat[fat$DIU.Result == "DIU" & fat$Major.Isoform.Switching == "YES",])
dim(fat[fat$DIU.Result == "DIU" & fat$Major.Isoform.Switching == "NO",])
dim(fat[fat$DIU.Result != "DIU" & fat$Major.Isoform.Switching == "YES",])
dim(fat[fat$DIU.Result != "DIU" & fat$Major.Isoform.Switching == "NO",])
dim(fat)

dim(muscle[muscle$DIU.Result == "DIU" & muscle$Major.Isoform.Switching == "YES",])
dim(muscle[muscle$DIU.Result == "DIU" & muscle$Major.Isoform.Switching == "NO",])
dim(muscle[muscle$DIU.Result != "DIU" & muscle$Major.Isoform.Switching == "YES",])
dim(muscle[muscle$DIU.Result != "DIU" & muscle$Major.Isoform.Switching == "NO",])
dim(muscle)
```

```{r}
fat = read.csv("fat_tappAS_DIUGene_Proteins-0.1filter.csv", header = T)
liver = read.csv("liver_tappAS_DIUGene_Proteins-0.1filter.csv", header = T)
muscle = read.csv("muscle_tappAS_DIUGene_Proteins-0.1filter.csv", header = T)
aa = merge(fat, liver, by.x = "X.Gene", by.y = "X.Gene", all.x = T, all.y = T)
bb = merge(aa, muscle, by.x = "X.Gene", by.y = "X.Gene", all.x = T, all.y = T)
write.table(bb, file = "all_tappAS_DIUGene_Proteins.csv", quote = F, sep = ",")


length(which(fat$DIU.Result == "DIU" & fat$Major.Isoform.Switching == "YES"))/dim(fat)[1]
length(which(fat$DIU.Result == "DIU" & fat$Major.Isoform.Switching == "NO"))/dim(fat)[1]
length(which(fat$DIU.Result == "Not DIU" & fat$Major.Isoform.Switching == "YES"))/dim(fat)[1]
length(which(fat$DIU.Result == "Not DIU" & fat$Major.Isoform.Switching == "NO"))/dim(fat)[1]

length(which(liver$DIU.Result == "DIU" & liver$Major.Isoform.Switching == "YES"))/dim(liver)[1]
length(which(liver$DIU.Result == "DIU" & liver$Major.Isoform.Switching == "NO"))/dim(liver)[1]
length(which(liver$DIU.Result == "Not DIU" & liver$Major.Isoform.Switching == "YES"))/dim(liver)[1]
length(which(liver$DIU.Result == "Not DIU" & liver$Major.Isoform.Switching == "NO"))/dim(liver)[1]

length(which(muscle$DIU.Result == "DIU" & muscle$Major.Isoform.Switching == "YES"))/dim(muscle)[1]
length(which(muscle$DIU.Result == "DIU" & muscle$Major.Isoform.Switching == "NO"))/dim(muscle)[1]
length(which(muscle$DIU.Result == "Not DIU" & muscle$Major.Isoform.Switching == "YES"))/dim(muscle)[1]
length(which(muscle$DIU.Result == "Not DIU" & muscle$Major.Isoform.Switching == "NO"))/dim(muscle)[1]

dim(fat)
dim(liver)
dim(muscle)

```

```{r}
counts = read.table("new_merge_estcounts_shortreads_ids.txt")
samples = read.csv("sample_map.csv", header = T)
cpm.counts = cpm(counts)


example = as.data.frame(cpm.counts["PB.6860.2",])
example.combo = merge(example, samples, by.x = 0, by.y = "SampleID")
colnames(example.combo)[2] <- "cpm"

example = as.data.frame(cpm.counts["PB.6860.1",])
example.combo = merge(example, samples, by.x = 0, by.y = "SampleID")
colnames(example.combo)[2] <- "cpm"


# Box plot
e <- ggplot(example.combo, aes(x = Tissue, y = cpm))
e + geom_boxplot()

# Boxplot with points
e2 <- e + geom_boxplot(
  aes(fill = Season),
  position = position_dodge(0.9) 
  ) +
  scale_fill_manual(values = c("#999999", "#E69F00"))
e2

# Side by side Boxplot with points
e + geom_dotplot(
  binaxis = "y", stackdir = "center",
  fill = "lightgray"
  ) + 
  stat_summary(
    fun.data = "mean_sdl", fun.args = list(mult=1), 
    geom = "pointrange", color = "red"
    )

# points
e + geom_boxplot(
  aes(color = Season), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  ) +
  geom_dotplot(
    aes(fill = Season, color = Season), trim = FALSE,
    binaxis='y', stackdir='center', dotsize = 0.8,
    position = position_dodge(0.8)
  )+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))



e + geom_jitter(
  aes(shape = Season, color = Season), 
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  size = 1.2
  ) +
  stat_summary(
    aes(color = Season),
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4,
    position = position_dodge(0.8)
    )+
  scale_color_manual(values =  c("#00AFBB", "#E7B800"))

# Points (this is the one i like best)
e + geom_jitter(
  aes(shape = Season, color = Season), 
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  size = 1.5
  ) + ggtitle("Isoform PB.6860.1")

```

```{r}
setwd("~/Documents/Projects-KelleyLab/BrownBears/Isoseq/new_merge/")
counts.all = read.table("new_merge_estcounts_shortreads_ids.txt")
samples = read.csv("sample_map.csv", header = T)
counts = counts.all[,-which(colnames(counts.all) %in% c("all_CF3N"))]

cpm.counts = cpm(counts)



example1 = as.data.frame(cpm.counts["PB.6860.1",])
example1.combo = merge(example1, samples, by.x = 0, by.y = "SampleID")
colnames(example1.combo)[2] <- "cpm"

example2 = as.data.frame(cpm.counts["PB.6860.2",])
example2.combo = merge(example2, samples, by.x = 0, by.y = "SampleID")
colnames(example2.combo)[2] <- "cpm"

lim = max(cbind(example1,example2))

require(gridExtra)
plot1 <- d + geom_jitter(
  aes(shape = Season, color = Season), 
  position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
  size = 2
  ) + ggtitle("Isoform PB.6860.1") + theme(legend.position = "none") + scale_y_continuous(limits=c(0, lim))

plot2 <- e + geom_jitter(
  aes(shape = Season, color = Season), 
  position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
  size = 2
  ) + ggtitle("Isoform PB.6860.2") + theme(legend.position = c(0.75, 0.85)) + scale_y_continuous(limits=c(0, lim))

grid.arrange(plot1, plot2, ncol=2)

####
total = example1 + example2
total.combo = merge(total, samples, by.x = 0, by.y = "SampleID")
colnames(total.combo)[2] <- "cpm"

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

total.fat = total.combo[total.combo$Tissue == "Fat",]

df2 <- data_summary(total.fat, varname="cpm", 
                    groupnames=c("Season"))
# Convert dose to a factor variable
df2$Season=as.factor(df2$Season)
head(df2)

p<- ggplot(df2, aes(x=Season, y=cpm)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=0, ymax=cpm+sd), width=.2,
                 position=position_dodge(0.05))
print(p)

```



```{r}
example1 = as.data.frame(cpm.counts["PB.6860.1",])
example1.combo = merge(example1, samples, by.x = 0, by.y = "SampleID")
colnames(example1.combo)[2] <- "cpm"

d <- ggplot(example1.combo, aes(x = Tissue, y = cpm))

example2 = as.data.frame(cpm.counts["PB.6860.2",])
example2.combo = merge(example2, samples, by.x = 0, by.y = "SampleID")
colnames(example2.combo)[2] <- "cpm"

e <- ggplot(example2.combo, aes(x = Tissue, y = cpm))

lim = max(cbind(example1,example2))

require(gridExtra)
plot1 <- d + geom_jitter(
  aes(shape = Season, color = Season), 
  position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
  size = 2
  ) + ggtitle("Isoform PB.6860.1") + theme(legend.position = "none") + scale_y_continuous(limits=c(0, lim))

plot2 <- e + geom_jitter(
  aes(shape = Season, color = Season), 
  position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
  size = 2
  ) + ggtitle("Isoform PB.6860.2") + theme(legend.position = c(0.75, 0.85)) + scale_y_continuous(limits=c(0, lim))

grid.arrange(plot1, plot2, ncol=2)

####
total = counts["PB.6860.1",]+ counts["PB.6860.2",]
total = as.data.frame(total)
total.combo = merge(t(total), samples, by.x = 0, by.y = "SampleID")
colnames(total.combo)[2] <- "cpm"

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

total.fat = total.combo[total.combo$Tissue == "Fat",]

df2 <- data_summary(total.fat, varname="cpm", 
                    groupnames=c("Season"))
# Convert dose to a factor variable
df2$Season=as.factor(df2$Season)
head(df2)

p<- ggplot(df2, aes(x=Season, y=cpm)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=cpm-sd, ymax=cpm+sd), width=.2,
                 position=position_dodge(0.05))
print(p)
```

```{r}
## CRY2

example1 = as.data.frame(cpm.counts["PB.16373.1",])
example1.combo = merge(example1, samples, by.x = 0, by.y = "SampleID")
colnames(example1.combo)[2] <- "cpm"

d <- ggplot(example1.combo, aes(x = Tissue, y = cpm))

example2 = as.data.frame(cpm.counts["PB.16373.10",])
example2.combo = merge(example2, samples, by.x = 0, by.y = "SampleID")
colnames(example2.combo)[2] <- "cpm"

e <- ggplot(example2.combo, aes(x = Tissue, y = cpm))

lim = max(cbind(example1,example2))

require(gridExtra)
plot1 <- d + geom_jitter(
  aes(shape = Season, color = Season), 
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  size = 2
  ) + ggtitle("Isoform PB.16373.1") + theme(legend.position = c(0.75, 0.85)) + scale_y_continuous(limits=c(0, lim))

plot2 <- e + geom_jitter(
  aes(shape = Season, color = Season), 
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  size = 2
  ) + ggtitle("Isoform PB.16373.10") + theme(legend.position = "none")  + scale_y_continuous(limits=c(0, lim))

grid.arrange(plot1, plot2, ncol=2)

```


```{r}

# PRRX1 
prrx1 = as.data.frame(cpm.counts["PB.14470.1",]+cpm.counts["PB.14470.2",]+cpm.counts["PB.14470.3",]+cpm.counts["PB.14470.4",])
colnames(prrx1)[1] <- "cpm_prrx1"
COL6A3 = as.data.frame(cpm.counts["PB.2902.1",]+cpm.counts["PB.2902.2",]+cpm.counts["PB.2902.3",]+cpm.counts["PB.2902.4",]+cpm.counts["PB.2902.5",]+cpm.counts["PB.2902.6",]+cpm.counts["PB.2902.7",])
colnames(COL6A3)[1] <- "cpm_col6a3"

aa = merge(prrx1,COL6A3, by.x = 0, by.y = 0)
aa.combo = merge(aa, samples, by.x = "Row.names", by.y = "SampleID")

aa.fat = aa.combo[aa.combo$Tissue == "Fat",]

sp <- ggplot(aa.fat, aes(x = cpm_prrx1, y = cpm_col6a3, colour = factor(Season)))+ 
  geom_point(size=2.5)
sp
```
