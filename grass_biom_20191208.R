## this is a notebook to work through the microbial data for publication - I have started by copying what I did in the
# r script 'grass_97_biom_clean', and I will make any changes I see fit to better fit the narrative that I am trying to tell-
# this will also include reorganization of the script...

## currrent work 3.14.2019 - line 866, although I haven't dealt with the meta data yet...
## 3.15.2019 - I've started in with the ordinations and permanovas - line 1131
## 3.21.2019 made it through variable selection, right before variance partitioning - line 2383
## 9.18.2019 The paper is being written, I went through to finalize variable selection for use in variance partitioning
#            

setwd('/Users/grahambailes/grass_endophyte_community/manuscript/git_repo/grass-endophyte-community/')

library(phyloseq); packageVersion('phyloseq')
library(ggplot2); packageVersion('ggplot2')
library('DESeq2')
library('MASS')
library('plyr')
library('dplyr')
library('tidyr')
library('scales')
library('biom')
library('grid')
library('vegan')
library('ade4')
library('repr')
library('perturb')
library('geosphere')
library('adespatial')
library('SoDA')
library('png')
library('devtools')
#install_github("ggbiplot", "vqv") ## don't need to install from devtools every time
#install_github("vqv/ggbiplot")  ## ''
#install_github('fawda123/ggord') ## ''
library('ggord')
#devtools::install_github("gavinsimpson/ggvegan") ## ''
library('ggvegan')
library('ggbiplot') # for env_pca things
library('multcompView') #for betadisper Tukey's hsd
library('ggrepel')
library('car') # type III Anovas
library('corrplot')
library('bipartite')
library('lme4')
library('nlme')
library('cowplot')
library('ggpubr')
library('ggbeeswarm')

options(scipen=999) # turn off scientific notation

## here is a living list of previous and current environment images for this script.  
# load the most recent to have access to all objects created in the script

#save.image(file = "biom_2019_05_15.RData")
#load(file = 'biom_2019_05_15.RData')

#save.image(file = "biom_2019_08_06.RData")
#load(file = 'biom_2019_08_06.RData')

#save.image(file = "biom_2019_09_24.RData")
#load(file = 'biom_2019_09_24.RData')

## I decided to use ordistep model selection instead of the correlation.  
#save.image(file = "biom_2019_10_16.RData")
#load(file = 'biom_2019_10_16.RData')

# save.image(file = 'biom_2019_12_08.RData') 
# load(file = 'biom_2019_12_08.RData')

## this is the start of re-working statistics and analyses as per Bitty's and Dan's recomendations
# save.image(file = 'biom_2020_01_15.RData')
# load(file = 'biom_2020_01_15.RData')

################import my biom table created via my bioinformatic pipeline
grass_biom <-import_biom('grass_97_wmeta.biom')
grass_biom

###################################################################################################
##################################################################################################
################################################# Workflow

## this will hopefully be a usefull inclusion as I work toward getting these analyses ready for publication
# I would like to lay out my questions/hypotheses, and then create a roadmap of analyses to answer the questions


## this workbook is broken up into three sections

# 1) final data massage using controls to remove contaminants (line )
#    - investigae mock community - do I need to remove reads because of tag-switching
#    - calculate diversity metrics - 
#    - DEseq variance stabilization
#    -

# 2) Microbial composition and diversity (line )
#    - composition at several taxonomic levels by host (line)
#    - Species accumulation and diversity (line)
#    - species richness among sites (line)
#    - NMDS visualizations and PERMANOVA of microbial community composition among hosts, region, and sites (line)

# 3) Analysis of microbial communities in the context of environmental predictors (line )
#    - Examine and trim PRISM climate data (line )
#    - Examine and trim meta data, including biotic and edaphic predictors (line )
#    - creation of spatial variables (line )
#        - Mantel test to examine spatial autocorrelation of communities (Line )
#        - create PCNM variables.  These will examine more complex spatial structures, 
#          and can be used in ordination techniques (as opposed to raw lat/long data) (line )
#    - Reduction of multicolinearity of predictor variables (line )
#    - Variation Partitioning
#    - 
#    - 
# 4) examination of functional diversity - sym, path, sap (line )
#    - Anova of richness among sites
#    - 
#    - 

###################################################################################################
##################################################################################################
################################################# clean up biom table using controls

################## Nguyen, 2015 reconmends actions for negative and positive controls
## for the negative, they recomend subtracting the read abundance of those present in the negative from every other sample
## we are also interested in trimming out low-read samples, and in determining what to do about tag switching

############################################### low read samples
## lets start to look at the read abundance in our samples
grass_biom_eco <- subset_samples(grass_biom, site != "Control")
grass_biom_hiread <- prune_samples(sample_sums(grass_biom_eco)>2000, grass_biom_eco)
grass_biom_hiread

# check to see if it worked
sample_sums(grass_biom_hiread)[sample_sums(grass_biom_hiread) < 2000]
sample_sums(grass_biom_hiread)[sample_sums(grass_biom_hiread) > 0] ## display reads/sample

##Yep, lets see if any samples were removed
grass_biom_hiread #155 samples -down from 162.  this seems right, as there are 7 control samples (for both dan and I)
sample_names(grass_biom_hiread)

##################################################### remove contamination found in control samples
# what reads do we have in our controls?
grass_controls <- prune_samples(grass_biom@sam_data$samplename == 'negative'|
                                  grass_biom@sam_data$samplename == 'mock_g'|
                                  grass_biom@sam_data$samplename == 'mock_its'|
                                  grass_biom@sam_data$samplename == 'Pos_X', grass_biom)

grass_controls

# get maximums
grassContamination <- apply(otu_table(grass_controls),1,max)
sum(grassContamination>0) # 173 otus

sort(grassContamination[grassContamination>0])

## remove contaminaton
grassCont_otu <- otu_table(grass_biom_hiread) - grassContamination

##################################################### minimum abundances of observations
## I will assume that the mock community members I saw in my negative control represent contamination, given
# the high number of reads (mc 16 - 4807 reads, mc 6 - 1437 reads).  For tag-switching in the mock community 
# we see a max of 141 reads in the genomic community (next lowest is 62 reads.  This work is included in a different notebook), 
# and 38 reads in the its community.  If I were to follow Dan's notebook, I would remove 150 reads from every observation of every otu (or peraps 70 reads if I'm feeling conservative)
# from every sample, probably leading to a large loss of information.  Devin, the bioinformatics tech in the McGuire lab,
# suggested not doing this, and instead writing in a disclaimer saying that we observed a rate of 1-2% tag-switching
# from our mock community, then removing mock-reads from our biom table - this seems attractive, yet doesn't do as much (or really anything at all)
# to solve the tag-switching issue.  


## lets do a quick check to see what would happen if we were to remove 150 reads from all observations...

grassThresh_otu <- grassCont_otu - 150 ## subtract minimum threshold - we'll see how many reads we loose here
grassThresh_otu[grassThresh_otu<0] <-0 # bring negatives up to zero
grass.min.ab <- grass_biom_eco # create test biom for this
otu_table(grass.min.ab) <- grassThresh_otu

## We'll keep the 'minus contaminants/controls' biom for now too -
grassCont_otu[grassCont_otu<0]<-0 # we'll keep this biom going for comparison - bring negatives to zero
grass.cont.rem <- grass_biom_eco
otu_table(grass.cont.rem) <- grassCont_otu


# sanity check
dim(otu_table(grass_biom_hiread)) #3880  155
dim(otu_table(grass.cont.rem)) # 3880 155
dim(otu_table(grass.min.ab)) # 3880  155

all(rownames(otu_table(grass.min.ab)) %in% rownames(otu_table(grass_biom_hiread))) #TRUE
all(colnames(otu_table(grass.min.ab)) %in% colnames(otu_table(grass_biom_hiread)))#TRUE

all(rownames(otu_table(grass_biom_hiread)) %in% rownames(otu_table(grass.min.ab)))#TRUE
all(colnames(otu_table(grass_biom_hiread)) %in% colnames(otu_table(grass.min.ab)))#TRUE

## before we move on, lets see what we loose in removal of contaminants and minimum threshold
# reads - 1,795,344 lost
sum(sample_sums(grass_biom)) # 9388226
sum(sample_sums(grass.cont.rem)) # 8887690 for removal of controls
sum(sample_sums(grass.min.ab)) # 7592882

# lets look at observations
sum(otu_table(grass_biom)>0) # 36375 for raw biom table
sum(otu_table(grass.cont.rem)>0) # 33806 removal of controls only
sum(otu_table(grass.min.ab)>0) # 5083 removal of controls and 150 reads from each sample

# Given the loss in observations and reads, I have chosen to work with the less conservative biom, 

################################################# remove mock community members

## create file containing all instances of matched blast hits for the mock community member - we'll use this to remove all
# OTUs matching the mock members...

## this isn't quite right yet... the problem is that not all OTUs have info to species level, so columns don't match...
# I'll need to figure out which file to use to blast...
# sed '/>OTU/ s/;size=.*//g' otus_97_uclust_relable.fasta > otus_97.fasta
# blastn -query otus_97_uclust_relabel.fasta -db mcsanger.fasta -out mcblast_allOTUs.csv -outfmt 10 -max_target_seqs 1

options(max.print = 999999)
reads <- taxa_sums(grass.cont.rem)[taxa_sums(grass.cont.rem) > 0]
reads <- sort(reads, decreasing = T)
sink('all_97_names.txt')
names(reads)
sink()

## import 'all_97_names.txt' into atom, and use regex to remove line numbers (\[.*]), and insert commas
# create python script to integrate otu names with sequences ' mock_list_all_97'

#sed '/^>/ s/;size=.*//' mcseq_all_97.txt | sed '/^>/ s/;size=.*//' mcseq_all_97.txt > mockseqs_all_97.fasta
#blastn -query mockseqs_all_97.fasta -db mcsanger.fasta -out mcblast_all_97.txt -num_descriptions 3 -num_alignments 3
#blastn -query mockseqs_all_97.fasta -db mcsanger.fasta -out mcblast_all_97.csv -outfmt 10 -max_target_seqs 1

## can't get these final three to work yet - maybe its not important...
#sed -i '' '1 i\
#qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore
#' mcblast_all_97.csv -i
#sed -i '' 's/_ITS[1,4],/,/g' mcblast_all_97.csv
#sed -i '' 's/Sample//g' mcblast_all_97.csv


## load into R
blast <- read.csv('/Users/grahambailes/grass_endophyte_community/alt_biom/mock_community/mcblast_all_97.csv', header=T)
head(blast);dim(blast)

## keep only strong matches
goodblast <- blast[blast$pident > 95 & blast$length > 100,]

nrow(blast) 
nrow(goodblast)

## which of our rownames (= OTU names) are not in this list of strong matches?
pcotus <- !(rownames(otu_table(grass_biom)) %in% goodblast$qseqid) 
pcotus
## keep only these:
grass_biom_rem_mc <- prune_taxa(pcotus, grass.cont.rem)
grass_biom_rem_mc
sum(sample_sums(grass_biom_rem_mc)) #8721596
sum(otu_table(grass.cont.rem)>0) #33806

## and here is the minimum abundance biom
grass_biom_min_ab <- prune_taxa(pcotus, grass.min.ab)
grass_biom_min_ab
sum(sample_sums(grass_biom_min_ab)) #7419381
sum(otu_table(grass_biom_min_ab)>0) # 5034

##################################################################################################
################################################ diversity metrics

## Dan suggested that we need to perform our diversity statistics bevore we we variance stabilization - instead we will use rarifacation
# to eaven out our read depth on our samples
# at this moment (2020/01/15), I'm not sure whether or not this includes everything (such as # otus, ), or just analyses
# such as microbial richness.  we'll start with that I suppose...

## 
# examine rarefaction curve
rarecurve(t(otu_table(grass_biom_eco)), step=50, cex=0.5, label = F) # biom table minus controls
rarecurve_mc <- rarecurve(t(otu_table(grass_biom_rem_mc)), step=50, cex=0.5, label = F) 
tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/git_repo/grass-endophyte-community/rarecurve.tiff", width = 2100, height = 1500, units = "px", res = 300)
rarecurve(t(otu_table(grass_biom_rem_mc)), step=50, cex=0.5, label = F) 
dev.off()

min(sample_sums(grass_biom_eco)) # 7101
mean(sample_sums(grass_biom_eco)) # 58566.17
max(sample_sums(grass_biom_eco)) # 111141

min(sample_sums(grass_biom_rem_mc)) # 7007
mean(sample_sums(grass_biom_rem_mc)) # 56268.36
max(sample_sums(grass_biom_rem_mc)) # 110609

sort((sample_sums(grass_biom_eco))) 
sort((sample_sums(grass_biom_rem_mc)))

# looks like sample 16 has the lowest read number... looks like it is the 12th sample from french flat (F. roemeri) 
# I guess we'll rarify to 7100 reads for our diversity metrics.

# I wonder if the difference in sampling depth had an effect on the richness metrics I calculated on the var stabilized samples.
# for instance, of the 12 (# at each site) samples with the lowest read number, six were FF F. roemeri.

grass_biom_eco_rare <- rarefy_even_depth(grass_biom_rem_mc, rngseed=1, sample.size=0.9*min(sample_sums(grass_biom_rem_mc)), replace=F)
rarecurve(t(otu_table(grass_biom_eco_rare)), step=50, cex=0.5, label = F)

# boxplots of richness
festuca_biom_rare_div <- subset_samples(grass_biom_eco_rare, host == 'F. roemeri')
festuca_biom_rare_div@sam_data$site <- factor(festuca_biom_rare_div@sam_data$site, levels = c('French_flat', 'Roxy_Ann', 'Upper_Table','Hazel_Dell', 'Horse_Rock', 'Upper_Weir', 'Whidbey'))
f_richness_rare_box <-plot_richness(festuca_biom_rare_div, x='site', measures='Observed', color = 'site' )
f_richness_rare_box <- f_richness_rare_box + scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699')) + geom_boxplot() + theme(legend.position="none", plot.subtitle = element_blank())
f_richness_rare_box <- f_richness_rare_box +scale_x_discrete(breaks=c('French_flat', 'Roxy_Ann', 'Upper_Table','Hazel_Dell', 'Horse_Rock', 'Upper_Weir', 'Whidbey'),
                                                             labels=c('French flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey')) +  theme(panel.grid.minor=element_blank())
f_richness_rare_box 

plot_richness(festuca_biom_rare_div, x="site", measures=c("Observed")) + geom_boxplot()

# extract values for Anova
f_richness_rare_values <- f_richness_rare_box$data$value
f_richness_rare_sites <- f_richness_rare_box$data$site

f_richness_rare_sites <- gsub( 'French_flat', 'French Flat', f_richness_rare_sites)
f_richness_rare_sites <- gsub( 'Roxy_Ann', 'Roxy Ann', f_richness_rare_sites)
f_richness_rare_sites <- gsub( 'Upper_Table', 'Upper Table', f_richness_rare_sites)
f_richness_rare_sites <- gsub( 'Lower_Table', 'Lower Table', f_richness_rare_sites)
f_richness_rare_sites <- gsub( 'Hazel_Dell', 'Hazel Dell', f_richness_rare_sites)
f_richness_rare_sites <- gsub( 'Horse_Rock', 'Horse Rock', f_richness_rare_sites)
f_richness_rare_sites <- gsub( 'Upper_Weir', 'Upper Weir', f_richness_rare_sites)

f_richness_rare_stats <- as.data.frame(f_richness_rare_sites)
f_richness_rare_stats$observations <- f_richness_rare_values

f_richness_stats_non_s_rare <- f_richness_rare_stats[-c(1:12),] ## without FF

## statistics - ANOVA
f_richness_rare_anova <- Anova(lm(observations~f_richness_rare_sites, data = f_richness_rare_stats), type = 3)
f_richness_rare_anova
#                         Sum Sq Df  F value  Pr(>F)    
#  (Intercept)           200208  1 230.9433 <2e-16 ***
#  f_richness_rare_sites   8407  6   1.6162  0.154    
#  Residuals              66752 77 
confint(f_richness_rare_anova)


# and ANOVA without french flat
f_richness_non_s_rare_anova <- Anova(lm(observations~f_richness_rare_sites, data = f_richness_stats_non_s_rare),type=3)
f_richness_non_s_rare_anova ## without ff, results are not significant
#                         Sum Sq Df  F value Pr(>F)    
#  (Intercept)           230187  1 257.9221 <2e-16 ***
#  f_richness_rare_sites   8302  5   1.8604 0.1133    
#  Residuals              58903 66 3778  

## tukey's post-hoc comparisons
library(emmeans)
f_richness_rare_tukey <- emmeans(lm(observations~f_richness_rare_sites, data = f_richness_rare_stats), pairwise ~ f_richness_rare_sites)  ##tukey test
plot(f_richness_rare_tukey, comparisons = T)
CLD(f_richness_rare_tukey$emmeans)


## richness plot showing mean and se
f_richness_rare_plot <- ggplot(data = f_richness_rare_stats, aes(x = f_richness_rare_sites, y = observations)) + 
  geom_boxplot(aes(fill = f_richness_rare_sites), outlier.shape = NA, alpha = 0.9)+ 
  scale_fill_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'))+ 
  geom_beeswarm(dodge.width = 0.75, size = 0.6, groupOnX = T)+
  
  theme(legend.position = 'none') + theme(axis.text.x= element_text(size = 12, angle=45, hjust = 1),axis.ticks = element_blank())

f_richness_rare_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/git_repo/grass-endophyte-community/f_richness_rarefied.tiff", width = 2100, height = 1500, units = "px", res = 300)
f_richness_rare_plot
dev.off()

# richness without FF
f_rich_summary_non_s <- ddply(f_richness_stats_non_s, c('f_richness_sites'), summarise,
                              N    = length(values), mean = mean(values),
                              sd   = sd(values),se   = sd / sqrt(N))

f_rich_summary_non_s$site <- c('Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey')
f_rich_summary_non_s$site <- factor(f_rich_summary_non_s$site, levels =  c('Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'))

## richness plot showing mean and se
f_richness_non_s_plot <- ggplot(f_rich_summary_non_s, aes(x=site, y=mean, colour=site)) + 
  scale_colour_manual(values=c('#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) + theme_bw() +
  geom_point(size = 1.2) + theme(legend.position="none", axis.text = element_text(size = 12, angle = 0)) + labs(x='Site', y='# OTUs') 

f_richness_non_s_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/f_richness_non_s.tiff", width = 2100, height = 1500, units = "px", res = 300)
f_richness_non_s_plot
dev.off()

## Danthonia
danthonia_biom_rare_div <- subset_samples(grass_biom_eco_rare, host == 'D. californica')
danthonia_biom_rare_div@sam_data$site <- factor(danthonia_biom_rare_div@sam_data$site, levels = c('French_flat', 'Whetstone', 'Lower_Table','Hazel_Dell', 'Horse_Rock', 'Whidbey'))
d_richness_rare_box <-plot_richness(danthonia_biom_rare_div, x='site', measures='Observed', color = 'site' )
d_richness_rare_box <- d_richness_rare_box + scale_colour_manual(values=c('#F8766D', '#93AA00', '#00BA38', '#00B9E3', '#619CFF', '#FF61C3')) + geom_boxplot() + theme(legend.position="none")
d_richness_rare_box <- d_richness_rare_box +scale_x_discrete(breaks=c('French_flat', 'Whetstone', 'Lower_Table','Hazel_Dell', 'Horse_Rock', 'Whidbey'),
                                                             labels=c('French Flat', 'Whetstone', 'Lower Table','Hazel Dell', 'Horse Rock', 'Whidbey')) + theme(panel.grid.minor=element_blank(),
                                                                                                                                                                panel.grid.major=element_blank())
d_richness_rare_box

d_richness_rare_values <- d_richness_rare_box$data$value
d_richness_rare_sites <- d_richness_rare_box$data$site

d_richness_rare_sites <- gsub( 'French_flat', 'French Flat', d_richness_rare_sites)
d_richness_rare_sites <- gsub( 'Lower_Table', 'Lower Table', d_richness_rare_sites)
d_richness_rare_sites <- gsub( 'Hazel_Dell', 'Hazel Dell', d_richness_rare_sites)
d_richness_rare_sites <- gsub( 'Horse_Rock', 'Horse Rock', d_richness_rare_sites)
d_richness_rare_sites <- gsub( 'Upper_Weir', 'Upper Weir', d_richness_rare_sites)

d_richness_rare_stats <- as.data.frame(d_richness_rare_sites)
d_richness_rare_stats$observations <- d_richness_rare_values
d_richness_stats_non_s_rare <- d_richness_rare_stats[-c(1:12),] ## without FF

d_richness_rare_anova <- Anova(lm(observations~d_richness_rare_sites, data = d_richness_rare_stats), type = 3)
(d_richness_rare_anova)   
#                    Sum Sq Df  F value Pr(>F)    
#  (Intercept)           118803  1 209.5124 < 2.2e-16 ***
#  d_richness_rare_sites  16196  5   5.7126 0.0001986 ***
#  Residuals              36858 65 

d_richness_anova_rare_non_s <- Anova(lm(observations~d_richness_rare_sites, data = d_richness_stats_non_s_rare), type = 3)
d_richness_anova_rare_non_s  # without ff, slightly more significant

#                   Sum Sq Df  F value Pr(>F)    
#  (Intercept)           148964  1 281.2124 < 2.2e-16 ***
#  d_richness_rare_sites  11374  4   5.3678  0.001031 ** 
#  Residuals              28605 54 

## data summary for figure 

d_richness_rare_tukey <- emmeans(lm(observations~d_richness_rare_sites, data = d_richness_rare_stats), pairwise ~ d_richness_rare_sites)  ##tukey test
plot(d_richness_rare_tukey, comparisons = T)
CLD(d_richness_rare_tukey$emmeans)


d_richness_rare_stats$site <- c(rep(c('French Flat', 'Whetstone', 'Lower Table','Hazel Dell'),each = 12),rep('Horse Rock',11), rep('Whidbey',12))
d_richness_rare_stats$site <- factor(d_richness_rare_stats$site, levels =  c('French Flat', 'Whetstone', 'Lower Table','Hazel Dell', 'Horse Rock', 'Whidbey'))

## richness plot showing mean and se
d_richness_rare_lettering <- c('A', 'AB', 'A', 'AB', 'B', 'B')


d_richness_rare_plot <- ggplot(data = d_richness_rare_stats, aes(x =site, y = observations)) + 
  geom_boxplot(aes(fill = d_richness_rare_sites), outlier.shape = NA, alpha = 0.9)+ 
  scale_fill_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'))+ 
  geom_beeswarm(dodge.width = 0.75, size = 0.6, groupOnX = T)+
  stat_summary(geom = 'text', label = d_richness_rare_lettering, fun.y =mean, vjust = -6) + 
  theme(legend.position = 'none') + theme(axis.text.x= element_text(size = 12, angle=45, hjust = 1),axis.ticks = element_blank()) +
  theme(axis.title.x = element_blank())

d_richness_rare_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Manuscript/git_repo/grass-endophyte-community/d_richness_rarefied.tiff", width = 2100, height = 1500, units = "px", res = 300)
d_richness_rare_plot
dev.off()

##################################################################################################
################################################# variance stabilization

## Most (every) study I am familiar with utilize some type of normalization for read variance among samples
## I would like to first examine the reads among all samples.  I'm not clear on what constitutes
# the need for variance stabilization.  Further, I don't know how to deal with the question of var.stab vs rarification


## Roo has written a function to stream-line the process of using deseq variance stabilization, 
## applied to its use with OTU tables.  DESeq2 was originally written for RNA-seq data, with dense matrices. 
## Species matrices are usually pretty sparse, so deseq needs a few workarounds to handle our kind of data. 
## Thanks Roo!

source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")

library(DESeq2)
deseq_varstab <- source('DESeq_varstab.txt')
save(file = 'Deseq_verstab.txt', mode = 'function')

## once the function is loaded into R, I can call it on my biom table.  
## It requires two inputs - a phyloseq object and a design variable.  According to Roo, this is an artifact from the 
## package's RNAseq past.  In theory, it shouldn't affect our douwnstream analysis.
set.seed(1)
grass_biom_vs_cont_rem <- DESeq_varstab(grass_biom_rem_mc, ~site)
grass_biom_vs_min_ab <- DESeq_varstab(grass_biom_min_ab, ~ site)
## converting counts to integer mode
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## Warning message:
##  In DESeqDataSet(se, design = design, ignoreRank) :
##  some variables in design formula are characters, converting to factors

## Okay, lets see if our sample number or taxa numbers have changed

grass_biom_hiread## OTU Table: [ 3881 taxa and 162 samples ]
grass_biom_vs_cont_rem ## OTU Table: [ 3713 taxa and 155 samples ]
grass_biom_vs_min_ab ## OTU Table: [ 764 taxa and 155 samples ]

## looks like a loss of ~135 taxa for the cont.rem biom, and 3109 for the min.ab biom

## okay, I'll need to write this variance stabalized biom to a biom file for the funguild assignment.
## I'll use the 'write_biom' function for this
new_vs_table <- cbind(grass_biom_vs_cont_rem@otu_table, grass_biom_vs_cont_rem@tax_table)
new_vs_table <-unite(as.data.frame(new_vs_table), taxonomy, c(Rank1,Rank2,Rank3,Rank4,Rank5,Rank6,Rank7), remove = F, sep = ';')
write.csv(new_vs_table, '/Users/grahambailes/grass_endophyte_community/alt_biom/grass_biom_vs_fun.csv')
# write_biom(grass_biom_vs, '/Users/grahambailes/grass_endophyte_community/alt_biom/grass_biom_vs.biom')
# biom_file <- system.file("grass_biom_vs", "/Users/grahambailes/grass_endophyte_community/alt_biom/grass_biom_vs.biom", package = "biomformat")

## first, I'll look at the variation in read number among my samples.  I've decided to go the rout of normalizing
## using DESeq2, per recomentation by roo, as well as the McMurdie, Holmes 2014 paper.
par(mfrow = c(2,2))
options(repr.plot.width = 10, repr.plot.height = 4)
barplot(sample_sums(grass_biom_hiread), axisnames=FALSE)
mtext(text='Festuca', side = 1, at=50)
mtext(text='Danthonia', side = 1, at=140)

## variance stabalized biome table
options(repr.plot.width = 10, repr.plot.height = 4)
barplot(sample_sums(grass_biom_vs_cont_rem), axisnames=FALSE)
mtext(text='Festuca', side = 1, at=50)
mtext(text='Danthonia', side = 1, at=140)

## variance stabalized biome table - 
options(repr.plot.width = 10, repr.plot.height = 4)
barplot(sample_sums(grass_biom_vs_min_ab), axisnames=FALSE)
mtext(text='Festuca', side = 1, at=50)
mtext(text='Danthonia', side = 1, at=140)


# rename taxa rankings
colnames(tax_table(grass_biom_vs_hiread)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
colnames(tax_table(grass_biom_vs_cont_rem)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
colnames(tax_table(grass_biom_vs_min_ab)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

##################################################################################################################################################
################################################# Subset otu matricies by host  ####################################

## At this point it may be usefull to create two new subsetted biom tables for each host.
## Many of my analyses will be used to examine communities of each individual host.  
#  My analyses will be performed using presence-absence data, so I'll convert to PA.
# I'll be using Jacard distances for my dissimilarity matricies, and apparently hellinger standardizations
# for any linear analyses I perform - I guess that Hellinger standardizations are better than raw data for ecological
# analyses because they are asymetric (ignore double zeros, or basing significance when two samples don't have a perticular species)
# 

grass_mat <- t(otu_table(grass_biom_vs_cont_rem)@.Data) # transpose to display as samples are rows
grass_mat <- grass_mat[,colSums(grass_mat)>0] # select taxa with non-zeroabundance
grass_mat_PA <- grass_mat
grass_mat_PA[grass_mat_PA > 0] <- 1 # convert to presence-absence
grass_mat_PA[1:10,1:10]
grass_mat_hel <- decostand(grass_mat_PA, 'hellinger') ## apply hellinger standardization
## looking at this in 2019, it looks like hellinger standardization is used 
grass_mat_hel[1:10,1:10]

festuca_biom <- subset_samples(grass_biom_vs_cont_rem, host == 'F. roemeri') ##subset to host
festuca_biom_pa <- festuca_biom
festuca_biom_pa[festuca_biom_pa@otu_table@.Data > 0] <- 1
festuca_mat <- t(otu_table(festuca_biom)@.Data) ##transpose table
festuca_mat <- festuca_mat[,colSums(festuca_mat) > 0] ##remove any empty columns
f_mat_PA <- festuca_mat ##create new object to be transformed incase I want to use abundance data
f_mat_PA[f_mat_PA > 0] <- 1 ## convert to presence/absence
danthonia_biom <- subset_samples(grass_biom_vs_cont_rem, host == 'D. californica')
danthonia_mat <- t(otu_table(danthonia_biom)@.Data) 
danthonia_mat <- danthonia_mat[,colSums(danthonia_mat) > 0] ## remove any empty columns
d_mat_PA <- danthonia_mat 
d_mat_PA[d_mat_PA > 0] <- 1 ## convert to presence/absence


## look at new matricies
festuca_mat[1:10, 1:5]
f_mat_PA[1:10,1:5]
danthonia_mat[1:10, 1:5]
d_mat_PA[1:10,1:5]

## create standardized community distance matricies from P/A matricies.
# 'The  chord  and  Hellinger  transformations  appear  to  be  the  best  for  general  use.
# Legendre  &  Gallagher  (2001)  showed  that  the  values  of  the  corresponding  distances are 
# monotonically increasing across a simulated ecological gradient and are maximallyrelated (R2) to the 
# spatial distances along the geographic gradient'
## for use with RDA or other tests with linear assumptions - 
f_mat_hel <- decostand(f_mat_PA, 'hellinger') ## standardization 
d_mat_hel <- decostand(d_mat_PA, 'hellinger') ## standardization

## distance matricies using jaccard distance for presence/absence data
f_mat_dist <- vegdist(f_mat_PA, method ="jaccard", binary = T)
d_mat_dist <- vegdist(d_mat_PA, method ="jaccard", binary = T)


############################################################################################################
############################################################################################################
############################################################ microbial data analysis

## Okay, we have our biom tables cleaned up and straitened out - now we will move on to data analysis.
# the first section will involve analyzing 

############################################################################################################
############################################################################################################
################################################ composition and diversity

###########################################################################################################
################################################ Taxa composition

## here we're looking at the distribution of different levels of taxa among host species.
# we'll make figures for the distributions of phyla, class, order, and family 

## convert otu tables to presence absence, insert back into phyloseq object, and remove empty taxa 
festucaPA <- subset_samples(grass_biom_vs_cont_rem, host == 'F. roemeri')
danthoniaPA <-  subset_samples(grass_biom_vs_cont_rem, host == 'D. californica')
festucaOTU <- otu_table(festucaPA) ## extract otu table matrix
festucaOTU[festucaOTU > 0] <- 1 ## convert to P/A
festucaOTU -> otu_table(festucaPA) ## put it back in
danthoniaOTU <- otu_table(danthoniaPA) ## same with wood
danthoniaOTU[danthoniaOTU > 0] <- 1
danthoniaOTU -> otu_table(danthoniaPA)
festucaPA <- prune_taxa(rowSums(otu_table(festucaPA)) > 0, festucaPA) ## get rid of empty taxa - loss of 887 taxa
danthoniaPA <- prune_taxa(rowSums(otu_table(danthoniaPA)) > 0, danthoniaPA) # loss of 1411 taxa

library(VennDiagram)
aa <- names(taxa_sums(festucaPA))
bb <- names(taxa_sums(danthoniaPA))

venn.plot <- venn.diagram(
  x=list(
    'F. roemeri' =aa,
    'D. californica' =bb),
  fill = c("#00BFC4", "#F8766D"),
  alpha = c(0.6, 0.8),
  cex = 2,
  cat.fontface = 4,
  cat.dist = 0.05,
  cat.just = list(c(0, -0.5), c(1, -0.5)),
  filename = NULL)

grid.draw(venn.plot)

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/git_repo/grass-endophyte-community/OTU_overlap.tiff", width = 1800, height = 1700, units = "px", res = 300)
grid.draw(venn.plot)
dev.off()

venn.diagram(list(Festuca = aa, Danthonia = bb),
             fill = c("#00BFC4", "#F8766D"),
             alpha = c(0.3, 0.7),
             cex = 2,
             cat.fontface = 4,
             fontfamily =3,imagetype = 'png',
             filename='overallvenn.png')

overlap = calculate.overlap(
  x=list(
    'F. roemeri' =aa,
    'D. californica' =bb
  )
) 

length(overlap[[3]])

venn.plot

## look at phyla level taxonomy among hosty species

festucaphyplot <- table(tax_table(festucaPA)[,"Phylum"])/nrow(otu_table(festucaPA))
danthoniaphyplot <- table(tax_table(danthoniaPA)[,"Phylum"])/nrow(otu_table(danthoniaPA))
bothphyplot <- as.table(rbind(festucaphyplot, danthoniaphyplot))
bothphyplot <- as.data.frame(bothphyplot)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
phyplot <- ggplot(bothphyplot, aes(x =Var1, y=Freq, fill=Var2)) + geom_bar(stat = 'identity', position = 'fill')+
  scale_fill_brewer(palette = 'Dark2')
phyplot

## or by phyla
bothphyplot <- as.table(rbind(festucaphyplot, danthoniaphyplot))
colnames(bothphyplot) <- gsub('p__*','',colnames(bothphyplot))
options(repr.plot.width = 10,repr.plot.height = 6)
par(mar = c(12,3,2,2))
phyplot2 <-barplot(bothphyplot,las=2, cex.names=1.3, ylim = c(0,1), beside=TRUE,col = c("#00BFC4", "#F8766D"))
legend('topright', box.lwd = 0,pt.cex = 0.5, xjust = 1,yjust = 1,
       fill = c("#00BFC4", "#F8766D"),cex = 1,
       legend= c(expression(italic('F. roemeri')),expression(italic('D. californica'))))
phyplot2

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/grass_phyla.tiff", width = 2100, height = 1500, units = "px", res = 300)
par(mar=c(12,3,2,2))
barplot(bothphyplot,las=2, cex.names=1.3, ylim = c(0,1), beside=TRUE,col = c("#00BFC4", "#F8766D"))
legend('topright', box.lwd = 0,pt.cex = 0.4, xjust = 1,yjust = 1,
       fill = c("#00BFC4", "#F8766D"),cex = 0.75,
       legend= c(expression(italic('F. roemeri')),expression(italic('D. californica'))))
dev.off()

festucaclasstable <- table(tax_table(festucaPA)[,"Class"])/nrow(otu_table(festucaPA))
danthoniaclasstable <- table(tax_table(danthoniaPA)[,"Class"])/nrow(otu_table(danthoniaPA))
festucaclassvector <- as.vector(festucaclasstable); names(festucaclassvector) <- names(festucaclasstable)
danthoniaclassvector <- as.vector(danthoniaclasstable); names(danthoniaclassvector) <- names(danthoniaclasstable)
## which classes that are in festuca are not observed in danthonia? and vice-versa:
notindanthoniaclassvector <- names(festucaclassvector)[!(names(festucaclassvector) %in% names(danthoniaclassvector))]
notinfestucaclassvector <- names(danthoniaclassvector)[!(names(danthoniaclassvector) %in% names(festucaclassvector))]
## match up the membership and order of these vectors so we can rbind them 
notindanthonia <- vector(length=length(notindanthoniaclassvector))
notindanthonia[] <- 0; names(notindanthonia) <- notindanthoniaclassvector
fulldanthoniaclassvector <- c(danthoniaclassvector,notindanthonia)
notinfestuca <- vector(length=length(notinfestucaclassvector))
notinfestuca[] <- 0; names(notinfestuca) <- notinfestucaclassvector
fullfestucaclassvector <- c(festucaclassvector,notinfestuca)
fulldanthoniaclassvector <- fulldanthoniaclassvector[names(fullfestucaclassvector)] ## match order
## put them in a matrix together
bothclass <- rbind(fullfestucaclassvector, fulldanthoniaclassvector)
bothclass_df <- as.data.frame(t(bothclass))
colnames(bothclass_df) <- c('F. roemeri', 'D. californica')
rownames(bothclass_df) <- gsub('c__*','',rownames(bothclass_df))
rownames(bothclass_df) <- gsub('*_cls_Incertae_sedis','_sp',rownames(bothclass_df))
bothclass_df <- bothclass_df[order(-bothclass_df$`F. roemeri`),]
bothclass_sorted <- t(bothclass_df)
barplot(bothclass_sorted,las=2, cex.names=1.3, ylim = c(0,0.5), beside=T,col = c("#00BFC4", "#F8766D"))
legend('topright', box.lwd = 0,pt.cex = 0.4, xjust = 1,yjust = 1,
       fill = c("#00BFC4", "#F8766D"),cex = 0.75,
       legend= c(expression(italic('F. roemeri')),expression(italic('D. californica'))))
#barplot(bothclasstable, 
tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/grass_class.tiff", width = 2100, height = 1500, units = "px", res = 300)
par(mar=c(12,3,4,2))
barplot(bothclass_sorted,las=2, ylim = c(0,0.5), beside=TRUE, col = c("#00BFC4", "#F8766D"))
legend('topright', box.lwd = 0,pt.cex = 0.4, xjust = 1,yjust = 1,
       fill = c("#00BFC4", "#F8766D"),cex = 0.75,
       legend= c(expression(italic('F. roemeri')),expression(italic('D. californica'))))
dev.off()


festucaordertable <- table(tax_table(festucaPA)[,"Order"])/nrow(otu_table(festucaPA))
danthoniaordertable <- table(tax_table(danthoniaPA)[,"Order"])/nrow(otu_table(danthoniaPA))
festucaordervector <- as.vector(festucaordertable); names(festucaordervector) <- names(festucaordertable)
danthoniaordervector <- as.vector(danthoniaordertable); names(danthoniaordervector) <- names(danthoniaordertable)
## which orderes that are in festuca are not observed in danthonia? and vice-versa:
notindanthoniaordervector <- names(festucaordervector)[!(names(festucaordervector) %in% names(danthoniaordervector))]
notinfestucaordervector <- names(danthoniaordervector)[!(names(danthoniaordervector) %in% names(festucaordervector))]
## match up the membership and order of these vectors so we can rbind them 
notindanthonia <- vector(length=length(notindanthoniaordervector))
notindanthonia[] <- 0; names(notindanthonia) <- notindanthoniaordervector
fulldanthoniaordervector <- c(danthoniaordervector,notindanthonia)
notinfestuca <- vector(length=length(notinfestucaordervector))
notinfestuca[] <- 0; names(notinfestuca) <- notinfestucaordervector
fullfestucaordervector <- c(festucaordervector,notinfestuca)
fulldanthoniaordervector <- fulldanthoniaordervector[names(fullfestucaordervector)] ## match order

bothorder <- (rbind(fullfestucaordervector, fulldanthoniaordervector))
bothorder_df <- as.data.frame(t(bothorder))
colnames(bothorder_df) <- c('F. roemeri', 'D. californica')
rownames(bothorder_df) <- gsub('o__*','',rownames(bothorder_df))
rownames(bothorder_df) <- gsub('*_ord_Incertae_sedis','_sp',rownames(bothorder_df))
bothorder_df <-bothorder_df[which(bothorder_df$`F. roemeri` >0.009 | bothorder_df$`D. californica`>0.009),]
bothorder_df <- bothorder_df[order(-bothorder_df$`F. roemeri`),]
bothorder_sorted <- t(bothorder_df)
barplot(bothorder_sorted,las=2, ylim = c(0,0.25), beside=TRUE, col = c("#00BFC4", "#F8766D"))
legend('topright', box.lwd = 0,pt.cex = 0.4, xjust = 1,yjust = 1,
       fill = c("#00BFC4", "#F8766D"),cex = 0.75,
       legend= c(expression(italic('F. roemeri')),expression(italic('D. californica'))))

festucafamilytable <- table(tax_table(festucaPA)[,"Family"])/nrow(otu_table(festucaPA))
danthoniafamilytable <- table(tax_table(danthoniaPA)[,"Family"])/nrow(otu_table(danthoniaPA))
festucafamilyvector <- as.vector(festucafamilytable); names(festucafamilyvector) <- names(festucafamilytable)
danthoniafamilyvector <- as.vector(danthoniafamilytable); names(danthoniafamilyvector) <- names(danthoniafamilytable)
## which familyes that are in festuca are not observed in danthonia? and vice-versa:
notindanthoniafamilyvector <- names(festucafamilyvector)[!(names(festucafamilyvector) %in% names(danthoniafamilyvector))]
notinfestucafamilyvector <- names(danthoniafamilyvector)[!(names(danthoniafamilyvector) %in% names(festucafamilyvector))]
## match up the membership and order of these vectors so we can rbind them 
notindanthonia <- vector(length=length(notindanthoniafamilyvector))
notindanthonia[] <- 0; names(notindanthonia) <- notindanthoniafamilyvector
fulldanthoniafamilyvector <- c(danthoniafamilyvector,notindanthonia)
notinfestuca <- vector(length=length(notinfestucafamilyvector))
notinfestuca[] <- 0; names(notinfestuca) <- notinfestucafamilyvector
fullfestucafamilyvector <- c(festucafamilyvector,notinfestuca)
fulldanthoniafamilyvector <- fulldanthoniafamilyvector[names(fullfestucafamilyvector)] ## match order
## put them in a matrix together
bothfamily <- rbind(fullfestucafamilyvector, fulldanthoniafamilyvector)
bothfamily_df <- as.data.frame(t(bothfamily))
colnames(bothfamily_df) <- c('F. roemeri', 'D. californica')
rownames(bothfamily_df) <- gsub('f__*','',rownames(bothfamily_df))
rownames(bothfamily_df) <- gsub('*_fam_Incertae_sedis','_sp',rownames(bothfamily_df))
bothfamily_df <- bothfamily_df[order(-bothfamily_df$`F. roemeri`),]
bothfamily_df <-bothfamily_df[which(bothfamily_df$`F. roemeri` >0.009 | bothfamily_df$`D. californica`>0.009),]
bothfamily_df <- bothfamily_df[order(-bothfamily_df$`F. roemeri`),]
bothfamily_sorted <- t(bothfamily_df)
barplot(bothfamily_sorted,las=2, ylim = c(0,0.1), beside=TRUE, col = c("#00BFC4", "#F8766D"))
legend('topright',
       fill = c("#00BFC4", "#F8766D"),cex = 1,
       legend= c(expression(italic('F. roemeri')),expression(italic('D. californica'))))

#barplot(bothfamilytable, 
tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/grass_family.tiff", width = 2100, height = 1500, units = "px", res = 300)
par(mar=c(10,3,4,2))
barplot(bothfamily_sorted,las=2, ylim = c(0,0.4), beside=TRUE, col = c("#00BFC4", "#F8766D"))
legend('topright',
       fill = c("#00BFC4", "#F8766D"),cex = 1,
       legend= c(expression(italic('F. roemeri')),expression(italic('D. californica'))))
dev.off()

## order figure
tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/grass_order.tiff", width = 2100, height = 1500, units = "px", res = 300)
par(mar=c(10,3,4,2))
barplot(bothorder_sorted,las=2, ylim = c(0,0.5), beside=TRUE, col = c("#00BFC4", "#F8766D"))
legend('topright',
       fill = c("#00BFC4", "#F8766D"),cex = 1,
       legend= c(expression(italic('F. roemeri')),expression(italic('D. californica'))))
dev.off()

#########################################################################################
#############################################species accumulation 

festuca_mat_1 <- as.data.frame(festuca_mat)
rownames(festuca_mat_1) <- f_environ$samplename
festuca_mat_1 <-  decostand(festuca_mat_1, 'hellinger')

danthonia_mat_1 <- as.data.frame(danthonia_mat)
rownames(danthonia_mat_1) <- d_environ$samplename
danthonia_mat_1 <- decostand(danthonia_mat_1, 'hellinger')

##species accumulation curve D. californica and f_roemeri on same figure
f_accum <-specaccum(festuca_mat_1, method = 'random')
d_accum <- specaccum(danthonia_mat_1, method = 'random')
fd_specaccum <-plot(f_accum, add=FALSE, ci.type = 'polygon', ci.col = 'grey', ylab = 'OTUs recovered', xlab = 'Samples', ylim = c(0,3000))
fd_specaccum <-plot(d_accum, ci.type = 'polygon', ci.col = 'grey40', ylab = 'OTUs recovered', xlab = 'Samples', add=TRUE)
fd_specaccum <-legend(x = 'bottomright', box.lwd = 0,pt.cex = 0.4, xjust = 1,yjust = 1,
                      legend = c(expression(italic('F. roemeri')), expression(italic('D. californica'))),fill = c('grey', 'grey40'), bty = 'n')

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/spec_accum.tiff", width = 3200, height = 3200, units = "px", res = 800)
fd_specaccum <-plot(f_accum, add=FALSE, ci.type = 'polygon', ci.col = 'grey', ylab = 'OTUs recovered', xlab = 'Samples', ylim = c(0,3000))
fd_specaccum <-plot(d_accum, ci.type = 'polygon', ci.col = 'grey40', ylab = 'OTUs recovered', xlab = 'Samples', add=TRUE)
fd_specaccum <-legend(x = 'bottomright', legend = c(expression(italic('F. roemeri')), expression(italic('D. californica'))),fill = c('grey', 'grey40'), bty = 'n', xjust=1)
dev.off()

f_specpool <- specpool(festuca_mat_1)
#      Species  chao   chao.se   jack1 jack1.se    jack2     boot  boot.se  n
# All   2826 3497.009 59.3622 3691.571 108.3336 4001.812 3238.327 63.39326 84

d_specpool <- specpool(danthonia_mat_1)
#      Species   chao   chao.se    jack1 jack1.se    jack2     boot  boot.se  n
# All    2302 3857.941 132.0162 3350.028 144.9751 4040.139 2747.699 66.95441 71

f_pool_accum <- poolaccum(festuca_mat_1)
d_pool_accum <- poolaccum(danthonia_mat_1)

plot(f_pool_accum)
plot(d_pool_accum)

taxa_sums(festuca_biom_pa)
sample_sums(festuca_biom_pa)


################################################################################################ 
################################################################################################ 
######################################## NMDS of microbial community data

## for our ordinations, we will use distance measures calculated using Jacard distances on poresence/absence data
# we'll look at all OTUs in respect to host, region, and site, then we will break it down to examine communities of 
# each host 

################################################################################################
################################################ Grass NMDS visualization

grass_biom_vs_cont_rem@sam_data$site <- gsub( 'French_flat', 'French Flat', grass_biom_vs_cont_rem@sam_data$site)
grass_biom_vs_cont_rem@sam_data$site <- gsub( 'Roxy_Ann', 'Roxy Ann', grass_biom_vs_cont_rem@sam_data$site)
grass_biom_vs_cont_rem@sam_data$site <- gsub( 'Upper_Table', 'Upper Table', grass_biom_vs_cont_rem@sam_data$site)
grass_biom_vs_cont_rem@sam_data$site <- gsub( 'Lower_Table', 'Lower Table', grass_biom_vs_cont_rem@sam_data$site)
grass_biom_vs_cont_rem@sam_data$site <- gsub( 'Hazel_Dell', 'Hazel Dell', grass_biom_vs_cont_rem@sam_data$site)
grass_biom_vs_cont_rem@sam_data$site <- gsub( 'Horse_Rock', 'Horse Rock', grass_biom_vs_cont_rem@sam_data$site)
grass_biom_vs_cont_rem@sam_data$site <- gsub( 'Upper_Weir', 'Upper Weir', grass_biom_vs_cont_rem@sam_data$site)

grass_biom_vs_cont_rem@sam_data$site <- factor(grass_biom_vs_cont_rem@sam_data$site, levels = c('French Flat', 'Roxy Ann', 'Whetstone', 'Lower Table', 'Upper Table',
                                                                                                'Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'))
levels(grass_biom_vs_cont_rem@sam_data$site) <- c('French Flat', 'Roxy Ann', 'Whetstone', 'Lower Table', 'Upper Table',
                                                  'Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey')


################################################ presence/absence data - communities from both hosts

## create a distance matrix for all grass samples
# this can be used when a function doesn't automatically compute a distance matrix
grass_dist_pa <- vegdist(grass_mat_PA, method = 'jaccard')

set.seed(1)
## we'll run the nmds on our presence/absence matrix, using jaccard distance
grass_ord <- metaMDS(grass_mat_PA, distance = 'jaccard', try = 20, trymax = 100)
grass_ord # stress is 0.219013 - not bad but not great...
stressplot(grass_ord)

## lets graph our 
fd_ord <- ordinate(grass_biom_vs_cont_rem, 'NMDS', 'jaccard')
fd_host <- plot_ordination(grass_biom_vs_cont_rem, fd_ord, color = "host")
fd_host <- fd_host + theme_bw() 
fd_host <- fd_host + theme(legend.text = element_text(face = "italic"))
plot(fd_host)

fd_host <- plot_ordination(grass_biom_vs_cont_rem, grass_ord, color = "host")
fd_host <- fd_host + theme_bw() + stat_ellipse() 
fd_host <- fd_host + theme(legend.text = element_text(face = "italic"))+
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 14))
plot(fd_host)

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/git_repo/grass-endophyte-community/grass_host_nmds.tiff", width = 2100, height = 1500, units = "px", res = 300)
fd_host
dev.off()

fd_host + facet_wrap(~site, 3) + theme(legend.position = 'none') ## figure showing seperation of host communities by site

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/git_repo/grass-endophyte-community/grass_host_nmds_facet.tiff", width = 2100, height = 1500, units = "px", res = 300)
fd_host + facet_wrap(~site, 3) + theme(legend.position = 'none')
dev.off()

adonis(grass_dist_pa ~ host, data = environmental, permutations = 10000, method = 'jaccard')

#Call:
#  adonis(formula = grass_dist_pa ~ host, data = environmental) 

# Permutation: free
# Number of permutations: 10000

#Terms added sequentially (first to last)

#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  host        1     1.953 1.95346  5.6684 0.03573  0.001 ***
#  Residuals 153    52.727 0.34462         0.96427           
#  Total     154    54.680                 1.00000   

## significan, but a pretty poor model fit - lets use the betadisper to determine group dispersion
# if we have uneven dispersion, this result should be taken with a grain of salt...

host_mod <- betadisper(grass_dist_pa, group = environmental$host, type = c('median', 'centroid'))
permutest(host_mod)

#             Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups      1 0.000763 0.00076306 1.1288    999  0.281
# Residuals 153 0.103429 0.00067600    

plot(host_mod)
boxplot(host_mod)

## taken together, these results suggest that we don't have uneven dispersion, and thus can trust our 
# permanova, despite the low R2.

## region plot
## first, I want to change the names to the actual ecoregions...
grass_biom_vs_cont_rem@sam_data$region <- gsub('SOR','Klamath Mtns',grass_biom_vs_cont_rem@sam_data$region)
grass_biom_vs_cont_rem@sam_data$region <- gsub('COR','Willamette Valley',grass_biom_vs_cont_rem@sam_data$region)
grass_biom_vs_cont_rem@sam_data$region <- gsub('WA','Puget Lowlands',grass_biom_vs_cont_rem@sam_data$region)

grass_biom_vs_cont_rem@sam_data$region <- factor(grass_biom_vs_cont_rem@sam_data$region, levels = c('Klamath Mtns', 'Willamette Valley', 'Puget Lowlands'))
fd_region <- plot_ordination(grass_biom_vs_cont_rem, grass_ord, color = "region")
fd_region <- fd_region + theme_bw() + stat_ellipse() + scale_colour_manual(values=c('#009900','#663399','#006699'), guide = guide_legend(reverse=TRUE)) +
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 14))
plot(fd_region)

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/grass_region_nmds.tiff", width = 2100, height = 1500, units = "px", res = 300)
fd_region
dev.off()

## region stats

adonis(grass_dist_pa ~ region, data = environmental, permutations = 10000, method = 'jaccard')

#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# region      2     3.48 1.74008  5.1658 0.06365  0.001 ***
# Residuals 152    51.52 0.33684         0.93635           
# Total     154    54.68                 1.00000 

region_mod <- betadisper(grass_dist_pa, group = environmental$region, type = c('median', 'centroid'))
permutest(region_mod) 

#             Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      2 0.002044 0.0010219 1.0176    999  0.355
# Residuals 152 0.152645 0.0010042        

plot(region_mod)
boxplot(region_mod)

## again, it looks like we don't have unequal dispersion here -

# I need to group the treatments that are not different each other together.

generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
site_mod$group <- factor(site_mod$group, levels = c('French_flat', 'Roxy_Ann', 'Whetstone', 'Lower_Table', 'Upper_Table',
                                                    'Hazel_Dell', 'Horse_Rock', 'Upper_Weir', 'Whidbey'))
LABELS=generate_label_df(site_mod_hsd , 'group')
my_colors <- c('grey20', 'grey35', 'grey50', 'grey65', 'grey80', 'grey95')
my_colors <- c('firebrick3', 'gold3', 'darkolivegreen3', 'forestgreen', 'turquoise', 'cadetblue3', 'deepskyblue4', 'darkorchid2', 'magenta1')
# I want to write the letter over each box. Over is how high I want to write it.
par(mar=c(8,8,1,1))
site_disp <- boxplot(site_mod, ylim=c(0.4,0.7),  col=my_colors, las=2) ## create figure colors represent tukey's designations
over=0.05*max( site_disp$stats[nrow(site_disp$stats),] ) ## text position
text( c(1:nlevels(site_mod$group)) , site_disp$stats[nrow(site_disp$stats),]+over , LABELS[,1]) ## add lables


##########################################################################  
##################################### festuca visualization

festuca_biom@sam_data$site <- gsub( 'French_flat', 'French Flat', festuca_biom@sam_data$site)
festuca_biom@sam_data$site <- gsub( 'Roxy_Ann', 'Roxy Ann', festuca_biom@sam_data$site)
festuca_biom@sam_data$site <- gsub( 'Upper_Table', 'Upper Table', festuca_biom@sam_data$site)
festuca_biom@sam_data$site <- gsub( 'Hazel_Dell', 'Hazel Dell', festuca_biom@sam_data$site)
festuca_biom@sam_data$site <- gsub( 'Horse_Rock', 'Horse Rock', festuca_biom@sam_data$site)
festuca_biom@sam_data$site <- gsub( 'Upper_Weir', 'Upper Weir', festuca_biom@sam_data$site)

festuca_biom@sam_data$site <- factor(festuca_biom@sam_data$site, levels = c('French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'))
festuca_biom@sam_data$region <- factor(festuca_biom@sam_data$region, levels = c('SOR', 'COR', 'WA'))
####################################### presence absence data - Festuca
## create ordination 
festuca_ord <- metaMDS(f_mat_PA, distance = 'jaccard', try = 20, trymax = 100) ## NMDS ordination
stressplot(festuca_ord)
# stress = 0.2109064
# non-metric R2 = 0.956, linear fit R@ = 0.792

festuca_dist_pa <- vegdist(f_mat_PA, method = 'jaccard') ## distance matrix

## region PERMANOVA and dispersion test 
f_region <- plot_ordination(festuca_biom, festuca_ord, color = 'region')
f_region <- f_region + theme_bw() + stat_ellipse(size=0.8) + scale_colour_manual(values=c('#009900','#663399','#006699'),
                                                                                 guide = guide_legend(reverse=TRUE))+
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 14))
plot(f_region)

f_region_adonis <- adonis(festuca_dist_pa ~ region, data = f_environ, permutations = 10000, by = 'jaccard')
f_region_adonis
#           Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)    
#  region     2    2.0711 1.03553  3.1891 0.073  0.001 ***
#  Residuals 81   26.3013 0.32471         0.927           
#  Total     83   28.3724                 1.000           

f_region_mod <- betadisper(festuca_dist_pa, group = f_environ$region, type = c('median', 'centroid'))
permutest(region_mod) 
#             Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      2 0.002044 0.0010219 1.0176    999  0.389
# Residuals 152 0.152645 0.0010042     

## site PERMANOVA and dispersion test
f_site <- plot_ordination(festuca_biom, festuca_ord, color = "site")
f_site <- f_site + theme_bw() + stat_ellipse(size=0.8) + scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'),
                                                                             labels=c('French Flat', 'Roxy Ann', 'Upper Table', "Hazel Dell", 'Horse Rock', 'Upper Weir', 'Whidbey'),
                                                                             guide = guide_legend(reverse=TRUE))
f_site <- f_site + theme(legend.text=element_text(size=13), legend.title = element_text(size = 14))
plot(f_site)

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/git_repo/grass-endophyte-community/f_site_nmds.tiff", width = 2100, height = 1500, units = "px", res = 300)
f_site
dev.off()

## permanova and permutation tests. also tukey's hsd for dispersion
f_site_adonis <-adonis(festuca_dist_pa ~ site, data = f_environ, permutations = 10000, by = 'jaccard')
f_site_adonis
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  site       6    5.6871 0.94785  3.2173 0.20045  0.001 ***
#  Residuals 77   22.6853 0.29461         0.79955           
#  Total     83   28.3724                 1.00000           

## check dispersion
f_site_mod <- betadisper(festuca_dist_pa, group = f_environ$site, type = c('median', 'centroid'))
permutest(f_site_mod)
#             Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#  Groups     6 0.010281 0.0017134 1.5478    999   0.17
#  Residuals 77 0.085238 0.0011070                     

##################################################################################################### 
###################################Danthonia visualization

danthonia_biom@sam_data$site <- factor(danthonia_biom@sam_data$site, levels = c('French_flat', 'Whetstone', 'Lower_Table','Hazel_Dell', 'Horse_Rock', 'Whidbey'))
danthonia_biom@sam_data$region <- factor(danthonia_biom@sam_data$region, levels = c('SOR', 'COR', 'WA'))
danthonia_dist_pa <- vegdist(d_mat_PA, method = 'jaccard') ## distance matrix

################################### Presence/absence - danthonia

danthonia_ord <- metaMDS(d_mat_PA, distance = 'jaccard', try = 20, trymax = 100) ## NMDS ordination
danthonia_ord # Stress:     0.1857302 
stressplot(danthonia_ord)
# non-metric R2: 0.966, linear R2: 0.842

d_region <- plot_ordination(danthonia_biom, danthonia_ord, color = 'region')
d_region <- d_region + theme_bw() + stat_ellipse(size = 0.8) + scale_colour_manual(values=c('#009900','#663399','#006699'),
                                                                                   guide = guide_legend(reverse=TRUE))
d_region <- d_region + theme(legend.text=element_text(size=13), legend.title = element_text(size = 14))
plot(d_region)

## Region PERMANOVA and dispersion test
d_region_adonis <- adonis(danthonia_dist_pa ~ region, data = d_environ)
d_region_adonis
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  region     2    3.9658 1.98290  6.6133 0.16284  0.001 ***
#  Residuals 68   20.3888 0.29983         0.83716           
#  Total     70   24.3546                 1.00000     

d_region_mod <- betadisper(danthonia_dist_pa, group = d_environ$region, type = c('median', 'centroid'))
permutest(d_region_mod)

#             Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#  Groups     2 0.081127 0.040563 38.503    999  0.001 ***
#  Residuals 68 0.071639 0.001054                         

plot(d_region_mod)
boxplot(d_region_mod)
## okay, so we see both a significant PERMANOVA and non-homogenous dispersion.  given the figure and statistics,
# it looks like we are seeing both non-homogeneous dispersion as well as an actuall difference in centroids
# 

d_region_mod_hsd <- TukeyHSD(d_region_mod)
plot(d_region_mod_hsd, las=2)


d_site <- plot_ordination(danthonia_biom, danthonia_ord, color = "site")
d_site <- d_site + theme_bw() + stat_ellipse(size = 0.8)+ scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                                                                              labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                                                                              guide = guide_legend(reverse=TRUE))
d_site <- d_site + theme(legend.text=element_text(size=13), legend.title = element_text(size = 14))
plot(d_site)

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/git_repo/grass-endophyte-community/d_site_nmds.tiff", width = 2100, height = 1500, units = "px", res = 300)
d_site
dev.off()

d_site_adonis <- adonis(danthonia_dist_pa ~ site, data = d_environ, permutations = 10000, by = 'jaccard')
d_site_adonis

#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  site       5     6.435 1.28699  4.6683 0.26422  0.001 ***
#  Residuals 65    17.920 0.27569         0.73578           
#  Total     70    24.355                 1.00000         

d_site_mod <- betadisper(danthonia_dist_pa, group = d_environ$site, type = c('median', 'centroid'))
permutest(d_site_mod)

#             Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)   
#  Groups     5 0.043665 0.0087330 4.6759    999  0.002 **
#  Residuals 65 0.121399 0.0018677   

plot(d_site_mod)
boxplot(d_site_mod)

## okay, so we see both a significant PERMANOVA and non-homogenous dispersion.  given the figure and statistics,
# it looks like we are seeing both non-homogeneous dispersion as well as an actuall difference in centroids


d_site_mod_hsd <- TukeyHSD(d_site_mod)
d_site_mod_hsd
plot(d_site_mod_hsd, las=2)

## looks like Whidbey is significantly different in dispersion from French Flat and Hazel Dell
d_site_df <- cbind(d_site_mod$group, d_site_mod$distances)


d_site_df <- cbind(melt(d_site_mod$group), melt(d_site_mod$distances)) ## combine data
colnames(d_site_df) <- c('site', 'value')
d_site_summary <- ddply(d_site_df, c('site'), summarise,
                        N    = length(value),
                        mean = mean(value),
                        sd   = sd(value),
                        se   = sd / sqrt(N))

## something's not right here, but I don't know if I care right now...
levels(d_site_summary$site) <- c('French Flat', 'Whetstone', 'L. Table Rock','Hazel Dell', 'Horse Rock', 'Whidbey')
d_site_summary$site <- factor(d_site_summary$site, levels = c('French Flat', 'Whetstone', 'L. Table Rock','Hazel Dell', 'Horse Rock', 'Whidbey'))
d_site_plot <- ggplot(d_site_summary, aes(x=site, y=mean)) + 
  scale_color_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'))+
  geom_errorbar(data = d_site_summary, aes(ymin=mean-se, ymax=mean+se, color = my_colors),width=.2) + theme_bw()+
  geom_point( size = 1.2) + theme(legend.position="none", axis.text = element_text(size=12)) + labs(x='Site', y='distance from centroid') #+
#geom_text(aes(label=f_rich_hsd_lab$Letters, colour = NULL, y=f_rich_summary$mean+(f_rich_summary$se/1)), vjust=-1.5)
d_site_plot

d_site_mod$group <- factor(d_site_mod$group, levels = c('French Flat', 'Whetstone', 'L.Table Rock','Hazel Dell', 'Horse Rock', 'Whidbey'))
d_LABELS=generate_label_df(d_site_mod_hsd , 'group')
my_colors <- c('grey20', 'grey35', 'grey50', 'grey65', 'grey80', 'grey95')
my_colors <- c('firebrick3', 'gold3', 'darkolivegreen3', 'forestgreen', 'turquoise', 'cadetblue3', 'deepskyblue4', 'darkorchid2', 'magenta1')
my_colors <- c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699')

# I want to write the letter over each box. Over is how high I want to write it.
par(mar=c(8,8,1,1))
d_site_disp <- boxplot(d_site_mod, ylim=c(0.25,0.65),  col=my_colors, las=2) ## create figure colors represent tukey's designations
over=0.05*max( d_site_disp$stats[nrow(d_site_disp$stats),] ) ## text position
text( c(1:nlevels(d_site_mod$group)) , d_site_disp$stats[nrow(d_site_disp$stats),]+over , d_LABELS[,1]) ## add lables

## making multi-pannel figures
NMDS_multi <- ggarrange(fd_host, fd_region, f_site, d_site, ncol = 2, nrow = 2, labels = c('A','B','C','D'))
NMDS_multi


plot_grid(fd_host, fd_region, f_site, d_site, align = 'hv')

# note I changed the legend size to 10 for the multi pannel figure...

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/nmds_multi_pannel.tiff", width = 2000, height = 1000, units = "px", res = 200)
NMDS_multi
dev.off()

## I'm going to skip over the envfit section here too - it is cool, but I don't expect to use that for my paper...

####################################################################################################
####################################################################################################
################################################## Ecological data analysis 

## here is where we will begin to examine our communities in the context of environmental predictors.
# we'll begin by examining and organizing our climatic and meta data.

#####################################################################################################
################################################### climatic data (Prism data)
## the climate data we are using is from the PRISM climate group at OSU - we chose to use the 30 year normals 
## (baseline datasets describing average monthly and annual conditions over the most recent three full decades -> 1981 - 2010)
## Lauren graciously extracted data for each of our sites
## I'll first take the prism data and subset sites by species, then by data type - I've also decided to include some additional calculated data-
## these are bioclimatic predictors described by the USGS, which highlight climate conditions best related to species physiology, and thus may be more predictive of communities

## ************************************************** note, this data was not used in any analysis!!
## I used the 2014-2015 seasonal data for my anaysis, but have left this here incase I changed my mind and decided to use normals

prism <- read.csv('PRISM_data.csv', stringsAsFactors = F)
## describe sites
f_prism <- prism[which(prism$SPECIES == 'Festuca roemeri'),]
d_prism <- prism[which(prism$SPECIES == 'Danthonia californica'),]

f_prism$SITE <- c("French Flat", "Roxy Ann", "Upper Table", "Hazel Dell", "Horse Rock",
                  "Upper Weir","Whidbey")
levels(f_prism$SITE) <- c("French Flat", "Roxy Ann", "Upper Table", "Hazel Dell", "Horse Rock",
                          "Upper Weir","Whidbey")

d_prism$SITE <- c('French Flat', 'Lower Table', 'Whetstone', 'Hazel Dell', 'Horse Rock','Whidbey')
d_prism$SITE <- factor(d_prism$SITE, levels = c('French Flat', 'Lower Table', 'Whetstone', 'Hazel Dell', 'Horse Rock','Whidbey'))

## bioclimate predictors
f_bioclim <- f_prism[,c(2,22:24,64:66,80:82,96:98)]
d_bioclim <- d_prism[,c(2,22:24,64:66,80:82,96:98)]
## festuca site precipitation
f_prism_precip <- (f_prism[,9:20])
rownames(f_prism_precip) <- f_prism$SITE
colnames(f_prism_precip) <-c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                             'September', 'October', 'November', 'December')

## festuca site temperature
f_prism_temp <- (f_prism[,38:49])
rownames(f_prism_temp) <- f_prism$SITE
colnames(f_prism_temp) <-c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                           'September', 'October', 'November', 'December')

## festuca dewpoint temperature
f_prism_dpt <- (f_prism[,67:78])
rownames(f_prism_dpt) <- f_prism$SITE
colnames(f_prism_dpt) <-c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                          'September', 'October', 'November', 'December')

## festuca dewpoint temperature
f_prism_RH <- (f_prism[,83:94])
rownames(f_prism_RH) <- f_prism$SITE
colnames(f_prism_RH) <-c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                         'September', 'October', 'November', 'December')

## data danthonia site precipitation
d_prism_precip <- (d_prism[,9:20])
rownames(d_prism_precip) <- d_prism$SITE
colnames(d_prism_precip) <-c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                             'September', 'October', 'November', 'December')

## danthonia site temperature
d_prism_temp <- (d_prism[,38:49])
rownames(d_prism_temp) <- d_prism$SITE
colnames(d_prism_temp) <-c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                           'September', 'October', 'November', 'December')

## danthonia dewpoint temperature
d_prism_dpt <- (d_prism[,67:78])
rownames(d_prism_dpt) <- d_prism$SITE
colnames(d_prism_dpt) <-c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                          'September', 'October', 'November', 'December')

## danthonia dewpoint temperature
d_prism_RH <- (d_prism[,83:94])
rownames(d_prism_RH) <- d_prism$SITE
colnames(d_prism_RH) <-c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                         'September', 'October', 'November', 'December')

library(ggplot2)
library(reshape2)
library(tibble)

#################################################### climate figures
## precipitation
f_precip <- as.data.frame(t(f_prism_precip)) ## new data frame for precipitation plot
f_precip <- rownames_to_column(f_precip, var="month")
f_precip <- melt(f_precip, id='month')
f_precip$month <- factor(f_precip$month, levels = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                                                    'September', 'October', 'November', 'December'))

## festuca precipitation plot for thesis
f_precip_plot <- ggplot(f_precip,aes(x=month,y=value,colour=variable,group=variable)) + geom_line( size=0.6) + theme_bw() +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'), name = 'site', guide = guide_legend(reverse=TRUE)) +
  labs(y='precipitation (cm)')+theme(legend.text=element_text(size=13), legend.title = element_text(size = 14), 
                                     axis.text.x = element_text(angle = 90, size = 12), axis.text.y = element_text(size = 12))

f_precip_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/f_precip.tiff", width = 2100, height = 1500, units = "px", res = 300)
f_precip_plot
dev.off()

## temperature
f_temp  <- as.data.frame(t(f_prism_temp)) ## new data frame for precipitation plot
f_temp <- rownames_to_column(f_temp, var="month")
f_temp <- melt(f_temp, id='month')
f_temp$month <-  c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                   'September', 'October', 'November', 'December')
f_temp$month <- factor(f_temp$month, levels = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                                                'September', 'October', 'November', 'December'))

f_temp_plot <- ggplot(f_temp, aes(x = month, y = value, colour = variable, group = variable)) + geom_line( size=0.6) + theme_bw() +
  scale_colour_manual( values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'), name = 'site',guide = guide_legend(reverse=TRUE))  +
  labs(y=(expression(~degree~C)))+ #+ stat_summary(geom="ribbon", fun.ymin=min, fun.ymax=max, aes(fill=type), alpha=0.3)  
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 14), 
        axis.text.x = element_text(angle = 90, size = 13),axis.text.y = element_text(size = 13), axis.title.y = element_text(angle = 0,vjust=0.5))

f_temp_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/f_temp.tiff", width = 2100, height = 1500, units = "px", res = 300)
f_temp_plot
dev.off()

## dew-point temperature
f_dpt  <- as.data.frame(t(f_prism_dpt)) ## new data frame for precipitation plot
f_dpt <- rownames_to_column(f_dpt, var="month")
f_dpt <- melt(f_dpt, id='month')
f_dpt$month <-  c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                  'September', 'October', 'November', 'December')
f_dpt$month <- factor(f_dpt$month, levels = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                                              'September', 'October', 'November', 'December'))

f_dpt_plot <- ggplot(f_dpt, aes(x = month, y = value, colour = variable, group = variable)) + geom_line( size=0.6) + theme_bw() +
  scale_colour_manual( values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'), name = 'site',guide = guide_legend(reverse=F))  +
  labs(y=(expression(~degree~C)))+ #+ stat_summary(geom="ribbon", fun.ymin=min, fun.ymax=max, aes(fill=type), alpha=0.3)  
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 14), 
        axis.text.x = element_text(angle = 90, size = 13),axis.text.y = element_text(size = 13), axis.title.y = element_text(angle = 0,vjust=0.5))

f_dpt_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/f_dpt.tiff", width = 2100, height = 1500, units = "px", res = 300)
f_dpt_plot
dev.off()

## Relative humidity
f_rh  <- as.data.frame(t(f_prism_RH)) ## new data frame for precipitation plot
f_rh <- rownames_to_column(f_rh, var="month")
f_rh <- melt(f_rh, id='month')
f_rh$month <-  c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                 'September', 'October', 'November', 'December')
f_rh$month <- factor(f_rh$month, levels = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                                            'September', 'October', 'November', 'December'))

f_rh_plot <- ggplot(f_rh, aes(x = month, y = value, colour = variable, group = variable)) + geom_line( size=0.6) + theme_bw() +
  scale_colour_manual( values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'), name = 'site',guide = guide_legend(reverse=TRUE))  +
  labs(y='%')+ #+ stat_summary(geom="ribbon", fun.ymin=min, fun.ymax=max, aes(fill=type), alpha=0.3)  
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 14), 
        axis.text.x = element_text(angle = 90, size = 13),axis.text.y = element_text(size = 13), axis.title.y = element_text(angle = 0,vjust=0.5))

f_rh_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/f_rh.tiff", width = 2100, height = 1500, units = "px", res = 300)
f_rh_plot
dev.off()


f_climate_multi <- ggarrange(f_dpt_plot, f_temp_plot, f_precip_plot, f_rh_plot, ncol = 2, nrow = 2, labels = c('A','B','C','D'), common.legend = TRUE)
f_climate_multi

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/f__climate_multi.tiff", width = 2100, height = 1500, units = "px", res = 300)
f_climate_multi
dev.off()

######################################## Danthonia

## precipitation
d_precip <- as.data.frame(t(d_prism_precip)) ## new data frame for precipitation plot
d_precip <- rownames_to_column(d_precip, var="month")
d_precip <- melt(d_precip, id='month')
d_precip$month <- factor(d_precip$month, levels = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                                                    'September', 'October', 'November', 'December'))
## Danthonia precipitation plot for thesis
d_precip_plot <- ggplot(d_precip,aes(x=month,y=value,colour=variable,group=variable)) + geom_line( size=0.6) + theme_bw() +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'), name = 'site', guide = guide_legend(reverse=TRUE)) +
  labs(y='precipitation (cm)') + theme(legend.text=element_text(size=13), legend.title = element_text(size = 14), 
                                       axis.text.x = element_text(angle = 90, size = 13),axis.text.y = element_text(size = 13))

d_precip_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/d_precip.tiff", width = 2100, height = 1500, units = "px", res = 300)
d_precip_plot
dev.off()

## temperature
d_temp  <- as.data.frame(t(d_prism_temp)) ## new data frame for precipitation plot
d_temp <- rownames_to_column(d_temp, var="month")
d_temp <- melt(d_temp, id='month')
d_temp$month <-  c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                   'September', 'October', 'November', 'December')
d_temp$month <- factor(d_temp$month, levels = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                                                'September', 'October', 'November', 'December'))

d_temp_plot <- ggplot(d_temp, aes(x = month, y = value, colour = variable, group = variable)) + geom_line( size=0.6) + theme_bw() +
  scale_colour_manual( values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'), name = 'site',guide = guide_legend(reverse=TRUE))  +
  labs(y=(expression(~degree~C)))+ #+ stat_summary(geom="ribbon", fun.ymin=min, fun.ymax=max, aes(fill=type), alpha=0.3)  
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 14), 
        axis.text.x = element_text(angle = 90, size = 13),axis.text.y = element_text(size = 13), axis.title.y = element_text(angle = 0,vjust=0.5))

d_temp_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/d_temp.tiff", width = 2100, height = 1500, units = "px", res = 300)
d_mean_temp_plot
dev.off()

## dew-point temperature
d_dpt  <- as.data.frame(t(d_prism_dpt)) ## new data frame for precipitation plot
d_dpt <- rownames_to_column(d_dpt, var="month")
d_dpt <- melt(d_dpt, id='month')
d_dpt$month <-  c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                  'September', 'October', 'November', 'December')
d_dpt$month <- factor(d_dpt$month, levels = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                                              'September', 'October', 'November', 'December'))

d_dpt_plot <- ggplot(d_dpt, aes(x = month, y = value, colour = variable, group = variable)) + geom_line( size=0.6) + theme_bw() +
  scale_colour_manual( values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'), name = 'site',guide = guide_legend(reverse=TRUE))  +
  labs(y=(expression(~degree~C)))+ #+ stat_summary(geom="ribbon", fun.ymin=min, fun.ymax=max, aes(fill=type), alpha=0.3)  
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 14), 
        axis.text.x = element_text(angle = 90, size = 13),axis.text.y = element_text(size = 13), axis.title.y = element_text(angle = 0,vjust=0.5))

d_dpt_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/d_dpt.tiff", width = 2100, height = 1500, units = "px", res = 300)
d_dpt_plot
dev.off()

## Relative humidity
d_rh  <- as.data.frame(t(d_prism_RH)) ## new data frame for precipitation plot
d_rh <- rownames_to_column(d_rh, var="month")
d_rh <- melt(d_rh, id='month')
d_rh$month <-  c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                 'September', 'October', 'November', 'December')
d_rh$month <- factor(d_rh$month, levels = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                                            'September', 'October', 'November', 'December'))

d_rh_plot <- ggplot(d_rh, aes(x = month, y = value, colour = variable, group = variable)) + geom_line( size=0.6) + theme_bw() +
  scale_colour_manual( values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#3399FF', '#006699'), name = 'site',guide = guide_legend(reverse=TRUE))  +
  labs(y='%')+ #+ stat_summary(geom="ribbon", fun.ymin=min, fun.ymax=max, aes(fill=type), alpha=0.3)  
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 14), 
        axis.text.x = element_text(angle = 90, size = 13),axis.text.y = element_text(size = 13), axis.title.y = element_text(angle = 0,vjust=0.5))

d_rh_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/d_rh.tiff", width = 2100, height = 1500, units = "px", res = 300)
d_rh_plot
dev.off()

###########################################################################################################
################################################ examine meta data

environmental <- read.csv('meta_data.csv')
environmental <- environmental[1:155,] # remove control data
rownames(environmental) <- environmental[,1] # rownames as sample id

# subset by host
f_environ <- environmental[environmental$host == 'F. roemeri',]
d_environ <- environmental[environmental$host == 'D. californica',]

# subset by 'class' of metadata -
# I'm actually going to use data modified from PRISM as the climatic data here -
f_environ$site <- gsub('_',' ',f_environ$site) ## before I merge, I need to make the site names match...
f_environ$site <- gsub('f', 'F', f_environ$site)
f_environ$site <- gsub('Horse Rock', 'Horse Rock Ridge', f_environ$site)
f_environ$site <- gsub('Upper Table', 'U. Table Rock', f_environ$site)

f_environ_climate <- merge.data.frame(f_environ,f_bioclim, by.x = c('site'), by.y = c('SITE'))
rownames(f_environ_climate) <- rownames(f_environ)
f_environ_climate <- f_environ_climate[,c(19,57:68)]
sapply(f_environ_climate,class) ## all numeric

f_environ_geography <- f_environ[,c(6,7,8,9,17,18,19)]
sapply(f_environ_geography,class) ## 

f_environ_edaphic <- f_environ[,20:32]
sapply(f_environ_edaphic,class) ## all numeric except soil texture

f_environ_bio <- f_environ[,c(11,12,13,15,16)]
sapply(f_environ_bio,class) ## numeric or integer

## create two different geographic data frames one which contains all unique points, and one that defines 
## lat/longitude at the site level
f_geo_data <- f_environ[,c(17,18,53,54,55,56)]
f_geo_data_u <- data.frame(f_environ$site, f_environ$samplename,f_environ$unique_longitude, f_environ$unique_latitude, f_environ$utm.E.unique, f_environ$utm.n_unique)
f_geo_data_s <- data.frame(f_environ$site, f_environ$samplename, f_environ$longitude, f_environ$latitude)
colnames(f_geo_data_u) <- c('site', 'sample name', 'longitude', 'latitude','UTM E', 'UTM N')
colnames(f_geo_data_s) <- c('site', 'sample name', 'longitude', 'latitude')
unique(f_geo_data_u)
unique(f_geo_data_s)

## danthonia
d_environ$site <- gsub('_',' ',d_environ$site) ## before I merge, I need to make the site names match...
d_environ$site <- gsub('f', 'F', d_environ$site)
d_environ$site <- gsub('Horse Rock', 'Horse Rock Ridge', d_environ$site)
d_environ$site <- gsub('Lower Table', 'L. Table Rock', d_environ$site)

d_environ_climate <- merge.data.frame(d_environ,d_bioclim, by.x = c('site'), by.y = c('SITE'))
rownames(d_environ_climate) <- rownames(d_environ)
d_environ_climate <- d_environ_climate[,c(19,57:68)]
sapply(d_environ_climate,class) ## all numeric

d_environ_geography <- d_environ[,c(6,7,8,9,17,18,19)]
sapply(d_environ_geography,class) ## 

d_environ_edaphic <- d_environ[,20:32]
sapply(d_environ_edaphic,class) ## all numeric except soil texture

d_environ_bio <- d_environ[,c(11,12,13,15,16)]
sapply(d_environ_bio,class) ## numeric or integer


d_geo_data <- d_environ[,c(53,54,55,56)]
d_geo_data_u <- data.frame(d_environ$site, d_environ$samplename, d_geo_data$unique_longitude, d_geo_data$unique_latitude, d_geo_data$utm.E.unique, d_geo_data$utm.n_unique)
d_geo_data_s <- data.frame(d_environ$site, d_environ$samplename, d_environ$longitude, d_environ$latitude)
colnames(d_geo_data_u) <- c('site', 'sample name', 'longitude', 'latitude','UTM E', 'UTM N')
colnames(d_geo_data_s) <- c('site', 'sample name', 'longitude', 'latitude')
unique(d_geo_data_u)
unique(d_geo_data_s)



## right now I'm looking to assign a unique long/lat or UTM coordinate in order to calculate the distance 
## between each point.  The following function seems to be effective in this, but I have so far only used it to find 
## distances among sites.  With this data, I plan to use a mantel (or partial mantel) to examine the effect of distance on 
## community dissimilarity.  Not sure how this will fit into my partitioning of variance, but it's a step forward atleast...

##########################################################################################################
##########################################################################################################
###################################################### Mantel Test for spatial correlation

## Distance matricies 


# matrix using lat/long - haversine distances were used to convert to 2d cartesian coordinates (meters)
f_geo_dist <- distm(f_geo_data_u[,3:4], fun = distHaversine) 
d_geo_dist <- distm(d_geo_data_u[,3:4], fun = distHaversine)

# matrix using UTM -these coordinates are already in meters, so a euclidian distance matrix will suffice
f_geo_UTM_d <- dist(f_geo_data_u[,5:6])
d_geo_UTM_d <- dist(d_geo_data_u[,5:6])

## matrix created using site definitions of lat/long
f_geo_dist_site <- distm(f_geo_data_s[,3:4], fun = distHaversine)
d_geo_dist_site <- distm(d_geo_data_s[,3:4], fun = distHaversine)

grass_geo_data <- rbind(f_geo_data_u, d_geo_data_u)
grass_geo_dist <- dist(grass_geo_data[,5:6])
## we have a matrix of community disimilarity (f_dist_mat/d_dist_mat), one of spatial distance (f_geo_dist/d_geo_dist)
## these can be used in a mantel test to examine spatial correlation

## mantel test

f_dist_matrix <- as.matrix(f_mat_dist)
f_mantel <- mantel(f_mat_dist, f_geo_dist, permutations = 9999)
f_mantelR <- mantel.rtest(f_mat_dist, f_geo_UTM_d, nrepet = 9999)

#Mantel statistic r: 0.1597
# Significance: 1e-04

# Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0359 0.0482 0.0593 0.0730   
# Permutation: free
# Number of permutations: 999

## mantel correlogram - visualize community dissimilarity with distance
f_mantel_corr <-mantel.correlog(f_dist_matrix,f_geo_dist)
summary(f_mantel_corr)
plot(f_mantel_corr)

host <- rep('F. roemeri',13)
f_mantel_mat <- as.matrix(f_mantel_corr$mantel.res)
f_mantel_mat <- as.data.frame(cbind(f_mantel_mat[,c(1,3,5)],host))
colnames(f_mantel_mat) <- c('distance', 'Mantel_correlation', 'p_value', 'host')
rownames(f_mantel_mat) <- c('f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8','f9', 'f10','f11','f12','f13')
# from the figure, we can see that there is a steep drop in community similarity until communities are more than
# ~75 km.  there is another steep drop in similarity around 180km or so.  This may correspond to the Kalamath mountain divide - 
# the for example, French Flat is about 213,5 km from Hazel dell, and Upper table is about 172.9km from Hazel Dell.  Hazel dell is 
#the closest site accross the Kla mtn divide with the south

d_dist_matrix <- as.matrix(d_mat_dist)
d_mantel <- mantel(d_mat_dist, d_geo_UTM_d, permutations = 99999)
d_mantel
## Mantel statistic based on Pearson's product-moment correlation 

## Call:
## mantel(xdis = d_mat_dist, ydis = d_geo_dist) 

## Mantel statistic r: 0.4515
## Significance: 1e-05 

## Upper quantiles of permutations (null model):
## 90%    95%  97.5%    99% 
## 0.0472 0.0625 0.0757 0.0925
## Permutation: free
## Number of permutations: 9999

## correlogram
d_mantel_corr1 <-mantel.correlog(d_dist_matrix,d_geo_dist, n.class = 10)
plot(d_mantel_corr1)
host <- rep('D. californica',10)
d_mantel_mat <- as.matrix(d_mantel_corr$mantel.res)
d_mantel_mat <- as.data.frame(cbind(d_mantel_mat[,c(1,3,5)], host), na.rm =T)
colnames(d_mantel_mat) <- c('distance', 'Mantel_correlation', 'p_value', 'host')
rownames(d_mantel_mat) <- c('d1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8','d9', 'd10')

## combine data into one data frame, force data into numeric
significant <- c(1,1,1,1,1,1,1,0,1,0)
mantel_mat <- rbind(d_mantel_mat[1:4,],f_mantel_mat[1:6,] )
mantel_mat <- cbind(mantel_mat, significant)
mantel_mat$distance <- as.numeric(as.character(mantel_mat$distance))
mantel_mat$Mantel_correlation <- as.numeric(as.character(mantel_mat$Mantel_correlation))

mantel_corr <- ggplot(mantel_mat, aes(x = distance, y = Mantel_correlation, group = host, color = host))  + 
  geom_line() + geom_point(aes(distance, Mantel_correlation, shape=factor(significant)), size = 3) +
  scale_shape_manual(values =c(1,19), guide='none') + theme(legend.position = c(0.95,0.96), legend.box = NULL, legend.background = element_rect(fill=NULL)) + 
  theme(legend.text = element_text(face = "italic", size = 10), text = element_text(size = 10)) 
mantel_corr <- mantel_corr + geom_hline(yintercept=0, linetype = 'dashed', color = 'black') + 
  ylab('Mantel r') + xlab('Distance (km)')+  theme_bw() + scale_x_continuous( limits = c(10000, 300000), breaks = c(50000,100000,150000,200000,250000,300000), labels = c('50', '100', '150', '200', '250', '300'))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=11,face="bold"),legend.position = c(.7,.8),text = element_text(size=10)) 

mantel_corr 


tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/git_repo/grass-endophyte-community/mantel.tiff", width = 1800, height = 1500, units = "px", res = 300)
plot(mantel_corr)
dev.off()

## from the data, there does seem to be some significant correlation between spatial distance and species matrix.
## both mantel tests are significant (festuca: p=0.001, r=0.1725; Danthonia: p=0.001, r=0.4459).
## we see positive spatial auto-correlation at close distances, with negative spatial auto-correlation at further distances.  
# Remember that our comparisons are using dissimilarity matricies
# I think this means that we are seeing more similar species at closer distance, and more dissimilar communities at further distances.

## I plan to use PCNM to refine the analysis, give more information to the system, and provide means to test the 
## explained variance in a variance partitioning model.

##########################################################################################################
##########################################################################################################
############################################ Trend surface analysis

## We will use trend surface analysis to examine any large-scale spatial patterns in out microbial data.  
# trend surface analysis is often used as a preliminary analysis before using pcnm/db-MEM analyses to tease out more
# fine sale structuring, but Dan has pointed out that given our sampling scheme, it may be more appropriate to 
# stop after the trend surface.  This anaysis will still yeild variables that we can use in out variation partitioning model

## center our data - not sure if it matters that we are using utm data and not lat lon data...
f_coord.c=scale(f_coord, scale=FALSE) #centring
plot(f_coord.c)
X=f_coord.c[,1]# give the x coordinate of each sample point
Y=f_coord.c[,2]# give the y coordinate of each sample point


# Play around by plot some first, second and third-degree functions of x and y, just to see how they graph
par(mfrow=c(3,3))

sr.value(f_coord, (X), clegend = NULL)
sr.value(f_coord, (Y), clegend = NULL)
sr.value(f_coord, (X^2), clegend = NULL)
sr.value(f_coord, (Y^2), clegend = NULL)
sr.value(f_coord, (X*Y))
sr.value(f_coord, (X^2*Y), clegend = NULL)
sr.value(f_coord, (Y^2*X), clegend = NULL)
sr.value(f_coord, (X^3), clegend = NULL)
sr.value(f_coord, (Y^3), clegend = NULL)
par(mfrow = c(1,1))

# Now let’s run the trend surface analysis of the actual data. 
# Firstwe will use a function poly() to construct third-degree orthogonal polynomials of the geographic coordinates. 
# Then we will run a RDA with the polynomials as explanatory variables

f_coord_poly_ortho=poly(as.matrix(f_coord.c), degree=3)  #raw=F default
colnames(f_coord_poly_ortho)=c("X","X2","X3","Y","XY","X2Y","Y2","XY2","Y3")

#RDA using all 9 polynomial terms
f_trend_rda_ortho <- rda(f_mat_hel~., data=as.data.frame(f_coord_poly_ortho))
#Computation of the adjusted R2.
f_R2adj_poly <- RsquareAdj(f_trend_rda_ortho)$adj.r.squared # 0.124299

# Forwardselection using Blanchet et al.(2008) double stopping criterion
f_coord_trend_fwd <- forward.sel(f_mat_hel, f_coord_poly_ortho, adjR2thresh=f_R2adj_poly)
#    variables order         R2      R2Cum   AdjR2Cum        F pvalue
# 1        Y2     7 0.06260545 0.06260545 0.05117381 5.476506  0.001
# 2        Y3     9 0.02867349 0.09127894 0.06884139 2.555848  0.001
# 3         X     1 0.02875166 0.12003060 0.08703175 2.613878  0.001
# 4        X2     2 0.02478729 0.14481789 0.10151753 2.289800  0.001
# 5        XY     5 0.02000668 0.16482458 0.11128769 1.868495  0.001
# 6       X2Y     6 0.01669935 0.18152393 0.11774657 1.571029  0.001

#New RDA using the 6 terms retained
f_coord_trend_rda <- rda(f_mat_hel~.,data=as.data.frame(f_coord_poly_ortho)[,f_coord_trend_fwd[,2]])
summary(f_coord_trend_rda)

RsquareAdj(rda(f_mat_hel~.,data=as.data.frame(f_coord_poly_ortho)[,f_coord_trend_fwd[,2]]))
# adjusted r2: 0.1177466

#overall test and test of canonical axes
anova.cca(f_coord_trend_rda, step=1000)

#           Df Variance      F Pr(>F)    
#  Model     6  0.12725 2.8462  0.001 ***
#  Residual 77  0.57377  

anova.cca(f_coord_trend_rda, step=1000, by="axis")

#           Df Variance      F Pr(>F)    
#  RDA1      1  0.05458 7.3246  0.001 ***
#  RDA2      1  0.03123 4.1908  0.001 ***
#  RDA3      1  0.02041 2.7392  0.001 ***
#  RDA4      1  0.00780 1.0471  0.853    
#  RDA5      1  0.00731 0.9804  0.879    
#  RDA6      1  0.00593 0.7952  0.992    
#  Residual 77  0.57377        

## looks like our first three RDAs are significant
# we'll visualize the trends

f_coord_trend_fit <- scores(f_coord_trend_rda, choices=c(1,2,3),display="lc", scaling=1)
par(mfrow=c(1,3))
sr.value(f_coord, f_coord_trend_fit[,1], clegend = NULL)
sr.value(f_coord, f_coord_trend_fit[,2])
sr.value(f_coord, f_coord_trend_fit[,3], clegend = NULL)
par(mfrow=c(1,1)) # return to default par

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/git_repo/grass-endophyte-community/festuca_trend_surface.tiff", width = 1800, height = 1500, units = "px", res = 300)
par(mfrow=c(1,3))
sr.value(f_coord, f_coord_trend_fit[,1], clegend = NULL)
sr.value(f_coord, f_coord_trend_fit[,2])
sr.value(f_coord, f_coord_trend_fit[,3], clegend = NULL)
dev.off()

############################################## danthonia

d_coord.c=scale(d_coord, scale=FALSE) #centring
plot(d_coord.c)
X=d_coord.c[,1]# give the x coordinate of each sample point
Y=d_coord.c[,2]# give the y coordinate of each sample point

# Now let’s run the trend surface analysis of the actual data. 
# Firstwe will use a function poly() to construct third-degree orthogonal polynomials of the geographic coordinates. 
# Then we will run a RDA with the polynomials as explanatory variables

d_coord_poly_ortho=poly(as.matrix(d_coord.c), degree=3)  #raw=F default
colnames(d_coord_poly_ortho)=c("X","X2","X3","Y","XY","X2Y","Y2","XY2","Y3")

#RDA using all 9 polynomial terms
d_trend_rda_ortho <- rda(d_mat_hel~., data=as.data.frame(d_coord_poly_ortho))
#Computation of the adjusted R2.
d_R2adj_poly <- RsquareAdj(d_trend_rda_ortho)$adj.r.squared # 0.179888

#Forwardselection using Blanchet et al.(2008) double stopping criterion
d_coord_trend_fwd <- forward.sel(d_mat_hel, d_coord_poly_ortho, adjR2thresh=d_R2adj_poly)
#    variables order         R2      R2Cum   AdjR2Cum        F pvalue
# 1         Y     4 0.08169125 0.08169125 0.06838243 6.138128  0.001
# 2        Y2     7 0.06266013 0.14435138 0.11918524 4.979718  0.001
# 3       X2Y     6 0.03961335 0.18396473 0.14742584 3.252426  0.001
# 4        X3     3 0.02701006 0.21097480 0.16315509 2.259325  0.001
# 5         X     1 0.02546297 0.23643776 0.17770220 2.167594  0.001


#New RDA using the 6 terms retained
d_coord_trend_rda <- rda(d_mat_hel~.,data=as.data.frame(d_coord_poly_ortho)[,d_coord_trend_fwd[,2]])
summary(d_coord_trend_rda)

RsquareAdj(rda(d_mat_hel~.,data=as.data.frame(d_coord_poly_ortho)[,d_coord_trend_fwd[,2]]))
# adjusted R2: 0.1777022

#overall test and test of canonical axes
anova.cca(d_coord_trend_rda, step=1000)

#           Df Variance      F Pr(>F)    
#  Model     5  0.16885 4.0255  0.001 ***
#  Residual 65  0.54528    

anova.cca(d_coord_trend_rda, step=1000, by="axis")

#           Df Variance      F Pr(>F)    
#  RDA1      1  0.07009 8.3556  0.001 ***
#  RDA2      1  0.03680 4.3868  0.001 ***
#  RDA3      1  0.02795 3.3312  0.001 ***
#  RDA4      1  0.01897 2.2608  0.001 ***
#  RDA5      1  0.01504 1.7928  0.001 ***
#  Residual 65  0.54528       

## looks like our first three RDAs are significant
# we'll visualize the trends

d_coord_trend_fit <- scores(d_coord_trend_rda, choices=c(1,2,3,4,5),display="lc", scaling=1)
par(mfrow=c(2,3))
sr.value(d_coord, d_coord_trend_fit[,1], clegend = NULL)
sr.value(d_coord, d_coord_trend_fit[,2])
sr.value(d_coord, d_coord_trend_fit[,3], clegend = NULL)
sr.value(d_coord, d_coord_trend_fit[,4], clegend = NULL)
sr.value(d_coord, d_coord_trend_fit[,5], clegend = NULL)
par(mfrow=c(1,1)) # return to default par

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/git_repo/grass-endophyte-community/danthonia_trend_surface.tiff", width = 1800, height = 1500, units = "px", res = 300)
par(mfrow=c(2,3))
sr.value(d_coord, d_coord_trend_fit[,1], clegend = NULL)
sr.value(d_coord, d_coord_trend_fit[,2])
sr.value(d_coord, d_coord_trend_fit[,3], clegend = NULL)
sr.value(d_coord, d_coord_trend_fit[,4], clegend = NULL)
sr.value(d_coord, d_coord_trend_fit[,5], clegend = NULL)
dev.off()

##########################################################################################################
##########################################################################################################
############################################ PCNM - principal coordinates of neighbourhood matricies

## PCNM is a method that produces orthogonal (linearly independent) spatial variable over a wide range of 
# spatial scales, and can work with many different sampling schemes, including irregular designs.  In irregular
# designs (such as mine)

## I've already created distance matricies for the site level: f_geo_dist_site and d_geo_dist_site
## now to assign truncation threshold - we'll follow the advice of Legende and Legende and look at a minimum spanning tree
## of euclidian distances
## numerical ecology pg 248

############################################################################################################
#################################################### source code for integral functions

##PCNM variables calculated from the package - now this in its self is complicated - I can't install the package directly
## since it is not supported for R 3.3.1, so I had to build a function from the source code on github

# load required functions - there isn't a working build of the packages PCNM and AEM that I can find for 
# r 3.5.2, so I have compiled the functions into a single .R file
source('PCNM_AEM_functions.R')

## load source code for sr.value graphing function
source('/Users/grahambailes/grass_endophyte_community/statistics/r_functions/sr.value.R')


########################################################################################################
########################################################################################################
###################################################### PCNM


################################################ Festuca - non-detrended data:
########################using non-detrended data seems to provide more of a midscale, site specific 
######################## spatial structure... In my anayses I used the detrended data, not this!!!!


## PCNM (or Moran's Eigenvector Mapping) basis workflow
# - construct distance matrix among sites
# - truncate matrix to include distnaces only among 'close' neighbours (the minimum distance which conects the furthest two sites)
# - compute a PCoA on the truncated distance matrix
# - retain eigenvectors that model positive spatial correlation
# - use these eigenvectors as spatial explanatory variables in downstream analyses

## Our minimum spanning distance is nearly 300km, which suggests that we will not be able to 
# investigate any spatial structuring that is less that that distance.  Essentially we can investigate
# the differences between regions.  Borlets try to add in an additional point to the 

# add points at portland (524648.9, 5038721), and Roseberg (470915.8, 4781253.2)
sub_points <- matrix(c(524648.9, 5038721,470915.8, 4781253.2), ncol = 2, byrow = T)
colnames(sub_points) <-c('UTM E', 'UTM N')
f_sub_geo <- rbind(f_geo_data_u[,5:6], sub_points)
f_sub_cord <- as.data.frame(scale(f_sub_geo, center = T, scale = T))

ptd <- dist(f_sub_cord) # construct a distance matrix using unique UTMs for each plant
span.ptd <- spantree(ptd)
plot.spantree(span.ptd, f_sub_geo)
dmin <- max(span.ptd$dist) ## 156457.3, down from 290359.3 m
## truncate our distance matrix using this distance.
# this should be the minimum distance capable of connecting all points.  Borcard and Legendre suggest  
ptd[ptd > dmin] <- 4*dmin
ptd.PCoA <- cmdscale(ptd, k=nrow(f_sub_geo)-1, eig = TRUE)
nb.ev <- length(which(ptd.PCoA$eig > 0.0000001)) ## keep only the positive eigenvalues
## here is where we can get rid of the two sub points we inserted
ptd.PCNM <- data.frame(ptd.PCoA$points[1:nrow(f_geo_data_u[5:6]), 1:nb.ev]) 

f_pcnm_rda <- rda(f_mat_hel, ptd.PCNM) ## rda of our festuca community by the PCNMs
f_pcnm_sigtest <- anova.cca(f_pcnm_rda) ## significance test
f_pcnm_sigtest
#Model: rda(X = f_mat_hel, Y = ptd.PCNM)
#Df               Variance      F Pr(>F)    
#    Model    37  0.34226 1.1987  0.001 ***
#    Residual 46  0.35499               

## okay, we need to use the ordiR2step function to determine which pcnms to use
f_pcnm_mod0 <- rda(f_mat_hel ~ 1, ptd.PCNM)
f_pcnm_mod1 <- rda(f_mat_hel ~ ., ptd.PCNM)
f_step_select <- ordiR2step(f_pcnm_mod0, f_pcnm_mod1, perm.max=999)

#Call: rda(formula = f_mat_hel ~ X4 + X2 + X3, data = ptd.PCNM)

#Inertia Proportion Rank
#Total         0.69725    1.00000     
#Constrained   0.07475    0.10720    3
#Unconstrained 0.62250    0.89280   80
#Inertia is variance 

#Eigenvalues for constrained axes:
#  RDA1    RDA2    RDA3 
#0.04211 0.02451 0.00812 

#Eigenvalues for unconstrained axes:
#  PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8 
#0.03699 0.02560 0.02185 0.02079 0.01609 0.01503 0.01403 0.01368 
#(Showed only 8 of all 80 unconstrained eigenvalues)

## significance?

f_step_select$anova

#                    R2.adj Df     AIC      F Pr(>F)   
#  + X4            0.030077  1 -30.881 3.5738  0.002 **
#  + X2            0.055596  1 -32.152 3.2158  0.002 **
#  + X3            0.073721  1 -32.823 2.5850  0.002 **
#  <All variables> 0.081360                            

## looks like they are all significant.  About 8% of system variance explained, which is cool

f_sigPCNM <- attributes(f_step_select$terms)$term.labels
f_PCNM <- ptd.PCNM[,f_sigPCNM] 

## lets look at the pcnm variables in conjunction with our species matrix
f_pcnm_rda <- rda(f_mat_hel ~., f_PCNM) ## rda of our festuca community by the PCNMs
f_pcnm_sigtest <- anova.cca(f_pcnm_rda)
f_pcnm_sigtest

# Model: rda(formula = f_mat_hel ~ X3 + X4 + X7 + X6 + X2 + X1, data = f_PCNM)
#           Df Variance      F Pr(>F)    
# Model     6  0.12678 2.8332  0.001 ***
# Residual 77  0.57424   

axes.test <- anova.cca(f_pcnm_rda, by='axis') ## takes a minute
axes.test
#M odel: rda(formula = f_mat_hel ~ X3 + X4 + X7 + X6 + X2 + X1, data = f_PCNM)
#          Df Variance      F Pr(>F)    
# RDA1      1  0.05431 7.2830  0.001 ***
# RDA2      1  0.03105 4.1640  0.001 ***
# RDA3      1  0.02039 2.7346  0.001 ***
# RDA4      1  0.00779 1.0450  0.891    
# RDA5      1  0.00731 0.9796  0.881    
# RDA6      1  0.00592 0.7934  0.991    
# Residual 77  0.57424  

## lets look at a RDA biplot to visualize the PCNM variables (blue arrows), 
options(repr.plot.height = 5)
plot(f_pcnm_rda, 
     display = c('sp','wa','bp'),
     scaling = 2,
     sub='Scaled to species',
     main='F.roemeri PCNMs'
)

options(repr.plot.height = 5)
plot(f_pcnm_rda, 
     display = c('sp','wa','bp'),
     scaling = 1,
     sub='Scaled to ssite',
     main='F.roemeri PCNMs'
)

##plot pcnms
sr.value(f_coord, f_PCNM[,1],clegend = NULL)
sr.value(f_coord, f_PCNM[,2],clegend = NULL)
sr.value(f_coord, f_PCNM[,3],clegend = NULL)
sr.value(f_coord, f_PCNM[,4],clegend = NULL)
sr.value(f_coord, f_PCNM[,5],clegend = NULL)

## run linear models to examine correlation among pcnm variables and 
# as a note - this is calling variables which are created further down on the script...
f_pcnm_test_vars <- cbind(f_host_model,f_env_model, f_space_model[,1])
f_pcnm_V3 <-lm(f_PCNM[,3] ~ ., data = f_pcnm_test_vars)
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          212656.52   72646.81   2.927  0.00451 ** 
#  plant_area               37.93      59.78   0.634  0.52770    
# stems                 -1011.41     852.45  -1.186  0.23913    
# density               -7962.85    1640.66  -4.853 6.32e-06 ***
# bulk_density         -98479.41   39994.45  -2.462  0.01607 *  
# percent_sand           1319.35     567.50   2.325  0.02275 *  
# WINTER.PRECIPITATION    -13.13     137.82  -0.095  0.92437    
# `f_space_model[, 1]`   -147.05      43.80  -3.357  0.00123 ** 

## Residual standard error: 65200 on 76 degrees of freedom
## Multiple R-squared:  0.4046,	Adjusted R-squared:  0.3498 
## F-statistic: 7.378 on 7 and 76 DF,  p-value: 9.514e-07

f_pcnm_V2 <-lm(f_PCNM[,2] ~ ., data = f_pcnm_test_vars)
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           306577.62   73358.98   4.179 7.74e-05 ***
# plant_area                48.47      60.36   0.803   0.4245    
# stems                   -518.71     860.80  -0.603   0.5486    
# density                -7215.42    1656.75  -4.355 4.10e-05 ***
# bulk_density         -220011.06   40386.52  -5.448 6.06e-07 ***
# percent_sand             898.56     573.06   1.568   0.1210    
# WINTER.PRECIPITATION     -34.07     139.17  -0.245   0.8073    
# `f_space_model[, 1]`     -99.65      44.23  -2.253   0.0272 *  

# Residual standard error: 65840 on 76 degrees of freedom
# Multiple R-squared:  0.5284,	Adjusted R-squared:  0.485 
# F-statistic: 12.16 on 7 and 76 DF,  p-value: 2.541e-10

f_pcnm_V1 <-lm(f_PCNM[,1] ~ ., data = f_pcnm_test_vars)

#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          212656.52   72646.81   2.927  0.00451 ** 
# plant_area               37.93      59.78   0.634  0.52770    
# stems                 -1011.41     852.45  -1.186  0.23913    
# density               -7962.85    1640.66  -4.853 6.32e-06 ***
# bulk_density         -98479.41   39994.45  -2.462  0.01607 *  
# percent_sand           1319.35     567.50   2.325  0.02275 *  
# WINTER.PRECIPITATION    -13.13     137.82  -0.095  0.92437    
# `f_space_model[, 1]`   -147.05      43.80  -3.357  0.00123 ** 

# Residual standard error: 65200 on 76 degrees of freedom
# Multiple R-squared:  0.4046,	Adjusted R-squared:  0.3498 
# F-statistic: 7.378 on 7 and 76 DF,  p-value: 9.514e-07

################################################ Danthonia- non-detrended data:
########################using non-detrended data seems to provide more of a midscale, site specific 
######################## spatial structure...


## PCNM (or Moran's Eigenvector Mapping) basis workflow
# - construct distance matrix among sites
# - truncate matrix to include distnaces only among 'close' neighbours (the minimum distance which conects the furthest two sites)
# - compute a PCoA on the truncated distance matrix
# - retain eigenvectors that model positive spatial correlation
# - use these eigenvectors as spatial explanatory variables in downstream analyses

## Our minimum spanning distance is nearly 300km, which suggests that we will not be able to 
# investigate any spatial structuring that is less that that distance.  Essentially we can investigate
# the differences between regions.  Borlets try to add in an additional point to the 

# add points at portland (524648.9, 5038721), and Roseberg (470915.8, 4781253.2), along with 
# olympia, since we don't have samples from there either 
d_sub_points <- matrix(c(524648.9, 5038721,470915.8, 4781253.2,522099.8, 5195176), ncol = 2, byrow = T)
colnames(d_sub_points) <-c('UTM E', 'UTM N')
d_sub_geo <- rbind(d_geo_data_u[,5:6], d_sub_points)
d_sub_cord <- as.data.frame(scale(d_sub_geo, center = T, scale = T))

ptd <- dist(d_sub_geo) # construct a distance matrix using unique UTMs for each plant
span.ptd <- spantree(ptd)
plot.spantree(span.ptd, d_sub_geo)
dmin <- max(span.ptd$dist) ## 156475.8, down from 290359.3 m
## truncate our distance matrix using this distance.
# this should be the minimum distance capable of connecting all points.  Borcard and Legendre suggest  
ptd[ptd > dmin] <- 4*dmin
ptd.PCoA <- cmdscale(ptd, k=nrow(d_sub_geo)-1, eig = TRUE) # 38 positive eigenvalues
nb.ev <- length(which(ptd.PCoA$eig > 0.0000001)) ## keep only the positive eigenvalues
## here is where we can get rid of the two sub points we inserted
ptd.PCNM <- data.frame(ptd.PCoA$points[1:nrow(d_geo_data_u[5:6]), 1:nb.ev]) 

d_pcnm_rda <- rda(d_mat_hel, ptd.PCNM) ## rda of our festuca community by the PCNMs
d_pcnm_sigtest <- anova.cca(d_pcnm_rda) ## significance test
d_pcnm_sigtest
#Model: rda(X = d_mat_hel, Y = ptd.PCNM)
#Df               Variance      F Pr(>F)    
#    Model    37  0.44425  1.3862  0.001 ***
#    Residual 46   0.26987               

## okay, we need to use the ordiR2step function to determine which pcnms to use
d_pcnm_mod0 <- rda(d_mat_hel ~ 1, ptd.PCNM)
d_pcnm_mod1 <- rda(d_mat_hel ~ ., ptd.PCNM)
d_step_select <- ordiR2step(d_pcnm_mod0, d_pcnm_mod1, perm.max=999)
d_step_select
#Call: rda(formula = d_mat_hel ~ X1 + X2 + X4 + X5, data = ptd.PCNM)

#               Inertia Proportion Rank
# Total          0.7141     1.0000     
# Constrained    0.1495     0.2094    4
# Unconstrained  0.5646     0.7906   66
# Inertia is variance 

#Eigenvalues for constrained axes:
#  RDA1    RDA2    RDA3    RDA4 
#0.06852 0.03636 0.02571 0.01893 

#Eigenvalues for unconstrained axes:
#  PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
#0.028361 0.019140 0.016910 0.016042 0.015123 0.014108 0.013829 0.013203 
#(Showing 8 of 66 unconstrained eigenvalues)

## significance?

d_step_select$anova

#                      R2.adj Df     AIC      F Pr(>F)   
#   + X1            0.075166  1 -27.483 6.6893  0.002 **
#   + X2            0.118683  1 -29.941 4.4070  0.002 **
#   + X4            0.145269  1 -31.168 3.1151  0.002 **
#   + X5            0.161461  1 -31.593 2.2937  0.002 **
#   <All variables> 0.173321                             

## looks like they are all significant.  About 17% of system variance explained, which is cool

d_sigPCNM <- attributes(d_step_select$terms)$term.labels
d_PCNM <- ptd.PCNM[,d_sigPCNM] 

## lets look at the pcnm variables in conjunction with our species matrix
d_pcnm_rda <- rda(d_mat_hel ~., d_PCNM) ## rda of our festuca community by the PCNMs
d_pcnm_sigtest <- anova.cca(d_pcnm_rda)
d_pcnm_sigtest

# Model: rda(formula = d_mat_hel ~ X1 + X2 + X4 + X5, data = d_PCNM)
#           Df Variance      F Pr(>F)    
# Model     4  0.14952 4.3696  0.001 ***
# Residual 66  0.56460           

axes.test <- anova.cca(d_pcnm_rda, by='axis') ## takes a minute
axes.test

# Model: rda(formula = d_mat_hel ~ X1 + X2 + X4 + X5, data = d_PCNM)
#           Df Variance      F Pr(>F)    
#  RDA1      1  0.06852 8.0097  0.001 ***
#  RDA2      1  0.03636 4.2499  0.001 ***
#  RDA3      1  0.02571 3.0059  0.001 ***
#  RDA4      1  0.01893 2.2129  0.001 ***
#  Residual 66  0.56460    

## lets look at a RDA biplot to visualize the PCNM variables (blue arrows), 
plot.new()
options(repr.plot.height = 5)
plot(d_pcnm_rda, 
     display = c('sp','wa','bp'),
     scaling = 2,
     sub='Scaled to species',
     main='F.roemeri PCNMs'
)

options(repr.plot.height = 5)
plot(d_pcnm_rda, 
     display = c('sp','wa','bp'),
     scaling = 1,
     sub='Scaled to ssite',
     main='F.roemeri PCNMs'
)

##plot pcnms
plot.new()
dev.off()
sr.value(d_sub_cord[1:71,], d_PCNM[,1],clegend = NULL, xlim = c(-3,3))
sr.value(d_sub_cord[1:71,], d_PCNM[,2],clegend = NULL)
sr.value(d_sub_cord[1:71,], d_PCNM[,3],clegend = NULL)
sr.value(d_sub_cord[1:71,], d_PCNM[,4],clegend = NULL)


## run linear models to examine correlation among pcnm variables and 
# as a note - this is calling variables which are created further down on the script...
d_pcnm_test_vars <- cbind(d_host_model,d_env_model, d_space_model[,1])
d_pcnm_V3 <-lm(d_PCNM[,3] ~ ., data = d_pcnm_test_vars)
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          212656.52   72646.81   2.927  0.00451 ** 
#  plant_area               37.93      59.78   0.634  0.52770    
# stems                 -1011.41     852.45  -1.186  0.23913    
# density               -7962.85    1640.66  -4.853 6.32e-06 ***
# bulk_density         -98479.41   39994.45  -2.462  0.01607 *  
# percent_sand           1319.35     567.50   2.325  0.02275 *  
# WINTER.PRECIPITATION    -13.13     137.82  -0.095  0.92437    
# `d_space_model[, 1]`   -147.05      43.80  -3.357  0.00123 ** 

## Residual standard error: 65200 on 76 degrees of freedom
## Multiple R-squared:  0.4046,	Adjusted R-squared:  0.3498 
## F-statistic: 7.378 on 7 and 76 DF,  p-value: 9.514e-07

d_pcnm_V2 <-lm(d_PCNM[,2] ~ ., data = d_pcnm_test_vars)
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           306577.62   73358.98   4.179 7.74e-05 ***
# plant_area                48.47      60.36   0.803   0.4245    
# stems                   -518.71     860.80  -0.603   0.5486    
# density                -7215.42    1656.75  -4.355 4.10e-05 ***
# bulk_density         -220011.06   40386.52  -5.448 6.06e-07 ***
# percent_sand             898.56     573.06   1.568   0.1210    
# WINTER.PRECIPITATION     -34.07     139.17  -0.245   0.8073    
# `d_space_model[, 1]`     -99.65      44.23  -2.253   0.0272 *  

# Residual standard error: 65840 on 76 degrees of freedom
# Multiple R-squared:  0.5284,	Adjusted R-squared:  0.485 
# F-statistic: 12.16 on 7 and 76 DF,  p-value: 2.541e-10

d_pcnm_V1 <-lm(d_PCNM[,1] ~ ., data = d_pcnm_test_vars)

#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          212656.52   72646.81   2.927  0.00451 ** 
# plant_area               37.93      59.78   0.634  0.52770    
# stems                 -1011.41     852.45  -1.186  0.23913    
# density               -7962.85    1640.66  -4.853 6.32e-06 ***
# bulk_density         -98479.41   39994.45  -2.462  0.01607 *  
# percent_sand           1319.35     567.50   2.325  0.02275 *  
# WINTER.PRECIPITATION    -13.13     137.82  -0.095  0.92437    
# `d_space_model[, 1]`   -147.05      43.80  -3.357  0.00123 ** 

# Residual standard error: 65200 on 76 degrees of freedom
# Multiple R-squared:  0.4046,	Adjusted R-squared:  0.3498 
# F-statistic: 7.378 on 7 and 76 DF,  p-value: 9.514e-07

################################################################################################
################################################ detrended data:

## this is the data which was used in the variation partitioning analyses -

## im having an issue with the 'AEM'package, which is relied upon for the PCNM function...
# I downloaded the source binary, and installed it, but 
# install.packages('/Users/grahambailes/Downloads/AEM_0.6.tar.gz', repos = NULL, type = 'source')

################################################ Festuca - detrended data:

## create data for PCNM!  I've already got my species matrix - f_mat_hel and d_mat_hel
## I've got my geopraphic data - f_geo_data_s.  I'll center this

## this is the matrix we'll use to detrend data
f_coord=as.data.frame(scale(f_geo_data_u[,5:6], center=TRUE, scale=FALSE))

# add points at portland (524648.9, 5038721), and Roseberg (470915.8, 4781253.2)
# we'll use this matrix to build our PCNM variables, 
sub_points <- matrix(c(524648.9, 5038721,470915.8, 4781253.2), ncol = 2, byrow = T)
colnames(sub_points) <-c('UTM E', 'UTM N')
f_sub_geo <- rbind(f_geo_data_u[,5:6], sub_points)
f_sub_cord <- as.data.frame(scale(f_sub_geo, center = T, scale = F))

## detrended data is essentially examining the residuals of a linear model using transformed species data against
## geographic data (unique UTM)

## test for linear trend between coordiantes and species matrix
anova(rda(f_mat_hel, f_coord))

##Model: rda(X = f_mat_hel, Y = f_geo_data_u[, 3:4])
#             Df Variance    F   Pr(>F)    
#  Model     2  0.02892 1.7428  0.001 ***
#  Residual 81  0.67210   

## We are seeing a significant 'trend' between our species matrix and geographic data - 
# this suggests that detrending the data will be benificial to building our PCNM variables

## detrend data - residuals from a linear model comparing geographic distance with species matrix
f_mat_detr <- resid(lm(as.matrix(f_mat_hel) ~., data = f_coord))

## create PCNM variables

## we'll test the impact of adding several data points to the distance matrix - Borcard and Legendre suggest
# that for irregular sampling schemes, the addition of several supplimentary points can help decrease the 
# truncated distance, and help compute more 'fine-scale' PCNM variables.  you apparently loose the 'orthogonality' of
# the PCNM variables, but they say that only adding a small fraction of points shouldn't change things much...

## add points at portland (524648.9, 5038721), and Roseberg (470915.8, 4781253.2)


f_dm <- dist(f_coord) ## distance matrix from geographic distances
f_sub_dist <- dist(f_sub_geo) # addition of two locations to decrease truncation distance

f_pcnm1 <- PCNM(f_dm) ## pcnm variable from coordinate distance matrix
f_pcnm2 <-f_pcnm1$vectors ## extract vectors from pcnm object.
plot(f_coord) ## plot of pcnm variables on coordinate space

plot.spantree(f_pcnm1$spanning, f_coord) ## span tree to visualize conectivity of sites
f_dmin <-f_pcnm1$thresh ## pcnm threshold = 290359.3
f_posEV <-length(f_pcnm1$values)
f_pcnm1$Moran_I ## see which pcnm varaibles are positive, there are 48 variables generated, here are the first 10

#            Moran p.value Positive
# 1   0.6846921822   0.001     TRUE
# 2   0.0007468302   0.116     TRUE
# 3  -0.0040728175   0.160     TRUE
# 4  -0.0171401359   0.502    FALSE
# 5  -0.0196811204   0.260    FALSE
# 6  -0.0208433153   0.010    FALSE
# 7  -0.0208433153   0.028    FALSE
# 8  -0.0208433153   0.019    FALSE
# 9  -0.0208433153   0.030    FALSE
# 10 -0.0208433153   0.021    FALSE



f_select=which(f_pcnm1$Moran_I$Positive==TRUE) ## select only positive pcnm variables - there are three
f_pcnm3=f_pcnm2[-c(85,86),f_select] # remove the rows we subed in -

## analysis with PCNM variables
f_out2=rda(f_mat_hel, f_pcnm3) ## rda with 
anova.cca(f_out2)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(X = f_mat_hel, Y = f_pcnm3)
#            Df Variance     F Pr(>F)    
#  Model     3  0.05791 2.401  0.001 ***
#   Residual 80  0.64311      
## forward selection of PCNM variables - this is probably unnecessary because I only have two
f_adj_R2 <-RsquareAdj(f_out2)$adj.r.squared # 0.04330077
f_out3 <- forward.sel(f_mat_hel, as.matrix(f_pcnm3)) #forward selection using our adjusted r2 as a threshold
f_out3 ## output
#   variables order         R2      R2Cum   AdjR2Cum        F pvalue
# 1        V2     2 0.03510437 0.03510437 0.02333735 2.983284  0.001
# 2        V3     3 0.02914696 0.06425132 0.04114642 2.523010  0.001
# 3        V1     1 0.01835042 0.08260174 0.04819931 1.600214  0.005

## not a ton of explained variance here -

f_pcnm4 <- f_pcnm3[,c(sort(f_out3[,2]))]
f_pcnm4 <- as.data.frame(f_pcnm4)

f_out4 <-rda(f_mat_hel~., data = f_pcnm4) ## re-do analysis with reduced set of variables - in this case, it is the same  as without forward selection
summary(f_out4)
anova.cca(f_out4, by = 'term') ## test significance of RDA - looks like the first two are significant

#Model: rda(formula = f_mat_hel ~ V1 + V2 + V3, data = f_pcnm4)
#            Df Variance      F Pr(>F)    
#  V1        1  0.01286 1.6002  0.006 ** 
#  V2        1  0.02461 3.0612  0.001 ***
#  V3        1  0.02043 2.5417  0.001 ***
#  Residual 80  0.64311     

## plots of RDA axes
f_out6 <- scores(f_out4, choices = c(1,2,3), display = 'lc', scaling = 1)
colnames(f_out6) <- c('F_PCNM1', 'F_PCNM2', 'F_PCNM#3')
par(mfrow = c(1,2))
s.value(f_coord, f_out6[,1],clegend = .2)
s.value(f_coord, f_out6[,2],clegend = .2)
s.value(f_coord, f_out6[,3], clegend = .2)
sr.value(f_coord, f_out6[,1], clegend = .2)
sr.value(f_coord, f_out6[,2], clegend = .2)
sr.value(f_coord, f_out6[,3], clegend = .2)


tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/f_pcnm1.tiff", width = 3200, height = 3200, units = "px", res = 800)
sr.value(f_coord, f_out6[,1], clegend = .2)
tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/f_pcnm2.tiff", width = 3200, height = 3200, units = "px", res = 800)
sr.value(f_coord, f_out6[,2], clegend = .2)
dev.off()

## with a map-

library('png')
eco_map <- readPNG('/Users/grahambailes/grass_endophyte_community/Thesis/figures_maps/ecoregions.png')

## interpretation of results:
## extract species scores from RDA axes of interest.  
## species with high and low scores deive positive and negative patterns.  
## find the species scores for 
f_rda.scores <- scores(f_out4, choices = c(1), display = 'species')
fix(f_rda.scores)

## plot the first rda axis and the abundances of several of the most common species
sort(colSums(f_mat_hel)) # sort for abundance of otus and select several of most abundant
par(mfrow=c(1,2))
s.value(f_coord, f_out6[,1])
s.value(f_coord, f_mat_hel[,'OTU12:5grass'])
s.value(f_coord, f_mat_hel[,'OTU17:5grass'])
s.value(f_coord, f_mat_hel[,'OTU18:5grass'])

## with map -
eco_map <- readPNG('/Users/grahambailes/grass_endophyte_community/Thesis/figures_maps/ecoregions.png')

################################################################################################
############################################### Danthonia

d_coord=as.data.frame(scale(d_geo_data_u[,5:6], center=TRUE, scale=FALSE))

## detrended data is essentially examining the residuals of a linear model using transformed species data against
## geographic data (unique UTM)

## test for linear trend between coordiantes and species matrix
anova(rda(d_mat_hel, d_coord))

##Model: rda(X = d_mat_hel, Y = d_geo_data_u[, 3:4])
#            Df Variance      F Pr(>F)    
#  Model     2  0.08257 4.4454  0.001 ***
#  Residual 68  0.63155    

## yes, but little variance explained by geopraphic distance.  The PCNM will explore more complex trends

## detrend data - residuals from a linear model comparing geographic distance with species matrix
d_mat_detr <- resid(lm(as.matrix(d_mat_hel) ~., data = d_coord))

## create PCNM variables

d_dm1 <- dist(d_coord) ## distance matrix from geographic distances
d_pcnm1 <- PCNM(d_dm1) ## pcnm variable from coordinate distance matrix
d_pcnm2 <-d_pcnm1$vectors ## extract vectors from pcnm object.
plot(d_coord) ## plot of pcnm variables on coordinate space

plot.spantree(d_pcnm1$spanning, d_coord) ## span tree to visualize conectivity of sites
d_dmin <-d_pcnm1$thresh ## pcnm threshold = 435011.9
d_posEV <-length(d_pcnm1$values)
d_pcnm1$Moran_I ## see which pcnm varaibles are positive, 31 in total
#         Moran p.value Positive
#1   0.378991939   0.001     TRUE
#2  -0.002002293   0.026     TRUE
#3  -0.012509726   0.244     TRUE
#4  -0.019724448   0.168    FALSE
#5  -0.020104740   0.003    FALSE
#6  -0.020104740   0.006    FALSE
#7  -0.020104740   0.047    FALSE
#8  -0.020104740   0.032    FALSE
#9  -0.020104740   0.032    FALSE
#10 -0.020104740   0.041    FALSE

d_select=which(d_pcnm1$Moran_I$Positive==TRUE) ## select only positive pcnm variables
d_pcnm3=d_pcnm2[,d_select]

## analysis with PCNM variables
d_out2=rda(d_mat_hel, d_pcnm3) ## rda with 
anova.cca(d_out2)

##Model: rda(X = d_mat_hel, Y = d_pcnm3)
##           Df Variance      F Pr(>F)    
#  Model     3  0.11129 4.1229  0.001 ***
#  Residual 67  0.60283                  

## forward selection of PCNM variables - this is probably unnecessary because I only have two
d_adj_R2 <-RsquareAdj(d_out2)$adj.r.squared
d_out3 <- forward.sel(d_mat_hel, as.matrix(d_pcnm3[,1:3]), alpha = .1) #forward selection using our adjusted r2 as a threshold
d_out3 ## output
#   variables order         R2      R2Cum   AdjR2Cum        F pvalue
# 1        V3     3 0.07747115 0.07747115 0.06410116 5.794409  0.001
# 2        V1     1 0.06323587 0.14070702 0.11543369 5.004160  0.001
# 3        V2     2 0.01513211 0.15583912 0.11804088 1.201017  0.075



d_pcnm4 <- as.data.frame(d_pcnm3)

d_out4 <-rda(d_mat_hel~., data = d_pcnm4) ## re-do analysis with reduced set of variables - in this case, it is the same  as without forward selection
summary(d_out4)
anova.cca(d_out4, by='term') ## test significance of RDA - looks like it is
##Permutation test for rda under reduced model
##Marginal tests for axes
##Permutation: free
##Number of permutations: 999

##Model: rda(formula = d_mat_hel ~ V1 + V2 + V3, data = d_pcnm4)
#            Df Variance      F Pr(>F)    
#  V1        1  0.04516 5.0190  0.001 ***
#  V2        1  0.01081 1.2010  0.074 .  
#  V3        1  0.05532 6.1488  0.001 ***
#  Residual 67  0.60283  778    


## plots of RDA axes
d_out6 <- scores(d_out4, choices = c(1:3), display = 'lc', scaling = 1)
colnames(d_out6) <- c('D_PCNM1','D_PCNM2', 'd_PCNM3')
par(mfrow = c(1,2))
sr.value(d_coord, d_out6[,1],clegend = .5, pch = 20)
sr.value(d_coord, d_out6[,2],clegend = .5, pch = 20)
sr.value(d_coord, d_out6[,3],clegend = .5, pch = 20)

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/d_pcnm1.tiff", width = 3200, height = 3200, units = "px", res = 800)
sr.value(d_coord, d_out6[,1], clegend = .2)
tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/d_pcnm2.tiff", width = 3200, height = 3200, units = "px", res = 800)
sr.value(d_coord, d_out6[,2], clegend = .2)
dev.off()

## interpretation of results:
## extract species scores from RDA axes of interest.  
## species with high and low scores deive positive and negative patterns.  
## find the species scores for 
d_rda.scores <- scores(d_out4, choices = c(1), display = 'species')

## Okay, so we have these spatial variables - lets see if we can understand what they mean. 
## Dan suggested runing linear models to test the relationship among the PCNM variables with environmental ones
d_pcnm_test_vars <- cbind(d_host_model,d_env_model, d_space_model[,1])
d_pcnm_V1 <-lm(d_out6[,1] ~ ., data = d_pcnm_test_vars)
d_pcnm_V2 <-lm(d_out6[,2] ~ ., data = d_pcnm_test_vars)

######################################################################################################
######################################################################################################
############################################################### Variable Reduction 

## I have measured over 50 predictor variables of microbial communities, including host traits, climate, edaphic, and space.
# unfortunately, there is high levels of correlation among many of these variables, which presents a problem for 
# many statistical analyses.  

## March 2019:  I have decided to use bioclimate predictors for my climate variables (via USGS) - hopefully this will
# reduce the overall number of variables, as well as potentially describe more variation as the predictors have
# been observed to influence communities (albeit plant communities, but fuck it!)

## workflow:  Examine correlation among predictor variables within larger groupings (ie climate, edaphic, biotic, spatial)
# and remove variables that are highly correlated.  I'm doing this systematicly so as to still include variables
# which seem meaningful to my project.

# bundle spatial variables
f_spatial <- data.frame(f_out6[,-3], f_environ$latitude)
colnames(f_spatial) <- c('PCNM1','PCNM2', 'latitude')
fr_env_all <- data.frame(f_environ_climate, f_environ_edaphic[,-1], f_environ_bio[,-1],f_spatial) ## all environmental variables and annual climate
write.csv(fr_env_all, file = 'festuca_meta_data.csv') # for use in other places

d_spatial <- data.frame(d_out6, d_environ$latitude)
colnames(d_spatial) <- c('PCNM1','PCNM2', 'latitude')
dc_env_all <- data.frame(d_environ_climate, d_environ_edaphic[,-1], d_environ_bio[,-1],d_spatial) ## all environmental variables and annual climate
write.csv(dc_env_all, file = 'danthonia_meta_data.csv') # for use in other places

# *****************************************************
## okay, I am starting to lean toward using the year-previous climate data instead of the 30-year normals...
# also, I'll be using the plant-level plant data (as opposed to the site-level data), 
#at least for the functional-diversity analyses.  I can address the heirarchy using mixed effect 
# modeling (using site as a random effect).  For the variance partitioning, I may need to use the site-level data
# if I can't find a way to model site as random...

f_meta_1415 <- read.csv('f_meta_1415.csv')
d_meta_1415 <- read.csv('d_meta_1415.csv')

# Lets look at the correlations overall -

# we'll remove id and edaphic data.  Also MAT, MAP, latitude, and the # stems - looks like latitude and temperature
# seasonality are strongly correlated...
f_env_1415_cor <- cor(f_meta_1415[,-c(1:3,5,10:20,21,23:26,30,33,34,36,38:40)], use = 'na.or.complete', method = 'spearman')
par(cex=0.7, mar = c(2,2,0,2))
corrplot.mixed(cor(f_meta_1415[,-c(1:3,5,10:23,30, 38:40)], use = 'na.or.complete', method = 'spearman'),
               lower.col = 'black', number.cex = 0.7, tl.pos = 'lt')
corrplot.mixed(f_env_1415_cor, lower.col = 'black', number.cex = 0.7, tl.pos = 'lt')
f_env_1415_cor[which(f_env_1415_cor >=.70)]
f_env_1415_cor[which(f_env_1415_cor<=(-.70))]
## I'll be keeping the plant data, rh coldest, temp coldest, temp seasonality, precip wettest and seasonality,
# elevation and PCNM1/PCNM3.
# This will remove correlations above 0.75 -If I need to I can remove more...

f_env_1415_eval <- f_meta_1415[,-c(1:3,5,10:20,21,23:26,30,33,34,36, 38:40)]


## danthonia
d_env_1415_cor <- cor(d_meta_1415[,-c(1:3,5,6,10:20,21,23,25,26,30,33,34,38:40)], use = 'na.or.complete', method = 'pearson')
corrplot.mixed(d_env_1415_cor, lower.col = 'black', number.cex = 0.7, tl.pos = 'lt')
d_env_1415_cor[which(d_env_1415_cor >=.70)]
d_env_1415_cor[which(d_env_1415_cor<=(-.70))]
## looks like mean RH,precip, dpt, and temp of warmest quarter need to go, as well as PCNM2
# This will remove correlations above 0.75 - 

d_env_1415_eval <- d_meta_1415[,-c(1:3,5,6,10:20,21,23,25,26,30,31,33,34,36,38:40)]


######################################### variable reduction 2020_01_25
## we've created different spatial variables (trend surfaces), and will try to use the ordi-step method for 
# variable selection.  Dan suggested that we also try to include host species into a 'global' var par model to
# back up our claim that host is important from the nmds analyses.

f_spatial_trend <- f_coord_trend_fit
colnames(f_spatial_trend) <- c('sp1', 'sp2', 'sp3')

d_spatial_trend <- d_coord_trend_fit 
colnames(d_spatial_trend) <- c('sp1', 'sp2', 'sp3','sp4', 'sp5')

f_env_1415_cor <- cor(data.frame(f_meta_1415[,-c(1:3,5,10:20,21,23:26,30,33,34:40)],f_spatial_trend), use = 'na.or.complete', method = 'spearman')
par(cex=0.7, mar = c(2,2,0,2))
corrplot.mixed(cor(f_meta_1415[,-c(1:3,5,10:23,30, 38:40)], use = 'na.or.complete', method = 'spearman'),
               lower.col = 'black', number.cex = 0.7, tl.pos = 'lt')
corrplot.mixed(f_env_1415_cor, lower.col = 'black', number.cex = 0.7, tl.pos = 'lt')
f_env_1415_cor[which(f_env_1415_cor >=.70)]
f_env_1415_cor[which(f_env_1415_cor<=(-.70))]
## I'll be keeping the plant data, rh coldest, temp coldest, temp seasonality, precip wettest and seasonality,
# elevation and PCNM1/PCNM3.
# This will remove correlations above 0.75 -If I need to I can remove more...

f_env_1415_eval <- data.frame(f_meta_1415[,-c(1:3,5,10:20,21,23:26,30,33,34:40)])

## for some reason, only the spatial variables are being chosen with ordiR2step - I may remove them from the 
# variable selection method, or try a different cutoff method...

f_bare_model <- rda(f_mat_hel ~1, f_env_1415_eval)
f_full_model <- rda(f_mat_hel ~., f_env_1415_eval)
f_mod_step <- ordiR2step(f_bare_model, f_full_model, perm.max = 1000) ## using 
f_mod_step$anova


## danthonia
d_env_1415_cor <- cor(data.frame(d_meta_1415[,-c(1:3,5,6,10:20,21,23,25,26,30,33,34:40)],d_spatial_trend), use = 'na.or.complete', method = 'pearson')
corrplot.mixed(d_env_1415_cor, lower.col = 'black', number.cex = 0.7, tl.pos = 'lt')
d_env_1415_cor[which(d_env_1415_cor >=.70)]
d_env_1415_cor[which(d_env_1415_cor<=(-.70))]
## looks like mean RH,precip, dpt, and temp of warmest quarter need to go, as well as PCNM2
# This will remove correlations above 0.75 - 

d_env_1415_eval <- d_meta_1415[,-c(1:3,5,6,10:20,21,23,25,26,30,31,33,34,36,38:40)]


d_env_1415_eval <- data.frame(d_spatial_trend[,1:5],d_meta_1415[,-c(1:3,5,10:20,21,23:26,30,33,34:40)])

d_bare_model <- rda(d_mat_hel ~1, d_env_1415_eval)
d_full_model <- rda(d_mat_hel ~., d_env_1415_eval)
d_mod_step <- ordiR2step(d_bare_model, d_full_model, perm.max = 1000) ## using 
d_mod_step$anova

## this is the result of variable selction without the trend surface variables included...
#                                    R2.adj Df     AIC      F Pr(>F)   
#  + mean.RH.coldest.quarter       0.040780  1 -31.360 4.5286  0.002 **
#  + precipitation.wettest.quarter 0.062532  1 -32.317 2.9027  0.002 **
#  + density                       0.080618  1 -32.997 2.5935  0.002 **
#  + temperature.seasonality       0.100981  1 -33.935 2.8120  0.002 **
#  + precipitation.seasonality     0.109987  1 -33.851 1.7994  0.002 **
#  + Elevation_m                   0.117737  1 -33.670 1.6852  0.002 **
#  + sum_damage                    0.120372  1 -33.019 1.2307  0.028 * 
#  <All variables>                 0.121253  
############################################################################################
############################################################################################
############################################ partitioning of Variance



############################################ Festuca - explained variance of total community 

## we are looking at model selection using stepwize AIC and/or VIF.  I'm not convinced that this is the most
# appropriate way to build by model, because it may haphazardly choose variables that I'm not interested in.
# 

#build a min and ful model, and test using ordistep
f_bare_model <- rda(f_mat_hel ~1, f_meta_1415[,-c(1:3,5,10:23,38:40)])
f_full_model <- rda(f_mat_hel ~., f_meta_1415[,-c(1:3,5,10:23,38:40)])
f_mod_step <- ordiR2step(f_bare_model, f_full_model, perm.max = 1000) ## using 
f_mod_step$anova

#                                     R2.adj Df     AIC      F Pr(>F)   
#  + PCNM1                         0.044595  1 -31.695 4.8741  0.002 **
#  + precipitation.seasonality     0.058823  1 -31.986 2.2396  0.002 **
#  + temperature.seasonality       0.075191  1 -32.503 2.4336  0.002 **
#  + density                       0.113122  1 -35.077 4.4215  0.002 **
#  + mean.dpt.coldest.quarter      0.118257  1 -34.635 1.4601  0.002 **
#  + sum_damage                    0.120941  1 -33.975 1.2382  0.014 * 
#  + precipitation.wettest.quarter 0.123630  1 -33.331 1.2363  0.020 * 
#  <All variables>                 0.126183                       



f_env_1415_model <-f_env_1415_eval[,c(7:11)] ## using correlation matrix
f_env_1415_model2 <-f_meta_1415[,c(28,29,31,34)] # using ordistep

f_host_1415_model <-f_env_1415_eval[,c(1,3:5)] ## correlation matrix
f_host_1415_model2 <-f_env_1415_eval[,c(4,5)] # ordistep

f_space_1415_model <- f_env_1415_eval[,c(12,13)] # correlation matrix
f_space_1415_model2 <- f_env_1415_eval[,12] # ordistep 


fr_var_par <- varpart(f_mat_dist, f_env_1415_model, f_space_1415_model, f_host_1415_model) #ordistep
fr_var_par2 <- varpart(f_mat_dist, f_env_1415_model2, f_space_1415_model2, f_host_1415_model2) # env = climate

plot(fr_var_par)
plot(fr_var_par2)

## The correlation model explains more varaince, but it doesn't seperate varaibles as well as the ordistep model
# I'm leaning toward the ordistep model because it 
## The output of these is a bit confusing... the output for X1, X2, and X2 are each 
# variable set tested alone, (in other words, their variance explained alone).  a, b, and c are 
# the unique variance explained by each variable set.

# Partition of squared Binary jaccard distance in dbRDA 

#Call: varpart(Y = f_mat_dist, X = f_env_1415_model2, f_space_1415_model2, f_host_1415_model2)

#No. of explanatory tables: 3 
#Total variation (SS): 28.372 
#No. of observations: 84 

#Explanatory tables:
# X1:  f_env_1415_model: temp seasonality, precip wettest quarter, precip seasonality, MAT
# X2:  f_space_1415_model: pcnm1
# X3:  f_host_1415_model :  damage, density

#Partition table:
#                       Df R.square Adj.R.square Testable
#[a+d+f+g] = X1         4  0.11189      0.06692     TRUE
#[b+d+e+g] = X2         1  0.06314      0.05172     TRUE
#[c+e+f+g] = X3         2  0.03873      0.01500     TRUE
#[a+b+d+e+f+g] = X1+X2  5  0.15229      0.09795     TRUE
#[a+c+d+e+f+g] = X1+X3  6  0.16472      0.09963     TRUE
#[b+c+d+e+f+g] = X2+X3  3  0.09676      0.06289     TRUE
#[a+b+c+d+e+f+g] = All  7  0.21326      0.14080     TRUE
#Individual fractions                                   
#[a] = X1 | X2+X3       4               0.07791     TRUE
#[b] = X2 | X1+X3       1               0.04117     TRUE
#[c] = X3 | X1+X2       2               0.04285     TRUE
#[d]                    0               0.00672    FALSE
#[e]                    0              -0.01014    FALSE
#[f]                    0              -0.03168    FALSE
#[g]                    0               0.01396    FALSE
#[h] = Residuals                        0.85920    FALSE
#Controlling 1 table X                                  
#[a+d] = X1 | X3        4               0.08464     TRUE
#[a+f] = X1 | X2        4               0.04623     TRUE
#[b+d] = X2 | X3        1               0.04789     TRUE
#[b+e] = X2 | X1        1               0.03103     TRUE
#[c+e] = X3 | X1        2               0.03271     TRUE
#[c+f] = X3 | X2        2               0.01117     TRUE

showvarparts(3)

## With such a low portion of unexplained variance, what are we missing?  Off the top of my head, I'd suggest that
## we did a pretty good job sampling environmental varaibles - 

##one point where we may have missed something was the weather.
## our climate variables were 30yr normals for the given sites, but it is thought that previous year weather data may be just as important,
## especially considering the anomolous weather we experienced that year.

##  spatial variables - I also think that my spatial analysis didn't adequately describe spatial variation patterns.  The PCNM analysis
## I used was only able to pick up large scale spatial patterns, given the conectivity (or lack thereof) of my sites.
## this may manifest in the form of dispersal, 




## stacked bar chart of varpar
f_varpar_stacked <- NULL
f_varpar_stacked$host <- rep('F. roemeri', 5)
f_varpar_stacked$variable <- c('climate', 'spatial structure', 'host effects', 'climate / space', 
                               'climate / space / host')
f_varpar_stacked$variance_explained <- c(0.07, 0.035,0.035,0.005,0.007)  ## ordistep
# I've adjusted values down a bit because I want the total variance explained to reflect the real value

f_varpar_stacked <- as.data.frame(f_varpar_stacked)


f_var_par <- ggplot(f_varpar_stacked, aes(x = host, y = variance_explained, fill = variable)) + coord_flip() +
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) +  labs(fill = NULL, x = NULL, y = 'variation explained') +
  scale_fill_brewer() + theme(legend.position = 'top', legend.text = element_text(size = 12)) + scale_y_continuous(breaks = (seq(0,0.20,0.02))) +
  theme(axis.text.y = element_blank())

f_var_par

############################################# significance testing of variation partitioning (festuca)
## this testing will be for conditional (unique) varaition explained by each filtering component.  

## total variation
fr_all_test <- capscale(f_mat_dist ~., data.frame(f_space_1415_model2, 
                                                  f_host_1415_model2, f_env_1415_model2), comm = f_mat_hel)
fr_summary_all_test <- summary(fr_all_test)
RsquareAdj(fr_all_test)
# 0.140821

anova.cca(fr_all_test, step = 200, perm.max = 1000)
#           Df SumOfSqs      F Pr(>F)    
#  Model     7   6.0508 2.9431  0.001 ***
#  Residual 76  22.3216                  


####################################### spatial variables
## full variation explained
RsquareAdj(capscale(f_mat_dist ~., data.frame(f_space_1415_model2) , comm = f_mat_hel))
# 0.05171844
anova.cca((capscale(f_mat_dist ~., data.frame(f_space_1415_model2) , comm = f_mat_hel)), step = 200, perm.max = 1000)
#           Df SumOfSqs      F Pr(>F)    
#  Model     1   1.7915 5.5267  0.001 ***
#  Residual 82  26.5809   

## conditional explained variation
fr_space_test <- capscale(f_mat_dist ~.+ Condition(f_env_1415_model2$temperature.seasonality) +
                            Condition(f_env_1415_model2$precipitation.wettest.quarter)+
                            Condition(f_env_1415_model2$precipitation.seasonality) + 
                            Condition(f_env_1415_model2$MAT) +Condition(f_host_1415_model2$sum_damage) + 
                            Condition(f_host_1415_model2$density), data.frame(f_space_1415_model2) , comm = f_mat_hel)
f_summary_space_test <- summary(fr_space_test)
fr_space_test_samples <- as.data.frame(f_summary_space_test$sites[,1:2])
fr_space_test_species <- as.data.frame(f_summary_space_test$species[,1:2])

RsquareAdj(fr_space_test)
# adj.r.squared 0.04116879

anova.cca(fr_space_test, step=200, perm.max=200)
## Model: capscale(formula = f_mat_dist ~ PCNM1, data = f_space_1415_model2, comm = f_mat_hel)
#            Df SumOfSqs      F Pr(>F)    
#  Model     1   1.3773 4.6895  0.001 ***
#  Residual 76  22.3216         

# space significance
anova(fr_space_test, by = 'margin', perm.max = 200)
##           Df SumOfSqs      F Pr(>F)    
#  PCNM1     1   1.3773 4.6895  0.001 ***
#  Residual 76  22.3216               


######################################## environmental variables
## full variation explained
RsquareAdj(capscale(f_mat_dist ~., f_env_1415_model2, comm = f_mat_hel))
#0.06691906
anova.cca((capscale(f_mat_dist ~., f_env_1415_model2, comm = f_mat_hel)),step = 200, perm.max = 1000)
#            Df SumOfSqs      F Pr(>F)    
#  Model     4   3.1745 2.4882  0.001 ***
#  Residual 79  25.1979  

## conditional explained variation
fr_env_test <- capscale(f_mat_dist ~. + Condition(f_space_1415_model2) + 
                          Condition(f_host_1415_model2$sum_damage) + Condition(f_host_1415_model2$density), f_env_1415_model2, comm = f_mat_hel)
f_summary_env_test <- summary(fr_env_test)
fr_env_test_samples <- as.data.frame(f_summary_env_test$sites[,1:2])
fr_env_test_species <- as.data.frame(f_summary_env_test$species[,1:2])
RsquareAdj(fr_env_test)
# adj.r.squared 0.07791192

anova(fr_env_test, step=200, perm.max=200)
## Model: capscale(formula = f_mat_dist ~ mean.temperature.coldest.quarter + temperature.seasonality + precipitation.wettest.quarter + precipitation.seasonality + Elevation_m, data = f_env_1415_model1, comm = f_mat_hel)
#           Df SumOfSqs      F Pr(>F)    
# Model     4   3.3055 2.8136  0.001 ***
# Residual 76  22.3216    

# test sig. on environmental variables
anova(fr_env_test, by = 'mar', perm.max = 500)

#                                   Df SumOfSqs      F Pr(>F)    
#  temperature.seasonality        1   0.9393 3.1980  0.001 ***
#  precipitation.wettest.quarter  1   0.3280 1.1167  0.179    
#  precipitation.seasonality      1   0.2375 0.8087  0.929    
#  MAT                            1   0.2850 0.9705  0.523    
#  Residual                      76  22.3216   



############################################# host effects
## full explained variation
RsquareAdj(capscale(f_mat_dist ~., f_host_1415_model2, comm = f_mat_hel))
# 0.01499724
anova.cca((capscale(f_mat_dist ~., f_host_1415_model2, comm = f_mat_hel)),step = 200, perm.max = 1000)
#           Df SumOfSqs      F Pr(>F)    
#  Model     2   1.0989 1.6319  0.001 ***
#  Residual 81  27.2735

## conditional variation explained
fr_host_test <- capscale(f_mat_dist ~.+ Condition(f_env_1415_model2$temperature.seasonality) +
                           Condition(f_env_1415_model2$precipitation.wettest.quarter)+
                           Condition(f_env_1415_model2$precipitation.seasonality) + 
                           Condition(f_env_1415_model2$MAT) +
                           Condition(f_space_1415_model2), f_host_1415_model2, comm = f_mat_hel)
f_summary_host_test <- summary(fr_host_test)
fr_host_test_samples <- as.data.frame(f_summary_host_test$sites[,1:2])
fr_host_test_species <- as.data.frame(f_summary_host_test$species[,1:2])
RsquareAdj(fr_host_test)
# adj.r.squared  0.04284974

anova(fr_host_test, step=200, perm.max=200)
# Model: capscale(formula = f_mat_dist ~ plant_area + seeds_produced + sum_damage + density, data = f_host_1415_model1, comm = f_mat_hel)
#           Df SumOfSqs      F Pr(>F)    
# Model     2   1.7299 2.945  0.001 ***
# Residual 76  22.3216                  

## significance of host effects
anova(fr_host_test, by = 'mar', perm.max = 500)

#                Df SumOfSqs      F Pr(>F)   
# sum_damage  1   0.3750 1.2768  0.047 *  
# density     1   1.3409 4.5655  0.001 ***
# Residual   76  22.3216                  



## festuca colors:  ('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699')
## danthonia colors:  ('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699')


## and create some figures  - 
## I'm wondering if it would make sense to only include significant predictors in these? - yes

fr_env_baseplot <-plot(fr_env_test)
mult <- attributes(fr_env_baseplot$biplot)$arrow.mul
fr_env_test_samples$site <- f_environ$site
fr_env_arrows <- as.data.frame(f_summary_env_test$biplot[,1:2])
fr_env_arrows$lables <- c('temperature seasonality', 'precip of wettest quarter','precipitation seasonality','MAT')
colnames(fr_env_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_env_arrows_sig <- fr_env_arrows[1,]
fr_env_test_plot <- ggplot(data = fr_env_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c(  'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  # theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## when multi-plot making...
  geom_segment(data = fr_env_arrows_sig,aes(x = 0, xend =mult*x_arrow, y = 0, yend =mult*y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_env_arrows_sig,aes(x=(mult+mult/10)*x_arrow, y =(mult+mult/10)*y_arrow, 
                                               label = lables), size = 3.5, color='black') +
  geom_point(data = fr_env_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

fr_env_test_plot

fr_space_baseplot <-plot(fr_space_test)
mult <- attributes(fr_space_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
fr_space_test_samples$site <- f_environ$site
colnames(fr_space_test_samples) <- c('CAP1', 'CAP2', 'site')
colnames(fr_space_test_species) <- c('CAP1', 'CAP2')
fr_space_arrows <- data.frame(1,  0, 'PCNM1')
colnames(fr_space_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_space_test_plot <- ggplot(data = fr_space_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c( 'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  #theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## when multi-plot making...
  geom_segment(data = fr_space_arrows,aes(x = 0, xend =mult*x_arrow, y = 0, yend =mult*y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_space_arrows,aes(x=(mult+mult/10)*x_arrow, y =(mult+mult/10)*y_arrow, 
                                             label = lables), size = 3.5, color='black') +
  geom_point(data = fr_space_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

fr_space_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_f_space.tiff", width = 3200, height = 2800, units = "px", res = 500)
fr_space_test_plot
dev.off()


fr_host_test_samples$site <- f_environ$site
fr_host_baseplot <-plot(fr_host_test)
mult <- attributes(fr_host_baseplot$biplot)$arrow.mul
fr_host_arrows <- as.data.frame(f_summary_host_test$biplot[,1:2])
fr_host_arrows$lables <- c('host damage', 'host density')
colnames(fr_host_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_host_arrows_sig <- fr_host_arrows
fr_host_test_plot <- ggplot(data = fr_host_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c( 'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 10), legend.title = element_blank())+
  geom_segment(data = fr_host_arrows_sig,aes(x = 0, xend =mult*x_arrow, y = 0, yend =mult*y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_host_arrows_sig,aes(x=(mult+mult/10)*x_arrow, y =(mult+mult/10)*y_arrow, 
                                                label = lables), size = 3.5, color='black') +
  geom_point(data = fr_host_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

fr_host_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_f_host.tiff", width = 3200, height = 2800, units = "px", res = 500)
fr_host_test_plot
dev.off()

## all together now - I'll temporarily remove the legend from two of the three plots so they can fit together nicely
f_biplot_multi <- ggarrange(fr_env_test_plot, fr_space_test_plot, fr_host_test_plot, nrow = 1, labels = c('A','B','C'), 
                            common.legend = T,legend = 'top')
# note I changed the legend size to 10 for the multi pannel figure...
f_biplot_multi

f_biplot_plus <-ggarrange(f_var_par, 
                          ggarrange(fr_env_test_plot, fr_space_test_plot, fr_host_test_plot, ncol = 3, labels = c('B','C','D'), 
                                    common.legend = T,legend = 'right'), nrow = 2, labels = 'A', legend = 'top')

f_biplot_plus

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/f_biplot_multi_pannel.tiff", width = 2000, height = 1000, units = "px", res = 300)
f_biplot_multi
dev.off()

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/f_biplot_plus_multi_pannel.tiff", width = 2500, height = 1700, units = "px", res = 300)
f_biplot_plus
dev.off()

###################################### Danthonia 


d_bare_env <- rda(d_mat_hel ~ 1, d_meta_1415[,-c(1:3,5,6,10:22,36,38:40)])
d_full_env <- rda(d_mat_hel ~ ., d_meta_1415[,-c(1:3,5,6,10:22,36,38:40)])
d_mod_env <- ordiR2step(d_bare_env, d_full_env, perm.max = 1000, direction = 'both')
d_mod_env$anova

## model for variance partitioning
#                                    R2.adj Df     AIC      F Pr(>F)   
#  + PCNM1                         0.083063  1 -28.091 7.3411  0.002 **
#  + precipitation.wettest.quarter 0.117677  1 -29.860 3.7069  0.002 **
#  + temperature.seasonality       0.142773  1 -30.961 2.9908  0.002 **
#  + mean.dpt.coldest.quarter      0.160128  1 -31.480 2.3844  0.002 **
#  <All variables>                 0.171434  

d_env_1415_model <- d_env_1415_eval[,c(6:9)] # correlation model
d_env_1415_model2 <- d_meta_1415[,c(24,28,29)] # OrdiR2step

## given that the ordiR2step model didn't include any host variables we will use our correlation matrix 
d_host_1415_model <- d_env_1415_eval[,c(1:4)] # all host variables
d_host_1415_model2 <- d_env_1415_eval[,c(3,4)]


d_space_1415_model <- d_env_1415_eval[,c(11,12)] # correlation model
d_space_1415_model2 <- d_env_1415_eval[,c(11)] # ordiR2step

dc_var_par <- varpart(d_mat_dist, d_host_1415_model, d_space_1415_model, d_env_1415_model) # correlation
dc_var_par2 <- varpart(d_mat_dist, d_env_1415_model2, d_space_1415_model2, d_host_1415_model) # ordistep suggested

plot(dc_var_par, bg = 4:7)
plot(dc_var_par2, bg = 4:7)

## ordistep model:
#Partition of squared Binary jaccard distance in dbRDA 

#Call: varpart(Y = d_mat_dist, X = d_host_1415_model, d_space_1415_model2, d_env_1415_model2)

#Explanatory tables:
#  X1:  d_env_1415_model2: dpt temp coldest quarter, temp coldest quarter, precipitation wettest quarter
# X2:  d_space_1415_model: pcnm 1
# X3:  d_host_1415_model2: plant area, seeds produced, damage, density

#No. of explanatory tables: 3 
#Total variation (SS): 24.355 
#No. of observations: 71 

#Partition table:
#                       Df R.square Adj.R.square Testable
#[a+d+f+g] = X1         3  0.17538      0.13845     TRUE
#[b+d+e+g] = X2         1  0.10893      0.09601     TRUE
#[c+e+f+g] = X3         4  0.12142      0.06817     TRUE
#[a+b+d+e+f+g] = X1+X2  4  0.23330      0.18684     TRUE
#[a+c+d+e+f+g] = X1+X3  7  0.26418      0.18242     TRUE
#[b+c+d+e+f+g] = X2+X3  5  0.18749      0.12499     TRUE
#[a+b+c+d+e+f+g] = All  8  0.28397      0.19158     TRUE
#Individual fractions                                   
#[a] = X1 | X2+X3       3               0.06659     TRUE
#[b] = X2 | X1+X3       1               0.00916     TRUE
#[c] = X3 | X1+X2       4               0.00474     TRUE
#[d]                    0               0.04766    FALSE
#[e]                    0               0.03922    FALSE
#[f]                    0               0.02424    FALSE
#[g]                    0              -0.00003    FALSE
#[h] = Residuals                        0.80842    FALSE
#Controlling 1 table X                                  
#[a+d] = X1 | X3        3               0.11425     TRUE
#[a+f] = X1 | X2        3               0.09083     TRUE
#[b+d] = X2 | X3        1               0.05682     TRUE
#[b+e] = X2 | X1        1               0.04838     TRUE
#[c+e] = X3 | X1        4               0.04396     TRUE
#[c+f] = X3 | X2        4               0.02898     TRUE

showvarparts(3)

## stacked bar chart of varpar
d_varpar_stacked <- NULL
d_varpar_stacked$host <- rep('D. californica', 6)
d_varpar_stacked$variable <- c('climate','spatial structure', 'host effects', 'climate / space','climate / host', 
                               'space / host')
d_varpar_stacked$variance_explained <- c(0.067, 0.009,0.005,0.047,0.024, 0.039)
d_varpar_stacked <- as.data.frame(d_varpar_stacked)


d_var_par <- ggplot(d_varpar_stacked, aes(x = host, y = variance_explained, fill = variable)) + coord_flip() +
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) +  labs(fill = NULL, x = NULL, y = 'variance explained') +
  scale_fill_brewer() + theme(legend.position = 'top', legend.text = element_text(size = 12)) + scale_y_continuous(breaks = (seq(0,0.22,0.02))) +
  theme(axis.text.y = element_blank()) 
d_var_par
#################################################################################################
############################################## RDAs for significance

## total variation
dc_all_test <- capscale(d_mat_dist ~., data.frame(d_space_1415_model2, 
                                                  d_host_1415_model, d_env_1415_model2), comm = d_mat_hel)
dc_summary_all_test <- summary(dc_all_test)
RsquareAdj(dc_all_test)
# 0.1915774

anova.cca(dc_all_test, step = 200, perm.max = 1000)
#           Df SumOfSqs      F Pr(>F)    
#  Model     8   6.9159 3.0735  0.001 ***
#  Residual 62  17.4386     

anova.cca(dc_all_test, step = 200, perm.max = 1000, by = 'mar')

## test of matricies - basic rundown - rda of model component, anova to test significance, plot of 
## rda triplot to show how individual components affect communities

############################################ spatial filters
## full explained variation
RsquareAdj(capscale(d_mat_dist ~., data.frame(d_space_1415_model2), comm = d_mat_hel))
# 0.09601174
anova.cca((capscale(d_mat_dist ~., data.frame(d_space_1415_model2), comm = d_mat_hel)), step = 200, perm.max = 1000)
#            Df SumOfSqs      F Pr(>F)    
#  Model     1   2.6528 8.4346  0.001 ***
#  Residual 69  21.7017  

## conditional variation explained
dc_space_test <- capscale(d_mat_dist ~.+ Condition(d_env_1415_model2$mean.dpt.coldest.quarter) +
                            Condition(d_env_1415_model2$temperature.seasonality)+
                            Condition(d_env_1415_model2$precipitation.wettest.quarter) +
                            Condition(d_host_1415_model$plant_area) +
                            Condition(d_host_1415_model$seeds_produced) +
                            Condition(d_host_1415_model$sum_damage) +
                            Condition(d_host_1415_model$density), data.frame(d_space_1415_model2), comm = d_mat_hel)
d_summary_space_test <- summary(dc_space_test)
dc_space_test_samples <- as.data.frame(d_summary_space_test$sites)[,1:2]
dc_space_test_species <- as.data.frame(d_summary_space_test$species[,1:2])
RsquareAdj(dc_space_test)
# 0.009159813 conditional


anova.cca(dc_space_test, step=200, perm.max=200)
# Model: capscale(formula = d_mat_dist ~ PCNM1 + PCNM2, data = d_space_1415_model, comm = d_mat_hel)
#            Df SumOfSqs      F Pr(>F)    
#  Model     1    0.482 1.7138  0.001 ***
#  Residual 62   17.439                         

anova.cca(dc_space_test, perm.max = 500, by = 'term')
##          Df SumOfSqs      F Pr(>F)    
#  PCNM1     1    0.482 1.7138  0.001 ***
#  Residual 62   17.439 7 

############################################ environmental filters
## full explained variation
RsquareAdj(capscale(d_mat_dist ~., d_env_1415_model2, comm = d_mat_hel))
# 0.138455
anova.cca(capscale(d_mat_dist ~., d_env_1415_model2, comm = d_mat_hel), step = 200, perm.max = 1000)
#            Df SumOfSqs      F Pr(>F)    
#  Model     3   4.2713 4.7498  0.001 ***
#  Residual 67  20.0833  

## conditional variation explained
dc_env_test <- capscale(d_mat_dist ~.+Condition(d_space_1415_model2) +
                          Condition(d_host_1415_model$plant_area) +
                          Condition(d_host_1415_model$seeds_produced) +
                          Condition(d_host_1415_model$sum_damage) +
                          Condition(d_host_1415_model$density), d_env_1415_model2, comm = d_mat_hel)
d_summary_env_test <- summary(dc_env_test)
dc_env_test_samples <- as.data.frame(d_summary_env_test$sites[,1:2])
dc_env_test_species <- as.data.frame(d_summary_env_test$species[,1:2])
RsquareAdj(dc_env_test)
# 0.06658561

anova.cca(dc_env_test, step=200, perm.max=1000)
# Model: capscale(formula = d_mat_dist ~ mean.dpt.coldest.quarter + mean.temperature.coldest.quarter + 
#precipitation.wettest.quarter + Elevation_m, data = d_env_1415_model, comm = d_mat_hel)

#          Df SumOfSqs      F Pr(>F)    
# Model     3   2.3496 2.7846  0.001 ***
# Residual 62  17.4386 

anova.cca(dc_env_test, step = 200, perm.max = 1000, by = 'mar')
##                               Df SumOfSqs      F Pr(>F)    
#  mean.dpt.coldest.quarter       1   0.4390 1.5609  0.006 **
#  temperature.seasonality        1   0.4430 1.5751  0.002 **
#  precipitation.wettest.quarter  1   0.4645 1.6514  0.002 **
#  Residual                      62  17.4386        


############################################### host effects
## full variation explained
RsquareAdj(capscale(d_mat_dist ~., d_host_1415_model, comm = d_mat_hel))
# 0.06816968
anova.cca((capscale(d_mat_dist ~., d_host_1415_model, comm = d_mat_hel)), step = 200, perm.max = 1000)
#            Df SumOfSqs      F Pr(>F)    
#  Model     4   2.9571 2.2802  0.001 ***
#  Residual 66  21.3975  

## conditional explained variation
dc_host_test <- capscale(d_mat_dist ~.+ Condition(d_env_1415_model2$mean.dpt.coldest.quarter) +
                           Condition(d_env_1415_model2$temperature.seasonality)+
                           Condition(d_env_1415_model2$precipitation.wettest.quarter) + 
                           Condition(d_space_1415_model2), d_host_1415_model, comm = d_mat_hel)
d_summary_host_test <- summary(dc_host_test)
dc_host_test_samples <- as.data.frame(d_summary_host_test$sites[,1:2])
dc_host_test_species <- as.data.frame(d_summary_host_test$species[,1:2])
RsquareAdj(dc_host_test)
# 0.047388 conditional

anova.cca(dc_host_test, step=200, perm.max=1000)

## Model: capscale(formula = d_mat_dist ~ plant_area + seeds_produced + sum_damage + density, data = d_host_1415_model, comm = d_mat_hel)
#            Df SumOfSqs      F Pr(>F)    
#  Model     4   1.2339 1.0967  0.045 *
#  Residual 62  17.4386 

anova.cca(dc_host_test, step = 200, perm.max = 500, by = 'mar')

#                Df SumOfSqs      F Pr(>F)
# plant_area      1   0.2318 0.8240  0.952   
# seeds_produced  1   0.2683 0.9541  0.635   
# sum_damage      1   0.2178 0.7742  0.989   
# density         1   0.4171 1.4828  0.002 **
# Residual       62  17.4386     

## festuca colors:  ('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699')
## danthonia colors:  ('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699')



dc_env_baseplot <-plot(dc_env_test)
mult <- attributes(dc_env_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
dc_env_test_samples$site <- d_environ$site
dc_env_arrows <- as.data.frame(d_summary_env_test$biplot[,1:2])
dc_env_arrows$lables <- c('dpt coldest quarter', 'temperature seasonality','precip wettest quarter')
colnames(dc_env_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_env_test_plot <- ggplot(data = dc_env_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = dc_env_arrows,aes(x = 0, xend =mult*x_arrow, y = 0, yend =mult*y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_env_arrows,aes(x=(mult)*x_arrow, y =(mult)*y_arrow, 
                                           label = lables), size = 3.5, color='black') +
  geom_point(data = dc_env_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)


dc_env_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_env.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_env_test_plot
dev.off()


dc_space_baseplot <-plot(dc_space_test)
mult <- attributes(dc_space_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
dc_space_test_samples$site <- d_environ$site
dc_space_arrows <- as.data.frame(d_summary_space_test$biplot[,1:2])
colnames(dc_space_test_samples) <- c('CAP1', 'CAP2', 'site')
colnames(dc_space_test_species) <- c('CAP1', 'CAP2')
dc_space_arrows <- data.frame(0.09993861,  0, 'PCNM1')
colnames(dc_space_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_space_test_plot <- ggplot(data = dc_space_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  #theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none')+ ## for multi-pannel plotting...
  geom_segment(data = dc_space_arrows,aes(x = 0, xend =mult*x_arrow, y = 0, yend =mult*y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_space_arrows,aes(x=mult*x_arrow, y =  mult*y_arrow, 
                                             label = lables), size = 3.5, color='black') +
  geom_point(data = dc_space_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

dc_space_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_space.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_space_test_plot
dev.off()

dc_host_baseplot <-plot(dc_host_test)
mult <- attributes(dc_host_baseplot$biplot)$arrow.mul
dc_host_test_samples$site <- d_environ$site
dc_host_arrows <- as.data.frame(d_summary_host_test$biplot[,1:2])
dc_host_arrows$lables <- c('plant size', 'reproductive output', 'damage', 'density')
colnames(dc_host_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_host_arrows_sig <- dc_host_arrows[4,]
dc_host_test_plot <- ggplot(data = dc_host_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 10), legend.title = element_blank())+
  geom_segment(data = dc_host_arrows_sig,aes(x = 0, xend =mult*x_arrow, y = 0, yend =mult*y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_host_arrows_sig,aes(x=mult*x_arrow, y =  mult*y_arrow, 
                                                label = lables), size = 3.5, color='black') +
  geom_point(data = dc_host_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

dc_host_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_host.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_host_test_plot
dev.off()


# all together now - I'll temporarily remove the legend from two of the three plots so they can fit together nicely
d_biplot_multi <- ggarrange(dc_env_test_plot, dc_space_test_plot, dc_host_test_plot, nrow = 1, labels = c('A','B','C'), 
                            common.legend = T,legend = 'top')
# note I changed the legend size to 10 for the multi pannel figure...
d_biplot_multi

## with varpar added
d_biplot_plus <-ggarrange(d_var_par, 
                          ggarrange(dc_env_test_plot, dc_space_test_plot, dc_host_test_plot, ncol = 3, labels = c('B','C','D'), 
                                    common.legend = T,legend = 'right'), nrow = 2, labels = 'A', legend = 'top')
d_biplot_plus

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/d_biplot_multi_pannel.tiff", width = 2000, height = 1000, units = "px", res = 300)
d_biplot_multi
dev.off()

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/d_biplot_plus_multi_pannel.tiff", width = 2500, height = 1700, units = "px", res = 300)
d_biplot_plus
dev.off()


##########################################################################################################
################################################# funguild  


## Okay, so we are going to import our funguild database and create a physoseq object from it.
# we need to assemble all of the components including the taxa table, otu table, and sample data.

# we'll start by looking at the ordinations of the trimmed-down funguild bioms

## all OTUs
guilds_all <- read.csv('grass_biom_vs_fun.guilds.txt', header = T, sep = '\t', row.names = 1)
guilds_assigned <- guilds_all[guilds_all$Notes !='Unassigned',]
guilds_hq <-guilds_assigned[guilds_assigned$Confidence.Ranking != 'Possible',] 
colnames(guilds_hq) <- gsub('X','',colnames(guilds_hq))

# taxa table
guilds_taxa <- guilds_hq$taxonomy # extract taxonomy 
guilds_taxa <- do.call(rbind,strsplit(as.character(guilds_taxa),';')) # split string into taxa classifications
guilds_taxa <- gsub('.*__','',guilds_taxa) # remove designation abreviations
colnames(guilds_taxa) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
rownames(guilds_taxa) <- rownames(guilds_hq)
guilds_taxa <- as.data.frame(guilds_taxa)

# meta data
f_meta <- cbind(f_environ_climate, f_spatial, f_environ_bio)
d_meta <- cbind(d_environ_climate, d_spatial, d_environ_bio)
guilds_meta <- rbind(f_meta, d_meta)  
guilds_meta <-cbind(guilds_meta, grass_biom@sam_data[1:155,c(5,13,14)])

# guilds
guilds_trophic <- guilds_hq$Trophic.Mod

# otu table
guilds_otu <-(guilds_hq[,1:155])
dim(grass_biom@sam_data)

# okay, lets assemble as a phyloseq object
guilds_biom <- phyloseq(otu_table(guilds_otu, taxa_are_rows = T), sample_data(guilds_meta), tax_table(guilds_taxa))

set.seed(1)
## we'll run the nmds on our presence/absence matrix, using jaccard distance
guild_mat_ord <- t(guilds_biom@otu_table)
guild_mat_ord[guild_mat_ord@.Data > 0] <- 1
guild_ord <- metaMDS(guild_mat_ord, distance = 'jaccard', try = 20, trymax = 100)
guild_ord # stress is 0.2100692 - not bad but not great...
stressplot(guild_ord)
plot(guild_ord, display = 'species')

## lets graph it - we can change the coloration and faceting to examine different aspects of the nmds
guild_host <- plot_ordination(guilds_biom, guild_ord, color = "host") + theme_bw() + stat_ellipse() +
  theme(legend.text = element_text(face = "italic")) #+ facet_wrap('site')
guild_host

## festuca
f_guild_biom <- subset_samples(guilds_biom, host == 'F. roemeri')
f_guild_biom@sam_data$site <- factor(f_guild_biom@sam_data$site, levels = c('French_flat', 'Roxy_Ann', 'Upper_Table','Hazel_Dell', 'Horse_Rock', 'Upper_Weir', 'Whidbey'))
f_guild_biom@sam_data$region <- factor(f_guild_biom@sam_data$region, levels = c('SOR', 'COR', 'WA'))
f_guild_mat <- t(f_guild_biom@otu_table)
f_guild_mat[f_guild_mat@.Data > 0] <- 1
f_guild_ord <- metaMDS(f_guild_mat, distance = 'jaccard', try = 20, trymax = 100)
f_guild_ord # stress = 0.206
stressplot(f_guild_ord)

f_guild_host <- plot_ordination(f_guild_biom, f_guild_ord, color = "site")
f_guild_host <- f_guild_host + theme_bw() + stat_ellipse() 
f_guild_host <- f_guild_host + theme(legend.text = element_text(face = "italic"))+
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 14), legend.position = c(0.89,0.9))
plot(f_guild_host)

## Dnthonia
d_guild_biom <- subset_samples(guilds_biom, host == 'D. californica')
d_guild_biom@sam_data$site <- factor(d_guild_biom@sam_data$site, levels = c('French_flat', 'Whetstone', 'Lower_Table','Hazel_Dell', 'Horse_Rock', 'Whidbey'))
d_guild_biom@sam_data$region <- factor(d_guild_biom@sam_data$region, levels = c('SOR', 'COR', 'WA'))
d_guild_mat <- t(d_guild_biom@otu_table)
d_guild_mat[d_guild_mat@.Data > 0] <- 1
d_guild_ord <- metaMDS(d_guild_mat, distance = 'jaccard', try = 20, trymax = 100)
d_guild_ord # stress = 0.173
stressplot(d_guild_ord)

d_guild_host <- plot_ordination(d_guild_biom, d_guild_ord, color = "site")
d_guild_host <- d_guild_host + theme_bw() + stat_ellipse() 
d_guild_host <- d_guild_host + theme(legend.text = element_text(face = "italic"))+
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 14), legend.position = c(0.89,0.9))
plot(d_guild_host)

## read in the trophic matricies I've worked on in my 'grass_funguild_2019_03_01.R' script
f_pathotroph_mat <- read.csv('f_pathotroph_mat.csv', row.names = 1)
f_saprotroph_mat <- read.csv('f_saprotroph_mat.csv', row.names = 1)
f_symbiotroph_mat <- read.csv('f_symbiotroph_mat.csv', row.names = 1)

# create distance matricies
f_path_dist <- vegdist(f_pathotroph_mat, method = 'jaccard', binar = T)
f_sap_dist <- vegdist(f_saprotroph_mat, method = 'jaccard', binar = T)
f_sym_dist <- vegdist(f_symbiotroph_mat, method = 'jaccard', binar = T)

## helenger 
f_path_hel <- decostand(f_pathotroph_mat, 'hellinger')
f_sap_hel <- decostand(f_saprotroph_mat, 'hellinger')
f_sym_hel <- decostand(f_symbiotroph_mat, 'hellinger')


# run variance partiotioning on the trophic distance matricies using the same models as the entire ensamble
f_path_varpart <- varpart(f_path_dist, f_env_1415_model, f_space_1415_model, f_host_1415_model)
f_path_varpart2 <- varpart(f_path_dist, f_env_1415_model2, f_space_1415_model2, f_host_1415_model2)

#No. of explanatory tables: 3 
#Total variation (SS): 26.825 
#No. of observations: 84 

#Partition table:
#  Df R.square Adj.R.square Testable
#[a+d+f+g] = X1         4  0.12417      0.07982     TRUE
#[b+d+e+g] = X2         1  0.07042      0.05908     TRUE
#[c+e+f+g] = X3         2  0.03928      0.01556     TRUE
#[a+b+d+e+f+g] = X1+X2  5  0.16970      0.11647     TRUE
#[a+c+d+e+f+g] = X1+X3  6  0.18127      0.11747     TRUE
#[b+c+d+e+f+g] = X2+X3  3  0.10174      0.06805     TRUE
#[a+b+c+d+e+f+g] = All  7  0.23394      0.16338     TRUE
#Individual fractions                                   
#[a] = X1 | X2+X3       4               0.09533     TRUE
#[b] = X2 | X1+X3       1               0.04591     TRUE
#[c] = X3 | X1+X2       2               0.04691     TRUE
#[d]                    0               0.00657    FALSE
#[e]                    0              -0.00927    FALSE
#[f]                    0              -0.03794    FALSE
#[g]                    0               0.01586    FALSE
#[h] = Residuals                        0.83662    FALSE
#Controlling 1 table X                                  
#[a+d] = X1 | X3        4               0.10191     TRUE
#[a+f] = X1 | X2        4               0.05739     TRUE
#[b+d] = X2 | X3        1               0.05249     TRUE
#[b+e] = X2 | X1        1               0.03665     TRUE
#[c+e] = X3 | X1        2               0.03765     TRUE
#[c+f] = X3 | X2        2               0.00897     TRUE


f_sap_varpart <- varpart(f_sap_dist, f_env_1415_model2, f_space_1415_model2, f_host_1415_model2)

#No. of explanatory tables: 3 
#Total variation (SS): 27.699 
#No. of observations: 84 

#Partition table:
#  Df R.square Adj.R.square Testable
#[a+d+f+g] = X1         5  0.16660      0.11317     TRUE
#[b+d+e+g] = X2         2  0.07663      0.05383     TRUE
#[c+e+f+g] = X3         4  0.06993      0.02284     TRUE
#[a+b+d+e+f+g] = X1+X2  7  0.21416      0.14178     TRUE
#[a+c+d+e+f+g] = X1+X3  9  0.23717      0.14439     TRUE
#[b+c+d+e+f+g] = X2+X3  6  0.13705      0.06980     TRUE
#[a+b+c+d+e+f+g] = All 11  0.26126      0.14839     TRUE
#Individual fractions                                   
#[a] = X1 | X2+X3       5               0.07859     TRUE
#[b] = X2 | X1+X3       2               0.00400     TRUE
#[c] = X3 | X1+X2       4               0.00661     TRUE
#[d]                    0               0.04296    FALSE
#[e]                    0               0.02460    FALSE
#[f]                    0               0.00936    FALSE
#[g]                    0              -0.01773    FALSE
#[h] = Residuals                        0.85161    FALSE
#Controlling 1 table X                                  
#[a+d] = X1 | X3        5               0.12155     TRUE
#[a+f] = X1 | X2        5               0.08795     TRUE
#[b+d] = X2 | X3        2               0.04696     TRUE
#[b+e] = X2 | X1        2               0.02861     TRUE
#[c+e] = X3 | X1        4               0.03122     TRUE
#[c+f] = X3 | X2        4               0.01597     TRUE


f_sym_varpart <- varpart(f_sym_dist,f_env_1415_model2[,c(1,3)], f_space_1415_model2, f_host_1415_model2[,2])

#Explanatory tables:
#X1:  f_env_1415_model2
#X2:  f_space_1415_model2
#X3:  f_host_1415_model2 

#No. of explanatory tables: 3 
#Total variation (SS): 26.798 
#No. of observations: 84 

#Partition table:
#  Df R.square Adj.R.square Testable
#[a+d+f+g] = X1         4  0.09136      0.04535     TRUE
#[b+d+e+g] = X2         1  0.04976      0.03818     TRUE
#[c+e+f+g] = X3         2  0.04127      0.01760     TRUE
#[a+b+d+e+f+g] = X1+X2  5  0.12753      0.07160     TRUE
#[a+c+d+e+f+g] = X1+X3  6  0.14419      0.07751     TRUE
#[b+c+d+e+f+g] = X2+X3  3  0.08098      0.04652     TRUE
#[a+b+c+d+e+f+g] = All  7  0.18038      0.10488     TRUE
#Individual fractions                                   
#[a] = X1 | X2+X3       4               0.05837     TRUE
#[b] = X2 | X1+X3       1               0.02738     TRUE
#[c] = X3 | X1+X2       2               0.03329     TRUE
#[d]                    0               0.00155    FALSE
#[e]                    0              -0.00113    FALSE
#[f]                    0              -0.02494    FALSE
#[g]                    0               0.01039    FALSE
#[h] = Residuals                        0.89512    FALSE
#Controlling 1 table X                                  
#[a+d] = X1 | X3        4               0.05991     TRUE
#[a+f] = X1 | X2        4               0.03342     TRUE
#[b+d] = X2 | X3        1               0.02892     TRUE
#[b+e] = X2 | X1        1               0.02624     TRUE
#[c+e] = X3 | X1        2               0.03215     TRUE
#[c+f] = X3 | X2        2               0.00834     TRUE

# figures
plot(f_path_varpart)
plot(f_sap_varpart)
plot(f_sym_varpart)


########## probably not going to create figures for these... they would only be suppliments anyway....
## stacked bar chart of varpar
## pathogens
f_path_varpar_stacked <- NULL
f_path_varpar_stacked$host <- rep('F. roemeri', 6)
f_path_varpar_stacked$variable <- c('climate', 'spatial distance', 'host effects',  'climate and space',
                                    'climate and host', 'space and host')
f_path_varpar_stacked$variance_explained <- c(0.07,0.01,0.01,0.07,0.03,0.03)
f_path_varpar_stacked <- as.data.frame(f_path_varpar_stacked)


f_path_var_par <- ggplot(f_path_varpar_stacked, aes(x = host, y = variance_explained, fill = variable)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) + coord_flip() + labs(fill = NULL, x = NULL, y = 'variance explained') +
  scale_fill_brewer() + theme(legend.position = 'top', legend.text = element_text(size = 14)) + scale_y_continuous(breaks = (seq(0,0.22,0.02))) +
  theme(axis.text.y = element_blank())
f_path_var_par

## saprotroph
f_sap_varpar_stacked <- NULL
f_sap_varpar_stacked$host <- rep('F. roemeri', 6)
f_sap_varpar_stacked$variable <- c('climate', 'spatial distance', 'host effects',  'climate and space',
                                   'climate and host', 'space and host')
f_sap_varpar_stacked$variance_explained <- c(0.07,0.0,0.01,0.05,0.02,0.03)
f_sap_varpar_stacked <- as.data.frame(f_sap_varpar_stacked)


f_sap_var_par <- ggplot(f_sap_varpar_stacked, aes(x = host, y = variance_explained, fill = variable)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) + coord_flip() + labs(fill = NULL, x = NULL, y = 'variance explained') +
  scale_fill_brewer() + theme(legend.position = 'top', legend.text = element_text(size = 14)) + scale_y_continuous(breaks = (seq(0,0.22,0.02))) +
  theme(axis.text.y = element_blank())
f_sap_var_par

f_sym_varpar_stacked <- NULL
f_sym_varpar_stacked$host <- rep('f. roemeri', 6)
f_sym_varpar_stacked$variable <- c('climate', 'spatial distance', 'host effects',  'climate and space',
                                   'climate and host', 'space and host')
f_sym_varpar_stacked$variance_explained <- c(0.05,0.01,0.00,0.03,0.02,0.02)
f_sym_varpar_stacked <- as.data.frame(f_sym_varpar_stacked)


f_sym_var_par <- ggplot(f_sym_varpar_stacked, aes(x = host, y = variance_explained, fill = variable)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) + coord_flip() + labs(fill = NULL, x = NULL, y = 'variance explained') +
  scale_fill_brewer() + theme(legend.position = 'top', legend.text = element_text(size = 14)) + scale_y_continuous(breaks = (seq(0,0.22,0.02))) +
  theme(axis.text.y = element_blank())
f_sym_var_par

############################################### danthonia
## danthonia trophic matricies
d_pathotroph_mat <- read.csv('d_pathotroph_mat.csv', row.names = 1)
d_saprotroph_mat <- read.csv('d_saprotroph_mat.csv', row.names = 1)
d_symbiotroph_mat <- read.csv('d_symbiotroph_mat.csv', row.names = 1)

# jaccard distance
d_path_dist <- vegdist(d_pathotroph_mat, method = 'jaccard', binar = T)
d_sap_dist <- vegdist(d_saprotroph_mat, method = 'jaccard', binar = T)
d_sym_dist <- vegdist(d_symbiotroph_mat, method = 'jaccard', binar = T)

d_path_hel <- decostand(d_pathotroph_mat, 'hellinger')
d_sap_hel <- decostand(d_saprotroph_mat, 'hellinger')
d_sym_hel <- decostand(d_symbiotroph_mat, 'hellinger')

d_bare_path <- rda(d_path_hel ~ 1, d_meta_1415[,-c(1:3,5,6,10:22,36,38:40)])
d_full_path <- rda(d_path_hel ~ ., d_meta_1415[,-c(1:3,5,6,10:22,36,38:40)])
d_mod_path <- ordiR2step(d_bare_path, d_full_path, perm.max = 1000, direction = 'both')
d_mod_path$anova


# variance partitioning
d_path_varpart <- varpart(d_path_dist, d_env_1415_model2, d_space_1415_model2, d_host_1415_model)
#Explanatory tables:
#X1:  d_env_1415_model2
#X2:  d_space_1415_model2
#X3:  d_host_1415_model2 

#No. of explanatory tables: 3 
#Total variation (SS): 22.082 
#No. of observations: 71 

#Partition table:
#                      Df R.square Adj.R.square Testable
#[a+d+f+g] = X1         3  0.19532      0.15929     TRUE
#[b+d+e+g] = X2         1  0.14865      0.13631     TRUE
#[c+e+f+g] = X3         4  0.12725      0.07436     TRUE
#[a+b+d+e+f+g] = X1+X2  4  0.25214      0.20681     TRUE
#[a+c+d+e+f+g] = X1+X3  7  0.28375      0.20417     TRUE
#[b+c+d+e+f+g] = X2+X3  5  0.21798      0.15783     TRUE
#[a+b+c+d+e+f+g] = All  8  0.29892      0.20846     TRUE
#Individual fractions                                   
#[a] = X1 | X2+X3       3               0.05064     TRUE
#[b] = X2 | X1+X3       1               0.00429     TRUE
#[c] = X3 | X1+X2       4               0.00165     TRUE
#[d]                    0               0.07918    FALSE
#[e]                    0               0.04323    FALSE
#[f]                    0               0.01986    FALSE
#[g]                    0               0.00961    FALSE
#[h] = Residuals                        0.79154    FALSE
#Controlling 1 table X                                  
#[a+d] = X1 | X3        3               0.12981     TRUE
#[a+f] = X1 | X2        3               0.07050     TRUE
#[b+d] = X2 | X3        1               0.08347     TRUE
#[b+e] = X2 | X1        1               0.04752     TRUE
#[c+e] = X3 | X1        4               0.04488     TRUE
#[c+f] = X3 | X2        4               0.02151     TRUE



d_sap_varpart <- varpart(d_sap_dist, d_env_1415_model2, d_space_1415_model2, d_host_1415_model)

#Explanatory tables:
#  X1:  d_env_1415_model2
#X2:  d_space_1415_model2
#X3:  d_host_1415_model2 

#No. of explanatory tables: 3 
#Total variation (SS): 23.552 
#No. of observations: 71 

#Partition table:
#                       Df R.square Adj.R.square Testable
#[a+d+f+g] = X1         3  0.19114      0.15492     TRUE
#[b+d+e+g] = X2         1  0.15384      0.14158     TRUE
#[c+e+f+g] = X3         4  0.14228      0.09030     TRUE
#[a+b+d+e+f+g] = X1+X2  4  0.26572      0.22122     TRUE
#[a+c+d+e+f+g] = X1+X3  7  0.29412      0.21568     TRUE
#[b+c+d+e+f+g] = X2+X3  5  0.23179      0.17269     TRUE
#[a+b+c+d+e+f+g] = All  8  0.31636      0.22815     TRUE
#Individual fractions                                   
#[a] = X1 | X2+X3       3               0.05545     TRUE
#[b] = X2 | X1+X3       1               0.01246     TRUE
#[c] = X3 | X1+X2       4               0.00693     TRUE
#[d]                    0               0.06993    FALSE
#[e]                    0               0.05384    FALSE
#[f]                    0               0.02418    FALSE
#[g]                    0               0.00535    FALSE
#[h] = Residuals                        0.77185    FALSE
#Controlling 1 table X                                  
#[a+d] = X1 | X3        3               0.12538     TRUE
#[a+f] = X1 | X2        3               0.07964     TRUE
#[b+d] = X2 | X3        1               0.08239     TRUE
#[b+e] = X2 | X1        1               0.06630     TRUE
#[c+e] = X3 | X1        4               0.06076     TRUE
#[c+f] = X3 | X2        4               0.03111     TRUE



d_sym_varpart <- varpart(d_sym_dist, d_env_1415_model2, d_space_1415_model2, d_host_1415_model)

#Explanatory tables:
#  X1:  d_env_1415_model2
#X2:  d_space_1415_model2
#X3:  d_host_1415_model2 

#No. of explanatory tables: 3 
#Total variation (SS): 19.64 
#No. of observations: 71 

#Partition table:
#                      Df R.square Adj.R.square Testable
#[a+d+f+g] = X1         3  0.11200      0.07223     TRUE
#[b+d+e+g] = X2         1  0.07747      0.06410     TRUE
#[c+e+f+g] = X3         4  0.10128      0.04681     TRUE
#[a+b+d+e+f+g] = X1+X2  4  0.15537      0.10418     TRUE
#[a+c+d+e+f+g] = X1+X3  7  0.19252      0.10280     TRUE
#[b+c+d+e+f+g] = X2+X3  5  0.13633      0.06989     TRUE
#[a+b+c+d+e+f+g] = All  8  0.20507      0.10250     TRUE
#Individual fractions                                   
#[a] = X1 | X2+X3       3               0.03261     TRUE
#[b] = X2 | X1+X3       1              -0.00030     TRUE
#[c] = X3 | X1+X2       4              -0.00168     TRUE
#[d]                    0               0.02338    FALSE
#[e]                    0               0.03225    FALSE
#[f]                    0               0.00747    FALSE
#[g]                    0               0.00878    FALSE
#[h] = Residuals                        0.89750    FALSE
#Controlling 1 table X                                  
#[a+d] = X1 | X3        3               0.05599     TRUE
#[a+f] = X1 | X2        3               0.04008     TRUE
#[b+d] = X2 | X3        1               0.02308     TRUE
#[b+e] = X2 | X1        1               0.03194     TRUE
#[c+e] = X3 | X1        4               0.03057     TRUE
#[c+f] = X3 | X2        4               0.00579     TRUE

# figures
plot(d_path_varpart)
plot(d_sap_varpart)
plot(d_sym_varpart)

## agin, I'm not going to mess with these figures right now...
## stacked bar chart of varpar
d_path_varpar_stacked <- NULL
d_path_varpar_stacked$host <- rep('D. californica', 5)
d_path_varpar_stacked$variable <- c('climate',  'climate and space',
                                    'climate and host', 'space and host','climate space and host')
d_path_varpar_stacked$variance_explained <- c(0.02,0.12,0.03,0.01,0.04)
d_path_varpar_stacked <- as.data.frame(d_path_varpar_stacked)


d_path_var_par <- ggplot(d_path_varpar_stacked, aes(x = host, y = variance_explained, fill = variable)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) + coord_flip() + labs(fill = NULL, x = NULL, y = 'variance explained') +
  scale_fill_brewer() + theme(legend.position = 'top', legend.text = element_text(size = 14)) + scale_y_continuous(breaks = (seq(0,0.22,0.02))) +
  theme(axis.text.y = element_blank())
d_path_var_par

d_sap_varpar_stacked <- NULL
d_sap_varpar_stacked$host <- rep('D. californica', 5)
d_sap_varpar_stacked$variable <- c('climate',  'climate and space',
                                   'climate and host', 'space and host','climate space and host')
d_sap_varpar_stacked$variance_explained <- c(0.03,0.12,0.04,0.03,0.03)
d_sap_varpar_stacked <- as.data.frame(d_sap_varpar_stacked)


d_sap_var_par <- ggplot(d_sap_varpar_stacked, aes(x = host, y = variance_explained, fill = variable)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) + coord_flip() + labs(fill = NULL, x = NULL, y = 'variance explained') +
  scale_fill_brewer() + theme(legend.position = 'top', legend.text = element_text(size = 14)) + scale_y_continuous(breaks = (seq(0,0.22,0.02))) +
  theme(axis.text.y = element_blank())
d_sap_var_par

d_sym_varpar_stacked <- NULL
d_sym_varpar_stacked$host <- rep('D. californica', 5)
d_sym_varpar_stacked$variable <- c('climate',  'climate and space',
                                   'climate and host', 'space and host','climate space and host')
d_sym_varpar_stacked$variance_explained <- c(0.02,0.06,0.01,0.01,0.03)
d_sym_varpar_stacked <- as.data.frame(d_sym_varpar_stacked)


d_sym_var_par <- ggplot(d_sym_varpar_stacked, aes(x = host, y = variance_explained, fill = variable)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) + coord_flip() + labs(fill = NULL, x = NULL, y = 'variance explained') +
  scale_fill_brewer() + theme(legend.position = 'top', legend.text = element_text(size = 14)) + scale_y_continuous(breaks = (seq(0,0.22,0.02))) +
  theme(axis.text.y = element_blank())
d_sym_var_par

#################################################################################################
############################################## RDAs for significance


##############################################
############################################## Festuca

## test of matricies - basic rundown - rda of model component, anova to test significance, plot of 
## rda triplot to show how individual components affect communities
f_path_space_test <- capscale(f_path_dist ~.+ Condition(f_env_1415_model2$temperature.seasonality) +
                                Condition(f_env_1415_model2$precipitation.wettest.quarter)+
                                Condition(f_env_1415_model2$precipitation.seasonality) + 
                                Condition(f_env_1415_model2$MAT) +
                                Condition(f_host_1415_model2$sum_damage) + 
                                Condition(f_host_1415_model2$density)
                              , data.frame(f_space_1415_model2), comm = f_path_hel)
f_path_summary_space_test <- summary(f_path_space_test)
f_path_space_test_samples <- as.data.frame(f_path_summary_space_test$sites)[,1:2]
f_path_space_test_species <- as.data.frame(f_path_summary_space_test$species[,1:2])

RsquareAdj(f_path_space_test)
# 0.04591483 conditional
# 0.05907894

anova.cca(f_path_space_test, step=200, perm.max=200)
# Model: capscale(formula = f_path_dist ~ PCNM1 + PCNM2, data = f_space_1415_model, comm = f_path_hel)
#           Df SumOfSqs      F Pr(>F)    
#  Model    1    1.413 5.2259  0.001 ***
# Residual 76   20.550                         

anova.cca(f_path_space_test, perm.max = 500, by = 'mar')
##          Df SumOfSqs      F Pr(>F)    
#  PCNM1     1    1.413 5.2259  0.001 ***
# Residual 76   20.550               


f_path_env_test <- capscale(f_path_dist ~.+ Condition(f_host_1415_model2$sum_damage) +
                              Condition(f_host_1415_model2$density) +
                              Condition(f_space_1415_model2), f_env_1415_model2, comm = f_path_hel)
f_path_summary_env_test <- summary(f_path_env_test)
f_path_env_test_samples <- as.data.frame(f_path_summary_env_test$sites)[,1:2]
f_path_env_test_species <- as.data.frame(f_path_summary_env_test$species[,1:2])

RsquareAdj(f_path_env_test)
# 0.09533244 conditional
# 0.07982232

anova.cca(f_path_env_test, step=200, perm.max=200)
#Model: capscale(formula = f_path_dist ~ mean.temperature.coldest.quarter + temperature.seasonality + precipitation.wettest.quarter + precipitation.seasonality + Elevation_m, data = f_env_1415_model1, comm = f_path_hel)
#           Df SumOfSqs      F Pr(>F)    
#  Model     4   3.5465 3.279  0.001 ***
#  Residual 76  20.5498                       

anova.cca(f_path_env_test, perm.max = 500, by = 'mar')
#Model: capscale(formula = f_path_dist ~ mean.temperature.coldest.quarter + temperature.seasonality + precipitation.wettest.quarter + precipitation.seasonality + Elevation_m, data = f_env_1415_model1, comm = f_path_hel)
#                                  Df SumOfSqs      F Pr(>F)    
# temperature.seasonality        1   0.9450 3.4950  0.001 ***
# precipitation.wettest.quarter  1   0.3330 1.2316  0.087 .  
# precipitation.seasonality      1   0.2398 0.8867  0.752    
# MAT                            1   0.2856 1.0561  0.337    
# Residual                      76  20.5498               


f_path_host_test <- capscale(f_path_dist ~.+ Condition(f_env_1415_model2$temperature.seasonality) +
                               Condition(f_env_1415_model2$precipitation.wettest.quarter)+
                               Condition(f_env_1415_model2$precipitation.seasonality) + 
                               Condition(f_env_1415_model2$MAT) +
                               Condition(f_space_1415_model2), f_host_1415_model2, comm = f_path_hel)
f_path_summary_host_test <- summary(f_path_host_test)
f_path_host_test_samples <- as.data.frame(f_path_summary_host_test$sites)[,1:2]
f_path_host_test_species <- as.data.frame(f_path_summary_host_test$species[,1:2])

RsquareAdj(f_path_host_test)
#  0.04691188 conditional
#  0.01556119

anova.cca(f_path_host_test, step=200, perm.max=200)
# Model: capscale(formula = f_path_dist ~ plant_area + seeds_produced + sum_damage + density, data = f_host_1415_model, comm = f_path_hel)
#           Df SumOfSqs      F Pr(>F)    
#  Model      2   1.7234 3.1869  0.001 ***
#   Residual 76  20.5498          

anova.cca(f_path_host_test, perm.max = 500, by = 'mar')
##                Df SumOfSqs      F Pr(>F)   
#  sum_damage  1   0.4075 1.5072  0.004 ** 
#  density     1   1.3444 4.9720  0.001 ***
#  Residual   76  20.5498    

## not making these at this moment...
fr_path_env_baseplot <-plot(f_path_env_test)
mult <- attributes(fr_path_env_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
f_path_env_test_samples$site <- f_environ$site

fr_path_env_arrows <- as.data.frame(f_path_summary_env_test$biplot[,1:2])
fr_path_env_arrows$lables <- c('temperature coldest quarter', 'temperature seasonality','precip wettest quarter','precip seasonality', 'elevation')
colnames(fr_path_env_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_path_env_test_plot <- ggplot(data = f_path_env_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c( 'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = fr_env_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_path_env_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                label = lables), size = 3, color='black') +
  geom_point(data = f_path_env_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

fr_path_env_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_f_path_env.tiff", width = 3200, height = 2800, units = "px", res = 500)
fr_path_env_test_plot
dev.off()

## space
fr_path_space_baseplot <-plot(f_path_space_test)
mult <- attributes(fr_path_space_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
f_path_space_test_samples$site <- f_environ$site

fr_path_space_arrows <- as.data.frame(f_path_summary_space_test$biplot[,1:2])
fr_path_space_arrows$lables <- c('PCNM1', 'PCNM2')
colnames(fr_path_space_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_path_space_test_plot <- ggplot(data = f_path_space_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c( 'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = fr_space_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_path_space_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                  label = lables), size = 3, color='black') +
  geom_point(data = f_path_space_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

fr_path_space_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_f_path_space.tiff", width = 3200, height = 2800, units = "px", res = 500)
fr_path_space_test_plot
dev.off()

## host 
fr_path_host_baseplot <-plot(f_path_host_test)
mult <- attributes(fr_path_host_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
f_path_host_test_samples$site <- f_environ$site

fr_path_host_arrows <- as.data.frame(f_path_summary_host_test$biplot[,1:2])
fr_path_host_arrows$lables <- c('plant size', 'reproduction', 'damage', 'density')
colnames(fr_path_host_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_path_host_test_plot <- ggplot(data = f_path_host_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c( 'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = fr_host_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_path_host_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                 label = lables), size = 3, color='black') +
  geom_point(data = f_path_host_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)


fr_path_host_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_f_path_host.tiff", width = 3200, height = 2800, units = "px", res = 500)
fr_path_host_test_plot
dev.off()



# all together now - I'll temporarily remove the legend from two of the three plots so they can fit together nicely
f_path_biplot_multi <- ggarrange(fr_path_env_test_plot, fr_path_space_test_plot, fr_path_host_test_plot, nrow = 1, labels = c('A','B','C'), 
                                 common.legend = T,legend = 'top')
# note I changed the legend size to 10 for the multi pannel figure...
f_path_biplot_multi

## with varpar added
f_path_biplot_plus <-ggarrange(f_path_var_par, 
                               ggarrange(fr_path_env_test_plot, fr_path_space_test_plot, fr_path_host_test_plot, ncol = 3, labels = c('B','C','D'), 
                                         common.legend = T,legend = 'top'), nrow = 2, labels = 'A')
f_path_biplot_plus

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/f_path_biplot_multi_pannel.tiff", width = 2000, height = 1000, units = "px", res = 300)
f_path_biplot_multi
dev.off()

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/f_path_biplot_plus_multi_pannel.tiff", width = 2500, height = 1700, units = "px", res = 300)
f_path_biplot_plus
dev.off()


#################################################### saprotrophs

## test of matricies - basic rundown - rda of model component, anova to test significance, plot of 
## rda triplot to show how individual components affect communities
f_sap_space_test <- capscale(f_sap_dist ~.+ Condition(f_env_1415_model2$temperature.seasonality) +
                               Condition(f_env_1415_model2$precipitation.wettest.quarter)+
                               Condition(f_env_1415_model2$precipitation.seasonality) + 
                               Condition(f_env_1415_model2$MAT) +
                               Condition(f_host_1415_model2$sum_damage) + 
                               Condition(f_host_1415_model2$density)
                             , data.frame(f_space_1415_model2), data.frame(f_space_1415_model2), comm = f_sap_hel)
f_sap_summary_space_test <- summary(f_sap_space_test)
f_sap_space_test_samples <- as.data.frame(f_sap_summary_space_test$sites)[,1:2]
f_sap_space_test_species <- as.data.frame(f_sap_summary_space_test$species[,1:2])

RsquareAdj(f_sap_space_test)
# 0.04871721 conditional
# 0.05334847

anova.cca(f_sap_space_test, step=200, perm.max=200)
# Model: capscale(formula = f_sap_dist ~ PCNM1 + PCNM2, data = f_space_1415_model, comm = f_sap_hel)
#           Df SumOfSqs      F Pr(>F)    
#  Model     1   1.5377 5.3797  0.001 ***
#  Residual 76  21.7232                     

anova.cca(f_sap_space_test, perm.max = 500, by = 'mar')
##           Df SumOfSqs      F Pr(>F)    
#  PCNM1     1   1.5377 5.3797  0.001 ***
#  Residual 76  21.7232             


f_sap_env_test <- capscale(f_sap_dist ~.+ Condition(f_host_1415_model2$sum_damage) +
                             Condition(f_host_1415_model2$density) +
                             Condition(f_space_1415_model2), f_env_1415_model2, comm = f_sap_hel)
f_sap_summary_env_test <- summary(f_sap_env_test)
f_sap_env_test_samples <- as.data.frame(f_sap_summary_env_test$sites)[,1:2]
f_sap_env_test_species <- as.data.frame(f_sap_summary_env_test$species[,1:2])

RsquareAdj(f_sap_env_test)
# 0.07649272 conditional
# 0.0674712

anova.cca(f_sap_env_test, step=200, perm.max=200)
# Model: capscale(formula = f_sap_dist ~ mean.temperature.coldest.quarter + temperature.seasonality + precipitation.wettest.quarter + precipitation.seasonality + Elevation_m, data = f_env_1415_model1, comm = f_sap_hel)
#           Df SumOfSqs      F Pr(>F)    
#  Model     4   3.1855 2.7862  0.001 ***
#  Residual 76  21.7232                         

anova.cca(f_sap_env_test, perm.max = 500, by = 'mar')
##                                   Df SumOfSqs      F Pr(>F)    
#  temperature.seasonality        1   0.9621 3.3659  0.001 ***
#  precipitation.wettest.quarter  1   0.3123 1.0925  0.258    
#  precipitation.seasonality      1   0.2465 0.8623  0.791    
#  MAT                            1   0.2691 0.9416  0.573    
#  Residual                      76  21.7232                  


f_sap_host_test <- capscale(f_sap_dist ~.+ Condition(f_env_1415_model2$temperature.seasonality) +
                              Condition(f_env_1415_model2$precipitation.wettest.quarter)+
                              Condition(f_env_1415_model2$precipitation.seasonality) + 
                              Condition(f_env_1415_model2$MAT) +
                              Condition(f_space_1415_model2), f_host_1415_model2, comm = f_sap_hel)
f_sap_summary_host_test <- summary(f_sap_host_test)
f_sap_host_test_samples <- as.data.frame(f_sap_summary_host_test$sites)[,1:2]
f_sap_host_test_species <- as.data.frame(f_sap_summary_host_test$species[,1:2])

RsquareAdj(f_sap_host_test)
# 0.05048057 conditional
# 0.01245613

anova.cca(f_sap_host_test, step=200, perm.max=200)
# Model: capscale(formula = f_sap_dist ~ plant_area + seeds_produced + sum_damage + density, data = f_host_1415_model, comm = f_sap_hel)
#            Df SumOfSqs     F Pr(>F)    
#  Model    2   1.8857 3.2986  0.001 ***
# Residual 76  21.7232               

anova.cca(f_sap_host_test, perm.max = 500, by = 'mar')
##                 Df SumOfSqs      F Pr(>F)   
#  sum_damage  1   0.3897 1.3634  0.033 *  
#  density     1   1.4860 5.1987  0.001 ***
#  Residual   76  21.7232  

## nope...
fr_sap_env_baseplot <-plot(f_sap_env_test)
mult <- attributes(fr_sap_env_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
f_sap_env_test_samples$site <- f_environ$site

fr_sap_env_arrows <- as.data.frame(f_sap_summary_env_test$biplot[,1:2])
fr_sap_env_arrows$lables <- c('temperature coldest quarter', 'temperature seasonality','precip wettest quarter','precip seasonality', 'elevation')
colnames(fr_sap_env_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_sap_env_test_plot <- ggplot(data = f_sap_env_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c( 'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = fr_env_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_sap_env_arrows,aes(x=x_arrow, y =  y_arrow, 
                                               label = lables), size = 3, color='black') +
  geom_point(data = f_sap_env_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

fr_sap_env_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_f_sap_env.tiff", width = 3200, height = 2800, units = "px", res = 500)
fr_sap_env_test_plot
dev.off()

## space
fr_sap_space_baseplot <-plot(f_sap_space_test)
mult <- attributes(fr_sap_space_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
f_sap_space_test_samples$site <- f_environ$site

fr_sap_space_arrows <- as.data.frame(f_sap_summary_space_test$biplot[,1:2])
fr_sap_space_arrows$lables <- c('PCNM1', 'PCNM2')
colnames(fr_sap_space_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_sap_space_test_plot <- ggplot(data = f_sap_space_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c( 'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = fr_space_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_sap_space_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                 label = lables), size = 3, color='black') +
  geom_point(data = f_sap_space_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

fr_sap_space_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_f_sap_space.tiff", width = 3200, height = 2800, units = "px", res = 500)
fr_sap_space_test_plot
dev.off()

## host 
fr_sap_host_baseplot <-plot(f_sap_host_test)
mult <- attributes(fr_sap_host_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
f_sap_host_test_samples$site <- f_environ$site

fr_sap_host_arrows <- as.data.frame(f_sap_summary_host_test$biplot[,1:2])
fr_sap_host_arrows$lables <- c('plant size', 'reproduction', 'damage', 'density')
colnames(fr_sap_host_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_sap_host_test_plot <- ggplot(data = f_sap_host_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c( 'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = fr_host_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_sap_host_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                label = lables), size = 3, color='black') +
  geom_point(data = f_sap_host_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)


fr_sap_host_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_f_sap_host.tiff", width = 3200, height = 2800, units = "px", res = 500)
fr_sap_host_test_plot
dev.off()



# all together now - I'll temporarily remove the legend from two of the three plots so they can fit together nicely
f_sap_biplot_multi <- ggarrange(fr_sap_env_test_plot, fr_sap_space_test_plot, fr_sap_host_test_plot, nrow = 1, labels = c('A','B','C'), 
                                common.legend = T,legend = 'top')
# note I changed the legend size to 10 for the multi pannel figure...
f_sap_biplot_multi

## with varpar added
f_sap_biplot_plus <-ggarrange(f_sap_var_par, 
                              ggarrange(fr_sap_env_test_plot, fr_sap_space_test_plot, fr_sap_host_test_plot, ncol = 3, labels = c('B','C','D'), 
                                        common.legend = T,legend = 'top'), nrow = 2, labels = 'A')
f_sap_biplot_plus

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/f_sap_biplot_multi_pannel.tiff", width = 2000, height = 1000, units = "px", res = 300)
f_sap_biplot_multi
dev.off()

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/f_sap_biplot_plus_multi_pannel.tiff", width = 2500, height = 1700, units = "px", res = 300)
f_sap_biplot_plus
dev.off()

#################################################### saprotrophs

## test of matricies - basic rundown - rda of model component, anova to test significance, plot of 
## rda triplot to show how individual components affect communities
f_sym_space_test <- capscale(f_sym_dist ~.+ Condition(f_env_1415_model2$temperature.seasonality) +
                               Condition(f_env_1415_model2$precipitation.wettest.quarter)+
                               Condition(f_env_1415_model2$precipitation.seasonality) + 
                               Condition(f_env_1415_model2$MAT) +
                               Condition(f_host_1415_model2$sum_damage) + 
                               Condition(f_host_1415_model2$density)
                             , data.frame(f_space_1415_model2), data.frame(f_space_1415_model2), comm = f_sym_hel)
f_sym_summary_space_test <- summary(f_sym_space_test)
f_sym_space_test_samples <- as.data.frame(f_sym_summary_space_test$sites)[,1:2]
f_sym_space_test_species <- as.data.frame(f_sym_summary_space_test$species[,1:2])

RsquareAdj(f_sym_space_test)
# 0.02742607 conditional
# 0.03819388

anova.cca(f_sym_space_test, step=200, perm.max=200)
# Model: capscale(formula = f_sym_dist ~ PCNM1 + PCNM2, data = f_space_1415_model, comm = f_sym_hel)
#           Df SumOfSqs      F Pr(>F)    
#  PCNM1     1   0.9707 3.3448  0.001 ***
#  Residual 76  22.0569                          

anova.cca(f_sym_space_test, perm.max = 500, by = 'mar')
#           Df SumOfSqs      F Pr(>F)    
#  PCNM1     1   0.9707 3.3448  0.001 ***
#  Residual 76  22.0569                 

f_sym_env_test <- capscale(f_sym_dist ~.+ Condition(f_host_1415_model2$sum_damage) +
                             Condition(f_host_1415_model2$density) +
                             Condition(f_space_1415_model2), f_env_1415_model2, comm = f_sym_hel)
f_sym_summary_env_test <- summary(f_sym_env_test)
f_sym_env_test_samples <- as.data.frame(f_sym_summary_env_test$sites)[,1:2]
f_sym_env_test_species <- as.data.frame(f_sym_summary_env_test$species[,1:2])

RsquareAdj(f_sym_env_test)
# 0.05856714 conditional
# 0.04552225

anova.cca(f_sym_env_test, step=200, perm.max=200)
# Model: capscale(formula = f_sym_dist ~ mean.temperature.coldest.quarter + temperature.seasonality + precipitation.wettest.quarter + precipitation.seasonality + Elevation_m, data = f_env_1415_model1, comm = f_sym_hel)
#           Df SumOfSqs      F Pr(>F)    
#  Model     4   2.6683 2.2985  0.001 ***
#  Residual 76  22.0569                                  

anova.cca(f_sym_env_test, perm.max = 500, by = 'mar')
##                                  Df SumOfSqs      F Pr(>F)    
# temperature.seasonality        1   0.8234 2.8372  0.001 ***
# precipitation.wettest.quarter  1   0.3441 1.1856  0.155    
# precipitation.seasonality      1   0.2219 0.7647  0.893    
# MAT                            1   0.2914 1.0039  0.448    
# Residual                      76  22.0569             


f_sym_host_test <- capscale(f_sym_dist ~.+ Condition(f_env_1415_model2$temperature.seasonality) +
                              Condition(f_env_1415_model2$precipitation.wettest.quarter)+
                              Condition(f_env_1415_model2$precipitation.seasonality) + 
                              Condition(f_env_1415_model2$MAT) +
                              Condition(f_space_1415_model2), f_host_1415_model2, comm = f_sym_hel)
f_sym_summary_host_test <- summary(f_sym_host_test)
f_sym_host_test_samples <- as.data.frame(f_sym_summary_host_test$sites)[,1:2]
f_sym_host_test_species <- as.data.frame(f_sym_summary_host_test$species[,1:2])

RsquareAdj(f_sym_host_test)
# 0.03337275 conditional
# 0.01769304

anova.cca(f_sym_host_test, step=200, perm.max=200)
# Model: capscale(formula = f_sym_dist ~ plant_area + seeds_produced + sum_damage + density, data = f_host_1415_model, comm = f_sym_hel)
#           Df SumOfSqs      F Pr(>F)    
#  Model     2   1.4182 2.4434  0.001 ***
#  Residual 76  22.0569                                   

anova.cca(f_sym_host_test, perm.max = 500, by = 'mar')
##               Df SumOfSqs      F Pr(>F)   
# sum_damage  1   0.2631 0.9066  0.647    
# density     1   1.0939 3.7692  0.001 ***
# Residual   76  22.0569   

# not now...
fr_sym_env_baseplot <-plot(f_sym_env_test)
mult <- attributes(fr_sym_env_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
f_sym_env_test_samples$site <- f_environ$site

fr_sym_env_arrows <- as.data.frame(f_sym_summary_env_test$biplot[,1:2])
fr_sym_env_arrows$lables <- c('temperature coldest quarter', 'temperature seasonality', 'precip wettest quarter','precip seasonality', 'elevation')
colnames(fr_sym_env_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_sym_env_test_plot <- ggplot(data = f_sym_env_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c( 'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = fr_env_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_sym_env_arrows,aes(x=x_arrow, y =  y_arrow, 
                                               label = lables), size = 3, color='black') +
  geom_point(data = f_sym_env_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

fr_sym_env_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_f_sym_env.tiff", width = 3200, height = 2800, units = "px", res = 500)
fr_sym_env_test_plot
dev.off()

## space
fr_sym_space_baseplot <-plot(f_sym_space_test)
mult <- attributes(fr_sym_space_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
f_sym_space_test_samples$site <- f_environ$site

fr_sym_space_arrows <- as.data.frame(f_sym_summary_space_test$biplot[,1:2])
fr_sym_space_arrows$lables <- c('PCNM1', 'PCNM2')
colnames(fr_sym_space_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_sym_space_test_plot <- ggplot(data = f_sym_space_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c( 'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = fr_space_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_sym_space_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                 label = lables), size = 3, color='black') +
  geom_point(data = f_sym_space_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

fr_sym_space_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_f_sym_space.tiff", width = 3200, height = 2800, units = "px", res = 500)
fr_sym_space_test_plot
dev.off()

## host 
fr_sym_host_baseplot <-plot(f_sym_host_test)
mult <- attributes(fr_sym_host_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
f_sym_host_test_samples$site <- f_environ$site

fr_sym_host_arrows <- as.data.frame(f_sym_summary_host_test$biplot[,1:2])
fr_sym_host_arrows$lables <- c('plant size', 'reproduction', 'damage', 'density')
colnames(fr_sym_host_arrows) <- c('x_arrow', 'y_arrow', 'lables')
fr_sym_host_test_plot <- ggplot(data = f_sym_host_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099','#3399FF', '#006699'),
                      labels=c( 'French Flat', 'Roxy Ann', 'Upper Table','Hazel Dell', 'Horse Rock', 'Upper Weir', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = fr_host_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = fr_sym_host_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                label = lables), size = 3, color='black') +
  geom_point(data = f_sym_host_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)


fr_sym_host_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_f_sym_host.tiff", width = 3200, height = 2800, units = "px", res = 500)
fr_sym_host_test_plot
dev.off()



# all together now - I'll temporarily remove the legend from two of the three plots so they can fit together nicely
f_sym_biplot_multi <- ggarrange(fr_sym_env_test_plot, fr_sym_space_test_plot, fr_sym_host_test_plot, nrow = 1, labels = c('A','B','C'), 
                                common.legend = T,legend = 'top')
# note I changed the legend size to 10 for the multi pannel figure...
f_sym_biplot_multi

## with varpar added
f_sym_biplot_plus <-ggarrange(f_sym_var_par, 
                              ggarrange(fr_sym_env_test_plot, fr_sym_space_test_plot, fr_sym_host_test_plot, ncol = 3, labels = c('B','C','D'), 
                                        common.legend = T,legend = 'top'), nrow = 2, labels = 'A')
f_sym_biplot_plus

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/f_sym_biplot_multi_pannel.tiff", width = 2000, height = 1000, units = "px", res = 300)
f_sym_biplot_multi
dev.off()

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/f_sym_biplot_plus_multi_pannel.tiff", width = 2500, height = 1700, units = "px", res = 300)
f_sym_biplot_plus
dev.off()



################################################
################################################ Danthonia

## test of matricies - basic rundown - rda of model component, anova to test significance, plot of 
## rda triplot to show how individual components affect communities
d_path_space_test <- capscale(d_path_dist ~.+Condition(d_env_1415_model2$mean.dpt.coldest.quarter) +
                                Condition(d_env_1415_model2$temperature.seasonality)+
                                Condition(d_env_1415_model2$precipitation.wettest.quarter) + 
                                Condition(d_host_1415_model$plant_area) +
                                Condition(d_host_1415_model$seeds_produced) +
                                Condition(d_host_1415_model$sum_damage) +
                                Condition(d_host_1415_model$density), data.frame(d_space_1415_model2), comm = d_path_hel)
d_path_summary_space_test <- summary(d_path_space_test)
d_path_space_test_samples <- as.data.frame(d_path_summary_space_test$sites)[,1:2]
d_path_space_test_species <- as.data.frame(d_path_summary_space_test$species[,1:2])

RsquareAdj(d_path_space_test)
# 0.00429226

anova.cca(d_path_space_test, step=200, perm.max=200)
# Model: capscale(formula = d_path_dist ~ PCNM1 + PCNM2, data = d_space_1415_model, comm = d_path_hel)
#          Df SumOfSqs      F Pr(>F)    
# Model     1    0.335 1.3416   0.02 *
# Residual 62   15.482                         

anova.cca(d_path_space_test, perm.max = 500, by = 'mar')
##           Df SumOfSqs       F Pr(>F)    
# PCNM1      1    0.335 1.3416  0.015 *
# Residual  62   15.482   


d_path_env_test <- capscale(d_path_dist ~.+Condition(d_host_1415_model$plant_area) +
                              Condition(d_host_1415_model$seeds_produced) +
                              Condition(d_host_1415_model$sum_damage) +
                              Condition(d_host_1415_model$density) + 
                              Condition(d_space_1415_model2), d_env_1415_model2, comm = d_path_hel)
d_path_summary_env_test <- summary(d_path_env_test)
d_path_env_test_samples <- as.data.frame(d_path_summary_env_test$sites)[,1:2]
d_path_env_test_species <- as.data.frame(d_path_summary_env_test$species[,1:2])

RsquareAdj(d_path_env_test)
# 0.05063638

anova.cca(d_path_env_test, step=200, perm.max=200)
#Model: capscale(formula = d_path_dist ~ mean.temperature.coldest.quarter + temperature.seasonality + precipitation.seasonality + Elevation_m, data = d_env_1415_model, comm = d_path_hel)
#          Df SumOfSqs      F Pr(>F)    
# Model     3   1.7874 2.3861  0.001 ***
# Residual 62  15.4815                    

anova.cca(d_path_env_test, perm.max = 500, by = 'mar')
##                                   Df SumOfSqs      F Pr(>F)    
#  mean.dpt.coldest.quarter       1   0.2979 1.1932  0.115  
#  temperature.seasonality        1   0.3300 1.3216  0.032 *
#  precipitation.wettest.quarter  1   0.3279 1.3133  0.031 *
#  Residual                      62  15.4815                  


d_path_host_test <- capscale(d_path_dist ~.+Condition(d_env_1415_model2$mean.dpt.coldest.quarter) +
                               Condition(d_env_1415_model2$temperature.seasonality)+
                               Condition(d_env_1415_model2$precipitation.wettest.quarter) + 
                               Condition(d_space_1415_model2), d_host_1415_model, comm = d_path_hel)
d_path_summary_host_test <- summary(d_path_host_test)
d_path_host_test_samples <- as.data.frame(d_path_summary_host_test$sites)[,1:2]
d_path_host_test_species <- as.data.frame(d_path_summary_host_test$species[,1:2])

RsquareAdj(d_path_host_test)
# 0.0.001649711

anova.cca(d_path_host_test, step=200, perm.max=200)
# Model: capscale(formula = d_path_dist ~ plant_area + seeds_produced + sum_damage + density, data = d_host_1415_model, comm = d_path_hel)
#           Df SumOfSqs      F Pr(>F)    
#  model     4   1.0332 1.0344  0.303
#  Residual 62  15.4815                   

anova.cca(d_path_host_test, perm.max = 500, by = 'mar')
##                Df SumOfSqs      F Pr(>F)    
# plant_area      1   0.2408 0.9642  0.568  
# seeds_produced  1   0.2199 0.8808  0.793  
# sum_damage      1   0.1914 0.7666  0.965  
# density         1   0.3211 1.2860  0.031 *
# Residual       62  15.4815      


# still no...
dc_path_env_baseplot <-plot(d_path_env_test)
mult <- attributes(dc_path_env_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
d_path_env_test_samples$site <- d_environ$site

dc_path_env_arrows <- as.data.frame(d_path_summary_env_test$biplot[,1:2])
dc_path_env_arrows$lables <- c('temperature coldest quarter', 'temperature seasonality','precip seasonality', 'elevation')
colnames(dc_path_env_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_path_env_test_plot <- ggplot(data = d_path_env_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = dc_env_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_path_env_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                label = lables), size = 3, color='black') +
  geom_point(data = d_path_env_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

dc_path_env_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_path_env.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_path_env_test_plot
dev.off()

## space
dc_path_space_baseplot <-plot(d_path_space_test)
mult <- attributes(dc_path_space_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
d_path_space_test_samples$site <- d_environ$site

dc_path_space_arrows <- as.data.frame(d_path_summary_space_test$biplot[,1:2])
dc_path_space_arrows$lables <- c('PCNM1', 'PCNM2')
colnames(dc_path_space_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_path_space_test_plot <- ggplot(data = d_path_space_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = dc_space_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_path_space_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                  label = lables), size = 3, color='black') +
  geom_point(data = d_path_space_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

dc_path_space_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_path_space.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_path_space_test_plot
dev.off()

## host 
dc_path_host_baseplot <-plot(d_path_host_test)
mult <- attributes(dc_path_host_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
d_path_host_test_samples$site <- d_environ$site

dc_path_host_arrows <- as.data.frame(d_path_summary_host_test$biplot[,1:2])
dc_path_host_arrows$lables <- c('plant size', 'reproduction', 'damage', 'density')
colnames(dc_path_host_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_path_host_test_plot <- ggplot(data = d_path_host_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = dc_host_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_path_host_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                 label = lables), size = 3, color='black') +
  geom_point(data = d_path_host_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)


dc_path_host_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_path_host.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_path_host_test_plot
dev.off()



# all together now - I'll temporarily remove the legend from two of the three plots so they can fit together nicely
d_path_biplot_multi <- ggarrange(dc_path_env_test_plot, dc_path_space_test_plot, dc_path_host_test_plot, nrow = 1, labels = c('A','B','C'), 
                                 common.legend = T,legend = 'top')
# note I changed the legend size to 10 for the multi pannel figure...
d_path_biplot_multi

## with varpar added
d_path_biplot_plus <-ggarrange(d_path_var_par, 
                               ggarrange(dc_path_env_test_plot, dc_path_space_test_plot, dc_path_host_test_plot, ncol = 3, labels = c('B','C','D'), 
                                         common.legend = T,legend = 'top'), nrow = 2, labels = 'A')
d_path_biplot_plus

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/d_path_biplot_multi_pannel.tiff", width = 2000, height = 1000, units = "px", res = 300)
d_path_biplot_multi
dev.off()

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/d_path_biplot_plus_multi_pannel.tiff", width = 2500, height = 1700, units = "px", res = 300)
d_path_biplot_plus
dev.off()


#################################################### saprotrophs

## test of matricies - basic rundown - rda of model component, anova to test significance, plot of 
## rda triplot to show how individual components affect communities
d_sap_space_test <- capscale(d_sap_dist ~.+Condition(d_env_1415_model2$mean.dpt.coldest.quarter) +
                               Condition(d_env_1415_model2$temperature.seasonality)+
                               Condition(d_env_1415_model2$precipitation.wettest.quarter) + 
                               Condition(d_host_1415_model$plant_area) +
                               Condition(d_host_1415_model$seeds_produced) +
                               Condition(d_host_1415_model$sum_damage) +
                               Condition(d_host_1415_model$density), data.frame(d_space_1415_model2), comm = d_sap_hel)
d_sap_summary_space_test <- summary(d_sap_space_test)
d_sap_space_test_samples <- as.data.frame(d_sap_summary_space_test$sites)[,1:2]
d_sap_space_test_species <- as.data.frame(d_sap_summary_space_test$species[,1:2])

RsquareAdj(d_sap_space_test)
# 0.0124616

anova.cca(d_sap_space_test, step=200, perm.max=200)
# Model: capscale(formula = d_sap_dist ~ PCNM1 + PCNM2, data = d_space_1415_model, comm = d_sap_hel)
#                   Df SumOfSqs      F Pr(>F)    
#  Model      1   0.5238 2.0171  0.001 ***
#  Residual  62  16.1012                    

anova.cca(d_sap_space_test, perm.max = 500, by = 'mar')
##           Df SumOfSqs       F Pr(>F)    
#  PCNM1     1   0.5238 2.0171  0.001 ***
#  Residual  62  16.1012      


d_sap_env_test <- capscale(d_sap_dist ~.+Condition(d_host_1415_model$plant_area) +
                             Condition(d_host_1415_model$seeds_produced) +
                             Condition(d_host_1415_model$sum_damage) +
                             Condition(d_host_1415_model$density) + 
                             Condition(d_space_1415_model2), d_env_1415_model2, comm = d_sap_hel)
d_sap_summary_env_test <- summary(d_sap_env_test)
d_sap_env_test_samples <- as.data.frame(d_sap_summary_env_test$sites)[,1:2]
d_sap_env_test_species <- as.data.frame(d_sap_summary_env_test$species[,1:2])

RsquareAdj(d_sap_env_test)
# 0.05545318

anova.cca(d_sap_env_test, step=200, perm.max=200)
#Model: capscale(formula = d_sap_dist ~ mean.temperature.coldest.quarter + temperature.seasonality + precipitation.seasonality + Elevation_m, data = d_env_1415_model, comm = d_sap_hel)
#           Df SumOfSqs      F Pr(>F)    
#  Model     3   1.9918 2.5566  0.001 ***
#  Residual 62  16.1012                   

anova.cca(d_sap_env_test, perm.max = 500, by = 'mar')
##                                  Df SumOfSqs      F Pr(>F)    
#  mean.dpt.coldest.quarter       1   0.4371 1.6831  0.001 ***
#  temperature.seasonality        1   0.4547 1.7509  0.002 ** 
#  precipitation.wettest.quarter  1   0.4566 1.7582  0.001 ***
#  Residual                      62  16.1012                 


d_sap_host_test <- capscale(d_sap_dist ~.+Condition(d_env_1415_model2$mean.dpt.coldest.quarter) +
                              Condition(d_env_1415_model2$temperature.seasonality)+
                              Condition(d_env_1415_model2$precipitation.wettest.quarter) + 
                              Condition(d_space_1415_model2), d_host_1415_model, comm = d_sap_hel)
d_sap_summary_host_test <- summary(d_sap_host_test)
d_sap_host_test_samples <- as.data.frame(d_sap_summary_host_test$sites)[,1:2]
d_sap_host_test_species <- as.data.frame(d_sap_summary_host_test$species[,1:2])

RsquareAdj(d_sap_host_test)
# 0.00692566

anova.cca(d_sap_host_test, step=200, perm.max=200)
# Model: capscale(formula = d_sap_dist ~ plant_area + seeds_produced + sum_damage + density, data = d_host_1415_model, comm = d_sap_hel)
#            Df SumOfSqs      F Pr(>F)    
#  Model     4   1.1926 1.1481  0.045 *
#  Residual 62  16.1012                  

anova.cca(d_sap_host_test, perm.max = 500, by = 'mar')
##                Df SumOfSqs      F Pr(>F)    
#  plant_area      1   0.2120 0.8164  0.918   
#  seeds_produced  1   0.2408 0.9274  0.661   
#  sum_damage      1   0.2081 0.8012  0.937   
#  density         1   0.4198 1.6167  0.003 **
#  Residual       62  16.1012

# nah...
dc_sap_env_baseplot <-plot(d_sap_env_test)
mult <- attributes(dc_sap_env_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
d_sap_env_test_samples$site <- d_environ$site

dc_sap_env_arrows <- as.data.frame(d_sap_summary_env_test$biplot[,1:2])
dc_sap_env_arrows$lables <- c('temperature coldest quarter', 'temperature seasonality','precip seasonality', 'elevation')
colnames(dc_sap_env_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_sap_env_test_plot <- ggplot(data = d_sap_env_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = dc_env_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_sap_env_arrows,aes(x=x_arrow, y =  y_arrow, 
                                               label = lables), size = 3, color='black') +
  geom_point(data = d_sap_env_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

dc_sap_env_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_sap_env.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_sap_env_test_plot
dev.off()

## space
dc_sap_space_baseplot <-plot(d_sap_space_test)
mult <- attributes(dc_sap_space_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
d_sap_space_test_samples$site <- d_environ$site

dc_sap_space_arrows <- as.data.frame(d_sap_summary_space_test$biplot[,1:2])
dc_sap_space_arrows$lables <- c('PCNM1', 'PCNM2')
colnames(dc_sap_space_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_sap_space_test_plot <- ggplot(data = d_sap_space_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = dc_space_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_sap_space_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                 label = lables), size = 3, color='black') +
  geom_point(data = d_sap_space_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

dc_sap_space_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_sap_space.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_sap_space_test_plot
dev.off()

## host 
dc_sap_host_baseplot <-plot(d_sap_host_test)
mult <- attributes(dc_sap_host_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
d_sap_host_test_samples$site <- d_environ$site

dc_sap_host_arrows <- as.data.frame(d_sap_summary_host_test$biplot[,1:2])
dc_sap_host_arrows$lables <- c('plant size', 'reproduction', 'damage', 'density')
colnames(dc_sap_host_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_sap_host_test_plot <- ggplot(data = d_sap_host_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = dc_host_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_sap_host_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                label = lables), size = 3, color='black') +
  geom_point(data = d_sap_host_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)


dc_sap_host_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_sap_host.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_sap_host_test_plot
dev.off()



# all together now - I'll temporarily remove the legend from two of the three plots so they can fit together nicely
d_sap_biplot_multi <- ggarrange(dc_sap_env_test_plot, dc_sap_space_test_plot, dc_sap_host_test_plot, nrow = 1, labels = c('A','B','C'), 
                                common.legend = T,legend = 'top')
# note I changed the legend size to 10 for the multi pannel figure...
d_sap_biplot_multi

## with varpar added
d_sap_biplot_plus <-ggarrange(d_sap_var_par, 
                              ggarrange(dc_sap_env_test_plot, dc_sap_space_test_plot, dc_sap_host_test_plot, ncol = 3, labels = c('B','C','D'), 
                                        common.legend = T,legend = 'top'), nrow = 2, labels = 'A')
d_sap_biplot_plus

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/d_sap_biplot_multi_pannel.tiff", width = 2000, height = 1000, units = "px", res = 300)
d_sap_biplot_multi
dev.off()

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/d_sap_biplot_plus_multi_pannel.tiff", width = 2500, height = 1700, units = "px", res = 300)
d_sap_biplot_plus
dev.off()

#################################################### saprotrophs

## test of matricies - basic rundown - rda of model component, anova to test significance, plot of 
## rda triplot to show how individual components affect communities
d_sym_space_test <- capscale(d_sym_dist ~.+Condition(d_env_1415_model2$mean.dpt.coldest.quarter) +
                               Condition(d_env_1415_model2$temperature.seasonality)+
                               Condition(d_env_1415_model2$precipitation.wettest.quarter) + 
                               Condition(d_host_1415_model$plant_area) +
                               Condition(d_host_1415_model$seeds_produced) +
                               Condition(d_host_1415_model$sum_damage) +
                               Condition(d_host_1415_model$density), data.frame(d_space_1415_model2), comm = d_sym_hel)
d_sym_summary_space_test <- summary(d_sym_space_test)
d_sym_space_test_samples <- as.data.frame(d_sym_summary_space_test$sites)[,1:2]
d_sym_space_test_species <- as.data.frame(d_sym_summary_space_test$species[,1:2])

RsquareAdj(d_sym_space_test)
# -8.228892e-05

anova.cca(d_sym_space_test, step=200, perm.max=1000)
# Model: capscale(formula = d_sym_dist ~ PCNM1 + PCNM2, data = d_space_1415_model, comm = d_sym_hel)
#            Df SumOfSqs      F Pr(>F)    
#  Model     1   0.2502 0.9866  0.448
#  Residual 62  15.7223                        

anova.cca(d_sym_space_test, perm.max = 500, by = 'mar')
##            Df SumOfSqs      F Pr(>F)    
#  PCNM1      1   0.2502 0.9866  0.509
# Residual   62  15.7223                   

d_sym_env_test <- capscale(d_sym_dist ~.+Condition(d_host_1415_model$plant_area) +
                             Condition(d_host_1415_model$seeds_produced) +
                             Condition(d_host_1415_model$sum_damage) +
                             Condition(d_host_1415_model$density)+
                             Condition(d_space_1415_model2), d_env_1415_model2, comm = d_sym_hel)
d_sym_summary_env_test <- summary(d_sym_env_test)
d_sym_env_test_samples <- as.data.frame(d_sym_summary_env_test$sites)[,1:2]
d_sym_env_test_species <- as.data.frame(d_sym_summary_env_test$species[,1:2])

RsquareAdj(d_sym_env_test)
# 0.0329098

anova.cca(d_sym_env_test, step=200, perm.max=200)
#Model: capscale(formula = d_sym_dist ~ mean.temperature.coldest.quarter + temperature.seasonality + precipitation.seasonality + Elevation_m, data = d_env_1415_model, comm = d_sym_hel)
#            Df SumOfSqs      F Pr(>F)    
#   Model     3   1.3551 1.7813  0.001 ***
#   Residual 62  15.7223                    

anova.cca(d_sym_env_test, perm.max = 500, by = 'mar')
##                                  Df SumOfSqs      F Pr(>F)    
#  mean.dpt.coldest.quarter       1   0.2320 0.9150  0.609
#  temperature.seasonality        1   0.2459 0.9696  0.510
#  precipitation.wettest.quarter  1   0.2615 1.0312  0.403
#  Residual                      62  15.7223                


d_sym_host_test <- capscale(d_sym_dist ~.+Condition(d_env_1415_model2$mean.dpt.coldest.quarter) +
                              Condition(d_env_1415_model2$temperature.seasonality)+
                              Condition(d_env_1415_model2$precipitation.wettest.quarter) + 
                              Condition(d_space_1415_model2), d_host_1415_model, comm = d_sym_hel)
d_sym_summary_host_test <- summary(d_sym_host_test)
d_sym_host_test_samples <- as.data.frame(d_sym_summary_host_test$sites)[,1:2]
d_sym_host_test_species <- as.data.frame(d_sym_summary_host_test$species[,1:2])

RsquareAdj(d_sym_host_test)
# -0.001217859

anova.cca(d_sym_host_test, step=200, perm.max=200)
# Model: capscale(formula = d_sym_dist ~ plant_area + seeds_produced + sum_damage + density, data = d_host_1415_model, comm = d_sym_hel)
#           Df SumOfSqs      F Pr(>F)    
#  Model     4    0.984 0.9701  0.589
#  Residual 62   15.722                  


anova.cca(d_sym_host_test, perm.max = 500, by = 'mar')
##                Df SumOfSqs      F Pr(>F)    
#  plant_area      1   0.2745 1.0825  0.309
#  seeds_produced  1   0.2393 0.9438  0.555
#  sum_damage      1   0.1752 0.6910  0.965
#  density         1   0.2796 1.1025  0.286
#  Residual       62  15.7223  


# not today...
dc_sym_env_baseplot <-plot(d_sym_env_test)
mult <- attributes(dc_sym_env_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
d_sym_env_test_samples$site <- d_environ$site

dc_sym_env_arrows <- as.data.frame(d_sym_summary_env_test$biplot[,1:2])
dc_sym_env_arrows$lables <- c('temperature coldest quarter', 'temperature seasonality','precip seasonality', 'elevation')
colnames(dc_sym_env_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_sym_env_test_plot <- ggplot(data = d_sym_env_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = dc_env_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_sym_env_arrows,aes(x=x_arrow, y =  y_arrow, 
                                               label = lables), size = 3, color='black') +
  geom_point(data = d_sym_env_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

dc_sym_env_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_sym_env.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_sym_env_test_plot
dev.off()

## space
dc_sym_space_baseplot <-plot(d_sym_space_test)
mult <- attributes(dc_sym_space_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
d_sym_space_test_samples$site <- d_environ$site

dc_sym_space_arrows <- as.data.frame(d_sym_summary_space_test$biplot[,1:2])
dc_sym_space_arrows$lables <- c('PCNM1', 'PCNM2')
colnames(dc_sym_space_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_sym_space_test_plot <- ggplot(data = d_sym_space_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = dc_space_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_sym_space_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                 label = lables), size = 3, color='black') +
  geom_point(data = d_sym_space_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)

dc_sym_space_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_sym_space.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_sym_space_test_plot
dev.off()

## host 
dc_sym_host_baseplot <-plot(d_sym_host_test)
mult <- attributes(dc_sym_host_baseplot$biplot)$arrow.mul
## making my own triplot from the RDA data - 2-21-17 now I need to create the arrows, and assign actual colors scheme
d_sym_host_test_samples$site <- d_environ$site

dc_sym_host_arrows <- as.data.frame(d_sym_summary_host_test$biplot[,1:2])
dc_sym_host_arrows$lables <- c('plant size', 'reproduction', 'damage', 'density')
colnames(dc_sym_host_arrows) <- c('x_arrow', 'y_arrow', 'lables')
dc_sym_host_test_plot <- ggplot(data = d_sym_host_test_samples, aes(x=CAP1, y=CAP2, color = site)) +
  scale_colour_manual(values=c('#00CC00', '#009900', '#003300', '#990066', '#660099', '#006699'),
                      labels=c('French Flat', 'Whetstone', 'Lower Table', "Hazel Dell", 'Horse Rock', 'Whidbey'),
                      guide = guide_legend(reverse=TRUE)) + geom_point(size = 1) + theme_bw()+
  theme(legend.text = element_text(size = 13), legend.title = element_text(size = 15))+
  theme(legend.position = 'none') + ## for multi-pannel plotting...
  geom_segment(data = dc_host_arrows,aes(x = 0, xend =x_arrow, y = 0, yend =y_arrow),
               arrow = arrow(length = unit(0.5, "cm")), colour = "grey") +
  geom_text_repel(data = dc_sym_host_arrows,aes(x=x_arrow, y =  y_arrow, 
                                                label = lables), size = 3, color='black') +
  geom_point(data = d_sym_host_test_species, aes(x=CAP1, y=CAP2), color='black', shape = 3)


dc_sym_host_test_plot

tiff(file = "/Users/grahambailes/grass_endophyte_community/Thesis/figures_for_thesis/rda_d_sym_host.tiff", width = 3200, height = 2800, units = "px", res = 500)
dc_sym_host_test_plot
dev.off()



# all together now - I'll temporarily remove the legend from two of the three plots so they can fit together nicely
d_sym_biplot_multi <- ggarrange(dc_sym_env_test_plot, dc_sym_space_test_plot, dc_sym_host_test_plot, nrow = 1, labels = c('A','B','C'), 
                                common.legend = T,legend = 'top')
# note I changed the legend size to 10 for the multi pannel figure...
d_sym_biplot_multi

## with varpar added
d_sym_biplot_plus <-ggarrange(d_sym_var_par, 
                              ggarrange(dc_sym_env_test_plot, dc_sym_space_test_plot, dc_sym_host_test_plot, ncol = 3, labels = c('B','C','D'), 
                                        common.legend = T,legend = 'top'), nrow = 2, labels = 'A')
d_sym_biplot_plus

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/d_sym_biplot_multi_pannel.tiff", width = 2000, height = 1000, units = "px", res = 300)
d_sym_biplot_multi
dev.off()

tiff(file = "/Users/grahambailes/grass_endophyte_community/manuscript/figures/d_sym_biplot_plus_multi_pannel.tiff", width = 2500, height = 1700, units = "px", res = 300)
d_sym_biplot_plus
dev.off()

## I am interested to look at this problem from a slightly different angle - variance partitioning relies
# upon linear relationships among variables, which may not be the case in nature.  

## I'm not sure how best to go about this, but I imagine that we should build a random forest model for our 

## total fuungal richness, pathogen richness, saprotroph richness, mutualist richness, and the composition
# of each of these.  I know it can be done, I'm just not sure at the moment how to use a dissimilarity matrix as a 
# response variable here... UPDATE - it looks like ther package randomforestSCR





##########################################################################################################
################################################# subset pathogen OTUs
## this hasn't ended up in my manuscript -

## lets look at the most common plant pathogens

# Dean et al 2012 - 'the top 10 fungal pathogens in molecular plant pathology'
#
grep('Magnaporthe', guilds_taxa$Genus) # none
grep('Botrytis', guilds_taxa$Genus)
grep('Puccinia', guilds_taxa$Genus)
grep('Fusarium', guilds_taxa$Genus) # none?
grep('Blumeria', guilds_taxa$Genus)
grep('Mycosphaerella', guilds_taxa$Genus)
grep('Colletotrichum', guilds_taxa$Genus)
grep('Ustilago', guilds_taxa$Genus)
grep('Melampsora', guilds_taxa$Genus)

#
grep('Claviceps', guilds_taxa$Genus)
grep('Pyrenophora', guilds_taxa$Genus)
grep('Sclerospora', guilds_taxa$Genus) # none
grep('Scolicotrichum', guilds_taxa$Genus) # none
grep('Selenophoma', guilds_taxa$Genus)
grep('Septoria', guilds_taxa$Genus)
grep('Sporisorium', guilds_taxa$Genus)
grep('Tilletia', guilds_taxa$Genus)
grep('Alternaria', guilds_taxa$Genus)
grep('Neotyphodium', guilds_taxa)
grep('Cladosporium', guilds_taxa)

f_path <- f
