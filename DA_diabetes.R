##libraries 
library(phyloseq)
library(ggplot2)
library(decontam)
library(microbiome)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(janitor)
library(stringr)
library(ape)
library(mixOmics)
library(ggfortify)
library(plotly)
library(Maaslin2)
library(ANCOMBC)
library(mia)
library(rlist)
library(vegan)
library(tidyr)
library(ALDEx2)
library(microViz)
library(MicrobiomeStat)
library(DESeq2)
library(GUniFrac)

######## For analysis of diabetes differences ##########
##this is written as for diabetes with BMI adjustment 
##for diabetes w/o BMI adj., just remove BMI from the fixed effects
##The exception for this is aldex for which both will be presented 

## For pathways - lines calculating approximate counts are not needed if RPK is used as input

#Read in dataset 
ps_obj <- ps_object ###phyloseq object containing data here 

filtered <- prune_taxa(taxa_sums(ps_obj) > 0, ps_obj) #filter out taxa that are zero for all 

taxo_level_names <- c("Phylum",
											"Class",
											"Order",
											"Family",
											"Genus",
											"Species") #creating list of taxonomy levels 


########## Maaslin2 ##########
meta <- data.frame(filtered@sam_data$Diabetes_status, filtered@sam_data$BMI)
row.names(meta) <- sample_names(filtered)
colnames(meta) <- c("Diabetes", "BMI")

### Method 1 - proportion ####

for (taxon in taxo_level_names) {
	taxo_level <- aggregate_taxa(filtered, level = taxon)
	
	#extract OTU table and convert %Relab to proportion 
	OTU.mix <- (taxo_level@otu_table)/100
	
	filename <- file.path('/path/to/output/folder/', 
												paste0(taxon, "file_name.csv")) #output file 
	
	#### RUN MAASLIN #########
	Maaslin2(OTU.mix,
					 meta,
					 filename,
					 min_abundance = 0.000000000000001,
					 min_prevalence = 0.1,
					 min_variance = 0.0,
					 normalization = "none",
					 transform = "AST",
					 analysis_method = "LM",
					 max_significance = 0.05, #qvalue 
					 fixed_effects = c("Diabetes", "BMI"),
					 correction = "BH",
					 standardize = TRUE,
					 cores = 4, 
					 plot_heatmap = TRUE, 
					 plot_scatter = TRUE,
					 heatmap_first_n = 20,
	)
}

### Method 2 - counts ###

for (taxon in taxo_level_names) {
	taxo_level <- aggregate_taxa(filtered, level = taxon)
	
	#calculate approximate counts 
	OTU.mix <- t((t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library))
	
	filename <- file.path('/path/to/output/folder/', 
												paste0(taxon, "file_name.csv")) #output file 
	
	#### RUN MAASLIN #########
	Maaslin2(OTU.mix,
					 meta,
					 filename,
					 min_abundance = 0.000000000000001,
					 min_prevalence = 0.1,
					 min_variance = 0.0,
					 normalization = "CLR",
					 transform = "none",
					 analysis_method = "LM",
					 max_significance = 0.05, #qvalue 
					 fixed_effects = c("Diabetes", "BMI"),
					 correction = "BH",
					 standardize = TRUE,
					 cores = 4, 
					 plot_heatmap = TRUE, 
					 plot_scatter = TRUE,
					 heatmap_first_n = 20,
	)
}

####### ANCOM-BC2 ########
n = 3
for (categ in taxa.cat) {
	####workaround to prevent error for non-unique taxonomy names 
	if (categ != "Species") {
		taxo_level <- aggregate_taxa(ps_object, level = categ)
		colnames(taxo_level@tax_table)[n] <- categ
		colnames(taxo_level@tax_table)[n-1] <- "ignore"
	} else {
		taxo_level <- aggregate_taxa(ps_object, level = categ)
	}
	###get "counts" by multiplying by library size 
	taxo_level@otu_table <- (t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library)
	
	##run ancombc 
	anc <- ancombc2(taxo_level, #phyloseq object 
									tax_level = categ, #NULL is lowest taxonomy level available 
									fix_formula = 'Diabetes_status + BMI',
									pseudo_sens = TRUE,
									prv_cut = 0.1, #bacteria present in less than 10% samples not included, #0.2 for infants 
									n_cl = 4, ##number of clusters to use
									pairwise = FALSE
	)
	
	res <- anc$res
	
	filename <- file.path('/path/to/output/folder/', 
												paste0(taxon, "file_name.csv")) #output file 
	write.csv(res, 
						filename)
	n = n + 1 
}

####### ALDEX-2 #########
n = 3
for (taxon in taxa.cat) {
	taxo_level <- aggregate_taxa(ps_object, level = taxon)
	
	###get "counts" by multiplying by library size 
	taxo_level@otu_table <- (t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library)
	
	#filter
	taxo_level <- core(taxo_level, detection = 0, prevalence = 0.1)
	
	#needs to be rounded aldex2
	taxo_level@otu_table <- round(taxo_level@otu_table)
	
	x <- aldex.clr(t(taxo_level@otu_table), taxo_level@sam_data$Diabetes_status,
								 mc.samples = 1000, #128 or more for t-test, 1000 for effect size, 16 anova  
								 verbose = T,
								 gamma = 0.25) #0.25 - 1 for most datasets 
	
	##t.test 
	x.tt <- aldex.ttest(x, hist.plot=F, paired.test = FALSE, verbose=FALSE)
	
	#effect size 
	x.effect <- aldex.effect(x, CI=T, verbose=F, include.sample.summary=F, 
													 paired.test=FALSE)
	
	
	x.all <- data.frame(x.tt,x.effect)
	
	filename <- file.path('/path/to/output/folder/', 
												paste0(taxon, "file_name.csv")) #output file 
	
	write.csv(x.all, 
						filename)
	n = n + 1 
}

### for BMI analysis ###

n = 3
for (taxon in taxa.cat) {
	taxo_level <- aggregate_taxa(ps_object, level = taxon)
	
	###get "counts" by multiplying by library size 
	taxo_level@otu_table <- (t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library)
	
	#filter
	taxo_level <- core(taxo_level, detection = 0, prevalence = 0.1)
	
	#needs to be rounded aldex2
	taxo_level@otu_table <- round(taxo_level@otu_table)
	
	##glm
	mm <- model.matrix(~ Diabetes_status + BMI, data.frame(taxo_level@sam_data))
	x.glm <- aldex.clr(t(taxo_level@otu_table), mm, mc.samples=128, denom="all", verbose=T)
	glm.test <- aldex.glm(x.glm, mm, fdr.method='holm')
	glm.eff<- aldex.glm.effect(x.glm)
	
	x.all <- data.frame(glm.test,glm.eff)
	
	filename <- file.path('/path/to/output/folder/', 
												paste0(taxon, "file_name.csv")) #output file 
	
	write.csv(x.all, 
						filename)
	n = n + 1 
}

######### LinDA ###########
#### Method 1 - proportions ###
n = 3
for (taxon in taxa.cat) {
	taxo_level <- aggregate_taxa(ps_object, level = taxon)
	otu.tab <- (taxo_level@otu_table)/100
	
	meta <- cbind.data.frame(Diabetes = taxo_level@sam_data$Diabetes_status,
													 BMI = taxo_level@sam_data$BMI)
	
	linda.obj <- linda(otu.tab, meta, formula = '~Diabetes + BMI', alpha = 0.05, feature.dat.type = "proportion",
										 prev.filter = 0.1)

filename_T2DM <- file.path('/path/to/output/folder/', 
													 paste0(taxon, "file_name.csv")) #output file

write.csv(linda.obj[["output"]][["DiabetesT2DM"]], 
					filename_T2DM)

n = n + 1 
}

#### Method 2 - counts ####

n = 3
for (taxon in taxa.cat) {
	taxo_level <- aggregate_taxa(ps_object, level = taxon)
	###get "counts" by multiplying by library size 
  otu.tab <- round(t((t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library)))
  otu.tab <- as.data.frame(otu.tab)
	
	meta <- cbind.data.frame(Diabetes = taxo_level@sam_data$Diabetes_status,
													 BMI = taxo_level@sam_data$BMI)
	
	linda.obj <- linda(otu.tab, meta, formula = '~ Diabetes + BMI', alpha = 0.05, feature.dat.type = c("count"),
										 prev.filter = 0.2, zero.handling = "imputation")
	
	filename_T2DM <- file.path('/path/to/output/folder/', 
														 paste0(taxon, "file_name.csv")) #output file
	
	write.csv(linda.obj[["output"]][["DiabetesT2DM"]], 
						filename_T2DM)
	
	n = n + 1 
}

#### DESeq2 #######
n = 3
for (taxon in taxa.cat) {
	taxo_level <- aggregate_taxa(ps_object, level = taxon)
	
	###get "counts" by multiplying by library size 
	taxo_level@otu_table <- (t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library)
	
	#filter
	taxo_level <- core(taxo_level, detection = 0, prevalence = 0.1)
	
	taxo1 = phyloseq_to_deseq2(taxo_level, ~ Diabetes_status)
	taxo2 = phyloseq_to_deseq2(taxo_level, ~ Diabetes_status + BMI)
	taxo1 = DESeq(taxo1, test="Wald", fitType="parametric", sfType = "poscounts")
	taxo2 = DESeq(taxo2, test="Wald", fitType="parametric", sfType = "poscounts")
	
	res1 = results(taxo1, cooksCutoff = FALSE)
	res2 = results(taxo2, cooksCutoff = FALSE, name = "Diabetes_status_T2DM_vs_Control")
	res3 = results(taxo2, cooksCutoff = FALSE, name = "BMI")
	resul1 <- data.frame(res1@rownames, res1@listData)
	resul2 <- data.frame(res2@rownames, res2@listData)
	resul3 <- data.frame(res3@rownames, res3@listData)
	
	filename1 <- file.path('/path/to/output/folder/', 
												 paste0(taxon, "file_name1.csv")) #output file
	filename2 <- file.path('/path/to/output/folder/', 
												 paste0(taxon, "file_name2.csv")) #output file
	filename3 <- file.path('/path/to/output/folder/', 
												 paste0(taxon, "file_name3.csv")) #output file
	
	write.csv(resul1, 
						filename1)
	write.csv(resul2, 
						filename2)
	write.csv(resul3, 
						filename3)
	
	
	n = n + 1 
}

####### ZicoSeq ##########

n = 3
for (taxon in taxa.cat) {
	taxo_level <- aggregate_taxa(ps_object, level = taxon)
	taxo_level@otu_table <- (taxo_level@otu_table)/100
	test <- as.matrix(taxo_level@otu_table)
	#prune all zeros 
	taxo_level <- prune_taxa(taxa_sums(taxo_level) > 0, taxo_level)
	
	ZicoSeq.obj.p <- ZicoSeq(meta.dat = data.frame(taxo_level@sam_data), feature.dat = as.matrix(taxo_level@otu_table), 
													 grp.name = 'Diabetes_status', adj.name = 'BMI', ## NULL for diabetes non-adjusted 
													 feature.dat.type = "proportion",
													 # Filter to remove rare taxa
													 prev.filter = 0.1, mean.abund.filter = 0,  max.abund.filter = 0, min.prop = 0, 
													 # Winsorization to replace outliers
													 is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
													 # Posterior sampling will be automatically disabled
													 is.post.sample = FALSE, post.sample.no = 25, 
													 # Use the square-root transformation
													 link.func = list(function (x) x^0.5, function (x) x^0.25), stats.combine.func = max,
													 # Permutation-based multiple testing correction
													 perm.no = 99,  strata = NULL, 
													 # Reference-based multiple stage normalization
													 ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
													 # Family-wise error rate control
													 is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)
	
	##ripped from zicoseq.plot function 
	grp.name <- ZicoSeq.obj.p$grp.name
	meta.dat <- ZicoSeq.obj.p$meta.dat
	
	abundance <- ZicoSeq.obj.p$feature.dat
	prevalence <- apply(ZicoSeq.obj.p$feature.dat, 1, function(x) mean(x > 
																																		 	0))
	
	coefs <- sapply(ZicoSeq.obj.p$coef.list, function(x) x[grep(grp.name, 
																															rownames(x)), ])
	colnames(coefs) <- paste0("coef_Func", 1:length(ZicoSeq.obj.p$coef.list))
	
	R2 <- rowMaxs(ZicoSeq.obj.p$R2)
	
	signs <- t(sign(coefs))[t(ZicoSeq.obj.p$R2 == R2)]
	signs[signs == 0] <- 1
	plot.data <- data.frame(fdr.pvals = ZicoSeq.obj.p[['p.adj.fdr']], 
													R2 = R2 * signs, taxa = rownames(ZicoSeq.obj.p$R2),
													prevalence = prevalence, abundance = abundance)
	
	
	filename <- file.path('/path/to/output/folder/', 
												 paste0(taxon, "file_name.csv")) 
	
	write.csv(plot.data, 
						filename)
	
	n = n + 1 
}

### Method 2 - counts ###
n = 3
for (taxon in taxa.cat) {
	taxo_level <- aggregate_taxa(ps_object, level = taxon)
	taxo_level@otu_table <- t((t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library))
	test <- as.matrix(taxo_level@otu_table)
	#prune all zeros 
	taxo_level <- prune_taxa(taxa_sums(taxo_level) > 0, taxo_level)
	
	ZicoSeq.obj.p <- ZicoSeq(meta.dat = data.frame(taxo_level@sam_data), feature.dat = as.matrix(taxo_level@otu_table), 
													 grp.name = 'Diabetes_status', adj.name = 'BMI', #NULL when not adjusted for BMI
													 feature.dat.type = "count",
													 # Filter to remove rare taxa
													 prev.filter = 0.1, mean.abund.filter = 0,  max.abund.filter = 0, min.prop = 0, 
													 # Winsorization to replace outliers
													 is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
													 # Posterior sampling will be automatically disabled
													 is.post.sample = TRUE, post.sample.no = 25, 
													 # Use the square-root transformation
													 link.func = list(function (x) x^0.5), stats.combine.func = max,
													 # Permutation-based multiple testing correction
													 perm.no = 99,  strata = NULL, 
													 # Reference-based multiple stage normalization
													 ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
													 # Family-wise error rate control
													 is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)
	
	##ripped from zicoseq.plot function 
	grp.name <- ZicoSeq.obj.p$grp.name
	meta.dat <- ZicoSeq.obj.p$meta.dat
	
	abundance <- ZicoSeq.obj.p$feature.dat
	prevalence <- apply(ZicoSeq.obj.p$feature.dat, 1, function(x) mean(x > 
																																		 	0))
	
	coefs <- sapply(ZicoSeq.obj.p$coef.list, function(x) x[grep(grp.name, 
																															rownames(x)), ])
	colnames(coefs) <- paste0("coef_Func", 1:length(ZicoSeq.obj.p$coef.list))
	
	R2 <- rowMaxs(ZicoSeq.obj.p$R2)
	
	signs <- t(sign(coefs))[t(ZicoSeq.obj.p$R2 == R2)]
	signs[signs == 0] <- 1
	plot.data <- data.frame(fdr.pvals = ZicoSeq.obj.p[['p.adj.fdr']], 
													R2 = R2 * signs, taxa = rownames(ZicoSeq.obj.p$R2),
													prevalence = prevalence, abundance = abundance)
	
	filename <- file.path('/path/to/output/folder/', 
												paste0(taxon, "file_name.csv")) 
	
	write.csv(plot.data, 
						filename)
	
	n = n + 1 
}
