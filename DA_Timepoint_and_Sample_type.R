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

######## For analysis of timepoint/sample_type differences ##########
##For sample type - just replace timepoint with sample type throughout script 

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

##### Maaslin2 ########
meta <- data.frame(filtered@sam_data$Timepoint, filtered@sam_data$ID) #extract timepoint and participant ID metadata 
row.names(meta) <- sample_names(filtered)

### Method 1 - proportion ###
for (taxon in taxo_level_names) { #loop for all taxonomy levels 
	taxo_level <- aggregate_taxa(filtered, level = taxon) #aggregate to desired taxonomy level
	
	OTU.mix <- (taxo_level@otu_table)/100 #extract OTU table from phyloseq object convert from %RelAb to proportion for AST 
	
	filename <- file.path('/path/to/output/folder/', 
												paste0(taxon, "file_name.csv")) #output file 
	
	#### RUN MAASLIN 
	Maaslin2(OTU.mix,
					 meta,
					 filename,
					 min_abundance = 0.000000000000001, #prevalence filter did not seem to work without this
					 min_prevalence = 0.1, 
					 min_variance = 0.0,
					 normalization = "none",
					 transform = "AST",
					 analysis_method = "LM",
					 max_significance = 0.05, #qvalue 
					 random_effects = c("ID"),
					 fixed_effects = c("Timepoint"),
					 correction = "BH",
					 standardize = TRUE,
					 cores = 4, 
					 plot_heatmap = TRUE, 
					 plot_scatter = TRUE,
					 heatmap_first_n = 20,
	)
}

## Method 2 - count ##
for (taxon in taxo_level_names) {
	taxo_level <- aggregate_taxa(filtered, level = taxon) #aggregate to desired taxonomy level
	
	OTU.mix <- t((t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library)) #multiple proportion by library size to get approx. count 
	
	filename <- file.path('/path/to/output/folder/', 
												paste0(taxon, "file_name.csv")) #output file 
	
	#### RUN MAASLIN 
	Maaslin2(OTU.mix,
					 meta,
					 filename,
					 min_abundance = 0.000000000000001, #prevalence filter did not seem to work without this
					 min_prevalence = 0.1,
					 min_variance = 0.0,
					 normalization = "CLR",
					 transform = "none",
					 analysis_method = "LM",
					 max_significance = 0.05, #qvalue 
					 random_effects = c("ID"),
					 fixed_effects = c("Timepoint"),
					 correction = "BH",
					 standardize = TRUE,
					 cores = 4, 
					 plot_heatmap = TRUE, 
					 plot_scatter = TRUE,
					 heatmap_first_n = 20,
	)
}

######## ANCOM-BC2 ###########
n = 3
for (categ in taxo_level_names) {
	####workaround to prevent error for non-unique taxonomy names 
	if (categ != "Species") {
		taxo_level <- aggregate_taxa(filtered, level = categ)
		colnames(taxo_level@tax_table)[n] <- categ
		colnames(taxo_level@tax_table)[n-1] <- "ignore"
	} else {
		taxo_level <- aggregate_taxa(filtered, level = categ)
	}
	###get "counts" by multiplying by library size 
	taxo_level@otu_table <- (t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library)
	
	##run ancombc 
	anc <- ancombc2(taxo_level, #phyloseq object 
									tax_level = categ, 
									fix_formula = 'Timepoint',
									rand_formula = '(1 | ID)',
									pseudo_sens = TRUE,
									prv_cut = 0.1, #bacteria present in less than 10% samples not included
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

########## ALDEx2 #############
n = 3
for (taxon in taxo_level_names) {
	taxo_level <- aggregate_taxa(filtered, level = taxon)
	
	###get "counts" by multiplying by library size 
	taxo_level@otu_table <- (t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library)
	
	#filter
	taxo_level <- core(taxo_level, detection = 0, prevalence = 0.1)
	
	#needs to be rounded aldex2
	taxo_level@otu_table <- round(taxo_level@otu_table)
	taxo_level <- ps_arrange(taxo_level, ID ,.target = "sample_data")
	
	#subset to only paired samples - this is only for timepoint, not necessary for sample type 
	taxo_level <- subset_samples(taxo_level, !(ID %in% c('mim004', 
																										 'mim046',
																										 'mim049',
																										 'mim054',
																										 'mim055',
																										 'mim064')))
	
	x <- aldex.clr(t(taxo_level@otu_table), taxo_level@sam_data$Timepoint,
								 mc.samples = 1000, #128 or more for t-test, 1000 for effect size, 16 anova  
								 verbose = T,
								 gamma = 0.25) #0.25 - 1 for most datasets 
	
	##t.test 
	x.tt <- aldex.ttest(x, hist.plot=F, paired.test = TRUE, verbose=FALSE)
	
	#effect size 
	x.effect <- aldex.effect(x, CI=T, verbose=F, include.sample.summary=F, 
													 paired.test=TRUE)
	
	#combine outputs 
	x.all <- data.frame(x.tt,x.effect)
	
	filename <- file.path('/path/to/output/folder/', 
												paste0(taxon, "file_name.csv")) #output file 
	
	write.csv(x.all, 
						filename)
	n = n + 1 
}

############ LinDA ###############
### Method 1 - proportion 
n = 3
for (taxon in taxo_level_names) {
	taxo_level <- aggregate_taxa(filtered, level = taxon)
	otu.tab <- (taxo_level@otu_table)/100 #convert to proportion 
	
	meta <- cbind.data.frame(ID = taxo_level@sam_data$ID,
													 Timepoint = taxo_level@sam_data$Timepoint)
	
	linda.obj <- linda(otu.tab, meta, formula = '~Timepoint + (1|ID)', alpha = 0.05, feature.dat.type = "proportion",
										 prev.filter = 0.1)
	#linda.plot(linda.obj, c('Timepoint', '(1|ID'),
	#					 alpha = 0.05,
	#					 legend = TRUE)
	
	filename <- file.path('/path/to/output/folder/', 
														 paste0(taxon, "file_name.csv")) #output file 
	
	write.csv(linda.obj[["output"]][["TimepointLate"]], 
						filename)
	
	n = n + 1 
}

### Method 2 - counts 
n = 3
for (taxon in taxa.cat) {
	taxo_level <- aggregate_taxa(filtered, level = taxon)
	###get "counts" by multiplying by library size 
	otu.tab <- round(t((t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library))) #convert to counts 
	otu.tab <- as.data.frame(otu.tab)
	
	meta <- cbind.data.frame(ID = taxo_level@sam_data$ID,
													 Timepoint = taxo_level@sam_data$Timepoint)
	
	linda.obj <- linda(otu.tab, meta, formula = '~Timepoint + (1|ID)', alpha = 0.05, feature.dat.type = c("count"),
										 prev.filter = 0.1, zero.handling = "imputation")
	
	filename <- file.path('/path/to/output/folder/', 
												paste0(taxon, "file_name.csv")) #output file 
	
	write.csv(linda.obj[["output"]][["TimepointLate"]], 
						filename)
	
	n = n + 1 
}

########## DESeq2 ##########
n = 3
for (taxon in taxo_level_names) {
	taxo_level <- aggregate_taxa(filtered, level = taxon)
	
	###get "counts" by multiplying by library size 
	taxo_level@otu_table <- (t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library)
	
	#filter
	taxo_level <- core(taxo_level, detection = 0, prevalence = 0.1)
	
	#arange? maybe not necessary for deseq2 
	taxo_level <- ps_arrange(taxo_level, ID ,.target = "sample_data")
	
	taxo1 = phyloseq_to_deseq2(taxo_level, ~ Timepoint + ID)
	taxo1 = DESeq(taxo1, test="Wald", fitType="parametric", sfType = "poscounts")
	
	res1 = results(taxo1, cooksCutoff = FALSE)
	resul1 <- data.frame(res1@rownames, res1@listData)
	
	filename <- file.path('/path/to/output/folder/', 
												paste0(taxon, "file_name.csv")) #output file 
	
	write.csv(resul1, 
						filename)
	
	n = n + 1 
}

########## ZicoSeq #############
### Method 1 - proportion 

n = 3
for (taxon in taxa.cat) {
	taxo_level <- aggregate_taxa(filtered, level = taxon)
	taxo_level@otu_table <- (taxo_level@otu_table)/100
	test <- as.matrix(taxo_level@otu_table)
	
	ZicoSeq.obj.p <- ZicoSeq(meta.dat = data.frame(taxo_level@sam_data), feature.dat = as.matrix(taxo_level@otu_table), 
													 grp.name = 'Timepoint', adj.name = NULL, feature.dat.type = "proportion",
													 # Filter to remove rare taxa
													 prev.filter = 0.1, mean.abund.filter = 0,  max.abund.filter = 0, min.prop = 0, 
													 # Winsorization to replace outliers
													 is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
													 # Posterior sampling will be automatically disabled
													 is.post.sample = FALSE, post.sample.no = 25, 
													 # Use the square-root transformation
													 link.func = list(function (x) x^0.5, function (x) x^0.25), stats.combine.func = max,
													 # Permutation-based multiple testing correction
													 perm.no = 99,  strata = c('ID'), 
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
												paste0(taxon, "file_name.csv")) #output file 
	
	write.csv(plot.data, 
						filename)
	
	n = n + 1 
}

### Method 2 - count ###
n = 3
for (taxon in taxa.cat) {
	taxo_level <- aggregate_taxa(filtered, level = taxon)
	taxo_level@otu_table <- t((t(taxo_level@otu_table)/100)*(taxo_level@sam_data$Library))
	test <- as.matrix(taxo_level@otu_table)
	
	ZicoSeq.obj.p <- ZicoSeq(meta.dat = data.frame(taxo_level@sam_data), feature.dat = as.matrix(taxo_level@otu_table), 
													 grp.name = 'Timepoint', adj.name = NULL, feature.dat.type = "count",
													 # Filter to remove rare taxa
													 prev.filter = 0.1, mean.abund.filter = 0,  max.abund.filter = 0, min.prop = 0, 
													 # Winsorization to replace outliers
													 is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
													 # Posterior sampling will be automatically disabled
													 is.post.sample = TRUE, post.sample.no = 25, 
													 # Use the square-root transformation
													 link.func = list(function (x) x^0.5), stats.combine.func = max,
													 # Permutation-based multiple testing correction
													 perm.no = 99,  strata = c('ID'), 
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
												paste0(taxon, "file_name.csv")) #output file 
	
	write.csv(plot.data, 
						filename)
	
	n = n + 1 
}






