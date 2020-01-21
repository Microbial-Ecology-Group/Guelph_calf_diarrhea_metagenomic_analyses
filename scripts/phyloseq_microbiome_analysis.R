## Phyloseq for unifrac distance comparison
microbiome.ps <- data.ps
#taxa_names(phylum_microbiome) <- fData(AMR_analytic_data[[1]])
otu_table(microbiome.ps) <- otu_table(MRcounts(microbiome, norm=TRUE), taxa_are_rows = TRUE)
phy_tree(microbiome.ps) <- phy_tree(OTU_table)

filtered_microbiome.ps <- tax_glom(microbiome.ps, "phylum")
#filtered_microbiome.ps <- filter_taxa(filtered_microbiome.ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
#filtered_microbiome.ps <- prune_samples(sample_sums(filtered_microbiome.ps)>=20, filtered_microbiome.ps)

#Make an PCoA ordination plot based on abundance based unifrac distances with the following commands
ord_weighted <- ordinate(filtered_microbiome.ps, method="PCoA", distance="unifrac", weighted=TRUE)
weighted.unifrac.scores <- as.data.table(ord_weighted$vectors, row.names = 1)
weighted.unifrac.scores$ID <- row.names(ord_weighted$vectors)
setkey(weighted.unifrac.scores, ID)
weighted.unifrac.scores <- weighted.unifrac.scores[microbiome_metadata]

ord_unweighted <- ordinate(filtered_microbiome.ps, method="PCoA", distance="unifrac", weighted=FALSE)
unweighted.unifrac.scores <- as.data.table(ord_unweighted$vectors, row.names = 1)
unweighted.unifrac.scores$ID <- row.names(ord_unweighted$vectors)
setkey(unweighted.unifrac.scores, ID)
unweighted.unifrac.scores <- unweighted.unifrac.scores[microbiome_metadata]


# as.integer(as.factor(sample.scores$farm_type))




# ##################### Only for testing for both #######################
# # Use phyloseq to print ordination
# p <- plot_ordination(filtered_microbiome.ps, ord, color="farm_type",title="Phyloseq's Weighted Unifrac")
# p <- p + geom_point(size=5) + theme_bw()
# pdf("phyloseq_unifrac.pdf",width=8)
# print(p)
# dev.off()


# 
# ### Testing various distance metrics
# dist_methods <- unlist(distanceMethodList)
# print(dist_methods)
# dist_methods = dist_methods[-which(dist_methods=="ANY")]
# 
# # # These require tree
# # dist_methods[(1:3)]
# # # Remove them from the vector
# # dist_methods <- dist_methods[-(1:3)]
# # # This is the user-defined method:
# # dist_methods["designdist"]
# # # Remove the user-defined distance
# # dist_methods = dist_methods[-which(dist_methods=="ANY")]
# 
# # Loop through each distance method, save each plot to a list, called plist.
# plist <- vector("list", length(dist_methods))
# names(plist) = dist_methods
# for( i in dist_methods ){
#   # Calculate distance matrix
#   iDist <- distance(filtered_microbiome.ps, method=i)
#   # Calculate ordination
#   iMDS  <- ordinate(filtered_microbiome.ps, "MDS", distance=iDist)
#   ## Make plot
#   # Don't carry over previous plot (if error, p will be blank)
#   p <- NULL
#   # Create plot, store as temp variable, p
#   p <- plot_ordination(filtered_microbiome.ps, iMDS, color="Group")
#   # Add title to each plot
#   p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
#   # Save the graphic to file.
#   plist[[i]] = p
# }
# 
# # Combine results and shade according to a metadata variable
# df = ldply(plist, function(x) x$data)
# names(df)[1] <- "distance"
# p = ggplot(df, aes(Axis.1, Axis.2, color="Group", shape="filtered_microbiome.ps"))
# p = p + geom_point(size=3, alpha=0.5)
# p = p + facet_wrap(~distance, scales="free")
# p = p + ggtitle("MDS on various distance metrics for Enterotype dataset")
# p