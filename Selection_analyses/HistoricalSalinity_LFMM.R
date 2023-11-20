library(lfmm)

#import snp data as 012 vcftools output format, add SNP position as colnames and indv as rownames
genotypes<-read.table("~/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/LFMM_analyses/HistoricalInput_exons_lfmm.012")

chr_pos<-read.table("~/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/LFMM_analyses/HistoricalInput_exons_lfmm.012.pos")

indv<-read.table("~/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/LFMM_analyses/HistoricalInput_exons_lfmm.012.indv")

#header2<-as.vector(header[,1])
#indv2<-as.vector(indv[,1])

#colnames(genotypes)<-header2
#rownames(genotypes)<-indv2
genotypes.fin<-genotypes[,-14583]

print(length(genotypes.fin))

#import sample salinity data, etc. 
env.data<-read.csv("~/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/LFMM_analyses/HistoricSamples_EnvData.csv",header=TRUE)

clim<-env.data[,c(13:19)]
pca.clim<-prcomp(clim)
print(pca.clim)
print(summary(pca.clim))

PC1<-pca.clim$x[,1]
PC2<-pca.clim$x[,2]
newclim.data<-cbind(env.data,PC1,PC2)

Env1<-newclim.data$Salinity

SAVS.data<-list('genotype'=genotypes.fin, 'chrpos'=chr_pos, 'phenotype'=Env1)

#mod.lfmm<-lfmm_ridge(Y = SAVS.data$genotype, X = SAVS.data$phenotype, K = 5)
#pv<-lfmm_test(Y = SAVS.data$genotype, X = SAVS.data$phenotype, lfmm=mod.lfmm, calibrate="gif")
#pvalues<-pv$calibrated.pvalue

### Generate GIF by averaging across 50 independent runs for k=1-8
#Define number of latent factors
 k <-4
z.table = NULL
gifs = NULL
for (i in 1:5){
 mod.lfmm <- lfmm_ridge(Y = SAVS.data$genotype, X = SAVS.data$phenotype, K = k)
   pv_COM <- lfmm_test(Y = SAVS.data$genotype, 
                     X = SAVS.data$phenotype, 
                       lfmm = mod.lfmm, 
                       calibrate = "gif")
   z.table <- cbind(z.table, pv_COM$score)
   gifs <- append(gifs,pv_COM$gif)
   rm(mod.lfmm)
   rm(pv_COM)
}
z.score.k6 = apply(z.table, MARGIN = 1, median) #combines z-scores
lambda.k6 = median(z.score.k6^2)/0.456
adjusted.p.values.k6 = pchisq(z.score.k6^2/lambda.k6, df = 1, lower = F) #re-adjust p-values
hist(adjusted.p.values.k6, col = 3, main="93ind.autosomes p-adjusted values, k=8")

# BH on adjusted values
# q is the false discovery rate
q = 0.05
L = length(adjusted.p.values.k6)
w = which(sort(adjusted.p.values.k6) < q * (1:L) / L)
candidates = order(adjusted.p.values.k6)[w]
SNP.COM <-c(1:(length(adjusted.p.values.k6)))
df_COM<- data.frame(SNP.COM,adjusted.p.values.k6,z.score.k6,chr_pos)
candidates.SNPs <- df_COM[candidates,]
print(nrow(candidates.SNPs))

# Write candidates
write.table(candidates.SNPs,file="~/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/LFMM_analyses/SAVS.CAL.Salinity.candidates.txt",row.names=FALSE,col.names = TRUE,sep = "\t")

# Write out all info
write.table(df_COM,file="~/Dropbox/CA_SAVS_temporalGenomics_analyses_ms/Bioinformatics_analyses/LFMM_analyses/SAVS.Salinity.LFMM_results.txt",row.names=FALSE,col.names = TRUE,sep = "\t")

#plot(-log10(pvalues), 
#      pch = 19, 
#      cex = .5, 
#      xlab = "SNP", ylab = "-Log P",
#      col = "grey")
#abline(h=2, col="red", lwd=0.75, lty=2)
#abline(h=1.3,col="blue", lwd=0.75, lty=2)      
# points(example.data$causal.set, 
#       -log10(pvalues)[example.data$causal.set], 
#        type = "h", 
#        col = "blue")
                  