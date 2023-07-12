library(PerformanceAnalytics)
library(ggplot2)

###############################################
# data prep

# set whether to use range size data from by-species or by-sediment sites
places <- 'by_sed' # 'by_species'

stage_dat <- read.csv('Data/DS1_stage_data.csv', stringsAsFactors=FALSE)
bins <- stage_dat$stage

if (places=='by_sed'){
  iucn <- read.csv('Data/DS3_IUCN_range_data_0.5degree_raster.csv', stringsAsFactors=FALSE)
  sims <- read.csv('Data/DS5_simulation_range_data_0.5degree_by_sed.csv', stringsAsFactors=FALSE)
} else {
  iucn <- read.csv('Data/DS2_IUCN_range_data_vector.csv', stringsAsFactors=FALSE)
  sims <- read.csv('Data/DS4_simulation_range_data_by_species.csv', stringsAsFactors=FALSE)
  # 1 Nevada endemic species is never sampled in the by-sediment site placement approach
  # but appears in by-species IUCN data:
  iucn <- iucn[!iucn$species=='Neotamias_palmeri',]
}
spp <- iucn$species

# Each simulation (i.e. parameter combination) has a unique ID
ids <- unique(sims$sim_id)
n_sim <- length(ids)

###############################################
# Summary statistics to report

# median true range size:
lapply(iucn[,-1], median)

# preserved range size
rep_cols <- grep('rep',colnames(sims))
Ms <- unique(sims$metric)
get_med <- function(m, dat){
  m_bool <- dat$metric==m
  dat <- dat[m_bool, rep_cols]
  rep_med <- apply(dat, 1, median, na.rm=TRUE)
  median(rep_med)
}
sapply(Ms, get_med, dat=sims)

# number of species preserved
get_n <- function(b, dat){
  b_bool <- dat$stage==b
  dat <- dat[b_bool,]
  ranges <- dat[,rep_cols]
  obs <- apply(ranges, 1, 
               function(x) any(! is.na(x)) 
               )
  b_spp <- dat$species[obs]
  length(unique(b_spp))
}
sort(sapply(bins, get_n, dat=sims))

###############################################
# calculate tau correlation between true and preserved geographic range size

get_tau <- function(sim,sim_spp,true){
  # Some species do not preserve so occur only in IUCN data
  # Note: by-species sim data already excludes NA species' ranges
  names(sim) <- sim_spp
  sim <- na.omit(sim)
  true_keep <- names(true) %in% names(sim)
  true <- true[true_keep]
  tau_test <- suppressWarnings(cor.test(sim, true, method='kend')) 
  c(est=tau_test$estimate, p=tau_test$p.value)
}

# Correlate fossil and true range size rank order
cor_dat <- data.frame(matrix(nrow=n_sim, ncol=8))
for (i in 1:n_sim){
  id <- ids[i]
  cur_sim <- sims$sim_id==id
  fsl <- sims[cur_sim,]
  
  # Separate simulation attributes from values
  atr_cols <- c('species','metric','stage','sites','sim_id')
  atr <- fsl[,atr_cols]
  fsl <- fsl[,! colnames(fsl) %in% atr_cols]
  cur_metric <- atr$metric[1]
  
  # Normalize preserved range sizes before taking a mean
  if (cur_metric=='gcell'){
    fsl <- apply(fsl, c(1,2), sqrt)
  } else {
    fsl <- apply(fsl, c(1,2), log)
    
    # prevent negative log range size values
    subst_neg <- function(v){
      v[v < 0] <- 1/exp(100)
      v
    }
    fsl <- apply(fsl, c(1,2), subst_neg)
  }
  
  true_range <- iucn[, cur_metric]
  names(true_range) <- spp
  tau <- apply(fsl, 2, get_tau, true=true_range, sim_spp=atr$species)
  
  med_tau <- quantile(tau['est.tau',], probs=c(0.05,0.5,0.95), na.rm=TRUE)
  # some replicates may give NA tau: standard deviation of zero
  
  med_p <- median(tau['p',], na.rm=TRUE)
  
  # average range size among species, but median among iterations
  range_itrs <- colMeans(fsl, na.rm = TRUE)
  med_range <- median(range_itrs)
  
  cur_atr <- atr[1, c('metric','stage','sites')]
  cor_dat[i,] <- cbind(cur_atr, med_range, t(med_tau), med_p)
}
colnames(cor_dat) <- c('metric','stage','sites','range','tau5','tau','tau95','p')

###############################################
# Output csv table for supplemental material

# Add stage data (dispersion, area, N spp recovered) to tau summary data
cols <- c('dispersion','area','speciesN')
cor_dat[,cols] <- 0
for (s in bins){
  stage_row <- bins==s
  d <- stage_dat[stage_row, cols]
  dat_row <- cor_dat$stage==s
  cor_dat[dat_row, cols] <- d
} 

# sort by sites then dispersion then metric
site_ordr <- order(cor_dat$sites)
cor_dat <- cor_dat[site_ordr,]
mst_ordr <- order(cor_dat$dispersion)
cor_dat <- cor_dat[mst_ordr,]
metric_ordr <- order(cor_dat$metric)
cor_dat <- cor_dat[metric_ordr,]

cor_nm <- paste0('Results/tau_table_',places,'.csv')
write.csv(cor_dat, cor_nm, row.names = FALSE)

summary(cor_dat$tau)
sd(cor_dat$tau)

###############################################
# 1st meta regression: selection of range metric

# round (preserve trailing zeroes) and combine into 1 string
round_spcl <- function(x,dig){
  fmt <- paste0("%.",dig,'f')
  sprintf(fmt, round(x, digits=dig))
}

pick_metric <- lm(tau ~ metric, data=cor_dat)
metric_mod_sumry <- summary(pick_metric)
metric_coef <- as.data.frame(metric_mod_sumry$coefficients)
colnames(metric_coef) <- c('Estimate','StandardError','tValue','pValue')
coef_m <- apply(metric_coef, 2, round_spcl, dig=3)
coef_df <- data.frame(coef_m)
coef_nms <- c('Intercept', unique(sims$metric)[-1])
row.names(coef_df) <- coef_nms
r2 <- round_spcl(metric_mod_sumry$r.squared, 3)
coef_nm <- paste0('Results/tau_meta_on_range_metric_',places,'_',r2,'R2.csv')
write.csv(coef_df, coef_nm, row.names=TRUE)

###############################################
# 2nd meta regression: tau predictions

chull <- cor_dat[cor_dat$metric=='chull',]
chull$log_sites <- log(chull$sites)
chull$log_mst <- log(chull$dispersion)
chull$log_area <- log(chull$area)

# investigate correlations between pairs of potential predictor variables
preds <- chull[,c('log_sites','range','log_mst','log_area','speciesN')]
# cor(preds, method='kendall')
pairs_fig_nm <- paste0('Results/predictor_corr_fig_chull_',places,'.pdf')
pdf(pairs_fig_nm, width=6, height=6)
chart.Correlation(preds, histogram=TRUE, pch=19, method='kendall')
dev.off()

# regression of tau on not-too-correlated predictors
meta2 <- lm(tau ~ log_sites + log_mst + speciesN, data=chull)
meta2sumry <- summary(meta2)
meta2coef <- as.data.frame(meta2sumry$coefficients)
colnames(meta2coef) <- c('Estimate','Standard error','t value','p value')
coef_rnd <- apply(meta2coef, 2, round_spcl, dig=3)
row.names(coef_rnd) <- row.names(meta2coef)
metaR2 <- round_spcl(meta2sumry$r.squared, 3)
meta2nm <- paste0('Results/tau_meta_2_',places,'_',metaR2,'R2.csv')
write.csv(coef_rnd, meta2nm, row.names=TRUE)

###############################################
# tau by stage plot (matching main fig. 2)

# sort stages in order of age
bin_ordr <- stage_dat$stage
cor_dat$stage <- factor(cor_dat$stage, levels=bin_ordr)

# use 2 site levels so the page is less busy
if (places=='by_sed'){
  n_sites <- c(20,403)
} else {
  n_sites <- c(5,15)
}

keep_rows <- cor_dat$sites %in% n_sites
cor_subset <- cor_dat[keep_rows,]

site_labl <- paste(n_sites,'sites')
names(site_labl) <- n_sites

p <- ggplot(data=cor_subset) + theme_bw() +
  scale_y_continuous(limits=c(0,1), expand=c(0,0), name='tau correlation') +
  geom_linerange(aes(x=stage, ymin=tau5, ymax=tau95, color=metric), position=position_dodge(width=0.5)) +
  geom_point(aes(x=stage, y=tau, color=metric), position=position_dodge(width=0.5)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  facet_wrap(~sites, dir='v', labeller=labeller(sites=site_labl)) +
  theme(axis.title = element_text(size=12), legend.position = 'top')

tseries_fig_nm <- paste0('Results/tau_by_stage_',places,'.pdf')
pdf(tseries_fig_nm, width=3.5, height=5)
p
dev.off()

###############################################
# 1. Probability of sampling vs. true range size
# 2. Degree of range truncation

# pick a representative simulation iteration to use as an example
# use a stage with median sediment coverage (mst & area)
# medium-small number of sites (so <100% of cells sampled)
# NB: picking another simulation leads to very similar results
if (places=='by_sed'){
  id <- 'Ki055chull' #'Wo055chull' #'Ki403chull'
} else {
  id <- 'Ki010chull'
}

set.seed(1)
iter <- sample(rep_cols, 1)
ranges <- sims[sims$sim_id==id, iter]

# convert to binary sampled/unsampled
blr_dat <- data.frame(sp=spp, true=log(iucn$chull), sampled=FALSE)
sim_spp <- sims$species[sims$sim_id==id]
# note - some species not sampled in a given iteration of by-sed, or chull=0
lost <- which(ranges==0 | is.na(ranges))
if (length(lost) > 0){
  sim_spp <- sim_spp[-lost]
ranges <- ranges[-lost]
}
found <- spp %in% sim_spp
blr_dat$sampled[found] <- TRUE
blr_dat$preserved[found] <- log(ranges)

# 1. Binary logistic regression of "discovery risk"
blr <- glm(sampled ~ true, data=blr_dat, family='binomial')
#summary(blr)
# convert to odds-ratio
lik <- cbind(coef(blr), confint(blr))
odds <- exp(lik)
odds <- apply(odds, 2, round_spcl, dig=4)
dimnames(odds) <- list(c('Intercept','Range size'), c('Estimate','2.5%','97.5%'))
odds_nm <- paste0('Results/discovery_odds_ratio_mod',id,'_',places,'.csv')
write.csv(odds, odds_nm)

# 2. Plot of preserved vs. true range size

# OLS model to annotate
fmla <- as.formula(preserved~true+0)
lm_eqn <- function(df,fmla){
  m <- lm(fmla, df)
  eq <- substitute(italic(y) == a %.% italic(x)*","~~italic(r)^2~"="~r2, # a + b
                   list(a = format(unname(coef(m)[1]), digits = 2),
                      #  b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq))
}

plot_dat <- subset(blr_dat, !is.na(blr_dat$preserved))
x <- min(plot_dat$true)*1.1
y <- max(plot_dat$preserved)
#y_max <- max(plot_dat$preserved)*1.1 # for sqrt-transformed data
#x_max <- max(plot_dat$true)*1.1
p <- ggplot(aes(x=true, y=preserved), data=plot_dat) + theme_bw() +
#  scale_y_continuous(limits=c(-4,y_max), expand=c(0,0)) + # for sqrt-transformed data
#  scale_x_continuous(limits=c(0,x_max), expand=c(0,0)) +
  geom_smooth(method='lm', formula=y~x+0) +
  geom_point() 
p1 <- p + 
  geom_text(x=x, y=y, label=lm_eqn(plot_dat,fmla), parse=TRUE, size=5)

scatter_nm <- paste0('Results/scatterplot_range_truncation_mod',id,'_',places,'.pdf')
pdf(scatter_nm, width=4, height=5)
p1
dev.off()

# how much are species' ranges truncated, on untransformed scale?
if (places=='by_sed'){
  prsrvd <- log(iucn$chull)*0.82
} else {
  prsrvd <- log(iucn$chull)*0.67
}
ratio <- exp(prsrvd)/iucn$chull
summary(ratio)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.05011 0.06558 0.07968 0.08860 0.10116 0.21592 # by-sed
#0.005323 0.008928 0.014458 0.021075 0.024711 0.276198 # by-sp
