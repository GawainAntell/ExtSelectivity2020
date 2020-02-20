# As long as the code has been run once, the intermediate output is saved as csv files.
# Look at the bottom of the script for a shortcut to read in these files and jump to line 225.																				   
library(beepr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)

###################################################
# Data prep

# Decide whether to use by-species or by-sediment ranges
places <- 'by_sed'  # 'by_species'

# Decide whether to run null (TRUE) or enforce selectivity on range size (FALSE)
h0 <- FALSE

# Decide whether to allow false extinctions (FALSE) or not (TRUE)
static <- FALSE

mus <- c(.25, .5, .75)
sp_loss <- c(0, 1/3, 2/3)

# Stage data
stage_dat <- read.csv('Data/DS1_stage_data.csv', stringsAsFactors=FALSE)
stage_dat$dispersion <- log(stage_dat$dispersion)
stage_dat$area <- log(stage_dat$area)
bins <- stage_dat$stage

# Calculate change in sed coverage, forcing the K to come after the L
mst_standing <- c(stage_dat$dispersion, stage_dat$dispersion[1])
stage_dat$mst_delta <- diff(mst_standing)
area_standing <- c(stage_dat$area, stage_dat$area[1])
stage_dat$area_delta <- diff(area_standing)

if (places=='by_sed'){
  iucn <- read.csv('Data/DS3_IUCN_range_data_0.5degree_raster.csv', stringsAsFactors=FALSE)
  sims <- read.csv('Data/DS5_simulation_range_data_0.5degree_by_sed.csv', stringsAsFactors=FALSE)
} else {
  iucn <- read.csv('Data/DS2_IUCN_range_data_vector.csv', stringsAsFactors=FALSE)
  sims <- read.csv('Data/DS4_simulation_range_data_by_species.csv', stringsAsFactors=FALSE)
  # 1 Nevada endemic species is never sampled by appears in by-sed IUCN data:
  iucn <- iucn[!iucn$species=='Neotamias_palmeri',]
}
spp <- iucn$species
iucn$chull <- log(iucn$chull)

chull_bool <- sims$metric=='chull'
sims <- sims[chull_bool,]

# Transform preserved range sizes (natural log)
# prevent negative log range size values
log_spcl <- function(v){
  v <- log(v)
  v[v < 0] <- 1/exp(100)
  v
}
range_cols <- grep('rep', colnames(sims))
sims[,range_cols] <- apply(sims[,range_cols], 2, log_spcl)

# Each simulation (i.e. parameter combination) has a unique ID
ids <- unique(sims$sim_id)
n_sim <- length(ids)
n_iter <- length(range_cols)

###################################################
# A separate meta analysis should be run for every value of mu.
# This is because all other parameters affect the independent var,
# but mu affects the dependent var (i.e. outcome).
fin <- data.frame(matrix(nrow=0, ncol=14)) 
for (loss in sp_loss){
  for (mu in mus){
    # Iterate probabilistic determination of survival
    # for every simulated dataset of fossil ranges
    # and for every level of selectivity.
    # Calculate accuracy of regression on fossil data.
    
    sim_dat <- data.frame(matrix(nrow=n_sim, ncol=12)) 
    
    for (i in 1:n_sim){
      cur_sim <- sims$sim_id == ids[i]
      fsl <- sims[cur_sim,]
      
      # Separate simulation attributes from values
      atr_cols <- c('species','metric','stage','sites','sim_id')
      atr <- fsl[,atr_cols]
      fsl <- fsl[,! colnames(fsl) %in% atr_cols]
      
      if (h0){
        p_surv <- rep(1-mu, length(spp)) 
      } else {
        # inverse logit to enforce a 'cutoff' of range size that increases extinction risk
        # transform so the given range size is at origin, will have 0 chance of survival
        true_range <- iucn$chull
        area_shift <- true_range - quantile(true_range, mu)
        p_surv <- exp(area_shift) / (1 + exp(area_shift)) 
      }
      names(p_surv) <- spp
      
      cur_bin <- atr$stage[1]
      if (static==FALSE){
        # species must be sampled in consecutive time bin to avoid "false extinction"
        next_bin <- which(bins == cur_bin) + 1 
        # Artificially loop stages, so Kinderhookian follows Lopingian
        if (next_bin > nrow(stage_dat)){
          next_bin <- 1
        }
        next_bin_nm <- bins[next_bin]
        next_bin_abrev <- substr(next_bin_nm, 1, 2)
        n_sites <- sprintf("%03d", atr$sites[1])
        next_id <- paste0(next_bin_abrev, n_sites, 'chull')
        next_sim <- sims$sim_id == next_id
        next_fsl <- sims[next_sim,]
        next_spp <- next_fsl$species #spp[next_spp_bool]
        unsampld <- setdiff(spp, next_spp)
      }
      
      # Run logistic regressions
      # save median regression accuracy among site replicates
      iter_dat <- data.frame(matrix(nrow=n_iter, ncol=8))
      for (j in 1:n_iter){
        # probabilistic survival assignment
        live <- sapply(p_surv, FUN=rbinom, size=1, n=1)
        
        range_fsl <- fsl[,j]
        sim_spp <- atr$species
        if (any(is.na(range_fsl))){
          lost_spp <- which(is.na(range_fsl)==TRUE)
          range_fsl <- range_fsl[-lost_spp]
          sim_spp <- sim_spp[-lost_spp]
        }
        live <- live[sim_spp]
        
        # experimentally manipulate species count
        n_sampled <- length(live)
        if (loss != 0){
          n_lost <- round(n_sampled*loss)
          pos_lost <- sample(1:n_sampled, n_lost)
          range_fsl <- range_fsl[-pos_lost]
          live <- live[-pos_lost]
        }
        
        # observed mu can differ from true value,
        # so save the proportion of observed survivors
        mu_obs <- 1 - sum(live) / length(live)
        
        if (static){
          n_false <- NA
        } else {
          surv <- names(which(live==1))
          false_death <- surv[surv %in% unsampld]
          live[false_death] <- 0
          n_false <- length(false_death)
        }
        
        # at low values of sites and mu, 
        # all species that go extinct may be unsampled (i.e. mu=0)
        if (mu_obs == 0){
          intrcpt_fsl <- beta_fsl <- NA
        } else {
          # define success as survival
          live <- relevel(factor(live), ref='0') 
          
          mod_fsl <- glm(live ~ range_fsl, family = binomial(link='logit')) 
          beta_fsl <- mod_fsl$coefficients['range_fsl']
          intrcpt_fsl <- mod_fsl$coefficients[1]
        }
        
        # average area observed among species
        fsl_range_avg <- mean(range_fsl)
        
        cur_atr <- atr[1,c('stage','sites')]
        iter_dat[j,] <- cbind(cur_atr, fsl_range_avg, n_sampled, n_false,
                              intrcpt_fsl, beta_fsl, mu_obs)
      }
      colnames(iter_dat) <- c('stage','sites','range','n_sampled','n_false_e',
                              'intercept','beta','mu_obs')
      beta_sd <- sd(iter_dat$beta, na.rm = TRUE)
      cols4sumry <- 3:ncol(iter_dat)
      col_meds <- apply(iter_dat[,cols4sumry], 2, median, na.rm=TRUE)
      out <-  cbind(cur_atr, t(col_meds), beta_sd)
      sim_dat[i, 1:ncol(out)] <- out
    } # end loop through simulations
    
    add_cols <- c('dispersion','mst_delta','area_delta')
    colnames(sim_dat) <- c('stage','sites','range','n_sampled','n_false_e',
                           'intercept','beta','mu_obs','beta_sd',add_cols)
    
    # Add stage-level variables to simulation summary data
    for (s in bins){
      stage_row <- bins==s
      vals <- stage_dat[stage_row, add_cols]
      dat_rows <- sim_dat$stage==s
      sim_dat[dat_rows, add_cols] <- vals
    }
    
    sim_dat$mu <- mu
    sim_dat$sp_loss <- loss

    fin <- rbind(fin, sim_dat)
  } # end loop through mu
} # end loop through sp loss thresholds
beep('ping')

# Save output with specific name
fl_nm <- paste0('Data/beta_tables_')		
if (h0){  fl_nm <- paste0(fl_nm,'null')
} else {  fl_nm <- paste0(fl_nm,'selective') }
if (static){  fl_nm <- paste0(fl_nm,'_static_')
} else {  fl_nm <- paste0(fl_nm,'_thruT_') }
fl_nm <- paste0(fl_nm, places,'.csv')

write.csv(fin, fl_nm, row.names=FALSE)

###################################################
# Meta regression on selectivity estimate accuracy

# check for collinearity
cors <- cor(fin[,c('sites','n_sampled','n_false_e','dispersion','mst_delta','area_delta','mu_obs')],
        method='kendall') 
cors > 0.6

# un-comment these lines to skip all lines above and read in files directly
  #for (places in c('by_sed','by_species')){
  #  for (h0 in c(FALSE,TRUE)){
  #    for (static in c(FALSE,TRUE)){
        
        # Collate outputs across mu values
  #      fl_nm <- paste0('Data/data_',places,'/beta_tables_')
  #      if (h0){  fl_nm <- paste0(fl_nm,'null')
  #      } else {  fl_nm <- paste0(fl_nm,'selective')}
  #      if (static){  fl_nm <- paste0(fl_nm,'_static_')
  #      } else {  fl_nm <- paste0(fl_nm,'_thruT_')}
  #      fl_nm <- paste0(fl_nm, places,'.csv')
      
      fin <- read.csv(fl_nm,stringsAsFactors = FALSE)
      
      fin$sites <- log(fin$sites)
      
      for (r in c('beta','beta_sd')){
        
        if (static){
          fmla <- as.formula(paste(r, '~ sites + n_sampled + dispersion + mu'))  # mst_delta 
        } else {
          fmla <- as.formula(paste(r, '~ sites + n_sampled + n_false_e + dispersion + mu'))  # mst_delta 
        }
        mod <- lm(fmla, data=fin)
        main_sum <- summary(mod)
        tbl <- as.data.frame(main_sum$coefficients)
        colnames(tbl) <- c('Estimate','Standard error','t value','p value')
        tbl_nm <- 'Results/beta_meta_'
        if (r=='beta'){ tbl_nm <- paste0(tbl_nm,'accuracy_') }
        if (r=='beta_sd'){ tbl_nm <- paste0(tbl_nm,'precision_') }
        
        if (h0){ tbl_nm <- paste0(tbl_nm,'null')
        } else { tbl_nm <- paste0(tbl_nm,'selective') }
        
        if (static){
          tbl_nm <- paste0(tbl_nm,'_static_') 
          # add a row for n false extinctions, so same row n as through-T models
          tbl <- rbind(tbl, n_false_e=rep(NA,4))
        } else { tbl_nm <- paste0(tbl_nm, '_thruT_')}
        
        r2 <- sprintf('%0.3f', main_sum$r.squared)
        mean_r <- sprintf('%0.3f', mean(fin[,r]))
        ci_r <- sprintf('%0.3f', quantile(fin[,r], c(0.05,0.95)))
        r_print <- paste0(mean_r, ' [', ci_r[1], ', ', ci_r[2], ']')
        tbl <- rbind(tbl, 
                     r2=c(r2,'','',''),
                     resp_mean_ci=c(r_print,'','',''))
        
        write.csv(tbl, paste0(tbl_nm, places,'.csv'), row.names=TRUE)
      }
#    }
#  }
#}

toss <- is.na(fin$beta)
cor(fin[!toss,c('beta','intercept')]) # ,method='kend'

###################################################
# table of meta regression estimates for accuracy

for (r in c('accuracy','precision')){
  for (places in c('by_species','by_sed')){  
    first <- paste0('beta_meta_',r)
    dat_fls <- list.files('Results/')
    tbl_fl_pos <- intersect(grep(first, dat_fls), grep(places, dat_fls))
    tbl_fls <- paste0('Results/', dat_fls[tbl_fl_pos])
    tbls <- lapply(tbl_fls, read.csv, stringsAsFactors=FALSE)
    var_nms <- c("(Intercept)","dispersion","n_sampled","sites","mu","n_false_e","r2","resp_mean_ci")
    l <- lapply(tbls, 
                # order to the same row names; save estimate+SE for variables in one string
                function(x) { 
                  row_ordr <- match(var_nms,x[,1])
                  x <- x[row_ordr,]
                  est <- sprintf('%0.3f', as.numeric(x$Estimate))
                  se <- sprintf('%0.3f', x$Standard.error)
                  str <- paste(est, se, sep=' +/- ')[1:6]
                  c(str, x$Estimate[7:8])
                }
                )
    m <- do.call(cbind, l)
    df <- data.frame(cbind(var_nms, m))
    grp_nms <- c('nullStatic','nullWF','selectStatic','selectWF')
    colnames(df) <- c('variable', grp_nms) 
    df_nm <- paste0('Results/combined_meta_table_',r,'_',places,'.csv')
    write.csv(df, df_nm, row.names = FALSE)
  }
}

###################################################
# plot the observed vs true values of mu (extinction proportion)

for (places in c('by_sed','by_species')){
    for (h0 in c(FALSE,TRUE)){
      for (static in c(FALSE,TRUE)){
  
# Collate outputs across mu values
fl_nm <- paste0('Data/beta_tables_')
if (h0){  fl_nm <- paste0(fl_nm,'null')
} else {  fl_nm <- paste0(fl_nm,'selective')}
if (static){  fl_nm <- paste0(fl_nm,'_static_')
} else {  fl_nm <- paste0(fl_nm,'_thruT_')}
fl_nm <- paste0(fl_nm, places,'.csv')
  
dat <- read.csv(fl_nm, stringsAsFactors = FALSE)
loss_bool <- dat$sp_loss==0
dat <- dat[loss_bool,]
dat$mu <- factor(dat$mu)

bars <- 
  ggplot(dat=dat, aes(x=mu, y=mu_obs)) + 
  theme_bw() +
  scale_x_discrete(name = 'True extinction rate') +
  scale_y_continuous(name = 'Observed extinction rate', limits=c(0,.85), expand=c(0,0)) +
  geom_boxplot(fill='grey') + # aes(fill=mu)
  theme(legend.position = 'none', 
        legend.text = element_text(size=11),
        legend.title = element_text(size=12),
        axis.text = element_text(size=11), 
        axis.title = element_text(size=12)
  ) 

fig_nm <- 'Results/mu_boxplots_'
if (h0){ fig_nm <- paste0(fig_nm,'null_')
} else { fig_nm <- paste0(fig_nm,'selective_') }
if (static){ fig_nm <- paste0(fig_nm,'static_') 
} else { fig_nm <- paste0(fig_nm,'thruT_') }
fig_nm <- paste0(fig_nm, places, '.pdf')

pdf(fig_nm, width=3, height=3)
  print(bars)
dev.off()

}}}

###################################################
# Make panels for beta plot that SD requested

# if skipping the top of the script, read in stage names in chron order:
stage_dat <- read.csv('Data/DS1_stage_data.csv', stringsAsFactors=FALSE)
chron <- stage_dat$stage

for (places in c('by_sed','by_species')){
  for (h0 in c(FALSE,TRUE)){
    for (static in c(FALSE,TRUE)){
      
      # Collate outputs across mu values
      fl_nm <- paste0('Data/beta_tables_')
      if (h0){  fl_nm <- paste0(fl_nm,'null')
      } else {  fl_nm <- paste0(fl_nm,'selective')}
      if (static){  fl_nm <- paste0(fl_nm,'_static_')
      } else {  fl_nm <- paste0(fl_nm,'_thruT_')}
      fl_nm <- paste0(fl_nm, places,'.csv')
      
      dat <- read.csv(fl_nm, stringsAsFactors=FALSE)
      loss_bool <- dat$sp_loss==0
      dat <- dat[loss_bool,]
      if (places=='by_species'){ s <- 10 }
      if (places=='by_sed'){ s <- 55 }
      s_bool <- dat$sites==s
      dat <- dat[s_bool,]
      dat$mu <- factor(dat$mu)
      dat$stage <- factor(dat$stage, levels=chron)
      
      if (h0){ ylim <- c(-0.05, 0.5) } else {
        ylim <- c(-0.05, 0.75)
      }
      
      p <- ggplot(data=dat) + theme_bw() +
        scale_y_continuous(limits=ylim, expand=c(0,0)) +
        geom_point(aes(colour=mu, x=stage, y=beta), position=position_dodge(width=0.3)) +
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25),
              axis.title = element_text(size=12), legend.position = 'top')
      
      SDfig_nm <- 'Results/beta_panels_SD_'
      if (h0){ SDfig_nm <- paste0(SDfig_nm,'null_')
      } else { SDfig_nm <- paste0(SDfig_nm,'selective_') }
      if (static){ SDfig_nm <- paste0(SDfig_nm,'static_') 
      } else { SDfig_nm <- paste0(SDfig_nm,'thruT_') }
      SDfig_nm <- paste0(SDfig_nm, places, '.pdf')
      
      pdf(SDfig_nm, width=3.5, height=4)
      print(p)
      dev.off()
      
}}}
