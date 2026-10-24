# =============================================================================
# 04_depression_cascade.R
# Exploratory depression care cascade with domain-specific income ranks.
# Repository name: R/04_depression_cascade.R   |   Development name: 16b_depression_cascade_domain_rank.R
# Code comments are in Portuguese. Cross-references in the comments use the
# development names; README.md maps them to the repository names.
# Paths: all input and output paths are relative to `root`, which defaults to
# the working directory (run from the repository root) or to the PNS_ROOT
# environment variable if it is set.
# =============================================================================

# 16b_depression_cascade_domain_rank.R
# Versão do 16_depression_cascade_exploratory.R com duas correções:
#
#  (1) PONTO 4 - RANKING DE RENDA NO DOMÍNIO. O script 16 usava `income_rank`, ridit
#      ponderado sobre os 90.846 adultos, em todos os modelos. Os RIIs co-primários usam
#      ridit dentro da coorte, então as escalas não eram comparáveis. Aqui cada modelo
#      recebe um ridit ponderado calculado só entre as pessoas do seu domínio de análise
#      (rank_in() antes do desenho; use_rank() troca a coluna no domínio). As fórmulas e a
#      extração de coeficientes continuam usando o nome `income_rank`, para mudar o mínimo.
#      O rank é definido na população do domínio, não nos casos completos do modelo
#      expandido (mesma convenção dos co-primários).
#
#  (2) RÓTULO DA CASCATA. A linha "Diagnóstico sem tratamento atual" estimava a média de
#      depression_current_treatment, que vale 1 para quem ESTÁ em tratamento. O 59,4% é a
#      proporção COM tratamento; sem tratamento = 40,6%. O rótulo foi corrigido.
#
# Todas as linhas alteradas em relação ao 16 estão marcadas com  # [16b]

options(timeout = 600)
options(survey.lonely.psu = 'adjust')
options(survey.adjust.domain.lonely = TRUE)
options(stringsAsFactors = FALSE)

library(PNSIBGE)
library(survey)
library(dplyr)
library(ggplot2)
cat('04_depression_cascade.R (16b_depression_cascade_domain_rank.R) | versão 26/09/2026\n')

root <- Sys.getenv('PNS_ROOT', unset = getwd())                                     # [16b]
data_dir <- file.path(root, 'data_pnsibge')
out_dir <- file.path(root, 'results', 'depression_cascade_domain_rank')        # [16b] não sobrescreve o 16
fig_dir <- file.path(root, 'figures')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

vars <- c(
  'V0001','V0024','UPA_PNS','V0006_PNS','V0015','V0025A','V0029','V00291','V00292','V00293',
  'C006','C008','C009','C011','V0026','VDD004A','VDF003','VDF004','M00203','M00302',
  'I00102','J009','J014','Q092','Q09201','Q09202','Q09301','Q094','Q09502',
  'Q09605','Q09606','Q09607','Q098','Q100','Q10101','Q10202','Q104','Q105','Q106','Q10701','Q109',
  'Q11006','Q11007','Q11008','Q11009','Q11010','Q111','Q11201','Q11405','Q11406','Q11407','Q11408','Q115',
  'N010','N011','N012','N013','N014','N015','N016','N017','N018'
)

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
pns <- PNSIBGE::get_pns(year=2019, selected=TRUE, vars=vars, labels=FALSE,
                         deflator=FALSE, design=FALSE, reload=FALSE, savedir=data_dir)
for (v in setdiff(vars, c('UPA_PNS','V0006_PNS','V0015','V0025A','V0029','V00291','V00292','V00293','M00302'))) {
  if (v %in% names(pns)) pns[[v]] <- num(pns[[v]])
}

region_exact <- c(`11`='Norte',`12`='Norte',`13`='Norte',`14`='Norte',`15`='Norte',`16`='Norte',`17`='Norte',
                  `21`='Nordeste',`22`='Nordeste',`23`='Nordeste',`24`='Nordeste',`25`='Nordeste',`26`='Nordeste',`27`='Nordeste',`28`='Nordeste',`29`='Nordeste',
                  `31`='Sudeste',`32`='Sudeste',`33`='Sudeste',`35`='Sudeste',
                  `41`='Sul',`42`='Sul',`43`='Sul',
                  `50`='Centro-Oeste',`51`='Centro-Oeste',`52`='Centro-Oeste',`53`='Centro-Oeste')
pns$region <- factor(unname(region_exact[as.character(pns$V0001)]),
                     levels=c('Norte','Nordeste','Sudeste','Sul','Centro-Oeste'))
pns$sex <- factor(pns$C006, levels=c(1,2), labels=c('Masculino','Feminino'))
pns$age <- pns$C008
pns$age_c <- pns$age - 45
pns$race3 <- factor(ifelse(pns$C009==1, 'Branca', ifelse(pns$C009 %in% 2:5, 'Não_branca', NA_character_)),
                    levels=c('Branca','Não_branca'))
pns$urban <- factor(pns$V0026, levels=c(1,2), labels=c('Urbana','Rural'))
pns$education3 <- cut(pns$VDD004A, breaks=c(0,3,5,7), labels=c('Baixa','Média','Alta'), right=TRUE)
pns$income3 <- cut(pns$VDF004, breaks=c(0,3,4,7), labels=c('Até_1_SM','Mais_1_a_2_SM','Mais_2_SM'), right=TRUE)
pns$income3 <- factor(pns$income3, levels=c('Até_1_SM','Mais_1_a_2_SM','Mais_2_SM'))
pns$income_reais <- pns$VDF003
pns$plan <- ifelse(pns$I00102 %in% c(1,2), ifelse(pns$I00102==1,1,0), NA_real_)

make_weighted_rank <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0 & x >= 0
  r <- rep(NA_real_, length(x))
  ord <- order(x[ok], na.last=NA)
  idx <- which(ok)[ord]
  ww <- w[idx]
  r[idx] <- (cumsum(ww)-0.5*ww)/sum(ww)
  r
}
pns$income_rank_pop <- make_weighted_rank(pns$income_reais, pns$V00291)   # [16b] ranking populacional, só para referência
pns$income_rank <- pns$income_rank_pop                                     # [16b] substituído por domínio em cada modelo

# PHQ-9: respostas 1--4 são convertidas para escores 0--3.
phq_vars <- paste0('N0', sprintf('%02d', 10:18))
phq_raw <- as.data.frame(pns[, phq_vars])
phq_score_items <- lapply(phq_raw, function(x) ifelse(x %in% 1:4, x-1, NA_real_))
phq_score_items <- as.data.frame(phq_score_items)
n_complete <- rowSums(!is.na(phq_score_items))
pns$phq9_total <- ifelse(n_complete==9, rowSums(phq_score_items), NA_real_)
pns$phq10 <- ifelse(!is.na(pns$phq9_total), as.numeric(pns$phq9_total >= 10), NA_real_)
pns$phq15 <- ifelse(!is.na(pns$phq9_total), as.numeric(pns$phq9_total >= 15), NA_real_)
pns$depression_limitation <- ifelse(pns$Q109 %in% 1:5, pns$Q109, NA_real_)
pns$high_depression_limitation <- ifelse(!is.na(pns$depression_limitation), as.numeric(pns$depression_limitation>=3), NA_real_)
pns$self_harm_thoughts <- ifelse(!is.na(pns$N018) & pns$N018 %in% 1:4, as.numeric(pns$N018 >= 2), NA_real_)
pns$self_harm_frequent <- ifelse(!is.na(pns$N018) & pns$N018 %in% 1:4, as.numeric(pns$N018 >= 3), NA_real_)
pns$self_harm_near_daily <- ifelse(!is.na(pns$N018) & pns$N018 %in% 1:4, as.numeric(pns$N018 == 4), NA_real_)

# Diagnóstico, tratamento e trajetória específica da depressão.
pns$depression_dx <- ifelse(pns$Q092==1, 1, ifelse(pns$Q092==2, 0, NA_real_))
pns$depression_current_treatment <- ifelse(pns$depression_dx==1 & rowSums(cbind(pns$Q09605==1, pns$Q09606==1, pns$Q09607==1), na.rm=TRUE)>0, 1,
                                            ifelse(pns$depression_dx==1 & rowSums(cbind(pns$Q09605 %in% 1:2, pns$Q09606 %in% 1:2, pns$Q09607 %in% 1:2), na.rm=TRUE)>0, 0, NA_real_))
pns$depression_regular_care <- ifelse(pns$Q094 %in% 1:3, as.numeric(pns$Q094==1), NA_real_)
pns$depression_no_regular_care <- ifelse(!is.na(pns$depression_regular_care), 1-pns$depression_regular_care, NA_real_)
pns$depression_recent_visit <- ifelse(pns$Q10101 %in% 1:6, as.numeric(pns$Q10101 %in% 1:4), NA_real_)
pns$referral_to_specialist <- ifelse(pns$Q106 %in% 1:3, as.numeric(pns$Q106==1), NA_real_)
pns$specialist_reached <- ifelse(pns$Q106==3, 1,
                                 ifelse(pns$Q106==1 & pns$Q10701 %in% 1:2, 1,
                                        ifelse(pns$Q106==1 & pns$Q10701==3, 0, NA_real_)))
pns$referral_failure <- ifelse(pns$Q106==1 & pns$Q10701 %in% 1:3, as.numeric(pns$Q10701==3), NA_real_)
pns$specialist_reached_all <- ifelse(pns$Q106==3 | (pns$Q106==1 & pns$Q10701 %in% 1:2), 1,
                                      ifelse(pns$Q106==2 | (pns$Q106==1 & pns$Q10701==3), 0, NA_real_))
pns$reason_group <- ifelse(pns$Q09502==1, 'remissao_percebida',
                           ifelse(pns$Q09502==3, 'sintomas_desanimo',
                                  ifelse(pns$Q09502 %in% 2:8, 'barreira_estrutural',
                                         ifelse(pns$Q09502==9, 'outro', NA_character_))))
pns$structural_reason <- ifelse(!is.na(pns$reason_group), as.numeric(pns$reason_group=='barreira_estrutural'), NA_real_)
pns$public_last_depression_visit <- ifelse(pns$Q105 %in% 1:3, ifelse(pns$Q105==1,1,0), NA_real_)
pns$poor_access_reason <- ifelse(pns$Q09502 %in% 1:9, pns$Q09502, NA_real_)
pns$structural_barrier <- ifelse(pns$poor_access_reason %in% 2:8, 1,
                                 ifelse(pns$poor_access_reason==1 | pns$poor_access_reason==9, 0, NA_real_))
pns$critical_depression_gap <- ifelse(pns$phq10==1 & pns$self_harm_thoughts==1 & pns$depression_dx==1 & pns$depression_current_treatment==0, 1,
                                      ifelse(pns$phq10==1 & pns$self_harm_thoughts==1 & pns$depression_dx==1 & pns$depression_current_treatment==1, 0, NA_real_))
pns$hidden_remission_gap <- ifelse(pns$depression_dx==1 & pns$Q094 %in% c(2,3) & pns$Q09502==1 & !is.na(pns$phq10) & pns$depression_current_treatment %in% 0:1,
                                    as.numeric(pns$phq10==1 & pns$depression_current_treatment==0), NA_real_)
pns$hidden_remission_severe <- ifelse(pns$depression_dx==1 & pns$Q094 %in% c(2,3) & pns$Q09502==1 & !is.na(pns$phq15) & pns$depression_current_treatment %in% 0:1,
                                      as.numeric(pns$phq15==1 & pns$depression_current_treatment==0), NA_real_)
pns$hidden_remission_selfharm <- ifelse(pns$depression_dx==1 & pns$Q094 %in% c(2,3) & pns$Q09502==1 & !is.na(pns$self_harm_thoughts) & pns$depression_current_treatment %in% 0:1,
                                        as.numeric(pns$self_harm_thoughts==1 & pns$depression_current_treatment==0), NA_real_)
dx_complete <- pns$depression_dx==1 & pns$Q094 %in% 1:3 & !is.na(pns$phq10) & !is.na(pns$phq15) & !is.na(pns$self_harm_thoughts) & pns$depression_current_treatment %in% 0:1
pns$hidden_remission_gap_all_dx <- ifelse(dx_complete, as.numeric(pns$Q094 %in% c(2,3) & pns$Q09502==1 & pns$phq10==1 & pns$depression_current_treatment==0), NA_real_)
pns$hidden_remission_severe_all_dx <- ifelse(dx_complete, as.numeric(pns$Q094 %in% c(2,3) & pns$Q09502==1 & pns$phq15==1 & pns$depression_current_treatment==0), NA_real_)
pns$hidden_remission_selfharm_all_dx <- ifelse(dx_complete, as.numeric(pns$Q094 %in% c(2,3) & pns$Q09502==1 & pns$self_harm_thoughts==1 & pns$depression_current_treatment==0), NA_real_)
pns$severe_no_treatment <- ifelse(pns$phq15==1 & pns$depression_dx==1 & pns$depression_current_treatment %in% 0:1,
                                  as.numeric(pns$depression_current_treatment==0), NA_real_)
pns$screen_positive_undiagnosed <- ifelse(pns$phq10==1 & pns$depression_dx %in% 0:1, as.numeric(pns$depression_dx==0), NA_real_)

# [16b] ---- Ridit de renda calculado dentro do domínio de cada modelo -----------------
rank_in <- function(mask) {
  r <- rep(NA_real_, nrow(pns))
  i <- which(mask %in% TRUE & !is.na(pns$income_reais) & !is.na(pns$V00291))
  r[i] <- make_weighted_rank(pns$income_reais[i], pns$V00291[i])
  r
}
hidden_vars <- c('hidden_remission_gap','hidden_remission_severe','hidden_remission_selfharm')
hidden_all_vars <- c('hidden_remission_gap_all_dx','hidden_remission_severe_all_dx','hidden_remission_selfharm_all_dx')
for (v in c(hidden_vars, hidden_all_vars)) pns[[paste0('rank_', v)]] <- rank_in(!is.na(pns[[v]]))
pns$rank_crit      <- rank_in(pns$phq10 %in% 1 & pns$self_harm_thoughts %in% 1 & pns$depression_dx %in% 1 & !is.na(pns$critical_depression_gap))
pns$rank_reached   <- rank_in(pns$depression_dx %in% 1 & pns$Q10101 %in% 1:4 & !is.na(pns$specialist_reached_all))
pns$rank_ref       <- rank_in(pns$depression_dx %in% 1 & pns$Q10101 %in% 1:4 & pns$Q106 %in% 1 & !is.na(pns$referral_failure))
pns$rank_reason    <- rank_in(pns$depression_dx %in% 1 & pns$Q094 %in% c(2,3) & !is.na(pns$reason_group))
pns$rank_remission <- rank_in(pns$depression_dx %in% 1 & pns$Q094 %in% c(2,3) & pns$reason_group %in% 'remissao_percebida' & !is.na(pns$phq10))
rank_cols <- grep('^rank_', names(pns), value = TRUE)
write.csv(data.frame(rank_var = rank_cols, n_domain = sapply(rank_cols, function(v) sum(!is.na(pns[[v]])))),
          file.path(out_dir, 'domain_rank_sizes.csv'), row.names = FALSE)

# Troca, dentro de um domínio, a coluna income_rank pelo ridit daquele domínio.
use_rank <- function(design, rank_var) {                                       # [16b]
  design$variables$income_rank <- design$variables[[rank_var]]
  design
}

# Criar desenho oficial depois das recodificações.
des <- PNSIBGE::pns_design(pns)

# Utilitários design-based.
prev <- function(design, var, domain=NULL) {
  f <- if (is.null(domain)) as.formula(paste0('~',var)) else as.formula(paste0('~',var,' + I(',var,'==0)'))
  if (is.null(domain)) {
    z <- svymean(as.formula(paste0('~',var)), design, na.rm=TRUE)
    out <- data.frame(variable=var, estimate=as.numeric(coef(z)[1]), SE=as.numeric(SE(z)[1]))
  } else {
    d <- subset(design, !!as.name(domain)==1)
    z <- svymean(as.formula(paste0('~',var)), d, na.rm=TRUE)
    out <- data.frame(variable=var, domain=domain, estimate=as.numeric(coef(z)[1]), SE=as.numeric(SE(z)[1]))
  }
  out$lower <- out$estimate - qt(.975, df=degf(design))*out$SE
  out$upper <- out$estimate + qt(.975, df=degf(design))*out$SE
  out
}

# Cascata populacional entre pessoas com PHQ-9 >=10.
stages <- data.frame(
  stage=c('PHQ-9 >=10','Diagnosticado com depressão',
          'Diagnóstico COM tratamento atual (depression_current_treatment = 1)',       # [16b] rótulo corrigido
          'Sem cuidado regular','PHQ-9 >=15 sem tratamento','Ideação de morte/autoagressão sem tratamento'),
  variable=c('phq10','depression_dx','depression_current_treatment','depression_no_regular_care','severe_no_treatment','critical_depression_gap'),
  domain=c(NA,'phq10==1','depression_dx==1','depression_dx==1','phq15==1 & depression_dx==1','phq10==1 & self_harm_thoughts==1 & depression_dx==1')
)
stage_rows <- list()
for (i in seq_len(nrow(stages))) {
  v <- stages$variable[i]
  dom <- stages$domain[i]
  if (is.na(dom)) {
    dsub <- des
  } else {
    dsub <- subset(des, eval(parse(text=dom)))
  }
  ok <- !is.na(dsub$variables[[v]])
  dsub <- subset(dsub, ok)
  z <- svymean(as.formula(paste0('~',v)), dsub, na.rm=TRUE)
  est <- as.numeric(coef(z)[1]); se <- as.numeric(SE(z)[1]); df <- degf(dsub)
  crit <- qt(.975, df=df)
  stage_rows[[i]] <- data.frame(stage=stages$stage[i], variable=v, denominator=dom, n_domain=sum(ok), estimate=est, SE=se, lower=est-crit*se, upper=est+crit*se, degf=df)
}
stage_tab <- dplyr::bind_rows(stage_rows)
write.csv(stage_tab, file.path(out_dir,'cascade_stages.csv'), row.names=FALSE)

# Fenótipos de lacuna oculta em toda a coorte diagnosticada.
hidden_rows <- lapply(hidden_vars, function(v) {
  dd <- subset(des, !is.na(get(v)))
  z <- svymean(as.formula(paste0('~',v)), dd, na.rm=TRUE)
  est <- as.numeric(coef(z)[1]); se <- as.numeric(SE(z)[1]); df <- degf(dd); cr <- qt(.975,df=df)
  data.frame(variable=v, n_domain=sum(!is.na(pns[[v]])), estimate=est, SE=se,
             lower=est-cr*se, upper=est+cr*se, degf=df)
})
write.csv(dplyr::bind_rows(hidden_rows), file.path(out_dir,'hidden_remission_gap_overall.csv'), row.names=FALSE)
hidden_all_rows <- lapply(hidden_all_vars, function(v) {
  dd <- subset(des, !is.na(get(v)))
  z <- svymean(as.formula(paste0('~',v)), dd, na.rm=TRUE)
  est <- as.numeric(coef(z)[1]); se <- as.numeric(SE(z)[1]); df <- degf(dd); cr <- qt(.975,df=df)
  data.frame(variable=v, n_domain=sum(!is.na(pns[[v]])), estimate=est, SE=se,
             lower=est-cr*se, upper=est+cr*se, degf=df)
})
write.csv(dplyr::bind_rows(hidden_all_rows), file.path(out_dir,'hidden_remission_gap_all_diagnosed.csv'), row.names=FALSE)
hidden_all_income_rows <- lapply(hidden_all_vars, function(v) {
  dd <- subset(des, !is.na(get(v)) & !is.na(income3))
  z <- svyby(as.formula(paste0('~',v)), ~income3, dd, svymean,
             na.rm=TRUE, vartype=c('se','ci'), keep.names=FALSE)
  names(z)[2] <- 'estimate'; z$variable <- v
  z$n_domain <- sum(!is.na(pns[[v]]) & !is.na(pns$income3)); z
})
write.csv(dplyr::bind_rows(hidden_all_income_rows), file.path(out_dir,'hidden_remission_gap_all_diagnosed_by_income.csv'), row.names=FALSE)
hidden_income_rows <- lapply(hidden_vars, function(v) {
  dd <- subset(des, !is.na(get(v)) & !is.na(income3))
  z <- svyby(as.formula(paste0('~',v)), ~income3, dd, svymean,
             na.rm=TRUE, vartype=c('se','ci'), keep.names=FALSE)
  names(z)[2] <- 'estimate'
  z$variable <- v
  z$n_domain <- sum(!is.na(pns[[v]]) & !is.na(pns$income3))
  z
})
write.csv(dplyr::bind_rows(hidden_income_rows), file.path(out_dir,'hidden_remission_gap_by_income.csv'), row.names=FALSE)
hidden_models <- list(minimum='income_rank + age_c + sex + region',
                      expanded='income_rank + age_c + sex + region + race3 + education3 + urban + plan')
hidden_pr <- do.call(rbind, lapply(hidden_vars, function(v) {
  dd <- use_rank(subset(des, !is.na(get(v)) & !is.na(income_rank)), paste0('rank_', v))   # [16b]
  do.call(rbind, lapply(names(hidden_models), function(nm) {
    formula_txt <- hidden_models[[nm]]
    m <- svyglm(as.formula(paste0(v,' ~ ',formula_txt)), dd,
                family=quasipoisson(link='log'), na.action=na.omit)
    b <- coef(m)['income_rank']; se <- sqrt(diag(vcov(m)))['income_rank']; df <- m$df.residual; cr <- qt(.975,df=df)
    data.frame(variable=v, model=nm, formula=formula_txt, rank_scope=paste0('rank_', v),          # [16b]
               n_domain=sum(m$prior.weights > 0), degf=df,  # [16b] n efetivo do modelo
               PR=exp(b), lower=exp(b-cr*se), upper=exp(b+cr*se), p=2*pt(abs(b/se),df=df,lower.tail=FALSE))
  }))
}))
write.csv(hidden_pr, file.path(out_dir,'hidden_remission_gap_PR_income.csv'), row.names=FALSE)
hidden_all_pr <- do.call(rbind, lapply(hidden_all_vars, function(v) {
  dd <- use_rank(subset(des, !is.na(get(v)) & !is.na(income_rank)), paste0('rank_', v))   # [16b]
  do.call(rbind, lapply(names(hidden_models), function(nm) {
    formula_txt <- hidden_models[[nm]]
    m <- svyglm(as.formula(paste0(v,' ~ ',formula_txt)), dd,
                family=quasipoisson(link='log'), na.action=na.omit)
    b <- coef(m)['income_rank']; se <- sqrt(diag(vcov(m)))['income_rank']; df <- m$df.residual; cr <- qt(.975,df=df)
    data.frame(variable=v, model=nm, formula=formula_txt, rank_scope=paste0('rank_', v),          # [16b]
               n_domain=sum(m$prior.weights > 0), degf=df,  # [16b] n efetivo do modelo
               PR=exp(b), lower=exp(b-cr*se), upper=exp(b+cr*se), p=2*pt(abs(b/se),df=df,lower.tail=FALSE))
  }))
}))
write.csv(hidden_all_pr, file.path(out_dir,'hidden_remission_gap_all_diagnosed_PR_income.csv'), row.names=FALSE)

# Prevalências por faixa de renda para os três fenótipos prioritários.
inc_rows <- list()
for (v in c('screen_positive_undiagnosed','severe_no_treatment','critical_depression_gap')) {
  dsub <- subset(des, !is.na(income3) & !is.na(get(v)))
  z <- svyby(as.formula(paste0('~',v)), ~income3, dsub, svymean, na.rm=TRUE, vartype=c('se','ci'), keep.names=FALSE)
  names(z)[2] <- 'estimate'
  z$variable <- v
  z$n_domain <- sum(!is.na(dsub$variables[[v]]) & !is.na(dsub$variables$income3))
  inc_rows[[v]] <- z
}
inc_tab <- dplyr::bind_rows(inc_rows)
write.csv(inc_tab, file.path(out_dir,'priority_phenotypes_by_income.csv'), row.names=FALSE)

# Modelo da lacuna crítica.
crit_domain <- use_rank(subset(des, phq10==1 & self_harm_thoughts==1 & depression_dx==1 & !is.na(critical_depression_gap) & !is.na(income_rank)),
                        'rank_crit')                                                              # [16b]
models <- list(
  minimum='income_rank + age_c + sex + region',
  expanded='income_rank + age_c + sex + region + race3 + education3 + urban + plan'
)
crit_pr <- do.call(rbind,lapply(names(models), function(nm) {
  m <- svyglm(as.formula(paste0('critical_depression_gap ~ ',models[[nm]])), crit_domain,
              family=quasipoisson(link='log'), na.action=na.omit)
  b <- coef(m)['income_rank']; se <- sqrt(diag(vcov(m)))['income_rank']; df <- m$df.residual; cr <- qt(.975, df=df)
  data.frame(model=nm, formula=models[[nm]], rank_scope='rank_crit', n_domain=sum(m$prior.weights > 0), degf=df, PR=exp(b), lower=exp(b-cr*se), upper=exp(b+cr*se), p=2*pt(abs(b/se),df=df,lower.tail=FALSE))  # [16b] n efetivo do modelo
}))
write.csv(crit_pr, file.path(out_dir,'critical_gap_PR_income.csv'), row.names=FALSE)

# Modelo de falha de encaminhamento.
ref_all_domain <- subset(des, depression_dx==1 & Q10101 %in% 1:4 & !is.na(Q106) & !is.na(income_rank))
ref_domain <- use_rank(subset(des, depression_dx==1 & Q10101 %in% 1:4 & Q106==1 & !is.na(referral_failure) & !is.na(income_rank)),
                       'rank_ref')                                                                # [16b]
# [13/10/2026] O modelo ajustado principal usa as mesmas covariáveis das outras linhas da Table 6 (bloco D).
# O ajuste adicional pelo setor do último atendimento (SUS, Q105) fica como sensibilidade (expanded_SUS).
ref_models <- list(minimum='income_rank + age_c + sex + region',
                   expanded='income_rank + age_c + sex + region + race3 + education3 + urban + plan',
                   expanded_SUS='income_rank + age_c + sex + region + race3 + education3 + urban + plan + public_last_depression_visit')
ref_pr <- do.call(rbind,lapply(names(ref_models), function(nm) {
  formula_txt <- ref_models[[nm]]
  m <- svyglm(as.formula(paste0('referral_failure ~ ',formula_txt)), ref_domain,
              family=quasipoisson(link='log'), na.action=na.omit)
  b <- coef(m)['income_rank']; se <- sqrt(diag(vcov(m)))['income_rank']; df <- m$df.residual; cr <- qt(.975, df=df)
  data.frame(model=nm, formula=formula_txt, rank_scope='rank_ref', n_domain=sum(m$prior.weights > 0), degf=df, PR=exp(b), lower=exp(b-cr*se), upper=exp(b+cr*se), p=2*pt(abs(b/se),df=df,lower.tail=FALSE))  # [16b] n efetivo do modelo
}))
write.csv(ref_pr, file.path(out_dir,'referral_failure_PR_income.csv'), row.names=FALSE)

# [16b] Chegou ao especialista em saúde mental (decisão de 26/09, substitui "foi encaminhado"):
#   1 = Q106 = 3 (todas as consultas para depressão já foram com profissional de saúde mental)
#       ou Q106 = 1 e Q107a = 1 ou 2 (encaminhado e foi a pelo menos uma consulta);
#   0 = Q106 = 2 (não encaminhado) ou Q106 = 1 e Q107a = 3 (encaminhado e não foi a nenhuma).
# "Foi encaminhado" (Q106 = 1) contava como 0 quem já era atendido por especialista (Q106 = 3),
# mais comum entre os mais ricos, e por isso não mede acesso ao especialista.
# Domínio: diagnosticados com atendimento por depressão nos últimos 3 anos (Q10101 = 1 a 4).
reached_domain <- use_rank(subset(des, depression_dx==1 & Q10101 %in% 1:4 & !is.na(specialist_reached_all) & !is.na(rank_reached)),
                           'rank_reached')
reached_models <- list(minimum='income_rank + age_c + sex + region',
                       expanded='income_rank + age_c + sex + region + race3 + education3 + urban + plan')
reached_pr <- do.call(rbind, lapply(names(reached_models), function(nm) {
  formula_txt <- reached_models[[nm]]
  m <- svyglm(as.formula(paste0('specialist_reached_all ~ ', formula_txt)), reached_domain,
              family=quasipoisson(link='log'), na.action=na.omit)
  b <- coef(m)['income_rank']; se <- sqrt(diag(vcov(m)))['income_rank']; df <- m$df.residual; cr <- qt(.975, df=df)
  data.frame(model=nm, formula=formula_txt, rank_scope='rank_reached',
             n_domain=sum(m$prior.weights > 0), events=sum(m$y[m$prior.weights > 0] == 1), degf=df,
             PR=exp(b), lower=exp(b-cr*se), upper=exp(b+cr*se), p=2*pt(abs(b/se), df=df, lower.tail=FALSE),
             max_prev_prevista=max(m$fitted.values[m$prior.weights > 0]))
}))
write.csv(reached_pr, file.path(out_dir,'specialist_reached_PR_income.csv'), row.names=FALSE)
reached_income <- svyby(~specialist_reached_all, ~income3, subset(reached_domain, !is.na(income3)), svymean,
                        na.rm=TRUE, vartype=c('se','ci'), keep.names=FALSE)
names(reached_income)[2] <- 'estimate'
write.csv(reached_income, file.path(out_dir,'specialist_reached_by_income.csv'), row.names=FALSE)

# Motivos de ausência de regularidade, entre diagnosticados sem cuidado regular.
bar_domain <- subset(des, depression_dx==1 & Q094 %in% c(2,3) & !is.na(poor_access_reason))
bar <- svytable(~poor_access_reason, bar_domain)
bar_tab <- data.frame(reason=names(bar), weighted_n=as.numeric(bar))
bar_tab$share <- bar_tab$weighted_n/sum(bar_tab$weighted_n)
write.csv(bar_tab, file.path(out_dir,'barriers_no_regular_care.csv'), row.names=FALSE)

reason_domain <- use_rank(subset(des, depression_dx==1 & Q094 %in% c(2,3) & !is.na(reason_group) & !is.na(income_rank)),
                          'rank_reason')                                                          # [16b]
reason_phq <- svyby(~phq10, ~reason_group, subset(reason_domain, !is.na(phq10)), svymean,
                    na.rm=TRUE, vartype=c('se','ci'), keep.names=FALSE)
names(reason_phq)[2] <- 'estimate'
reason_phq$n_domain <- sum(pns$depression_dx==1 & pns$Q094 %in% c(2,3) & !is.na(pns$reason_group) & !is.na(pns$phq10))
write.csv(reason_phq, file.path(out_dir,'reason_group_by_phq10.csv'), row.names=FALSE)

severity_vars <- c('phq10','phq15','self_harm_thoughts','self_harm_frequent','self_harm_near_daily','high_depression_limitation')
severity_rows <- lapply(severity_vars, function(v) {
  dd <- subset(reason_domain, !is.na(get(v)))
  z <- svyby(as.formula(paste0('~',v)), ~reason_group, dd, svymean,
             na.rm=TRUE, vartype=c('se','ci'), keep.names=FALSE)
  names(z)[2] <- 'estimate'
  z$indicator <- v
  z$n_domain <- sum(pns$depression_dx==1 & pns$Q094 %in% c(2,3) & !is.na(pns$reason_group) & !is.na(pns[[v]]))
  z
})
severity_tab <- dplyr::bind_rows(severity_rows)
write.csv(severity_tab, file.path(out_dir,'reason_group_severity.csv'), row.names=FALSE)

reason_given_symptoms <- lapply(c('phq10','phq15','self_harm_thoughts'), function(v) {
  dd <- subset(reason_domain, get(v)==1)
  tab <- svytable(~reason_group, dd)
  out <- data.frame(indicator=v, reason_group=names(tab), weighted_n=as.numeric(tab))
  out$share <- out$weighted_n/sum(out$weighted_n)
  out$n_domain <- sum(pns$depression_dx==1 & pns$Q094 %in% c(2,3) & !is.na(pns$reason_group) & pns[[v]]==1, na.rm=TRUE)
  out
})
write.csv(dplyr::bind_rows(reason_given_symptoms), file.path(out_dir,'reason_distribution_among_symptomatic.csv'), row.names=FALSE)

remission_domain <- use_rank(subset(reason_domain, reason_group=='remissao_percebida' & !is.na(phq10)),
                             'rank_remission')                                                    # [16b]
remission_overall <- svymean(~phq10, remission_domain, na.rm=TRUE)
remission_est <- data.frame(estimate=as.numeric(coef(remission_overall)[1]), SE=as.numeric(SE(remission_overall)[1]), degf=degf(remission_domain))
remission_est$lower <- remission_est$estimate-qt(.975,df=remission_est$degf)*remission_est$SE
remission_est$upper <- remission_est$estimate+qt(.975,df=remission_est$degf)*remission_est$SE
remission_est$n_domain <- sum(pns$depression_dx==1 & pns$Q094 %in% c(2,3) & pns$reason_group=='remissao_percebida' & !is.na(pns$phq10), na.rm=TRUE)
write.csv(remission_est, file.path(out_dir,'remission_symptom_discordance_overall.csv'), row.names=FALSE)
remission_income <- svyby(~phq10, ~income3, subset(remission_domain, !is.na(income3)), svymean,
                           na.rm=TRUE, vartype=c('se','ci'), keep.names=FALSE)
names(remission_income)[2] <- 'estimate'
remission_income$n_domain <- sum(pns$depression_dx==1 & pns$Q094 %in% c(2,3) & pns$reason_group=='remissao_percebida' & !is.na(pns$phq10) & !is.na(pns$income3), na.rm=TRUE)
write.csv(remission_income, file.path(out_dir,'remission_symptom_discordance_by_income.csv'), row.names=FALSE)
remission_models <- list(minimum='income_rank + age_c + sex + region',
                         expanded='income_rank + age_c + sex + region + race3 + education3 + urban + plan')
remission_pr <- do.call(rbind, lapply(names(remission_models), function(nm) {
  formula_txt <- remission_models[[nm]]
  m <- svyglm(as.formula(paste0('phq10 ~ ',formula_txt)), remission_domain,
              family=quasipoisson(link='log'), na.action=na.omit)
  b <- coef(m)['income_rank']; se <- sqrt(diag(vcov(m)))['income_rank']; df <- m$df.residual; cr <- qt(.975,df=df)
  data.frame(model=nm, formula=formula_txt, rank_scope='rank_remission',
             n_domain=sum(m$prior.weights > 0),  # [16b] n efetivo do modelo
             degf=df, PR=exp(b), lower=exp(b-cr*se), upper=exp(b+cr*se),
             p=2*pt(abs(b/se),df=df,lower.tail=FALSE))
}))
write.csv(remission_pr, file.path(out_dir,'remission_symptom_discordance_PR_income.csv'), row.names=FALSE)

reason_income <- svyby(~structural_reason, ~income3,
                        subset(reason_domain, !is.na(structural_reason) & !is.na(income3)),
                        svymean, na.rm=TRUE, vartype=c('se','ci'), keep.names=FALSE)
names(reason_income)[2] <- 'estimate'
reason_income$n_domain <- sum(!is.na(reason_domain$variables$structural_reason) & !is.na(reason_domain$variables$income3) & reason_domain$prob < Inf)  # [16b] só o domínio
write.csv(reason_income, file.path(out_dir,'structural_barrier_by_income.csv'), row.names=FALSE)

reason_models <- list(minimum='income_rank + age_c + sex + region',
                      expanded='income_rank + age_c + sex + region + race3 + education3 + urban + plan')
reason_pr <- do.call(rbind, lapply(names(reason_models), function(nm) {
  formula_txt <- reason_models[[nm]]
  m <- svyglm(as.formula(paste0('structural_reason ~ ',formula_txt)), reason_domain,
              family=quasipoisson(link='log'), na.action=na.omit)
  b <- coef(m)['income_rank']; se <- sqrt(diag(vcov(m)))['income_rank']; df <- m$df.residual; cr <- qt(.975,df=df)
  data.frame(model=nm, formula=formula_txt, rank_scope='rank_reason',
             n_domain=sum(m$prior.weights > 0),  # [16b] n efetivo do modelo
             degf=df, PR=exp(b), lower=exp(b-cr*se), upper=exp(b+cr*se),
             p=2*pt(abs(b/se),df=df,lower.tail=FALSE))
}))
write.csv(reason_pr, file.path(out_dir,'structural_barrier_PR_income.csv'), row.names=FALSE)

# Perfil descritivo do encaminhamento.
ref_profile <- data.frame(
  component=c('Encaminhados a especialista entre consultas recentes',
              'Chegaram a alguma consulta especializada entre consultas recentes',
              'Falharam completamente após encaminhamento'),
  variable=c('referral_to_specialist','specialist_reached_all','referral_failure'),
  denominator=c('all_recent','all_recent','referred')
)
ref_rows <- lapply(seq_len(nrow(ref_profile)), function(i) {
  v <- ref_profile$variable[i]
  dd <- if (ref_profile$denominator[i]=='all_recent') ref_all_domain else ref_domain
  dd <- subset(dd, !is.na(get(v)))
  z <- svymean(as.formula(paste0('~',v)), dd, na.rm=TRUE)
  est <- as.numeric(coef(z)[1]); se <- as.numeric(SE(z)[1]); cr <- qt(.975,df=degf(dd))
  data.frame(component=ref_profile$component[i], variable=v, denominator=ref_profile$denominator[i],
             n_domain=if (ref_profile$denominator[i]=='all_recent') sum(pns$depression_dx==1 & pns$Q10101 %in% 1:4 & !is.na(pns$Q106) & !is.na(pns[[v]]), na.rm=TRUE) else sum(pns$depression_dx==1 & pns$Q10101 %in% 1:4 & pns$Q106==1 & !is.na(pns$referral_failure), na.rm=TRUE), estimate=est, SE=se,
             lower=est-cr*se, upper=est+cr*se, degf=degf(dd))
})
write.csv(dplyr::bind_rows(ref_rows), file.path(out_dir,'referral_profile.csv'), row.names=FALSE)

# Auditoria simples de tamanhos e valores brutos.
write.csv(as.data.frame(table(pns$Q092, useNA='ifany')), file.path(out_dir,'audit_Q092.csv'), row.names=FALSE)
write.csv(as.data.frame(table(pns$phq9_total, useNA='ifany')), file.path(out_dir,'audit_phq9_total.csv'), row.names=FALSE)
write.csv(as.data.frame(table(pns$Q106, useNA='ifany')), file.path(out_dir,'audit_Q106.csv'), row.names=FALSE)
write.csv(as.data.frame(table(pns$Q10701, useNA='ifany')), file.path(out_dir,'audit_Q10701.csv'), row.names=FALSE)

capture.output(sessionInfo(), file=file.path(out_dir,'sessionInfo.txt'))
cat('Depression cascade (domain-rank) analysis complete\n')
cat('N selected:', nrow(pns), '\n')
cat('PHQ9 >=10:', sum(pns$phq10==1, na.rm=TRUE), '\n')
cat('Depression diagnosis:', sum(pns$depression_dx==1, na.rm=TRUE), '\n')
cat('Referral-failure domain:', sum(pns$depression_dx==1 & pns$Q10101 %in% 1:4 & pns$Q106==1 & !is.na(pns$referral_failure) & !is.na(pns$income_reais), na.rm=TRUE), '\n')
cat('Critical gap domain:', sum(!is.na(pns$rank_crit)), '\n')
cat('\nChegou ao especialista, por faixa de renda:\n'); print(reached_income, digits=3)
cat('\nModelos de acesso ao especialista (RII com rank do domínio):\n')
print(rbind(cbind(desfecho='chegou ao especialista',     reached_pr[, c('model','n_domain','PR','lower','upper','p')]),
            cbind(desfecho='falha após encaminhamento',  ref_pr[, c('model','n_domain','PR','lower','upper','p')])), digits=3)
cat('\nCom tratamento atual entre diagnosticados (era o "59,4% sem tratamento"):\n')
print(stage_tab[stage_tab$variable == 'depression_current_treatment', c('stage','n_domain','estimate','lower','upper')])
