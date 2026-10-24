# =============================================================================
# 08_referral_failure_blockD.R
# Exploratory depression cascade: referral failure re-estimated with the block-D covariates used in the other rows of Table 6 (minimal, block D, and block D plus public-sector care as a sensitivity model).
# Repository name: R/08_referral_failure_blockD.R   |   Development name: 33_falha_encaminhamento_blocoD.R
# Code comments are in Portuguese. Cross-references in the comments use the
# development names; README.md maps them to the repository names.
# Paths: all input and output paths are relative to `root`, which defaults to
# the working directory (run from the repository root) or to the PNS_ROOT
# environment variable if it is set.
# =============================================================================

# 33_falha_encaminhamento_blocoD.R
# Refaz SÓ o modelo ajustado da falha após encaminhamento (Table 6, linha "No specialist visit
# after referral"), agora com as mesmas covariáveis das outras linhas (bloco D):
#   rank + idade + sexo + região + raça/cor + escolaridade + urbano + plano de saúde.
# O modelo antigo, que também incluía se o último atendimento por depressão foi pelo SUS
# (public_last_depression_visit, Q105), passa a ser análise de sensibilidade.
#
# POR QUÊ (decisão de 13/10/2026)
#   Era a única linha da Table 6 com uma covariável a mais, não descrita nos Métodos, e o uso
#   do SUS pode ser mediador (renda -> setor do atendimento -> falha), não confundidor.
#
# O QUE FAZ, no domínio dos 2.020 encaminhados (diagnosticados, atendimento médico por
# depressão há menos de 3 anos, Q106 = 1, Q107a válida), com o rank do domínio (como no 16b):
#   M0 mínimo            rank + idade + sexo + região                (deve repetir RII 0,80)
#   M1 bloco D           + raça/cor + escolaridade + urbano + plano  (NOVO principal)
#   M2 bloco D + SUS     M1 + public_last_depression_visit           (deve repetir RII 1,27)
#   Para cada um: n, eventos, RII Poisson (IC, p, gl), previsões > 1, máximo previsto, e o
#   RII log-binomial com início válido (mesma rotina do 32), para a Tabela S10.
#   Também conta as respostas da Q105 no domínio (quantos "não sabe").
#
# Variáveis e domínio: idênticos ao 16b e ao 32 (a preparação abaixo é copiada do 32).
# Saída: results/referral_failure_blockD/33_falha_encaminhamento_modelos.csv
#        results/referral_failure_blockD/33_Q105_no_dominio.csv
#        results/referral_failure_blockD/sessionInfo.txt
options(timeout = 600)
options(survey.lonely.psu = 'adjust')
options(survey.adjust.domain.lonely = TRUE)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(PNSIBGE)
  library(survey)
  library(dplyr)
})
cat('33_falha_encaminhamento_blocoD.R | versão 13/10/2026\n')

root <- Sys.getenv('PNS_ROOT', unset = getwd())
data_dir <- file.path(root, 'data_pnsibge')
out_dir <- file.path(root, 'results', 'referral_failure_blockD')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Preparação idêntica ao 32 (linhas 42-167), que repete o 16b ----------------------
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

des <- PNSIBGE::pns_design(pns)

# ---------------------------------------------------------------------------
# Modelos da falha após encaminhamento
# ---------------------------------------------------------------------------
dom_ref <- 'depression_dx == 1 & Q10101 %in% 1:4 & Q106 == 1 & !is.na(referral_failure)'
modelos <- list(
  M0_minimo      = 'income_rank + age_c + sex + region',
  M1_blocoD      = 'income_rank + age_c + sex + region + race3 + education3 + urban + plan',
  M2_blocoD_SUS  = 'income_rank + age_c + sex + region + race3 + education3 + urban + plan + public_last_depression_visit'
)

fitm <- function(dd, y, rhs, family, start = NULL, ...)
  tryCatch(suppressWarnings(svyglm(as.formula(paste(y, '~', rhs)), design = dd, family = family,
                                   start = start, ...)), error = function(e) e)
rii <- function(m) {
  if (inherits(m, 'error') || !('income_rank' %in% names(coef(m))))
    return(c(est = NA, lower = NA, upper = NA, p = NA, df = NA))
  b <- coef(m)[['income_rank']]; se <- sqrt(diag(vcov(m)))[['income_rank']]
  df <- m$df.residual; cr <- qt(.975, df)
  c(est = exp(b), lower = exp(b - cr * se), upper = exp(b + cr * se), p = 2 * pt(abs(b / se), df, lower.tail = FALSE), df = df)
}

dd <- use_rank(subset(des, eval(parse(text = dom_ref)) & !is.na(rank_ref) & !is.na(income_rank_pop)), 'rank_ref')
cat('Domínio da falha após encaminhamento: n =', sum(is.finite(dd$prob)), '(esperado 2.020)\n')

linhas <- list()
for (nm in names(modelos)) {
  rhs <- modelos[[nm]]
  cat('  ', nm, '\n'); flush.console()
  mp <- fitm(dd, 'referral_failure', rhs, quasipoisson(link = 'log'))
  if (inherits(mp, 'error')) {
    linhas[[nm]] <- data.frame(modelo = nm, formula = rhs, erro = conditionMessage(mp)); next
  }
  keep <- mp$prior.weights > 0
  f <- mp$fitted.values[keep]
  st <- coef(mp); eta_max <- max(mp$linear.predictors[keep])
  if (eta_max > log(0.99)) st[['(Intercept)']] <- st[['(Intercept)']] - (eta_max - log(0.99))
  mb <- fitm(dd, 'referral_failure', rhs, quasibinomial(link = 'log'), start = st, control = glm.control(maxit = 200))
  rp <- rii(mp); rb <- rii(mb)
  linhas[[nm]] <- data.frame(
    modelo = nm, formula = rhs, n = sum(keep), eventos = sum(mp$y[keep] == 1),
    RII = rp[['est']], lower = rp[['lower']], upper = rp[['upper']], p = rp[['p']], gl_residual = rp[['df']],
    n_previsto_acima_1 = sum(f > 1), max_prev_poisson = max(f),
    RII_logbin = rb[['est']], lower_logbin = rb[['lower']], upper_logbin = rb[['upper']], p_logbin = rb[['p']],
    logbin_convergiu = if (inherits(mb, 'error')) FALSE else isTRUE(mb$converged),
    max_prev_logbin = if (inherits(mb, 'error')) NA else max(mb$fitted.values[mb$prior.weights > 0]),
    erro = if (inherits(mb, 'error')) paste('log-binomial:', conditionMessage(mb)) else NA_character_)
}
res <- bind_rows(linhas)
write.csv(res, file.path(out_dir, '33_falha_encaminhamento_modelos.csv'), row.names = FALSE)

# Respostas da Q105 ("Esse atendimento foi feito pelo SUS?") no domínio: 1 Sim, 2 Não, 3 Não sabe.
q105 <- as.data.frame(table(Q105 = dd$variables$Q105[is.finite(dd$prob)], useNA = 'ifany'))
write.csv(q105, file.path(out_dir, '33_Q105_no_dominio.csv'), row.names = FALSE)

cat('\n=== Falha após encaminhamento: RII por modelo ===\n')
print(res[, intersect(c('modelo','n','eventos','RII','lower','upper','p','n_previsto_acima_1',
                        'max_prev_poisson','RII_logbin','lower_logbin','upper_logbin','logbin_convergiu'),
                      names(res))], digits = 3, row.names = FALSE)
cat('\nConferência: M0 deve dar RII 0,80 (0,41-1,57) e M2 deve dar 1,27 (0,45-3,60), como no 16b.\n')
cat('\nQ105 no domínio:\n'); print(q105, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, 'sessionInfo.txt'))
cat('Fim. Arquivos em', out_dir, '\n')
