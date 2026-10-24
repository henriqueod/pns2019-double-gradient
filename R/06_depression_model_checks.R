# =============================================================================
# 06_depression_model_checks.R
# Prediction-range audit and log-binomial refits for the 38 exploratory depression models.
# Repository name: R/06_depression_model_checks.R   |   Development name: 32_checagem_previsoes_depressao.R
# Code comments are in Portuguese. Cross-references in the comments use the
# development names; README.md maps them to the repository names.
# Paths: all input and output paths are relative to `root`, which defaults to
# the working directory (run from the repository root) or to the PNS_ROOT
# environment variable if it is set.
# =============================================================================

# 32_checagem_previsoes_depressao.R
# Checagem de previsões > 1 nos modelos de Poisson modificado da cascata de depressão (16b)
# e da lacuna de detecção (18b), com reajuste por log-binomial.
#
# POR QUÊ
#   No 16b e no 18b, alguns modelos expandidos tiveram probabilidade prevista máxima acima de 1
#   (chegou ao especialista: 1,04; sem diagnóstico entre PHQ-9 >= 10, domínio D1: 1,02).
#   O Poisson modificado pode prever acima de 1 sem que isso vicie o RII, mas é preciso contar
#   quantas previsões passam de 1, ver quem são e conferir o RII com um log-binomial.
#   Mesma lógica da seção 9 do 30 (script 30b) para os co-primários.
#
# O QUE FAZ, para cada indicador, modelo mínimo e expandido, rank do domínio:
#   - n do modelo, eventos, nº e % de previsões > 1, previsão máxima
#   - RII do Poisson (igual ao 16b/18b) e RII do log-binomial com início válido
#     (coeficientes do Poisson com o intercepto deslocado até todas as previsões < 0,99)
#   - perfil das previsões > 1 (quintil do rank do domínio x plano x escolaridade)
#   Detecção: também com o rank populacional (como no 18b), só para comparação.
#
# Variáveis e domínios: idênticos ao 16b_depression_cascade_domain_rank.R e ao
# 18b_detection_gap_domain_rank.R (a preparação abaixo é copiada do 16b, linhas 37-161).
#
# Saída: results/depression_prediction_check/
#   32_checagem_previsoes_depressao.csv  (uma linha por indicador x modelo x rank)
#   32_perfil_previsoes_acima_de_1.csv   (só se houver previsões > 1)

options(timeout = 600)
options(survey.lonely.psu = 'adjust')
options(survey.adjust.domain.lonely = TRUE)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(PNSIBGE)
  library(survey)
  library(dplyr)
})
cat('06_depression_model_checks.R (32_checagem_previsoes_depressao.R) | versão 27/09/2026\n')

root <- Sys.getenv('PNS_ROOT', unset = getwd())
data_dir <- file.path(root, 'data_pnsibge')
out_dir <- file.path(root, 'results', 'depression_prediction_check')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Preparação idêntica ao 16b (linhas 37-161) -------------------------------------
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

# ---- Variáveis e ranks da lacuna de detecção (idênticos ao 18b) ----------------------
pns$oth_dx <- ifelse(pns$Q11006 == 1, 1, ifelse(pns$Q11006 == 2, 0, NA_real_))
pns$undiagnosed <- ifelse(pns$depression_dx %in% 0, 1, ifelse(pns$depression_dx %in% 1, 0, NA_real_))
pns$undiagnosed_any <- ifelse(pns$depression_dx %in% 1 | pns$oth_dx %in% 1, 0,
                        ifelse(pns$depression_dx %in% 0 & pns$oth_dx %in% 0, 1, NA_real_))
pns$dx_regular <- as.numeric(pns$depression_dx %in% 1 & pns$Q094 %in% 1)
dominios_det <- list(
  D1_todos_phq10          = pns$phq10 %in% 1,
  D2_phq10_sem_dx_regular = pns$phq10 %in% 1 & pns$dx_regular == 0,
  D3_phq15                = pns$phq15 %in% 1
)
for (dn in names(dominios_det)) for (y in c('undiagnosed', 'undiagnosed_any')) {
  pns[[paste0('dom_', dn, '_', y)]]  <- as.numeric(dominios_det[[dn]] & !is.na(pns[[y]]))
  pns[[paste0('rank_', dn, '_', y)]] <- rank_in(dominios_det[[dn]] & !is.na(pns[[y]]))
}

des <- PNSIBGE::pns_design(pns)

# ---------------------------------------------------------------------------
# Especificação dos modelos (domínio, desfecho, rank, covariáveis extras no expandido)
# ---------------------------------------------------------------------------
minimo    <- 'income_rank + age_c + sex + region'
expandido <- 'income_rank + age_c + sex + region + race3 + education3 + urban + plan'
spec <- list(
  list('remissao_percebida_PHQ10_sem_tratamento', 'hidden_remission_gap',      '!is.na(hidden_remission_gap)',      'rank_hidden_remission_gap', ''),
  list('remissao_percebida_PHQ15_sem_tratamento', 'hidden_remission_severe',   '!is.na(hidden_remission_severe)',   'rank_hidden_remission_severe', ''),
  list('remissao_percebida_ideacao_sem_tratamento','hidden_remission_selfharm','!is.na(hidden_remission_selfharm)', 'rank_hidden_remission_selfharm', ''),
  list('idem_PHQ10_todos_diagnosticados',  'hidden_remission_gap_all_dx',      '!is.na(hidden_remission_gap_all_dx)',      'rank_hidden_remission_gap_all_dx', ''),
  list('idem_PHQ15_todos_diagnosticados',  'hidden_remission_severe_all_dx',   '!is.na(hidden_remission_severe_all_dx)',   'rank_hidden_remission_severe_all_dx', ''),
  list('idem_ideacao_todos_diagnosticados','hidden_remission_selfharm_all_dx', '!is.na(hidden_remission_selfharm_all_dx)', 'rank_hidden_remission_selfharm_all_dx', ''),
  list('lacuna_critica', 'critical_depression_gap',
       'phq10 == 1 & self_harm_thoughts == 1 & depression_dx == 1 & !is.na(critical_depression_gap)', 'rank_crit', ''),
  list('falha_apos_encaminhamento', 'referral_failure',
       'depression_dx == 1 & Q10101 %in% 1:4 & Q106 == 1 & !is.na(referral_failure)', 'rank_ref', ''),  # [13/10/2026] bloco D, como as outras linhas; antes: ' + public_last_depression_visit'
  list('chegou_ao_especialista', 'specialist_reached_all',
       'depression_dx == 1 & Q10101 %in% 1:4 & !is.na(specialist_reached_all)', 'rank_reached', ''),
  list('barreira_estrutural_motivo', 'structural_reason',
       'depression_dx == 1 & Q094 %in% c(2,3) & !is.na(reason_group)', 'rank_reason', ''),
  list('remissao_percebida_com_PHQ10', 'phq10',
       "depression_dx == 1 & Q094 %in% c(2,3) & reason_group %in% 'remissao_percebida' & !is.na(phq10)", 'rank_remission', '')
)
for (dn in names(dominios_det)) for (y in c('undiagnosed', 'undiagnosed_any')) {
  spec[[length(spec) + 1]] <- list(paste0('deteccao_', dn, '_', y), y, paste0('dom_', dn, '_', y, ' == 1'),
                                   paste0('rank_', dn, '_', y), '')
  if (dn == 'D1_todos_phq10')
    spec[[length(spec) + 1]] <- list(paste0('deteccao_', dn, '_', y, '_rank_populacional'), y,
                                     paste0('dom_', dn, '_', y, ' == 1'), 'income_rank_pop', '')
}

# ---------------------------------------------------------------------------
# Checagem
# ---------------------------------------------------------------------------
fitm <- function(dd, y, rhs, family, start = NULL, ...)
  tryCatch(suppressWarnings(svyglm(as.formula(paste(y, '~', rhs)), design = dd, family = family,
                                   start = start, ...)), error = function(e) e)
rii <- function(m) {
  if (inherits(m, 'error') || !('income_rank' %in% names(coef(m))))
    return(c(est = NA, lower = NA, upper = NA, p = NA))
  b <- coef(m)[['income_rank']]; se <- sqrt(diag(vcov(m)))[['income_rank']]
  df <- m$df.residual; cr <- qt(.975, df)
  c(est = exp(b), lower = exp(b - cr * se), upper = exp(b + cr * se), p = 2 * pt(abs(b / se), df, lower.tail = FALSE))
}

linhas <- list(); perfis <- list()
for (s in spec) for (mod in c('minimo', 'expandido')) {
  lab <- s[[1]]; y <- s[[2]]; dom <- s[[3]]; rk <- s[[4]]; extra <- s[[5]]
  rhs <- if (mod == 'minimo') minimo else paste0(expandido, extra)
  cat('  ', lab, '|', mod, '\n'); flush.console()
  dd <- tryCatch(use_rank(subset(des, eval(parse(text = dom)) & !is.na(get(rk)) & !is.na(income_rank_pop)), rk),
                 error = function(e) e)
  if (inherits(dd, 'error')) { cat('     erro no domínio:', conditionMessage(dd), '\n'); next }
  mp <- fitm(dd, y, rhs, quasipoisson(link = 'log'))
  if (inherits(mp, 'error')) {
    linhas[[length(linhas) + 1]] <- data.frame(indicador = lab, modelo = mod, rank = rk, erro = conditionMessage(mp))
    next
  }
  keep <- mp$prior.weights > 0
  f <- mp$fitted.values[keep]; acima <- f > 1
  st <- coef(mp); eta_max <- max(mp$linear.predictors[keep])
  if (eta_max > log(0.99)) st[['(Intercept)']] <- st[['(Intercept)']] - (eta_max - log(0.99))
  mb <- fitm(dd, y, rhs, quasibinomial(link = 'log'), start = st, control = glm.control(maxit = 200))
  rp <- rii(mp); rb <- rii(mb)
  linhas[[length(linhas) + 1]] <- data.frame(
    indicador = lab, modelo = mod, rank = rk, n = sum(keep), eventos = sum(mp$y[keep] == 1),
    n_previsto_acima_1 = sum(acima), pct_previsto_acima_1 = round(100 * mean(acima), 2), max_prev_poisson = max(f),
    RII_poisson = rp[['est']], lower_poisson = rp[['lower']], upper_poisson = rp[['upper']], p_poisson = rp[['p']],
    RII_logbin = rb[['est']], lower_logbin = rb[['lower']], upper_logbin = rb[['upper']], p_logbin = rb[['p']],
    logbin_convergiu = if (inherits(mb, 'error')) FALSE else isTRUE(mb$converged),
    max_prev_logbin = if (inherits(mb, 'error')) NA else max(mb$fitted.values[mb$prior.weights > 0]),
    erro = if (inherits(mb, 'error')) paste('log-binomial:', conditionMessage(mb)) else NA_character_)
  if (any(acima)) {
    rows <- match(rownames(model.frame(mp)), rownames(dd$variables))[keep][acima]
    q <- cut(dd$variables$income_rank[rows], breaks = c(0, .2, .4, .6, .8, 1), labels = paste0('Q', 1:5), include.lowest = TRUE)
    tb <- as.data.frame(table(quintil_rank_dominio = q, plano = dd$variables$plan[rows],
                              escolaridade = dd$variables$education3[rows], useNA = 'ifany'))
    perfis[[length(perfis) + 1]] <- cbind(indicador = lab, modelo = mod, rank = rk, tb[tb$Freq > 0, ])
  }
}
res <- bind_rows(linhas)
write.csv(res, file.path(out_dir, '32_checagem_previsoes_depressao.csv'), row.names = FALSE)
if (length(perfis)) write.csv(bind_rows(perfis), file.path(out_dir, '32_perfil_previsoes_acima_de_1.csv'), row.names = FALSE)

cat('\nResumo (modelos com alguma previsão > 1):\n')
cols <- intersect(c('indicador','modelo','n','n_previsto_acima_1','pct_previsto_acima_1','max_prev_poisson',
                    'RII_poisson','lower_poisson','upper_poisson','RII_logbin','lower_logbin','upper_logbin',
                    'logbin_convergiu','max_prev_logbin'), names(res))
print(subset(res, !is.na(n_previsto_acima_1) & n_previsto_acima_1 > 0)[, cols], digits = 3)
cat('\nModelos checados:', nrow(res), '| com previsões > 1:', sum(res$n_previsto_acima_1 > 0, na.rm = TRUE),
    '| log-binomial sem convergir:', sum(!res$logbin_convergiu, na.rm = TRUE), '\n')
writeLines(capture.output(sessionInfo()), file.path(out_dir, 'sessionInfo.txt'))
cat('Fim. Arquivos em', out_dir, '\n')
