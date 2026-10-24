# =============================================================================
# 05_detection_gap.R
# Exploratory detection gap among PHQ-9 screen-positive respondents.
# Repository name: R/05_detection_gap.R   |   Development name: 18b_detection_gap_domain_rank.R
# Code comments are in Portuguese. Cross-references in the comments use the
# development names; README.md maps them to the repository names.
# Paths: all input and output paths are relative to `root`, which defaults to
# the working directory (run from the repository root) or to the PNS_ROOT
# environment variable if it is set.
# =============================================================================

# 18b_detection_gap_domain_rank.R
# Lacuna de detecção: entre quem rastreia positivo no PHQ-9 (>= 10), quem nunca recebeu
# diagnóstico profissional de depressão, e o gradiente de renda dessa lacuna.
#
# Decisões de 24-25/09:
#   - RANK DE RENDA NO DOMÍNIO: o RII principal usa um ridit ponderado calculado só entre as
#     pessoas de cada domínio de análise. O ridit sobre os 90.846 moradores selecionados
#     é reestimado só para comparação.
#   - PHQ-9 pontuado como no 16/16b: itens N010-N018, respostas 1-4 -> 0-3, só com os 9 itens válidos.
#
# DOMÍNIOS (todos com PHQ-9 >= 10, independentemente de diagnóstico)
#   D1  todos os rastreados positivos                                   (principal)
#   D2  rastreados positivos, excluindo quem tem diagnóstico de depressão E visita
#       regularmente por causa dela (Q094 = 1)
#   D3  PHQ-9 >= 15 (sintomas mais graves)                               (sensibilidade)
#
# DESFECHOS
#   undiagnosed      1 = nunca recebeu diagnóstico de depressão (Q092 = 2); 0 = recebeu (Q092 = 1)
#   undiagnosed_any  1 = sem diagnóstico de depressão E sem diagnóstico de outra doença mental
#                    (Q092 = 2 e Q11006 = 2), isto é, fora da coorte "qualquer doença mental";
#                    0 = com qualquer um dos dois diagnósticos
#
# MODELOS: Poisson modificado (quasipoisson, link log) do rank de renda
#   mínimo:    rank + idade + sexo + região
#   expandido: + raça/cor + escolaridade + urbano + plano de saúde

options(timeout = 600)
options(survey.lonely.psu = 'adjust')
options(survey.adjust.domain.lonely = TRUE)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(PNSIBGE)
  library(survey)
  library(dplyr)
})
cat('05_detection_gap.R (18b_detection_gap_domain_rank.R) | versão 26/09/2026\n')

root     <- Sys.getenv('PNS_ROOT', unset = getwd())
data_dir <- file.path(root, 'data_pnsibge')
out_dir  <- file.path(root, 'results', 'detection_gap_domain_rank')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wcsv <- function(x, name) write.csv(x, file.path(out_dir, name), row.names = FALSE)

vars <- c(
  'V0001','V0024','UPA_PNS','V0006_PNS','V0015','V0025A','V0029','V00291','V00292','V00293',
  'C006','C008','C009','V0026','VDD004A','VDF003','VDF004','M00203','I00102',
  'Q092','Q094','Q11006',
  'N010','N011','N012','N013','N014','N015','N016','N017','N018'
)
pns <- PNSIBGE::get_pns(year = 2019, selected = TRUE, vars = vars, labels = FALSE,
                        deflator = FALSE, design = FALSE, reload = FALSE, savedir = data_dir)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
for (v in setdiff(vars, c('UPA_PNS','V0006_PNS','V0015','V0025A','V0029','V00291','V00292','V00293'))) {
  if (v %in% names(pns)) pns[[v]] <- num(pns[[v]])
}

# ---------------------------------------------------------------------------
# Covariáveis (idênticas ao 01c / 16b)
# ---------------------------------------------------------------------------
region_exact <- c(`11`='Norte',`12`='Norte',`13`='Norte',`14`='Norte',`15`='Norte',`16`='Norte',`17`='Norte',
                  `21`='Nordeste',`22`='Nordeste',`23`='Nordeste',`24`='Nordeste',`25`='Nordeste',`26`='Nordeste',`27`='Nordeste',`28`='Nordeste',`29`='Nordeste',
                  `31`='Sudeste',`32`='Sudeste',`33`='Sudeste',`35`='Sudeste',
                  `41`='Sul',`42`='Sul',`43`='Sul',
                  `50`='Centro-Oeste',`51`='Centro-Oeste',`52`='Centro-Oeste',`53`='Centro-Oeste')
pns$region <- factor(unname(region_exact[as.character(pns$V0001)]),
                     levels = c('Norte','Nordeste','Sudeste','Sul','Centro-Oeste'))
pns$sex <- factor(pns$C006, levels = c(1,2), labels = c('Masculino','Feminino'))
pns$age <- pns$C008
pns$age_c <- pns$age - 45
pns$race3 <- factor(ifelse(pns$C009 == 1, 'Branca', ifelse(pns$C009 %in% 2:5, 'Não_branca', NA_character_)),
                    levels = c('Branca','Não_branca'))
pns$urban <- factor(pns$V0026, levels = c(1,2), labels = c('Urbana','Rural'))
pns$education3 <- cut(pns$VDD004A, breaks = c(0,3,5,7), labels = c('Baixa','Média','Alta'), right = TRUE)
pns$income3 <- cut(pns$VDF004, breaks = c(0,3,4,7), labels = c('Até_1_SM','Mais_1_a_2_SM','Mais_2_SM'), right = TRUE)
pns$income3 <- factor(pns$income3, levels = c('Até_1_SM','Mais_1_a_2_SM','Mais_2_SM'))
pns$income_reais <- pns$VDF003
pns$plan <- ifelse(pns$I00102 %in% c(1,2), ifelse(pns$I00102 == 1, 1, 0), NA_real_)

make_weighted_rank <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0 & x >= 0
  r <- rep(NA_real_, length(x))
  ord <- order(x[ok], na.last = NA)
  idx <- which(ok)[ord]
  ww <- w[idx]
  r[idx] <- (cumsum(ww) - 0.5 * ww) / sum(ww)
  r
}
rank_in <- function(mask) {
  r <- rep(NA_real_, nrow(pns))
  i <- which(mask %in% TRUE & !is.na(pns$income_reais) & !is.na(pns$V00291))
  r[i] <- make_weighted_rank(pns$income_reais[i], pns$V00291[i])
  r
}
pns$income_rank_pop <- make_weighted_rank(pns$income_reais, pns$V00291)   # 90.846 moradores selecionados

# ---------------------------------------------------------------------------
# PHQ-9, diagnósticos e desfechos
# ---------------------------------------------------------------------------
phq_vars <- paste0('N0', 10:18)
phq_items <- as.data.frame(lapply(pns[, phq_vars], function(x) ifelse(x %in% 1:4, x - 1, NA_real_)))
pns$phq9_total <- ifelse(rowSums(!is.na(phq_items)) == 9, rowSums(phq_items), NA_real_)
pns$phq10 <- ifelse(is.na(pns$phq9_total), NA_real_, as.numeric(pns$phq9_total >= 10))
pns$phq15 <- ifelse(is.na(pns$phq9_total), NA_real_, as.numeric(pns$phq9_total >= 15))

pns <- pns %>%
  mutate(
    dep_dx = case_when(Q092 == 1 ~ 1, Q092 == 2 ~ 0, TRUE ~ NA_real_),
    oth_dx = case_when(Q11006 == 1 ~ 1, Q11006 == 2 ~ 0, TRUE ~ NA_real_),
    undiagnosed = case_when(dep_dx == 0 ~ 1, dep_dx == 1 ~ 0, TRUE ~ NA_real_),
    # Fora da coorte "qualquer doença mental": precisa de Não válido nas duas perguntas.
    undiagnosed_any = case_when(dep_dx == 1 | oth_dx == 1 ~ 0,
                                dep_dx == 0 & oth_dx == 0 ~ 1,
                                TRUE ~ NA_real_),
    # Diagnosticado com depressão e em seguimento regular (Q094 = 1: "Sim").
    # Q094: 1 Sim | 2 Não, só quando tem algum problema | 3 Nunca vai.
    dx_regular = as.numeric(dep_dx %in% 1 & Q094 %in% 1)
  )

# ---------------------------------------------------------------------------
# Domínios e ridit de renda de cada domínio (calculado ANTES do desenho)
# ---------------------------------------------------------------------------
dominios <- list(
  D1_todos_phq10          = pns$phq10 %in% 1,
  D2_phq10_sem_dx_regular = pns$phq10 %in% 1 & pns$dx_regular == 0,
  D3_phq15                = pns$phq15 %in% 1
)
outcomes <- c('undiagnosed', 'undiagnosed_any')
for (dn in names(dominios)) for (y in outcomes) {
  pns[[paste0('dom_', dn, '_', y)]]  <- as.numeric(dominios[[dn]] & !is.na(pns[[y]]))
  pns[[paste0('rank_', dn, '_', y)]] <- rank_in(dominios[[dn]] & !is.na(pns[[y]]))
}

des <- PNSIBGE::pns_design(pns)

use_rank <- function(design, rank_var) {
  design$variables$income_rank <- design$variables[[rank_var]]
  design
}
in_dom <- function(dd) is.finite(dd$prob)
fit_pr <- function(dd, y, rhs) tryCatch(suppressWarnings(
  svyglm(as.formula(paste(y, '~', rhs)), design = dd, family = quasipoisson(link = 'log'))),
  error = function(e) e)

# ---------------------------------------------------------------------------
# 1. Descritiva: prevalência da lacuna por domínio, por faixa de SM e por quintil do domínio
# ---------------------------------------------------------------------------
desc_rows <- list(); faixa_rows <- list()
for (dn in names(dominios)) for (y in outcomes) {
  dd <- subset(des, get(paste0('dom_', dn, '_', y)) == 1)
  z <- svymean(as.formula(paste0('~', y)), dd, na.rm = TRUE)
  cr <- qt(.975, degf(dd))
  desc_rows[[paste(dn, y)]] <- data.frame(dominio = dn, desfecho = y, n = sum(in_dom(dd)),
    casos = sum(dd$variables[[y]][in_dom(dd)] == 1, na.rm = TRUE),
    prevalencia = as.numeric(coef(z)), lower = as.numeric(coef(z)) - cr * as.numeric(SE(z)),
    upper = as.numeric(coef(z)) + cr * as.numeric(SE(z)), degf = degf(dd))
  dd$variables$quintil_dominio <- cut(dd$variables[[paste0('rank_', dn, '_', y)]], breaks = c(0, .2, .4, .6, .8, 1),
                                      labels = paste0('Q', 1:5), include.lowest = TRUE)
  for (by in c('income3', 'quintil_dominio')) {
    dby <- subset(dd, !is.na(get(by)))
    zb <- svyby(as.formula(paste0('~', y)), as.formula(paste0('~', by)), dby, svymean,
                na.rm = TRUE, vartype = c('se', 'ci'), keep.names = FALSE)
    names(zb)[1:2] <- c('categoria', 'prevalencia')
    faixa_rows[[paste(dn, y, by)]] <- cbind(dominio = dn, desfecho = y, por = by, zb)
  }
}
descritiva <- bind_rows(desc_rows)
wcsv(descritiva, '01_lacuna_prevalencia_por_dominio.csv')
wcsv(bind_rows(faixa_rows), '02_lacuna_por_faixa_e_quintil.csv')

# ---------------------------------------------------------------------------
# 2. Gradiente de renda: rank do domínio (principal) e rank populacional (comparação)
# ---------------------------------------------------------------------------
modelos <- list(minimo    = 'income_rank + age_c + sex + region',
                expandido = 'income_rank + age_c + sex + region + race3 + education3 + urban + plan')
pr_rows <- list()
for (dn in names(dominios)) for (y in outcomes) for (escopo in c('dominio', 'populacao')) {
  rank_var <- if (escopo == 'dominio') paste0('rank_', dn, '_', y) else 'income_rank_pop'
  dd <- use_rank(subset(des, get(paste0('dom_', dn, '_', y)) == 1 & !is.na(get(rank_var))), rank_var)
  for (nm in names(modelos)) {
    m <- fit_pr(dd, y, modelos[[nm]])
    if (inherits(m, 'error')) {
      pr_rows[[paste(dn, y, escopo, nm)]] <- data.frame(dominio = dn, desfecho = y, rank = escopo, modelo = nm, erro = conditionMessage(m))
      next
    }
    b <- coef(m)[['income_rank']]; se <- sqrt(diag(vcov(m)))[['income_rank']]
    df <- m$df.residual; cr <- qt(.975, df)
    pr_rows[[paste(dn, y, escopo, nm)]] <- data.frame(
      dominio = dn, desfecho = y, rank = escopo, modelo = nm,
      n = sum(m$prior.weights > 0), casos = sum(m$y[m$prior.weights > 0] == 1),
      RII = exp(b), lower = exp(b - cr * se), upper = exp(b + cr * se),
      p = 2 * pt(abs(b / se), df, lower.tail = FALSE), df_residual = df,
      max_prev_prevista = max(m$fitted.values[m$prior.weights > 0]), erro = NA_character_)
  }
}
gradiente <- bind_rows(pr_rows)
wcsv(gradiente, '03_lacuna_RII_renda.csv')

# ---------------------------------------------------------------------------
# Resumo
# ---------------------------------------------------------------------------
cat('\n=== CONFERÊNCIA (esperado: 9.495 com PHQ-9 >= 10; 64% sem diagnóstico) ===\n')
cat('PHQ-9 completo:', sum(!is.na(pns$phq9_total)), '| PHQ-9 >= 10:', sum(pns$phq10 %in% 1),
    '| PHQ-9 >= 15:', sum(pns$phq15 %in% 1), '\n')
cat('\n=== PREVALÊNCIA DA LACUNA POR DOMÍNIO (ponderada) ===\n')
print(descritiva[, c('dominio','desfecho','n','casos','prevalencia','lower','upper')], digits = 3)
cat('\n=== GRADIENTE DE RENDA (RII; < 1 = lacuna maior entre os mais pobres) ===\n')
print(gradiente[, intersect(c('dominio','desfecho','rank','modelo','n','casos','RII','lower','upper','p','erro'), names(gradiente))], digits = 3)
cat('\nResultados em:', out_dir, '\n')
