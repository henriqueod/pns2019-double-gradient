# =============================================================================
# 03_income_bands.R
# Five per-capita income bands: prevalence, prevalence ratios and shape tests.
# Repository name: R/03_income_bands.R   |   Development name: 31_faixas_renda_finas.R
# Code comments are in Portuguese. Cross-references in the comments use the
# development names; README.md maps them to the repository names.
# Paths: all input and output paths are relative to `root`, which defaults to
# the working directory (run from the repository root) or to the PNS_ROOT
# environment variable if it is set.
# =============================================================================

# 31_faixas_renda_finas.R
# Refaz as análises por FAIXA DE RENDA com mais resolução nas pontas.
#
# POR QUÊ
#   Com a coorte de 10.956, as faixas extremas de renda têm casos suficientes para serem
#   abertas (três faixas, <=1, >1-2, >2 SM, esconderiam a forma do gradiente). A VDF004 do IBGE tem
#   7 classes de renda domiciliar per capita: 1 <=1/4 | 2 >1/4-1/2 | 3 >1/2-1 | 4 >1-2 |
#   5 >2-3 | 6 >3-5 | 7 >5 SM.
#
# FAIXAS NOVAS (padrão; ver FAIXAS abaixo para mudar)
#   <=1/2 SM (VDF004 1-2) | >1/2-1 (3) | >1-2 (4) | >2-5 (5-6) | >5 SM (7)
#   A extrema pobreza (<=1/4) e o topo (>5) não são separados por padrão porque as
#   células ficariam pequenas; a seção 0 imprime as contagens das 7 classes para decidir.
#
# O QUE NÃO MUDA: todos os RII/SII (usam o ridit contínuo), quintis e modelos em blocos.
# O QUE É REFEITO (só o que usa a faixa categórica):
#   A. Coorte (results/mental_care_design_anymh.rds, do 01c):
#      A0 contagens por classe VDF004 (7) e pelas faixas novas, com eventos -> decisão
#      A1 distribuição da faixa na coorte e por grupo diagnóstico (linha da Tabela 1)
#      A2 prevalência dos desfechos por faixa (Tabelas 3 e 6, figura do leque)
#      A3 PR por faixa, bloco B e bloco D (ref = faixa mais pobre)
#      A4 forma do gradiente com as faixas novas (associação e desvio de linearidade)
#      A5 grupo diagnóstico x faixa (Tabela 4), com n e eventos por célula
#   B. Cascata de depressão e lacuna de detecção (recria as variáveis do 16b a partir do
#      PNSIBGE): prevalência por faixa de cada indicador da Tabela 6, com n e eventos por
#      célula, e PR por faixa no modelo mínimo (faixa + idade + sexo + região).
#
# Saída: results/income_bands_fine/*.csv
# Células com menos de 20 eventos são marcadas (coluna flag_eventos_lt20).

options(timeout = 600)
options(survey.lonely.psu = 'adjust')
options(survey.adjust.domain.lonely = TRUE)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(PNSIBGE)
  library(survey)
  library(dplyr)
})
cat('03_income_bands.R (31_faixas_renda_finas.R) | versão 26/09/2026\n')

root <- Sys.getenv('PNS_ROOT', unset = getwd())
data_dir <- file.path(root, 'data_pnsibge')
out  <- file.path(root, 'results', 'income_bands_fine')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
wcsv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE)

# ---------------------------------------------------------------------------
# Definição das faixas (classes da VDF004 -> faixa). Para testar outro corte, mude aqui.
# ---------------------------------------------------------------------------
FAIXAS <- list(
  breaks = c(0, 2, 3, 4, 6, 7),                          # (0,2] (2,3] (3,4] (4,6] (6,7]
  labels = c('ate_meio_SM', 'meio_a_1_SM', '1_a_2_SM', '2_a_5_SM', 'mais_5_SM')
)
SM_2019 <- 998   # só usado se a VDF004 não estiver disponível (reconstrói a partir da VDF003)
vdf004_from_vdf003 <- function(x) as.numeric(cut(x / SM_2019, c(-Inf, .25, .5, 1, 2, 3, 5, Inf), right = TRUE))
make_band <- function(vdf004) factor(cut(vdf004, breaks = FAIXAS$breaks, labels = FAIXAS$labels, right = TRUE),
                                     levels = FAIXAS$labels)
vdf004_labels <- c('1_ate_1/4', '2_1/4_a_1/2', '3_1/2_a_1', '4_1_a_2', '5_2_a_3', '6_3_a_5', '7_mais_5')

in_dom <- function(dd) is.finite(dd$prob) & dd$prob > 0
n_dom  <- function(dd) sum(in_dom(dd))

# Prevalência ponderada por categoria, com n e eventos não ponderados por célula.
prev_cells <- function(dd, y, by, label = y) {
  dd <- subset(dd, !is.na(get(y)) & !is.na(get(by)))
  v <- dd$variables[in_dom(dd), , drop = FALSE]
  cnt <- v %>% group_by(cat = .data[[by]], .drop = FALSE) %>%
    summarise(n = n(), eventos = sum(.data[[y]] == 1), .groups = 'drop')
  z <- tryCatch(svyby(as.formula(paste0('~', y)), as.formula(paste0('~', by)), dd, svymean,
                      na.rm = TRUE, vartype = c('se', 'ci'), keep.names = FALSE, drop.empty.groups = FALSE),
                error = function(e) NULL)
  if (is.null(z)) {
    res <- data.frame(cat = cnt$cat, prevalencia = NA, se = NA, ci_l = NA, ci_u = NA)
  } else {
    res <- data.frame(cat = z[[by]], prevalencia = z[[y]], se = z$se, ci_l = z$ci_l, ci_u = z$ci_u)
  }
  res <- merge(cnt, res, by = 'cat', all = TRUE)
  res$cat <- factor(res$cat, levels = levels(v[[by]])); res <- res[order(res$cat), ]
  cbind(desfecho = label, por = by, res, flag_eventos_lt20 = res$eventos < 20)
}

# PR por faixa (referência = primeira faixa) num modelo Poisson modificado.
pr_band <- function(dd, y, rhs_extra, band = 'income_band', label = y, modelo = '') {
  dd <- subset(dd, !is.na(get(y)) & !is.na(get(band)))
  m <- tryCatch(suppressWarnings(svyglm(as.formula(paste(y, '~', band, '+', rhs_extra)), design = dd,
                                        family = quasipoisson(link = 'log'))), error = function(e) e)
  terms <- paste0(band, levels(dd$variables[[band]])[-1])
  if (inherits(m, 'error'))
    return(data.frame(desfecho = label, modelo = modelo, termo = terms, PR = NA, lower = NA, upper = NA,
                      p = NA, df = NA, nota = conditionMessage(m)))
  bind_rows(lapply(terms, function(tm) {
    if (!(tm %in% names(coef(m))) || is.na(coef(m)[[tm]]))
      return(data.frame(desfecho = label, modelo = modelo, termo = tm, PR = NA, lower = NA, upper = NA,
                        p = NA, df = m$df.residual, nota = 'termo ausente (célula vazia?)'))
    b <- coef(m)[[tm]]; se <- sqrt(diag(vcov(m)))[[tm]]; df <- m$df.residual; cr <- qt(.975, df)
    data.frame(desfecho = label, modelo = modelo, termo = tm, PR = exp(b), lower = exp(b - cr * se),
               upper = exp(b + cr * se), p = 2 * pt(abs(b / se), df, lower.tail = FALSE), df = df, nota = NA)
  }))
}
wald <- function(m, f) {
  if (inherits(m, 'error')) return(data.frame(F = NA, df_num = NA, df_den = NA, p = NA, nota = conditionMessage(m)))
  t <- tryCatch(regTermTest(m, f, method = 'Wald'), error = function(e) e)
  if (inherits(t, 'error')) return(data.frame(F = NA, df_num = NA, df_den = NA, p = NA, nota = conditionMessage(t)))
  data.frame(F = as.numeric(t$Ftest), df_num = as.numeric(t$df), df_den = as.numeric(t$ddf), p = as.numeric(t$p), nota = NA)
}
fitp <- function(dd, y, rhs) tryCatch(suppressWarnings(svyglm(as.formula(paste(y, '~', rhs)), design = dd,
                                       family = quasipoisson(link = 'log'))), error = function(e) e)

# ===========================================================================
# PARTE A — coorte "qualquer doença mental"
# ===========================================================================
cat('\n== PARTE A: coorte ==\n')
des <- readRDS(file.path(root, 'results', 'mental_care_design_anymh.rds'))
d <- des$variables
bin_and <- function(a, b) ifelse(is.na(a) | is.na(b), NA_real_, as.numeric(a == 1 & b == 1))
d$high_limitation <- ifelse(is.na(d$limitation), NA_real_, as.numeric(d$limitation >= 3))
d$no_regular      <- ifelse(is.na(d$regular_care), NA_real_, 1 - d$regular_care)
d$high_need_gap   <- bin_and(d$high_limitation, d$no_regular)
if (!('VDF004' %in% names(d)) || all(is.na(d$VDF004))) {
  cat('VDF004 ausente no desenho: reconstruída a partir da VDF003 (SM 2019 = R$', SM_2019, ')\n')
  d$VDF004 <- vdf004_from_vdf003(d$VDF003)
}
d$vdf004_cls  <- factor(d$VDF004, levels = 1:7, labels = vdf004_labels)
d$income_band <- make_band(d$VDF004)

# Mesma amostra analítica do 30 (todas as covariáveis dos blocos A-E completas).
cov_all <- c('income_rank_mental','age_c','sex','region','race3','education3','urban','plan',
             'family_support0','friend_support0','income3')
cc <- rep(TRUE, nrow(d)); for (v in cov_all) cc <- cc & !is.na(d[[v]])
d$analytic <- as.numeric(d$mental_dx %in% 1 & cc)
des$variables <- d
coorte <- subset(des, analytic == 1)
cat('Coorte analítica: n =', n_dom(coorte), '(esperado 10.956)\n')
outs <- c('psychotherapy', 'high_need_gap', 'regular_care', 'high_limitation')

# A0. Contagens para decidir as faixas: 7 classes da VDF004 e faixas novas.
a0 <- bind_rows(lapply(c('vdf004_cls', 'income_band'), function(by)
  bind_rows(lapply(outs, function(y) prev_cells(coorte, y, by)))))
wcsv(a0, 'A0_contagens_classes_e_faixas.csv')
a0g <- bind_rows(lapply(levels(d$dx_hier), function(lv) {
  gg <- subset(coorte, dx_hier == lv)
  cbind(grupo = lv, bind_rows(lapply(c('psychotherapy', 'high_need_gap'), function(y) prev_cells(gg, y, 'income_band'))))
}))
cat('\nEventos por faixa (coorte):\n')
print(subset(a0, por == 'income_band' & desfecho %in% c('psychotherapy', 'high_need_gap'))[, c('desfecho','cat','n','eventos','prevalencia')], digits = 3)
cat('\nEventos por classe VDF004 (coorte):\n')
print(subset(a0, por == 'vdf004_cls' & desfecho %in% c('psychotherapy', 'high_need_gap'))[, c('desfecho','cat','n','eventos','prevalencia')], digits = 3)

# A1. Distribuição ponderada da faixa (coorte e por grupo) - linha da Tabela 1.
dist_band <- function(dd, label) {
  z <- svymean(~income_band, subset(dd, !is.na(income_band)), na.rm = TRUE)
  zc <- svymean(~vdf004_cls, subset(dd, !is.na(vdf004_cls)), na.rm = TRUE)
  rbind(data.frame(grupo = label, variavel = 'income_band', nivel = sub('^income_band', '', names(coef(z))),
                   pct = 100 * as.numeric(coef(z)), se = 100 * as.numeric(SE(z))),
        data.frame(grupo = label, variavel = 'vdf004_cls', nivel = sub('^vdf004_cls', '', names(coef(zc))),
                   pct = 100 * as.numeric(coef(zc)), se = 100 * as.numeric(SE(zc))))
}
a1 <- bind_rows(dist_band(coorte, 'coorte'),
                bind_rows(lapply(levels(d$dx_hier), function(lv) dist_band(subset(coorte, dx_hier == lv), lv))))
wcsv(a1, 'A1_distribuicao_faixas_tabela1.csv')

# A2. Prevalência por faixa (já está em A0; gravada à parte para as tabelas).
wcsv(subset(a0, por == 'income_band'), 'A2_prevalencia_por_faixa.csv')

# A3. PR por faixa, blocos B e D.
rhsB <- 'age_c + sex + region'
rhsD <- 'age_c + sex + region + race3 + education3 + urban + plan'
a3 <- bind_rows(lapply(c('psychotherapy', 'high_need_gap'), function(y)
  rbind(pr_band(coorte, y, rhsB, modelo = 'B'), pr_band(coorte, y, rhsD, modelo = 'D'))))
wcsv(a3, 'A3_PR_faixas_blocos_B_D.csv')

# A4. Forma do gradiente com as faixas novas (bloco B): associação e desvio de linearidade.
a4 <- bind_rows(lapply(c('psychotherapy', 'high_need_gap'), function(y) {
  dd <- subset(coorte, !is.na(get(y)))
  m_assoc  <- fitp(dd, y, paste('income_band +', rhsB))
  m_desvio <- fitp(dd, y, paste('income_rank_mental + income_band +', rhsB))
  rbind(cbind(desfecho = y, teste = 'associação (sem o rank)', wald(m_assoc, ~income_band)),
        cbind(desfecho = y, teste = 'desvio de linearidade (condicional ao rank)', wald(m_desvio, ~income_band)))
}))
wcsv(a4, 'A4_forma_faixas_novas.csv')

# A5. Grupo diagnóstico x faixa.
wcsv(a0g, 'A5_grupo_x_faixa.csv')

cat('\nPR bloco B por faixa (ref = <=1/2 SM):\n'); print(subset(a3, modelo == 'B')[, c('desfecho','termo','PR','lower','upper','p')], digits = 3)
cat('\nForma:\n'); print(a4[, c('desfecho','teste','F','df_num','p')], digits = 3)

# ---- Preparação idêntica ao 16b_depression_cascade_domain_rank.R (linhas 37-136) ----
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

# ===========================================================================
# PARTE B — cascata de depressão e lacuna de detecção (variáveis idênticas ao 16b/18b)
# ===========================================================================
cat('\n== PARTE B: depressão e detecção ==\n')
pns$oth_dx <- ifelse(pns$Q11006 == 1, 1, ifelse(pns$Q11006 == 2, 0, NA_real_))
pns$any_dx_among_phq10 <- ifelse(pns$phq10 %in% 1 & (pns$depression_dx %in% 1 | pns$oth_dx %in% 1), 1,
                          ifelse(pns$phq10 %in% 1 & pns$depression_dx %in% 0 & pns$oth_dx %in% 0, 0, NA_real_))
pns$remission_phq10 <- ifelse(pns$depression_dx %in% 1 & pns$Q094 %in% c(2,3) & pns$reason_group %in% 'remissao_percebida',
                              pns$phq10, NA_real_)
if (!('VDF004' %in% names(pns)) || all(is.na(pns$VDF004))) pns$VDF004 <- vdf004_from_vdf003(pns$VDF003)
pns$vdf004_cls  <- factor(pns$VDF004, levels = 1:7, labels = vdf004_labels)
pns$income_band <- make_band(pns$VDF004)

desB <- PNSIBGE::pns_design(pns)

# Indicador, domínio (expressão), rótulo. Os domínios são os mesmos do 16b/18b.
ind <- list(
  list('depression_dx',              'phq10 %in% 1',                                              'diagnosticado_depressao_entre_PHQ10'),
  list('any_dx_among_phq10',         'phq10 %in% 1',                                              'qualquer_diagnostico_entre_PHQ10'),
  list('specialist_reached_all',     'depression_dx %in% 1 & Q10101 %in% 1:4',                    'chegou_ao_especialista'),
  list('structural_reason',          'depression_dx %in% 1 & Q094 %in% c(2,3)',                   'barreira_estrutural_motivo'),
  list('hidden_remission_gap',       'TRUE',                                                      'remissao_percebida_PHQ10_sem_tratamento'),
  list('hidden_remission_severe',    'TRUE',                                                      'remissao_percebida_PHQ15_sem_tratamento'),
  list('hidden_remission_selfharm',  'TRUE',                                                      'remissao_percebida_ideacao_sem_tratamento'),
  list('remission_phq10',            'TRUE',                                                      'remissao_percebida_com_PHQ10'),
  list('critical_depression_gap',    'TRUE',                                                      'lacuna_critica_ideacao_sem_tratamento'),
  list('referral_failure',           'depression_dx %in% 1 & Q10101 %in% 1:4',                    'falha_apos_encaminhamento'),
  list('depression_current_treatment','depression_dx %in% 1',                                     'depressao_em_tratamento_atual'),
  list('depression_no_regular_care', 'depression_dx %in% 1',                                      'depressao_sem_seguimento_regular')
)
b1 <- list(); b2 <- list()
for (it in ind) {
  y <- it[[1]]; dom <- it[[2]]; lab <- it[[3]]
  cat('  ', lab, '\n'); flush.console()
  dd <- subset(desB, eval(parse(text = dom)) & !is.na(get(y)))
  b1[[lab]] <- bind_rows(prev_cells(dd, y, 'income_band', lab), prev_cells(dd, y, 'vdf004_cls', lab))
  b2[[lab]] <- pr_band(dd, y, 'age_c + sex + region', label = lab, modelo = 'minimo')
}
b1 <- bind_rows(b1); b2 <- bind_rows(b2)
wcsv(b1, 'B1_depressao_prevalencia_por_faixa.csv')
wcsv(b2, 'B2_depressao_PR_faixas_modelo_minimo.csv')

cat('\nDepressão/detecção por faixa nova (prevalência e eventos):\n')
print(subset(b1, por == 'income_band')[, c('desfecho','cat','n','eventos','prevalencia','flag_eventos_lt20')], digits = 3)
cat('\nCélulas com < 20 eventos (faixas novas):', sum(subset(b1, por == 'income_band')$flag_eventos_lt20, na.rm = TRUE),
    '| (coorte):', sum(subset(a0, por == 'income_band')$flag_eventos_lt20, na.rm = TRUE), '\n')
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
cat('\nFim. Arquivos em', out, '\n')
