# =============================================================================
# 02_main_analysis.R
# Co-primary outcomes: Table 1, adjustment blocks A-E, shape tests, interactions, diagnostic groups, sensitivity analyses and prediction-range audit.
# Repository name: R/02_main_analysis.R   |   Development name: 30_anymh_analysis.R
# Code comments are in Portuguese. Cross-references in the comments use the
# development names; README.md maps them to the repository names.
# Paths: all input and output paths are relative to `root`, which defaults to
# the working directory (run from the repository root) or to the PNS_ROOT
# environment variable if it is set.
# =============================================================================

# 30_anymh_analysis.R
# Análise principal na coorte "qualquer doença mental" (Q092 OU Q11006),
# com comparação por grupo diagnóstico hierárquico.
#
# Entrada: results/mental_care_design_anymh.rds, gerado pelo 01c_prepare_any_mental.R.
# Saída:   results/anymh/*.csv e um resumo no console.
#
# O QUE ESTE SCRIPT FAZ
#   1. Descritiva: Tabela 1 ponderada (coorte e por grupo) e prevalência dos desfechos
#      por faixa de renda, por grupo e por grupo x faixa.
#   2. Duplo gradiente na coorte: blocos A-E para psicoterapia e high need gap
#      (RII por Poisson modificado, SII por modelo linear), mesma amostra em todos os blocos.
#   3. Forma do gradiente: associação categórica (faixas de SM e quintis, SEM o rank) e
#      desvio de linearidade (bloco categórico testado CONDICIONALMENTE ao rank linear).
#   4. Interações pré-especificadas: renda x esquizofrenia e renda x bipolar (Wald conjunto).
#   5. Comparação entre grupos (post hoc): renda x grupo hierárquico nos blocos B e D;
#      RII de cada grupo a partir do modelo com interação (mesma escala de rank).
#   6. Duplo gradiente dentro de cada grupo (blocos A-E), com eventos por parâmetro.
#   7. Sensibilidades dentro dos grupos: itens do próprio bloco do questionário e ridit
#      calculado dentro do grupo.
#   8. Sensibilidades na coorte: 18+, só autorrespondentes, definições alternativas do
#      high need gap, ajuste por grupo diagnóstico, família de modelo (log-binomial e logística).
#
# CONVENÇÕES
#   - Rank de renda: ridit ponderado dentro da coorte (income_rank_mental), o mesmo em
#     todos os modelos de coorte e de comparação entre grupos. O ridit dentro do grupo
#     entra só como sensibilidade (seção 7), porque muda a escala do RII.
#   - Graus de liberdade reportados: df residual do svyglm (degf + 1 - nº de coeficientes).
#   - Blocos: A renda; B + idade, sexo, região; C + raça/cor, escolaridade, urbano;
#     D + plano de saúde; E + apoio de família e de amigos.

options(survey.lonely.psu = 'adjust')
options(survey.adjust.domain.lonely = TRUE)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survey)
  library(dplyr)
})

cat('02_main_analysis.R (30_anymh_analysis.R) | versão 26/09/2026 | seções 1-9 (inclui checagem de previsões > 1 e log-binomial)\n')
root <- Sys.getenv('PNS_ROOT', unset = getwd())
des  <- readRDS(file.path(root, 'results', 'mental_care_design_anymh.rds'))
out  <- file.path(root, 'results', 'anymh')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
wcsv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE)

# ---------------------------------------------------------------------------
# 0. Variáveis derivadas (antes de qualquer subset)
# ---------------------------------------------------------------------------
make_weighted_rank <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0 & x >= 0
  r <- rep(NA_real_, length(x))
  ord <- order(x[ok], na.last = NA)
  idx <- which(ok)[ord]
  ww <- w[idx]
  r[idx] <- (cumsum(ww) - 0.5 * ww) / sum(ww)
  r
}
bin_and <- function(a, b) ifelse(is.na(a) | is.na(b), NA_real_, as.numeric(a == 1 & b == 1))

d <- des$variables
d$high_limitation <- ifelse(is.na(d$limitation), NA_real_, as.numeric(d$limitation >= 3))
d$no_regular      <- ifelse(is.na(d$regular_care), NA_real_, 1 - d$regular_care)
d$high_need_gap   <- bin_and(d$high_limitation, d$no_regular)

# Definições alternativas do high need gap (sensibilidade; mesma lógica da Tabela 5).
d$hng_lim4     <- bin_and(ifelse(is.na(d$limitation), NA_real_, as.numeric(d$limitation >= 4)), d$no_regular)
d$hng_strict   <- bin_and(d$high_limitation, ifelse(is.na(d$continuity3), NA_real_, as.numeric(d$continuity3 == 2)))
d$hng_episodic <- bin_and(d$high_limitation, ifelse(is.na(d$continuity3), NA_real_, as.numeric(d$continuity3 == 3)))
# Nas duas últimas, quem tem a outra forma de não continuidade conta como 0 (não é excluído).

# Itens do próprio bloco do questionário, para a sensibilidade dentro dos grupos:
# grupos "grave" e "outras" -> bloco Q110a (originais); "depressao" -> bloco de depressão.
g <- as.character(d$dx_hier)
use_q110 <- g %in% c('grave', 'outras'); use_dep <- g %in% 'depressao'
d$blk_regular <- ifelse(use_q110 & d$orig_Q11201 %in% 1:3, as.numeric(d$orig_Q11201 == 1),
                 ifelse(use_dep  & d$Q094        %in% 1:3, as.numeric(d$Q094 == 1), NA_real_))
d$blk_limit   <- ifelse(use_q110 & d$orig_Q115 %in% 1:5, d$orig_Q115,
                 ifelse(use_dep  & d$Q109      %in% 1:5, d$Q109, NA_real_))
d$blk_psychotherapy <- ifelse(use_q110 & d$orig_Q11405 %in% 1:2, as.numeric(d$orig_Q11405 == 1),
                       ifelse(use_dep  & d$Q09605      %in% 1:2, as.numeric(d$Q09605 == 1), NA_real_))
d$blk_high_need_gap <- bin_and(ifelse(is.na(d$blk_limit), NA_real_, as.numeric(d$blk_limit >= 3)),
                               ifelse(is.na(d$blk_regular), NA_real_, 1 - d$blk_regular))

# Quintis do ridit da coorte (ridit ponderado -> cortes fixos = quintis ponderados).
d$income5 <- cut(d$income_rank_mental, breaks = c(0, .2, .4, .6, .8, 1),
                 labels = paste0('Q', 1:5), include.lowest = TRUE)

# Ridit dentro de cada grupo hierárquico (só para a sensibilidade da seção 7).
d$rank_group <- NA_real_
for (lv in levels(d$dx_hier)) {
  i <- which(d$dx_hier %in% lv & !is.na(d$income_reais) & !is.na(d$V00291))
  d$rank_group[i] <- make_weighted_rank(d$income_reais[i], d$V00291[i])
}
des$variables <- d

# ---------------------------------------------------------------------------
# Amostra analítica comum (todas as covariáveis dos blocos A-E completas),
# para que a mudança entre blocos seja de especificação, não de denominador.
# ---------------------------------------------------------------------------
blocks <- list(
  A = 'income_rank_mental',
  B = 'income_rank_mental + age_c + sex + region',
  C = 'income_rank_mental + age_c + sex + region + race3 + education3 + urban',
  D = 'income_rank_mental + age_c + sex + region + race3 + education3 + urban + plan',
  E = 'income_rank_mental + age_c + sex + region + race3 + education3 + urban + plan + family_support0 + friend_support0'
)
cov_all <- c('income_rank_mental','age_c','sex','region','race3','education3','urban','plan',
             'family_support0','friend_support0','income3')
cc <- rep(TRUE, nrow(d)); for (v in cov_all) cc <- cc & !is.na(d[[v]])
des$variables$analytic <- as.numeric(d$mental_dx %in% 1 & cc)
coorte <- subset(des, analytic == 1)
outcomes <- c('psychotherapy', 'high_need_gap')

# ---------------------------------------------------------------------------
# Utilitários
# ---------------------------------------------------------------------------
in_dom <- function(dd) is.finite(dd$prob)                     # linhas do domínio
n_dom  <- function(dd) sum(in_dom(dd))
fit <- function(dd, y, rhs, family = quasipoisson(link = 'log'), start = NULL, ...) {
  tryCatch(suppressWarnings(
    svyglm(as.formula(paste(y, '~', rhs)), design = dd, family = family, start = start, ...)),
    error = function(e) e)
}
est_row <- function(m, term, transform = exp) {
  if (inherits(m, 'error') || !(term %in% names(coef(m))))
    return(data.frame(est = NA, lower = NA, upper = NA, p = NA, df = NA,
                      converged = NA, note = if (inherits(m, 'error')) conditionMessage(m) else 'termo ausente'))
  b <- coef(m)[[term]]; se <- sqrt(diag(vcov(m)))[[term]]; df <- m$df.residual; cr <- qt(.975, df)
  data.frame(est = transform(b), lower = transform(b - cr * se), upper = transform(b + cr * se),
             p = 2 * pt(abs(b / se), df, lower.tail = FALSE), df = df,
             converged = isTRUE(m$converged), note = NA_character_)
}
wald <- function(m, term_formula) {
  if (inherits(m, 'error')) return(data.frame(F = NA, df_num = NA, df_den = NA, p = NA, note = conditionMessage(m)))
  t <- tryCatch(regTermTest(m, term_formula, method = 'Wald'), error = function(e) e)
  if (inherits(t, 'error')) return(data.frame(F = NA, df_num = NA, df_den = NA, p = NA, note = conditionMessage(t)))
  data.frame(F = as.numeric(t$Ftest), df_num = as.numeric(t$df), df_den = as.numeric(t$ddf),
             p = as.numeric(t$p), note = NA_character_)
}
n_coef <- function(rhs, dd) {
  mf <- model.frame(as.formula(paste('~', rhs)), dd$variables[in_dom(dd), , drop = FALSE])
  ncol(model.matrix(as.formula(paste('~', rhs)), mf))
}

# Duplo gradiente: blocos A-E, RII e SII, numa amostra fixa.
run_blocks <- function(dd, y, label, rank_var = 'income_rank_mental', block_list = blocks) {
  dd <- subset(dd, !is.na(get(y)))
  n <- n_dom(dd); ev <- sum(dd$variables[[y]][in_dom(dd)] == 1)
  rows <- lapply(names(block_list), function(b) {
    rhs <- gsub('income_rank_mental', rank_var, block_list[[b]], fixed = TRUE)
    mp <- fit(dd, y, rhs)
    ml <- fit(dd, y, rhs, family = gaussian())
    k  <- n_coef(rhs, dd)
    rii <- est_row(mp, rank_var); sii <- est_row(ml, rank_var, transform = identity)
    maxfit <- if (inherits(mp, 'error')) NA else max(fitted(mp)[mp$prior.weights > 0])
    lin_range <- if (inherits(ml, 'error')) c(NA, NA) else range(fitted(ml)[ml$prior.weights > 0])
    data.frame(analise = label, desfecho = y, bloco = b, n = n, eventos = ev, n_coef = k,
               eventos_por_coef = round(ev / k, 1),
               RII = rii$est, RII_lower = rii$lower, RII_upper = rii$upper, RII_p = rii$p,
               SII = sii$est, SII_lower = sii$lower, SII_upper = sii$upper, SII_p = sii$p,
               df_residual = rii$df, convergiu = rii$converged,
               max_prev_poisson = maxfit, min_prev_linear = lin_range[1], max_prev_linear = lin_range[2],
               nota = rii$note)
  })
  res <- bind_rows(rows)
  ref <- log(res$RII[res$bloco == 'B'])
  res$atenuacao_vs_B_pct <- round(100 * (1 - log(res$RII) / ref), 1)
  res
}

# Tabela descritiva ponderada.
describe <- function(dd, label) {
  num_vars <- c('age')
  cat_vars <- c('sex','race3','urban','education3','plan','self_report','income3',
                'dx_hier','depression','schizophrenia','bipolar','psychotherapy','high_need_gap',
                'regular_care','high_limitation')
  rows <- list()
  for (v in num_vars) {
    z <- svymean(as.formula(paste0('~', v)), dd, na.rm = TRUE)
    rows[[v]] <- data.frame(grupo = label, variavel = v, nivel = 'média', estimativa = as.numeric(coef(z)),
                            se = as.numeric(SE(z)), n_valido = sum(in_dom(dd) & !is.na(dd$variables[[v]])))
  }
  for (v in cat_vars) {
    x <- dd$variables[[v]]
    if (all(is.na(x[in_dom(dd)]))) next
    dv <- dd; dv$variables$.v <- factor(x)
    dv <- subset(dv, !is.na(.v))
    if (nlevels(droplevels(dv$variables$.v[in_dom(dv)])) < 2) next
    z <- svymean(~.v, dv, na.rm = TRUE)
    rows[[v]] <- data.frame(grupo = label, variavel = v, nivel = sub('^\\.v', '', names(coef(z))),
                            estimativa = 100 * as.numeric(coef(z)), se = 100 * as.numeric(SE(z)),
                            n_valido = n_dom(dv))
  }
  bind_rows(rows)
}

prev_by <- function(dd, y, by) {
  dd <- subset(dd, !is.na(get(y)) & !is.na(get(by)))
  z <- svyby(as.formula(paste0('~', y)), as.formula(paste0('~', by)), dd, svymean,
             na.rm = TRUE, vartype = c('se', 'ci'), keep.names = FALSE, drop.empty.groups = FALSE)
  names(z)[names(z) == y] <- 'prevalencia'
  cbind(desfecho = y, z)
}

# ---------------------------------------------------------------------------
# 1. Descritiva
# ---------------------------------------------------------------------------
tab1 <- bind_rows(describe(coorte, 'coorte'),
                  bind_rows(lapply(levels(d$dx_hier), function(lv) describe(subset(coorte, dx_hier == lv), lv))))
wcsv(tab1, '01_tabela1_ponderada.csv')

prev <- bind_rows(lapply(outcomes, function(y) {
  z <- svymean(as.formula(paste0('~', y)), subset(coorte, !is.na(get(y))), na.rm = TRUE)
  data.frame(desfecho = y, prevalencia = as.numeric(coef(z)), se = as.numeric(SE(z)),
             n = n_dom(subset(coorte, !is.na(get(y)))))
}))
wcsv(prev, '02_prevalencia_coorte.csv')
wcsv(bind_rows(lapply(outcomes, prev_by, dd = coorte, by = 'income3')), '03_prevalencia_por_faixa_renda.csv')
wcsv(bind_rows(lapply(outcomes, prev_by, dd = coorte, by = 'dx_hier')), '04_prevalencia_por_grupo.csv')
dd_gx <- coorte; dd_gx$variables$grupo_faixa <- interaction(dd_gx$variables$dx_hier, dd_gx$variables$income3, sep = ' | ')
wcsv(bind_rows(lapply(outcomes, prev_by, dd = dd_gx, by = 'grupo_faixa')), '05_prevalencia_grupo_x_faixa.csv')

# ---------------------------------------------------------------------------
# 2. Duplo gradiente na coorte (blocos A-E)
# ---------------------------------------------------------------------------
blocos_coorte <- bind_rows(lapply(outcomes, function(y) run_blocks(coorte, y, 'coorte')))
wcsv(blocos_coorte, '06_blocos_coorte.csv')

# PR por faixa de renda (bloco B), referência = até 1 SM.
cat_pr <- bind_rows(lapply(outcomes, function(y) {
  dd <- subset(coorte, !is.na(get(y)))
  m <- fit(dd, y, 'income3 + age_c + sex + region')
  bind_rows(lapply(c('income3Mais_1_a_2_SM', 'income3Mais_2_SM'),
                   function(tm) cbind(desfecho = y, termo = tm, est_row(m, tm))))
}))
wcsv(cat_pr, '07_PR_faixas_renda_blocoB.csv')

# PR e prevalência por quintil do ridit da coorte (bloco B), referência = Q1: descreve a forma.
q_pr <- bind_rows(lapply(outcomes, function(y) {
  dd <- subset(coorte, !is.na(get(y)))
  m <- fit(dd, y, 'income5 + age_c + sex + region')
  bind_rows(lapply(paste0('income5Q', 2:5), function(tm) cbind(desfecho = y, termo = tm, est_row(m, tm))))
}))
wcsv(q_pr, '07b_PR_quintis_blocoB.csv')
wcsv(bind_rows(lapply(outcomes, prev_by, dd = coorte, by = 'income5')), '03b_prevalencia_por_quintil.csv')

# ---------------------------------------------------------------------------
# 3. Forma do gradiente
# ---------------------------------------------------------------------------
forma <- bind_rows(lapply(outcomes, function(y) {
  dd <- subset(coorte, !is.na(get(y)))
  bind_rows(lapply(c('income3', 'income5'), function(cv) {
    m_assoc <- fit(dd, y, paste0(cv, ' + age_c + sex + region'))
    m_desvio <- fit(dd, y, paste0('income_rank_mental + ', cv, ' + age_c + sex + region'))
    rbind(cbind(desfecho = y, resolucao = cv, teste = 'associação (sem o rank)', wald(m_assoc, as.formula(paste0('~', cv)))),
          cbind(desfecho = y, resolucao = cv, teste = 'desvio de linearidade (condicional ao rank)', wald(m_desvio, as.formula(paste0('~', cv)))))
  }))
}))
wcsv(forma, '08_forma_do_gradiente.csv')

# ---------------------------------------------------------------------------
# 4. Interações pré-especificadas (bloco B): renda x esquizofrenia, renda x bipolar
# ---------------------------------------------------------------------------
inter_pre <- bind_rows(lapply(outcomes, function(y) {
  dd <- subset(coorte, !is.na(get(y)))
  bind_rows(lapply(c('schizophrenia', 'bipolar'), function(dx) {
    ddx <- subset(dd, !is.na(get(dx)))
    m <- fit(ddx, y, paste0('income_rank_mental * ', dx, ' + age_c + sex + region'))
    cbind(desfecho = y, diagnostico = dx, n = n_dom(ddx),
          casos = sum(ddx$variables[[dx]][in_dom(ddx)] == 1),
          wald(m, as.formula(paste0('~income_rank_mental:', dx))))
  }))
}))
wcsv(inter_pre, '09_interacoes_pre_especificadas.csv')

# ---------------------------------------------------------------------------
# 5. Comparação entre grupos hierárquicos (post hoc), blocos B e D
# ---------------------------------------------------------------------------
cov_B <- 'age_c + sex + region'
cov_D <- 'age_c + sex + region + race3 + education3 + urban + plan'
grupo_inter <- list(); grupo_rii <- list()
for (y in outcomes) for (b in c('B', 'D')) {
  dd <- subset(coorte, !is.na(get(y)) & !is.na(dx_hier))
  covs <- if (b == 'B') cov_B else cov_D
  m <- fit(dd, y, paste0('income_rank_mental * dx_hier + ', covs))
  grupo_inter[[paste(y, b)]] <- cbind(desfecho = y, bloco = b, n = n_dom(dd),
                                      wald(m, ~income_rank_mental:dx_hier))
  # RII de cada grupo a partir do mesmo modelo (efeitos das covariáveis compartilhados).
  for (lv in levels(d$dx_hier)) {
    if (inherits(m, 'error')) next
    cf <- c(income_rank_mental = 1)
    if (lv != levels(d$dx_hier)[1]) cf[paste0('income_rank_mental:dx_hier', lv)] <- 1
    full <- setNames(rep(0, length(coef(m))), names(coef(m))); full[names(cf)] <- cf
    sc <- svycontrast(m, full)
    b_ <- as.numeric(coef(sc)); se_ <- as.numeric(SE(sc)); cr <- qt(.975, m$df.residual)
    grupo_rii[[paste(y, b, lv)]] <- data.frame(desfecho = y, bloco = b, grupo = lv,
      n_grupo = sum(in_dom(dd) & dd$variables$dx_hier %in% lv),
      eventos_grupo = sum(in_dom(dd) & dd$variables$dx_hier %in% lv & dd$variables[[y]] %in% 1),
      RII = exp(b_), lower = exp(b_ - cr * se_), upper = exp(b_ + cr * se_),
      p = 2 * pt(abs(b_ / se_), m$df.residual, lower.tail = FALSE))
  }
}
wcsv(bind_rows(grupo_inter), '10_interacao_renda_x_grupo.csv')
wcsv(bind_rows(grupo_rii), '11_RII_por_grupo_modelo_interacao.csv')

# ---------------------------------------------------------------------------
# 6. Duplo gradiente dentro de cada grupo (blocos A-E, modelos separados)
# ---------------------------------------------------------------------------
blocos_grupo <- bind_rows(lapply(levels(d$dx_hier), function(lv) {
  dg <- subset(coorte, dx_hier == lv)
  bind_rows(lapply(outcomes, function(y) run_blocks(dg, y, paste0('grupo: ', lv))))
}))
wcsv(blocos_grupo, '12_blocos_por_grupo.csv')

# ---------------------------------------------------------------------------
# 7. Sensibilidades dentro dos grupos (bloco B)
# ---------------------------------------------------------------------------
sens_grupo <- bind_rows(lapply(levels(d$dx_hier), function(lv) {
  dg <- subset(coorte, dx_hier == lv)
  bind_rows(
    # (a) itens do próprio bloco do questionário
    run_blocks(dg, 'blk_psychotherapy', paste0('grupo: ', lv, ' | itens do próprio bloco'), block_list = blocks['B']),
    run_blocks(dg, 'blk_high_need_gap', paste0('grupo: ', lv, ' | itens do próprio bloco'), block_list = blocks['B']),
    # (b) ridit calculado dentro do grupo
    run_blocks(subset(dg, !is.na(rank_group)), 'psychotherapy', paste0('grupo: ', lv, ' | ridit do grupo'),
               rank_var = 'rank_group', block_list = blocks['B']),
    run_blocks(subset(dg, !is.na(rank_group)), 'high_need_gap', paste0('grupo: ', lv, ' | ridit do grupo'),
               rank_var = 'rank_group', block_list = blocks['B'])
  )
}))
wcsv(sens_grupo, '13_sensibilidades_por_grupo.csv')

# ---------------------------------------------------------------------------
# 8. Sensibilidades na coorte (bloco B, salvo indicação)
# ---------------------------------------------------------------------------
sens <- bind_rows(
  bind_rows(lapply(outcomes, function(y) run_blocks(subset(coorte, age >= 18), y, 'idade 18+', block_list = blocks['B']))),
  bind_rows(lapply(outcomes, function(y) run_blocks(subset(coorte, self_report == 1), y, 'só autorrespondentes', block_list = blocks['B']))),
  run_blocks(coorte, 'hng_lim4',     'high need gap: limitação >= 4',                    block_list = blocks['B']),
  run_blocks(coorte, 'hng_strict',   'high need gap: só não visita (sem episódico)',     block_list = blocks['B']),
  run_blocks(coorte, 'hng_episodic', 'high need gap: só episódico',                      block_list = blocks['B']),
  # As duas variantes acima dependem da distinção "não visita" x "só quando tem problema", a
  # parte menos concordante entre os blocos (kappa 0,58 nos 2.279 com os dois diagnósticos,
  # contra 0,73-0,84 nos demais componentes). Por isso são refeitas também só em quem
  # respondeu um único bloco, onde nenhuma regra de harmonização entra.
  run_blocks(subset(coorte, dx_group != 'ambos'), 'hng_strict',   'high need gap: só não visita | só quem respondeu um bloco', block_list = blocks['B']),
  run_blocks(subset(coorte, dx_group != 'ambos'), 'hng_episodic', 'high need gap: só episódico | só quem respondeu um bloco',  block_list = blocks['B']),
  bind_rows(lapply(outcomes, function(y) run_blocks(subset(coorte, !is.na(dx_hier)), y, 'ajustado por grupo diagnóstico',
                   block_list = list(B = paste0(blocks$B, ' + dx_hier'), D = paste0(blocks$D, ' + dx_hier')))))
)
wcsv(sens, '14_sensibilidades_coorte.csv')

# Família de modelo (bloco B): Poisson modificado (principal), log-binomial e logística (auxiliar, OR).
familia <- bind_rows(lapply(outcomes, function(y) {
  dd <- subset(coorte, !is.na(get(y)))
  mp <- fit(dd, y, blocks$B)
  st <- if (inherits(mp, 'error')) NULL else coef(mp)
  mb <- fit(dd, y, blocks$B, family = quasibinomial(link = 'log'), start = st)
  ml <- fit(dd, y, blocks$B, family = quasibinomial(link = 'logit'))
  bind_rows(cbind(desfecho = y, modelo = 'Poisson modificado (RII)', est_row(mp, 'income_rank_mental')),
            cbind(desfecho = y, modelo = 'log-binomial (RII)',       est_row(mb, 'income_rank_mental')),
            cbind(desfecho = y, modelo = 'logística (OR, auxiliar)', est_row(ml, 'income_rank_mental')))
}))
wcsv(familia, '15_familia_de_modelo.csv')

# ---------------------------------------------------------------------------
# 9. Checagem dos modelos com e sem ajuste por grupo diagnóstico (blocos B e D)
#    (a) quantas prevalências previstas pelo Poisson passam de 1, e de que perfil;
#    (b) o mesmo modelo por log-binomial, que por construção não prevê acima de 1.
#    Se o RII log-binomial ficar próximo do Poisson, as previsões > 1 não distorcem a estimativa.
# ---------------------------------------------------------------------------
checar_modelo <- function(y, ajuste, b) {
  dd  <- subset(coorte, !is.na(get(y)) & !is.na(dx_hier))
  rhs <- if (ajuste == 'ajustado por grupo') paste0(blocks[[b]], ' + dx_hier') else blocks[[b]]
  mp  <- fit(dd, y, rhs)
  if (inherits(mp, 'error')) stop('Poisson falhou: ', conditionMessage(mp))
  keep <- mp$prior.weights > 0
  f    <- mp$fitted.values[keep]
  acima <- f > 1
  # Início válido para o log-binomial: coeficientes do Poisson com o intercepto deslocado
  # até que todas as previsões fiquem abaixo de 0,99 (sem isso o log-binomial nem começa).
  st <- coef(mp); eta_max <- max(mp$linear.predictors[keep])
  if (eta_max > log(0.99)) st[['(Intercept)']] <- st[['(Intercept)']] - (eta_max - log(0.99))
  mb <- fit(dd, y, rhs, family = quasibinomial(link = 'log'), start = st,
            control = glm.control(maxit = 200))
  rp <- est_row(mp, 'income_rank_mental'); rb <- est_row(mb, 'income_rank_mental')
  linha <- data.frame(
    desfecho = y, ajuste = ajuste, bloco = b, n = sum(keep),
    n_previsto_acima_1 = sum(acima), pct_previsto_acima_1 = round(100 * mean(acima), 2),
    max_prev_poisson = max(f),
    RII_poisson = rp$est, lower_poisson = rp$lower, upper_poisson = rp$upper, p_poisson = rp$p,
    RII_logbin = rb$est, lower_logbin = rb$lower, upper_logbin = rb$upper, p_logbin = rb$p,
    logbin_convergiu = if (inherits(mb, 'error')) FALSE else isTRUE(mb$converged),
    max_prev_logbin = if (inherits(mb, 'error')) NA else max(mb$fitted.values[mb$prior.weights > 0]),
    nota_logbin = rb$note, erro = NA_character_)
  perfil <- NULL
  if (any(acima)) {
    # Perfil lido das variáveis do desenho, casando as linhas usadas no modelo pelo nome.
    rows <- match(rownames(model.frame(mp)), rownames(dd$variables))[keep][acima]
    q <- cut(dd$variables$income_rank_mental[rows], breaks = c(0, .2, .4, .6, .8, 1),
             labels = paste0('Q', 1:5), include.lowest = TRUE)
    tb <- as.data.frame(table(grupo = dd$variables$dx_hier[rows], quintil_renda = q,
                              plano = dd$variables$plan[rows], useNA = 'ifany'))
    perfil <- cbind(desfecho = y, ajuste = ajuste, bloco = b, tb[tb$Freq > 0, ])
  }
  list(linha = linha, perfil = perfil)
}

check_rows <- list(); perfil_rows <- list()
for (y in outcomes) for (ajuste in c('sem ajuste por grupo', 'ajustado por grupo')) for (b in c('B', 'D')) {
  cat('checagem:', y, '|', ajuste, '| bloco', b, '...\n'); flush.console()
  r <- tryCatch(checar_modelo(y, ajuste, b), error = function(e)
    list(linha = data.frame(desfecho = y, ajuste = ajuste, bloco = b, erro = conditionMessage(e)), perfil = NULL))
  check_rows[[paste(y, ajuste, b)]] <- r$linha
  if (!is.null(r$perfil)) perfil_rows[[paste(y, ajuste, b)]] <- r$perfil
}
checagem <- bind_rows(check_rows)
wcsv(checagem, '16_checagem_previsoes_e_logbinomial.csv')
if (length(perfil_rows)) wcsv(bind_rows(perfil_rows), '17_perfil_previsoes_acima_de_1.csv')
cat('CSV 16 salvo em', file.path(out, '16_checagem_previsoes_e_logbinomial.csv'),
    if (length(perfil_rows)) '| CSV 17 salvo (houve previsões > 1)' else '| sem previsões > 1: CSV 17 não é gerado', '\n')
if (any(!is.na(checagem$erro))) { cat('ATENÇÃO, modelos com erro:\n'); print(checagem[!is.na(checagem$erro), c('desfecho','ajuste','bloco','erro')]) }

# ---------------------------------------------------------------------------
# Resumo
# ---------------------------------------------------------------------------
cat('\n=== AMOSTRA ANALÍTICA ===\n')
cat('Coorte (mental_dx = 1):', sum(d$mental_dx %in% 1), '| amostra comum dos blocos:', n_dom(coorte),
    '| degf do desenho no domínio:', degf(coorte), '\n')
print(table(grupo = coorte$variables$dx_hier[in_dom(coorte)], useNA = 'ifany'))
cat('\n=== PREVALÊNCIAS ===\n'); print(prev, digits = 3)
cat('\n=== DUPLO GRADIENTE NA COORTE ===\n')
print(blocos_coorte[, c('desfecho','bloco','n','eventos','RII','RII_lower','RII_upper','RII_p','SII','df_residual','atenuacao_vs_B_pct')], digits = 3)
cat('\n=== PR POR FAIXA DE SM E POR QUINTIL (bloco B) ===\n'); print(bind_rows(cat_pr, q_pr)[, c('desfecho','termo','est','lower','upper','p')], digits = 3)
cat('\n=== FORMA DO GRADIENTE ===\n'); print(forma[, c('desfecho','resolucao','teste','F','p')], digits = 3)
cat('\n=== INTERAÇÕES PRÉ-ESPECIFICADAS ===\n'); print(inter_pre[, c('desfecho','diagnostico','casos','F','p')], digits = 3)
cat('\n=== RENDA x GRUPO (post hoc) ===\n'); print(bind_rows(grupo_inter)[, c('desfecho','bloco','n','F','p')], digits = 3)
cat('\n=== RII POR GRUPO (modelo com interação) ===\n'); print(bind_rows(grupo_rii), digits = 3)
cat('\n=== CHECAGEM: PREVISÕES > 1 E LOG-BINOMIAL ===\n')
print(checagem[, intersect(c('desfecho','ajuste','bloco','n_previsto_acima_1','max_prev_poisson','RII_poisson','lower_poisson','upper_poisson',
                   'RII_logbin','lower_logbin','upper_logbin','p_logbin','logbin_convergiu','max_prev_logbin','erro'), names(checagem))], digits = 3)
if (length(perfil_rows)) { cat('\n=== PERFIL DE QUEM TEM PREVISÃO > 1 ===\n'); print(bind_rows(perfil_rows)) }
instaveis <- blocos_grupo[!is.na(blocos_grupo$eventos_por_coef) & blocos_grupo$eventos_por_coef < 10, c('analise','desfecho','bloco','eventos','n_coef','eventos_por_coef')]
if (nrow(instaveis)) { cat('\n=== ATENÇÃO: menos de 10 eventos por coeficiente ===\n'); print(instaveis) }
cat('\nResultados em:', out, '\n')
