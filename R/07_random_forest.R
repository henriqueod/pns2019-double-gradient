# =============================================================================
# 07_random_forest.R
# Exploratory random-forest check (Supplementary Section S6).
# Repository name: R/07_random_forest.R   |   Development name: 15b_forest_anymh.R
# Code comments are in Portuguese. Cross-references in the comments use the
# development names; README.md maps them to the repository names.
# Paths: all input and output paths are relative to `root`, which defaults to
# the working directory (run from the repository root) or to the PNS_ROOT
# environment variable if it is set.
# =============================================================================

# 15b_forest_anymh.R
# Floresta aleatória exploratória (Material Suplementar) na coorte "qualquer doença mental"
# (Q092 OU Q11006), a partir do desenho do 01c.
#   Desfecho: high need gap harmonizado (limitação >= 3 e sem seguimento regular em nenhum
#     bloco), lido do mental_care_design_anymh.rds.
#   Preditores incluem depressão (Q092); multiple_mental_dx conta depressão (definição do 01c).
#   Saídas: validação por partição, resumo (média, DP, mínimo e máximo), importância média
#     nas 5 partições e importância da floresta da primeira partição.
#
# Análise exploratória (TRIPOD+AI); não é ferramenta clínica nem de triagem. Plano de saúde, fonte usual, busca recente, medicação e psicoterapia
# ficam fora dos preditores, e nenhum componente do desfecho entra como preditor.
# Partição treino/teste por UPA (70/30), 5 repetições, pesos amostrais como case.weights.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ranger))
cat('07_random_forest.R (15b_forest_anymh.R) | versão 26/09/2026\n')

root <- Sys.getenv('PNS_ROOT', unset = getwd())
des  <- readRDS(file.path(root, 'results', 'mental_care_design_anymh.rds'))
out  <- file.path(root, 'results', 'forest_anymh')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- des$variables

d$high_limitation <- ifelse(is.na(d$limitation), NA_real_, as.numeric(d$limitation >= 3))
d$care_gap <- ifelse(is.na(d$regular_care) | is.na(d$high_limitation), NA_real_,
                     as.numeric(d$regular_care == 0 & d$high_limitation == 1))

vars <- c('care_gap', 'age', 'income_rank_mental', 'sex', 'race3', 'education3', 'urban', 'region',
          'depression', 'schizophrenia', 'bipolar', 'ocd', 'other_mental', 'multiple_mental_dx',
          'family_support0', 'friend_support0', 'self_report')

mental_rows <- which(d$mental_dx %in% 1)
psu_all <- unique(as.character(d$UPA_PNS[mental_rows])); psu_all <- psu_all[!is.na(psu_all)]

weighted_auc <- function(y, p, w) {
  pos <- which(y == 1); neg <- which(y == 0)
  cmp <- outer(p[pos], p[neg], '>') + .5 * outer(p[pos], p[neg], '==')
  sum(cmp * outer(w[pos], w[neg])) / (sum(w[pos]) * sum(w[neg]))
}
as_target <- function(x) factor(ifelse(x == 1, 'high_gap', 'other'), levels = c('other', 'high_gap'))

particao <- function(b) {
  set.seed(20260827 + b)
  trpsu <- sample(psu_all, round(.70 * length(psu_all)))
  split <- ifelse(as.character(d$UPA_PNS[mental_rows]) %in% trpsu, 'train', 'test')
  tr <- mental_rows[split == 'train']; te <- mental_rows[split == 'test']
  list(tr = tr[complete.cases(d[tr, vars])], te = te[complete.cases(d[te, vars])])
}

res <- list(); imps <- list()
for (b in 1:5) {
  cat('partição', b, 'de 5...\n'); flush.console()
  p_ <- particao(b)
  train <- d[p_$tr, vars]; test <- d[p_$te, vars]
  train$care_gap <- as_target(train$care_gap); test$care_gap <- as_target(test$care_gap)
  fit <- ranger(care_gap ~ ., data = train, num.trees = 500, mtry = 4, min.node.size = 100,
                probability = TRUE, importance = 'permutation', case.weights = d$V00291[p_$tr],
                respect.unordered.factors = 'partition', seed = 20260827 + b)
  p <- predict(fit, data = test)$predictions[, 'high_gap']
  y <- as.numeric(test$care_gap == 'high_gap'); w <- d$V00291[p_$te]
  q <- quantile(p, seq(0, 1, .2)); qid <- cut(p, breaks = unique(q), include.lowest = TRUE, labels = FALSE)
  top <- qid == max(qid, na.rm = TRUE)
  res[[b]] <- data.frame(rep = b, n_train = length(p_$tr), n_test = length(p_$te),
                         prevalence = weighted.mean(y, w), auc = weighted_auc(y, p, w),
                         brier = sum(w * (y - p)^2) / sum(w),
                         top20_lift = weighted.mean(y[top], w[top]) / weighted.mean(y, w))
  imps[[b]] <- data.frame(rep = b, variable = names(fit$variable.importance),
                          importance = as.numeric(fit$variable.importance))
}
res <- do.call(rbind, res); imps <- do.call(rbind, imps)

resumo <- data.frame(metric = c('auc', 'brier', 'top20_lift'),
                     mean = sapply(res[c('auc', 'brier', 'top20_lift')], mean),
                     sd   = sapply(res[c('auc', 'brier', 'top20_lift')], sd),
                     min  = sapply(res[c('auc', 'brier', 'top20_lift')], min),
                     max  = sapply(res[c('auc', 'brier', 'top20_lift')], max), row.names = NULL)
imp_media <- aggregate(importance ~ variable, data = imps, FUN = mean)
imp_media$sd <- aggregate(importance ~ variable, data = imps, FUN = sd)$importance
imp_media <- imp_media[order(-imp_media$importance), ]

# Importância da floresta ajustada na primeira partição, com 800 árvores (como no 15).
p1 <- particao(0)
train1 <- d[p1$tr, vars]; train1$care_gap <- as_target(train1$care_gap)
fit1 <- ranger(care_gap ~ ., data = train1, num.trees = 800, mtry = 4, min.node.size = 100,
               probability = TRUE, importance = 'permutation', case.weights = d$V00291[p1$tr],
               respect.unordered.factors = 'partition', seed = 20260827)
imp1 <- data.frame(variable = names(fit1$variable.importance), importance = as.numeric(fit1$variable.importance))
imp1 <- imp1[order(-imp1$importance), ]

write.csv(res,       file.path(out, '49b_validacao_floresta_anymh.csv'), row.names = FALSE)
write.csv(resumo,    file.path(out, '49b_resumo_validacao_anymh.csv'), row.names = FALSE)
write.csv(imp1,      file.path(out, '50b_importancia_floresta_anymh.csv'), row.names = FALSE)
write.csv(imp_media, file.path(out, '50c_importancia_media_5_particoes_anymh.csv'), row.names = FALSE)

cat('\n=== VALIDAÇÃO (5 partições por UPA) ===\n'); print(res, digits = 3)
cat('\n=== RESUMO ===\n'); print(resumo, digits = 3)
cat('\n=== IMPORTÂNCIA POR PERMUTAÇÃO (floresta da 1ª partição, 800 árvores) ===\n'); print(imp1, digits = 3, row.names = FALSE)
cat('\n=== IMPORTÂNCIA MÉDIA NAS 5 PARTIÇÕES ===\n'); print(imp_media, digits = 3, row.names = FALSE)
cat('\nResultados em:', out, '\n')
