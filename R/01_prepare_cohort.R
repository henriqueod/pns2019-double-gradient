# =============================================================================
# 01_prepare_cohort.R
# Builds the analytic cohort (Q092 OR Q11006), harmonizes care items across questionnaire blocks and saves the survey design.
# Repository name: R/01_prepare_cohort.R   |   Development name: 01c_prepare_any_mental.R
# Code comments are in Portuguese. Cross-references in the comments use the
# development names; README.md maps them to the repository names.
# Paths: all input and output paths are relative to `root`, which defaults to
# the working directory (run from the repository root) or to the PNS_ROOT
# environment variable if it is set.
# =============================================================================

# 01c_prepare_any_mental.R
# Coorte "qualquer doença mental" = depressão (Q092) OU outra doença mental (Q11006).
# POR QUE SÓ A UNIÃO NÃO BASTA
# Os itens dos desfechos co-primários (Q11201, Q11405, Q115) só são perguntados a quem
# responde Sim na Q110a ("Se Q110a = 2, passe Q116a"). Quem tem só depressão tem esses
# itens como "não aplicável". Sem harmonização, o complete-case descartaria em silêncio
# os ~6 mil respondentes só com depressão.
# Por isso os desfechos são harmonizados com os itens paralelos do bloco de depressão:
#   continuidade  Q11201 <-> Q094    (1 regular; categorias 2 e 3 em ordem inversa entre os blocos: ver seção 2)
#   psicoterapia  Q11405 <-> Q09605  (1 sim / 2 não)
#   medicamento   Q11406 <-> Q09606
#   integrativas  Q11407 <-> Q09607
#   limitação     Q115   <-> Q109    (1 a 5)
#   Q11408 (acompanhamento regular com profissional) não tem par: fica definido só no bloco Q110a.
#
# REGRA PARA QUEM TEM OS DOIS DIAGNÓSTICOS (nível da pessoa)
#   psicoterapia / medicamento / integrativas = 1 se Sim em qualquer bloco aplicável;
#   continuidade regular = 1 se regular em qualquer bloco aplicável
#     (quem é acompanhado regularmente por um transtorno está em seguimento de saúde mental);
#   limitação = maior grau entre os blocos aplicáveis.
#   Um 0 (ou o grau de limitação) só é atribuído quando TODOS os blocos aplicáveis têm
#   resposta válida; um único Sim basta para 1. Qualquer outra combinação vira NA.
#
# Os nomes das variáveis derivadas (mental_dx, income_rank_mental, regular_care,
# psychotherapy, limitation, ...) são os usados pelos scripts a jusante.

options(timeout = 600)
options(survey.lonely.psu = 'adjust')
options(survey.adjust.domain.lonely = TRUE)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(survey)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# Funções de recodificação (testáveis isoladamente)
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

# Item binário 1 = Sim / 2 = Não, harmonizado entre os blocos aplicáveis.
harm_bin <- function(app_o, x_o, app_d, x_d) {
  case_when(
    (app_o & x_o %in% 1) | (app_d & x_d %in% 1)                         ~ 1,
    (app_o | app_d) & (!app_o | x_o %in% 2) & (!app_d | x_d %in% 2)     ~ 0,
    TRUE                                                                ~ NA_real_
  )
}

# Continuidade em 3 níveis: 1 regular; 3 só quando tem problema; 2 não visita.
# Regular em qualquer bloco -> 1. Sem regularidade em nenhum, com "só quando tem
# problema" em algum -> 3. "Não" em todos os aplicáveis -> 2.
harm_cont <- function(app_o, x_o, app_d, x_d) {
  all_valid_nonreg <- (app_o | app_d) & (!app_o | x_o %in% 2:3) & (!app_d | x_d %in% 2:3)
  any_episodic     <- (app_o & x_o %in% 3) | (app_d & x_d %in% 3)
  case_when(
    (app_o & x_o %in% 1) | (app_d & x_d %in% 1) ~ 1,
    all_valid_nonreg & any_episodic             ~ 3,
    all_valid_nonreg                            ~ 2,
    TRUE                                        ~ NA_real_
  )
}

# Limitação 1-5: maior grau entre os blocos aplicáveis; NA se algum bloco aplicável
# estiver sem resposta válida (o máximo seria desconhecido).
harm_lim <- function(app_o, x_o, app_d, x_d) {
  vo <- ifelse(app_o & x_o %in% 1:5, x_o, NA_real_)
  vd <- ifelse(app_d & x_d %in% 1:5, x_d, NA_real_)
  mx <- suppressWarnings(pmax(vo, vd, na.rm = TRUE))
  case_when(
    !(app_o | app_d)                            ~ NA_real_,
    (app_o & is.na(vo)) | (app_d & is.na(vd))   ~ NA_real_,
    TRUE                                        ~ as.numeric(mx)
  )
}

recode_any_mental <- function(pns) {
  num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  num_vars <- c(
    'V0001','V0024','C006','C008','C009','C011','V0026','VDD004A','VDE001','VDE002','VDF003','VDF004',
    'M001','M00203','M00302','M01401','M01501','M01601','M01701','M01801','M01901',
    'Q092','Q094','Q09605','Q09606','Q09607','Q109',
    'Q11006','Q11007','Q11008','Q11009','Q11010','Q111','Q11201','Q11405','Q11406','Q11407','Q11408','Q115',
    'I00102','J009','J01002','J014','J01502','J01602','J01701','J01802','J019','J02002','J021','J022010',
    'J023','J024','J025','J026','J02702','J02901','J03001','J03102','J034','J035','J03602'
  )
  for (v in num_vars) if (v %in% names(pns)) pns[[v]] <- num(pns[[v]])

  # --- Covariáveis ---
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
  pns$age_grp <- cut(pns$age, breaks = c(14,29,44,59,Inf),
                     labels = c('15-29','30-44','45-59','60+'), right = TRUE)
  pns$race <- factor(ifelse(pns$C009 %in% 1:5, pns$C009, NA_real_), levels = 1:5,
                     labels = c('Branca','Preta','Amarela','Parda','Indígena'))
  pns$race3 <- factor(ifelse(pns$C009 == 1, 'Branca',
                             ifelse(pns$C009 %in% c(2,3,4,5), 'Não_branca', NA_character_)),
                      levels = c('Branca','Não_branca'))
  pns$urban <- factor(pns$V0026, levels = c(1,2), labels = c('Urbana','Rural'))
  pns$education3 <- cut(pns$VDD004A, breaks = c(0,3,5,7),
                        labels = c('Baixa','Média','Alta'), right = TRUE)
  # Faixas de renda domiciliar per capita em salários mínimos (NÃO são tercis).
  pns$income3 <- cut(pns$VDF004, breaks = c(0,3,4,7),
                     labels = c('Até_1_SM','Mais_1_a_2_SM','Mais_2_SM'), right = TRUE)
  pns$income3 <- factor(pns$income3, levels = c('Até_1_SM','Mais_1_a_2_SM','Mais_2_SM'))
  pns$income_reais <- pns$VDF003
  pns$income_rank <- make_weighted_rank(pns$income_reais, pns$V00291)   # população adulta total

  pns$self_report <- ifelse(pns$M00203 == 1, 1, ifelse(pns$M00203 == 2, 0, NA_real_))
  pns$proxy_type <- factor(ifelse(pns$M00302 == 1, 'Outro_morador',
                                  ifelse(pns$M00302 == 2, 'Não_morador', NA_character_)))
  pns$usual_source <- ifelse(pns$J009 %in% c(1,2), ifelse(pns$J009 == 1, 1, 0), NA_real_)
  pns$recent_care_seek <- ifelse(pns$J014 %in% c(1,2), ifelse(pns$J014 == 1, 1, 0), NA_real_)
  pns$plan <- ifelse(pns$I00102 %in% c(1,2), ifelse(pns$I00102 == 1, 1, 0), NA_real_)
  pns$family_support0 <- ifelse(pns$M01401 == 0, 1, ifelse(pns$M01401 %in% 1:3, 0, NA_real_))
  pns$friend_support0 <- ifelse(pns$M01501 == 0, 1, ifelse(pns$M01501 %in% 1:3, 0, NA_real_))
  pns$social_isolation <- ifelse(pns$M01401 %in% 0:3 & pns$M01501 %in% 0:3,
                                 ifelse(pns$M01401 == 0 & pns$M01501 == 0, 1, 0), NA_real_)

  # --- 1. Coorte: união de Q092 e Q11006 ---
  # Um Sim em qualquer pergunta basta (quem tem os dois entra uma única vez).
  # Um 0 exige Não válido nas duas. Ignorado/NA sem nenhum Sim -> NA (não vira 0).
  pns <- pns %>%
    mutate(
      dep_dx = case_when(Q092 == 1 ~ 1, Q092 == 2 ~ 0, TRUE ~ NA_real_),
      oth_dx = case_when(Q11006 == 1 ~ 1, Q11006 == 2 ~ 0, TRUE ~ NA_real_),
      transtorno_mental = case_when(
        dep_dx == 1 | oth_dx == 1 ~ 1,
        dep_dx == 0 & oth_dx == 0 ~ 0,
        TRUE                      ~ NA_real_
      ),
      dx_group = case_when(
        dep_dx == 1 & oth_dx == 1 ~ 'ambos',
        dep_dx == 1 & oth_dx == 0 ~ 'so_depressao',
        dep_dx == 0 & oth_dx == 1 ~ 'so_outra',
        transtorno_mental == 1    ~ 'indeterminado',   # Sim numa, sem resposta válida na outra
        TRUE                      ~ NA_character_
      ),
      app_d = dep_dx %in% 1,   # bloco de depressão aplicável
      app_o = oth_dx %in% 1    # bloco Q110a aplicável
    )

  # --- 2. Desfechos harmonizados ---
  # ATENÇÃO: Q11201 e Q094 têm as categorias 2 e 3 em ordem INVERSA no questionário/dicionário:
  #   Q11201 (Q110a):     1 Sim | 2 Não                               | 3 Não, só quando tem algum problema
  #   Q094   (depressão): 1 Sim | 2 Não, só quando tem algum problema | 3 Nunca vai
  # Q094 é recodificado para a escala do Q11201 antes de harmonizar (2 <-> 3). A continuidade
  # regular (código 1) não é afetada; só a distinção "não visita" x "só quando tem problema".
  pns$Q094_h <- ifelse(pns$Q094 %in% 2, 3, ifelse(pns$Q094 %in% 3, 2, pns$Q094))
  pns <- pns %>%
    mutate(
      continuity3      = harm_cont(app_o, Q11201, app_d, Q094_h),
      regular_care     = case_when(continuity3 == 1 ~ 1, continuity3 %in% 2:3 ~ 0, TRUE ~ NA_real_),
      unmet_continuity = 1 - regular_care,
      episodic_care    = case_when(continuity3 == 3 ~ 1, continuity3 %in% 1:2 ~ 0, TRUE ~ NA_real_),
      psychotherapy    = harm_bin(app_o, Q11405, app_d, Q09605),
      medication       = harm_bin(app_o, Q11406, app_d, Q09606),
      integrative      = harm_bin(app_o, Q11407, app_d, Q09607),
      limitation       = harm_lim(app_o, Q115,   app_d, Q109),
      regular_professional = ifelse(app_o & Q11408 %in% 1:2, as.numeric(Q11408 == 1), NA_real_)
    )

  # --- 3. Diagnósticos ---
  # Q11007-Q11010 são saltados quando Q11006 = 2. Para quem disse Não à Q110a, a ausência
  # desses diagnósticos é informação, não dado faltante: codifica-se 0. Sem isso, os
  # modelos de interação renda x esquizofrenia/bipolar perderiam em silêncio todos os
  # respondentes só com depressão.
  pns <- pns %>%
    mutate(
      mental_dx     = transtorno_mental,     # nome usado pelos scripts a jusante
      depression    = dep_dx,
      schizophrenia = case_when(Q11007 == 1 ~ 1, Q11007 == 2 ~ 0, oth_dx == 0 ~ 0, TRUE ~ NA_real_),
      bipolar       = case_when(Q11008 == 1 ~ 1, Q11008 == 2 ~ 0, oth_dx == 0 ~ 0, TRUE ~ NA_real_),
      ocd           = case_when(Q11009 == 1 ~ 1, Q11009 == 2 ~ 0, oth_dx == 0 ~ 0, TRUE ~ NA_real_),
      other_mental  = case_when(Q11010 == 1 ~ 1, Q11010 == 2 ~ 0, oth_dx == 0 ~ 0, TRUE ~ NA_real_)
    )
  # Contagem de diagnósticos agora inclui depressão (muda a definição de multimorbidade).
  pns$n_mental_diagnoses <- rowSums(cbind(pns$depression, pns$schizophrenia, pns$bipolar,
                                          pns$ocd, pns$other_mental), na.rm = TRUE)
  pns$multiple_mental_dx <- ifelse(pns$mental_dx %in% 1, as.numeric(pns$n_mental_diagnoses >= 2), NA_real_)

  # --- 3b. Hierarquia diagnóstica (grupos mutuamente exclusivos) ---
  # Transtorno mental grave (esquizofrenia OU bipolar) > depressão sem transtorno grave >
  # só outras doenças mentais (sem depressão e sem transtorno grave: TOC e "outro",
  # onde caem ansiedade, pânico e psicoses sem subitem próprio).
  # Um grupo só é atribuído quando a informação necessária para excluir os níveis acima
  # é válida; caso contrário, NA.
  pns <- pns %>%
    mutate(
      smi = case_when(
        schizophrenia == 1 | bipolar == 1 ~ 1,
        schizophrenia == 0 & bipolar == 0 ~ 0,
        TRUE                              ~ NA_real_
      ),
      dx_hier = case_when(
        !(mental_dx %in% 1)                         ~ NA_character_,
        smi == 1                                    ~ 'grave',
        smi == 0 & depression == 1                  ~ 'depressao',
        smi == 0 & depression == 0 & oth_dx == 1    ~ 'outras',
        TRUE                                        ~ NA_character_
      ),
      # Referência = depressão (maior grupo, estimativas mais estáveis).
      dx_hier = factor(dx_hier, levels = c('depressao', 'grave', 'outras'))
    )

  # --- 4. Rankings de renda ---
  # Principal: ridit ponderado dentro da coorte.
  pns$income_rank_mental <- NA_real_
  i <- which(pns$mental_dx %in% 1 & !is.na(pns$income_reais) & !is.na(pns$V00291))
  pns$income_rank_mental[i] <- make_weighted_rank(pns$income_reais[i], pns$V00291[i])
  pns$income_rank_schizophrenia <- NA_real_
  i <- which(pns$schizophrenia %in% 1 & !is.na(pns$income_reais) & !is.na(pns$V00291))
  pns$income_rank_schizophrenia[i] <- make_weighted_rank(pns$income_reais[i], pns$V00291[i])

  # --- 5. Itens de origem regravados com o valor harmonizado ---
  # Alguns scripts a jusante leem os itens brutos (ex.: sensibilidades Q11201 = 2 e
  # Q115 >= 4). Se esses itens continuassem só com as respostas do bloco Q110a, esses
  # scripts ignorariam o bloco de depressão. Por isso Q11201, Q115, Q11405, Q11406 e
  # Q11407 passam a conter o valor harmonizado, no MESMO código do questionário; os
  # originais ficam em orig_*. Q11408 não tem par e não é alterado. Q11006 também não:
  # qualquer script que filtre por Q11006 == 1 em vez de mental_dx == 1 precisa ser trocado.
  for (v in c('Q11201','Q115','Q11405','Q11406','Q11407')) pns[[paste0('orig_', v)]] <- pns[[v]]
  to12 <- function(x) ifelse(x %in% 1, 1, ifelse(x %in% 0, 2, NA_real_))
  pns$Q11201 <- pns$continuity3
  pns$Q115   <- pns$limitation
  pns$Q11405 <- to12(pns$psychotherapy)
  pns$Q11406 <- to12(pns$medication)
  pns$Q11407 <- to12(pns$integrative)

  pns
}

# Concordância entre os dois blocos entre quem tem os dois diagnósticos (descritiva,
# não ponderada): mostra se os itens paralelos medem a mesma coisa.
kappa_simple <- function(a, b, levels) {
  ok <- a %in% levels & b %in% levels
  t <- table(factor(a[ok], levels), factor(b[ok], levels))
  n <- sum(t)
  if (n == 0) return(c(n = 0, agreement = NA, kappa = NA))
  po <- sum(diag(t)) / n
  pe <- sum(rowSums(t) * colSums(t)) / n^2
  c(n = n, agreement = po, kappa = if (pe < 1) (po - pe) / (1 - pe) else NA)
}

audit_any_mental <- function(pns) {
  m <- pns[pns$mental_dx %in% 1, ]
  flow <- data.frame(
    item = c('selecionados', 'Q092 = 1 (depressão)', 'Q11006 = 1 (outra doença mental)',
             'ambos', 'so_depressao', 'so_outra', 'indeterminado',
             'coorte (transtorno_mental = 1)', 'transtorno_mental NA',
             'hierarquia: grave (esquizofrenia ou bipolar)', 'hierarquia: depressao (sem grave)',
             'hierarquia: outras (sem grave, sem depressão)', 'hierarquia: NA dentro da coorte',
             'grave: esquizofrenia', 'grave: bipolar', 'grave: ambos'),
    n = c(nrow(pns), sum(pns$dep_dx %in% 1), sum(pns$oth_dx %in% 1),
          sum(pns$dx_group %in% 'ambos'), sum(pns$dx_group %in% 'so_depressao'),
          sum(pns$dx_group %in% 'so_outra'), sum(pns$dx_group %in% 'indeterminado'),
          sum(pns$transtorno_mental %in% 1), sum(is.na(pns$transtorno_mental)),
          sum(pns$dx_hier %in% 'grave'), sum(pns$dx_hier %in% 'depressao'),
          sum(pns$dx_hier %in% 'outras'), sum(pns$mental_dx %in% 1 & is.na(pns$dx_hier)),
          sum(pns$schizophrenia %in% 1), sum(pns$bipolar %in% 1),
          sum(pns$schizophrenia %in% 1 & pns$bipolar %in% 1))
  )
  hier_x_depression <- as.data.frame(table(dx_hier = pns$dx_hier[pns$mental_dx %in% 1],
                                           depressao = pns$depression[pns$mental_dx %in% 1],
                                           useNA = 'ifany'))
  # Respostas não válidas nos itens de origem, por bloco aplicável (checa os saltos).
  item_na <- data.frame(
    bloco = c(rep('Q110a (app_o)', 3), rep('depressão (app_d)', 3)),
    item  = c('Q11201','Q11405','Q115','Q094','Q09605','Q109'),
    n_aplicavel = c(rep(sum(pns$app_o), 3), rep(sum(pns$app_d), 3)),
    n_sem_resposta_valida = c(
      sum(pns$app_o & !(pns$orig_Q11201 %in% 1:3)), sum(pns$app_o & !(pns$orig_Q11405 %in% 1:2)),
      sum(pns$app_o & !(pns$orig_Q115 %in% 1:5)),
      sum(pns$app_d & !(pns$Q094 %in% 1:3)),   sum(pns$app_d & !(pns$Q09605 %in% 1:2)),
      sum(pns$app_d & !(pns$Q109 %in% 1:5)))
  )
  # NA dos desfechos harmonizados por grupo diagnóstico (o que o complete-case vai perder).
  outcome_na <- m %>%
    group_by(dx_group) %>%
    summarise(n = n(),
              na_regular_care  = sum(is.na(regular_care)),
              na_psychotherapy = sum(is.na(psychotherapy)),
              na_limitation    = sum(is.na(limitation)),
              na_high_need_gap_inputs = sum(is.na(regular_care) | is.na(limitation)),
              .groups = 'drop')
  amb <- pns[pns$dx_group %in% 'ambos', ]
  concord <- rbind(
    continuidade_regular = kappa_simple(as.numeric(amb$orig_Q11201 == 1), as.numeric(amb$Q094 == 1), 0:1),
    continuidade_3niveis = kappa_simple(amb$orig_Q11201, amb$Q094_h, 1:3),   # Q094 já na escala do Q11201
    psicoterapia         = kappa_simple(amb$orig_Q11405, amb$Q09605, 1:2),
    limitacao_alta_ge3   = kappa_simple(as.numeric(amb$orig_Q115 >= 3), as.numeric(amb$Q109 >= 3), 0:1)
  )
  concord <- data.frame(item = rownames(concord), concord, row.names = NULL)
  list(flow = flow, hier_x_depression = hier_x_depression, item_na = item_na, outcome_na = outcome_na, concordance_ambos = concord,
       crosstab_continuidade = as.data.frame(table(Q11201 = amb$orig_Q11201, Q094_recodificado = amb$Q094_h, useNA = 'ifany')),
       crosstab_limitacao    = as.data.frame(table(Q115 = amb$orig_Q115, Q109 = amb$Q109, useNA = 'ifany')),
       crosstab_psicoterapia = as.data.frame(table(Q11405 = amb$orig_Q11405, Q09605 = amb$Q09605, useNA = 'ifany')))
}

# ---------------------------------------------------------------------------
# Execução principal
# ---------------------------------------------------------------------------
if (!isTRUE(getOption('anymh.test', FALSE))) {
  library(PNSIBGE)

  root     <- Sys.getenv('PNS_ROOT', unset = getwd())
  data_dir <- file.path(root, 'data_pnsibge')
  out_dir  <- file.path(root, 'results', 'any_mental')
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Desenho amostral salvo para os scripts 02 a 08.
  out_rds <- file.path(root, 'results', 'mental_care_design_anymh.rds')

  vars <- c(
    'V0001','V0024','UPA_PNS','V0006_PNS','V0015','V0025A','V0029','V00291','V00292','V00293',
    'C006','C008','C009','C011','V0026','VDD004A','VDE001','VDE002','VDF003','VDF004',
    'M001','M00203','M00302','M01401','M01501','M01601','M01701','M01801','M01901',
    'Q092','Q094','Q09605','Q09606','Q09607','Q109',
    'Q11006','Q11007','Q11008','Q11009','Q11010','Q111','Q11201','Q11405','Q11406','Q11407','Q11408','Q115',
    'I00102',
    'J009','J01002','J014','J01502','J01602','J01701','J01802','J019','J02002','J021','J022010',
    'J023','J024','J025','J026','J02702','J02901','J03001','J03102','J034','J035','J03602','J060'
  )
  pns <- PNSIBGE::get_pns(year = 2019, selected = TRUE, vars = vars, labels = FALSE,
                          deflator = FALSE, design = FALSE, reload = FALSE, savedir = data_dir)

  # Toda a recodificação acontece ANTES de declarar o desenho.
  pns <- recode_any_mental(pns)
  aud <- audit_any_mental(pns)
  for (nm in names(aud)) write.csv(aud[[nm]], file.path(out_dir, paste0('audit_', nm, '.csv')), row.names = FALSE)

  des <- PNSIBGE::pns_design(pns)

  # Composição ponderada da coorte por grupo diagnóstico (linha nova da Tabela 1).
  coorte <- subset(des, mental_dx == 1 & !is.na(dx_group))
  comp <- svymean(~factor(dx_group), coorte, na.rm = TRUE)
  comp_tab <- data.frame(grupo = sub('factor\\(dx_group\\)', '', names(coef(comp))),
                         pct = 100 * as.numeric(coef(comp)), se = 100 * as.numeric(SE(comp)))
  write.csv(comp_tab, file.path(out_dir, 'cohort_composition_weighted.csv'), row.names = FALSE)
  comp_h <- svymean(~dx_hier, subset(des, mental_dx == 1 & !is.na(dx_hier)), na.rm = TRUE)
  comp_hier <- data.frame(grupo = sub('dx_hier', '', names(coef(comp_h))),
                          pct = 100 * as.numeric(coef(comp_h)), se = 100 * as.numeric(SE(comp_h)))
  write.csv(comp_hier, file.path(out_dir, 'cohort_hierarchy_weighted.csv'), row.names = FALSE)
  prev_pop <- svymean(~transtorno_mental, subset(des, !is.na(transtorno_mental)))

  saveRDS(des, out_rds, compress = 'xz')

  cat('\n=== FLUXO (conferir: ambos 2.279; só depressão 6.053; coorte 10.958) ===\n'); print(aud$flow)
  cat('\n=== ITENS SEM RESPOSTA VÁLIDA NO BLOCO APLICÁVEL (checa saltos) ===\n');   print(aud$item_na)
  cat('\n=== NA DOS DESFECHOS HARMONIZADOS POR GRUPO ===\n');                        print(aud$outcome_na)
  cat('\n=== CONCORDÂNCIA ENTRE BLOCOS (quem tem os dois diagnósticos) ===\n');      print(aud$concordance_ambos, digits = 3)
  cat('\n=== COMPOSIÇÃO PONDERADA DA COORTE ===\n');                                 print(comp_tab, digits = 3)
  cat('\n=== HIERARQUIA DIAGNÓSTICA (ponderada) ===\n');                             print(comp_hier, digits = 3)
  cat('\n=== HIERARQUIA x DEPRESSÃO (não ponderada) ===\n');                         print(aud$hier_x_depression)
  cat('\nPrevalência populacional ponderada de qualquer diagnóstico:', round(100 * coef(prev_pop), 2), '%\n')
  cat('Graus de liberdade do desenho no domínio da coorte:', degf(subset(des, mental_dx == 1)), '\n')
  cat('Desenho salvo em:', out_rds, '\n')
}
