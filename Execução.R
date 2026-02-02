# ==============================================================================
# 0. INÍCIO
# ==============================================================================
inicio_execucao <- Sys.time()

library(nimble)
library(coda)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

cat("--- Carregando Valores Iniciais e Dados ---\n")
if(file.exists("ValoresIniciais.R")) {
  source("ValoresIniciais.R")
} else {
  stop("ERRO: Arquivo 'ValoresIniciais.R' não encontrado.")
}

# ==============================================================================
# 1. AMOSTRADOR FFBS (Ajustado para Epsilon Fixo)
# ==============================================================================

ffbs_sampler_Lambda_it <- nimbleFunction(
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    n_regions <- control$n_regions
    n_times   <- control$n_times
    p         <- control$p
    a0        <- control$a0
    b0        <- control$b0
    
    buf_size <- n_regions * (n_times + 1)
    at_buf <- nimNumeric(buf_size, 0)
    bt_buf <- nimNumeric(buf_size, 0)
    
    calcNodes   <- model$getDependencies(target, self = FALSE)
    targetNodes <- model$expandNodeNames(target)
    
    setupOutputs(at_buf, bt_buf)
  },
  
  run = function() {
    declare(i, integer())
    declare(t, integer())
    declare(k, integer())
    declare(prod_val, double())
    declare(g_it, double())
    declare(att_t, double())
    declare(btt_t, double())
    declare(shape_tmp, double())
    declare(rate_tmp, double())
    declare(lambda_futuro, double())
    declare(nu, double())
    declare(idx, integer())
    declare(idx_next, integer())
    
    for(i in 1:n_regions) {
      
      # 1. FORWARD FILTERING
      idx <- (i-1)*(n_times+1) + 1
      at_buf[idx] <<- a0
      bt_buf[idx] <<- b0
      
      for(t in 1:n_times) {
        idx      <- (i-1)*(n_times+1) + t
        idx_next <- idx + 1
        
        att_t <- model$w[i] * at_buf[idx]
        btt_t <- model$w[i] * bt_buf[idx]
        
        prod_val <- 0
        for(k in 1:p) prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
        
        # MUDANÇA AQUI: epsilon[i] não tem índice t
        g_it <- model$E[i, t] * model$epsilon[i] * exp(prod_val)
        
        at_buf[idx_next] <<- att_t + model$Y[i, t]
        bt_buf[idx_next] <<- btt_t + g_it
      }
      
      # 2. BACKWARD SAMPLING
      idx <- (i-1)*(n_times+1) + n_times + 1
      shape_tmp <- at_buf[idx]
      rate_tmp  <- bt_buf[idx]
      model$lambda[i, n_times] <<- rgamma(1, shape = shape_tmp, rate = rate_tmp)
      
      for(t_idx in 1:(n_times-1)) {
        t_back  <- n_times - t_idx
        idx_buf <- (i-1)*(n_times+1) + t_back + 1
        lambda_futuro <- model$lambda[i, t_back + 1]
        
        shape_tmp <- (1 - model$w[i]) * at_buf[idx_buf]
        rate_tmp  <- bt_buf[idx_buf]
        
        nu <- rgamma(1, shape = shape_tmp, rate = rate_tmp)
        model$lambda[i, t_back] <<- nu + model$w[i] * lambda_futuro
      }
    }
    
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
  },
  methods = list(reset = function() {})
)

# ==============================================================================
# 2. CÓDIGO DO MODELO (Gamma Fixo)
# ==============================================================================

code_static_gamma <- nimbleCode({
  
  # Priors Beta
  for (j in 1:p) {
    beta[j] ~ dnorm(mu_beta[j], sd = 10)
  }
  
  # Priors W (CORRIGIDO: for numérico)
  for (i in 1:n_regions) {
    w[i] ~ dunif(min = li_w, max = ls_w)
  }
  
  # Priors Gamma (Sem loop de tempo)
  gamma[1] ~ dunif(min = a_unif, max = b_unif)
  for(j in 2:K){
    gamma[j] ~ dunif(min = 0, max = (1 - sum(gamma[1:(j-1)])))
  }
  
  # Likelihood
  for (i in 1:n_regions) {
    
    # Epsilon Fixo (Fora do loop de tempo)
    epsilon[i] <- 1 - inprod(h[i, 1:K], gamma[1:K])
    
    for(t in 1:n_times){
      lambda[i, t] ~ dgamma(1, 1) # Dummy prior para FFBS
      
      # mu usa epsilon[i] fixo
      mu[i, t] <- lambda[i, t] * E[i, t] * epsilon[i] * exp(inprod(beta[1:p], x[i, t, 1:p]))
      
      Y[i, t] ~ dpois(mu[i, t])
      logLik_Y[i, t] <- log(dpois(Y[i, t], mu[i, t]))
    }
  }
})

# ==============================================================================
# 3. MCMC
# ==============================================================================

cat("--- Configurando MCMC ---\n")
model <- nimbleModel(code = code_static_gamma, constants = constants_nimble, 
                     data = data_nimble, inits = inits_list_nimble[[1]])

conf <- configureMCMC(model)
conf$removeSamplers("lambda")
conf$addSampler(target = "lambda", type = ffbs_sampler_Lambda_it,
                control = list(n_regions = constants_nimble$n_regions,
                               n_times = constants_nimble$n_times,
                               p = constants_nimble$p,
                               a0 = constants_nimble$a0, b0 = constants_nimble$b0))

# MUDANÇA: Adicionado "w" nos monitores
conf$addMonitors(c("beta", "gamma", "lambda", "epsilon", "logLik_Y", "w"))

cat("--- Compilando ---\n")
Cmodel <- compileNimble(model)
Rmcmc  <- buildMCMC(conf)
Cmcmc  <- compileNimble(Rmcmc, project = model)

cat("--- Rodando MCMC (20.000 iterações / 5.000 burn-in) ---\n")
samples_nimble <- runMCMC(Cmcmc, niter = 20000, nburnin = 5000, nchains = 2, 
                          inits = inits_list_nimble, samplesAsCodaMCMC = TRUE)

cat("--- MCMC Concluído ---\n")

# ==============================================================================
# 4. PÓS-PROCESSAMENTO E MÉTRICAS
# ==============================================================================

# Definir Clusters (Mapeamento: clAI=4 -> Cluster 1, clAI=1 -> Cluster 4)
clAI_data <- c(rep(4, 17), rep(3, 14), rep(2, 16), rep(1, 28))
cluster_mapping <- 5 - clAI_data
all_regions_info <- data.frame(Region = 1:75, Cluster = cluster_mapping)

samples_mat <- as.matrix(samples_nimble)
metrics_list <- list()

calc_metrics_relative <- function(posterior_vec, true_val) {
  mean_est <- mean(posterior_vec)
  bias_abs <- mean_est - true_val
  mse_abs  <- bias_abs^2 + var(posterior_vec)
  denom <- ifelse(abs(true_val) < 1e-12, 1e-12, true_val)
  rel_bias <- bias_abs / denom
  rel_mse  <- mse_abs / (denom^2)
  rel_rmse <- sqrt(rel_mse)
  ci       <- quantile(posterior_vec, probs = c(0.025, 0.975))
  coverage <- (true_val >= ci[1] && true_val <= ci[2])
  return(c(Mean = mean_est, True = true_val, RelBias = rel_bias, RelMSE = rel_mse, RelRMSE = rel_rmse, Cov95 = as.integer(coverage)))
}

cat("Calculando Métricas...\n")

# Betas
for(j in 1:constants_nimble$p) {
  p_name <- paste0("beta[", j, "]")
  if(p_name %in% colnames(samples_mat)) {
    metrics_list[[length(metrics_list) + 1]] <- c(Parameter = p_name, Type = "Beta", calc_metrics_relative(samples_mat[, p_name], beta_true[j]))
  }
}

# Gammas (Agora Fixo - Sem loop de tempo)
for(k in 1:constants_nimble$K) {
  p_name <- paste0("gamma[", k, "]")
  if(p_name %in% colnames(samples_mat)) {
    # gamma_true agora é um vetor, usamos gamma_true[k] direto
    metrics_list[[length(metrics_list) + 1]] <- c(Parameter = p_name, Type = "Gamma", calc_metrics_relative(samples_mat[, p_name], gamma_true[k]))
  }
}

# Lambda/Epsilon por Cluster
for(k in 1:4) {
  regs <- all_regions_info$Region[all_regions_info$Cluster == k]
  count <- 0; l_bias <- 0; l_mse <- 0; l_cov <- 0; e_bias <- 0; e_mse <- 0; e_cov <- 0
  
  for(i in regs) {
    # Métricas de Lambda (Temporal)
    for(t in 1:constants_nimble$n_times) {
      nm_l <- paste0("lambda[", i, ", ", t, "]")
      if(nm_l %in% colnames(samples_mat)) {
        res_l <- calc_metrics_relative(samples_mat[, nm_l], lambda_true[i, t])
        l_bias <- l_bias + res_l["RelBias"]; l_mse <- l_mse + res_l["RelMSE"]; l_cov <- l_cov + res_l["Cov95"]
        count <- count + 1 # Conta amostras totais de lambda (Regiões * Tempos)
      }
    }
    
    # Métricas de Epsilon (Fixo - Uma por região)
    nm_e <- paste0("epsilon[", i, "]")
    if(nm_e %in% colnames(samples_mat)) {
      res_e <- calc_metrics_relative(samples_mat[, nm_e], epsilon_true[i])
      # Como epsilon é fixo por região, somamos aqui e dividimos pelo numero de regioes no final
      e_bias <- e_bias + res_e["RelBias"]; e_mse <- e_mse + res_e["RelMSE"]; e_cov <- e_cov + res_e["Cov95"]
    }
  }
  
  if(count > 0) {
    # Para Lambda dividimos pelo count total (N_regs * T)
    metrics_list[[length(metrics_list)+1]] <- c(Parameter=paste0("Lambda_Cl_", k), Type="Lambda_Cl", Mean=NA, True=NA, RelBias=l_bias/count, RelMSE=l_mse/count, RelRMSE=sqrt(l_mse/count), Cov95=l_cov/count)
    
    # Para Epsilon dividimos pelo número de regiões no cluster (pois só tem 1 epsilon por região)
    n_regs_k <- length(regs)
    metrics_list[[length(metrics_list)+1]] <- c(Parameter=paste0("Epsilon_Cl_", k), Type="Epsilon_Cl", Mean=NA, True=NA, RelBias=e_bias/n_regs_k, RelMSE=e_mse/n_regs_k, RelRMSE=sqrt(e_mse/n_regs_k), Cov95=e_cov/n_regs_k)
  }
}

df_metrics <- do.call(rbind, metrics_list) %>% as.data.frame()
cols_num <- c("Mean", "True", "RelBias", "RelMSE", "RelRMSE", "Cov95")
df_metrics[cols_num] <- lapply(df_metrics[cols_num], as.numeric)
write.csv(df_metrics, "metricas_finais_relativas_clusters.csv", row.names = FALSE)
print(df_metrics)

# WAIC/LPML (CORRIGIDO PARA logLik_Y)
logProb_cols <- grep("logLik_Y", colnames(samples_mat))
if(length(logProb_cols) > 0) {
  logProb_vals <- samples_mat[, logProb_cols]
  cpo_inv <- colMeans(exp(-logProb_vals)); lpml <- sum(log(1/cpo_inv))
  lppd <- sum(log(colMeans(exp(logProb_vals)))); p_waic <- sum(apply(logProb_vals, 2, var))
  waic <- -2 * (lppd - p_waic)
  write.csv(data.frame(Metric=c("WAIC","LPML"), Value=c(waic,lpml)), "metricas_ajuste.csv", row.names=FALSE)
  cat(sprintf("WAIC: %.4f | LPML: %.4f\n", waic, lpml))
} else {
  cat("Aviso: colunas logLik_Y não encontradas para WAIC.\n")
}

# ==============================================================================
# 5. GERAÇÃO DE TODOS OS GRÁFICOS
# ==============================================================================

if(!dir.exists("plots_output")) dir.create("plots_output")
cat("Gerando gráficos em 'plots_output'...\n")

# A. PAINEL LAMBDA (Amostra de Regiões)
regions_of_interest <- c(1, 8, 15, 19, 22, 31, 34, 40, 46, 55, 65, 75)
roi_info <- all_regions_info[all_regions_info$Region %in% regions_of_interest, ]
roi_info <- roi_info[order(roi_info$Cluster, roi_info$Region), ]
roi_info$Label <- paste0("Região ", roi_info$Region, " (Cl ", roi_info$Cluster, ")")

df_lambda <- data.frame()

for(i in roi_info$Region) {
  lbl <- roi_info$Label[roi_info$Region == i]
  for(t in 1:constants_nimble$n_times) {
    # Lambda
    nm <- paste0("lambda[", i, ", ", t, "]")
    if(nm %in% colnames(samples_mat)) {
      vec <- samples_mat[, nm]
      df_lambda <- rbind(df_lambda, data.frame(Region=i, Label=lbl, Time=t, True=lambda_true[i,t], Est=mean(vec), Lower=quantile(vec,0.025), Upper=quantile(vec,0.975)))
    }
  }
}
df_lambda$Label <- factor(df_lambda$Label, levels=roi_info$Label)

p1 <- ggplot(df_lambda, aes(x=Time)) + geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="grey70", alpha=0.5) +
  geom_line(aes(y=Est), color="black") + geom_line(aes(y=True), color="red", linetype="dashed") +
  facet_wrap(~Label, ncol=3, scales="fixed") + theme_bw() + labs(title="Painel Lambda", y=expression(lambda))
ggsave("plots_output/painel_lambda.png", p1, width=12, height=10)

# B. TRACEPLOT GAMMA (Alterado para Traceplot)
df_gamma_trace <- data.frame()
for(k in 1:constants_nimble$K) {
  nm <- paste0("gamma[", k, "]")
  if(nm %in% colnames(samples_mat)) {
    df_gamma_trace <- rbind(df_gamma_trace, data.frame(Iter=1:nrow(samples_mat), Value=samples_mat[,nm], Param=nm, True=gamma_true[k]))
  }
}

p3 <- ggplot(df_gamma_trace, aes(x=Iter, y=Value)) + geom_line(alpha=0.6, size=0.3) +
  geom_hline(aes(yintercept=True), color="red", linetype="dashed") +
  facet_grid(Param~., scales="free_y") + theme_bw() + labs(title="Traceplots Gamma")
ggsave("plots_output/traceplots_gamma.png", p3, width=8, height=6)

# C. TRACEPLOTS BETA
df_beta <- data.frame()
for(j in 1:constants_nimble$p) {
  nm <- paste0("beta[", j, "]")
  if(nm %in% colnames(samples_mat)) {
    df_beta <- rbind(df_beta, data.frame(Iter=1:nrow(samples_mat), Value=samples_mat[,nm], Param=nm, True=beta_true[j]))
  }
}
p4 <- ggplot(df_beta, aes(x=Iter, y=Value)) + geom_line(alpha=0.6, size=0.3) +
  geom_hline(aes(yintercept=True), color="red", linetype="dashed") +
  facet_grid(Param~., scales="free_y") + theme_bw() + labs(title="Traceplots Beta")
ggsave("plots_output/traceplots_beta.png", p4, width=8, height=6)

# D. W (ESPACIAL) - MANTIDO
# 1. Extração e Cálculo de Métricas por Região
cat("Calculando métricas para os 75 w's...\n")

df_w_results <- data.frame(
  Region = 1:75,
  True = w_true,  # Do script ValoresIniciais.R
  Est_Mean = NA, Est_Lower = NA, Est_Upper = NA, Bias = NA, MSE = NA, Coverage95 = NA
)

for(i in 1:75) {
  param_name <- paste0("w[", i, "]")
  if(param_name %in% colnames(samples_mat)) {
    samp_vec <- samples_mat[, param_name]
    mean_val <- mean(samp_vec)
    ci       <- quantile(samp_vec, probs = c(0.025, 0.975))
    true_val <- w_true[i]
    
    df_w_results$Est_Mean[i]  <- mean_val
    df_w_results$Est_Lower[i] <- ci[1]
    df_w_results$Est_Upper[i] <- ci[2]
    df_w_results$Bias[i]      <- mean_val - true_val
    df_w_results$MSE[i]       <- (mean_val - true_val)^2 + var(samp_vec)
    df_w_results$Coverage95[i] <- as.integer(true_val >= ci[1] & true_val <= ci[2])
  }
}

# 2. Métricas Gerais (Resumo Global)
global_metrics <- data.frame(
  Metric = c("Average Bias", "Global RMSE", "Average MSE", "Coverage Rate (95%)"),
  Value = c(
    mean(df_w_results$Bias, na.rm=TRUE),
    sqrt(mean(df_w_results$MSE, na.rm=TRUE)),
    mean(df_w_results$MSE, na.rm=TRUE),
    mean(df_w_results$Coverage95, na.rm=TRUE)
  )
)
print("--- Métricas Globais para w ---")
print(global_metrics)
write.csv(global_metrics, "metricas_globais_w.csv", row.names = FALSE)

# 3. Gráficos de Diagnóstico
# Gráfico A: Estimado vs Verdadeiro
p_w1 <- ggplot(df_w_results, aes(x = True, y = Est_Mean)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_errorbar(aes(ymin = Est_Lower, ymax = Est_Upper), color = "grey70", width = 0.005) +
  geom_point(color = "blue", size = 2, alpha = 0.7) +
  labs(title = "Recuperação do w: Estimado vs Verdadeiro",
       subtitle = "Barras cinzas indicam IC 95%",
       x = "w Verdadeiro", y = "w Estimado (Média a Posteriori)") + theme_bw()

# Gráfico B: Viés por Região
p_w2 <- ggplot(df_w_results, aes(x = Region, y = Bias)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_segment(aes(x = Region, xend = Region, y = 0, yend = Bias), color = "grey50") +
  geom_point(aes(color = abs(Bias)), size = 2) +
  scale_color_gradient(low = "blue", high = "red", name = "|Viés|") +
  labs(title = "Viés do w por Região", x = "Índice da Região", y = "Viés (Média - True)") + theme_bw()

# Gráfico C: EQM (MSE) por Região
p_w3 <- ggplot(df_w_results, aes(x = Region, y = MSE)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  labs(title = "Erro Quadrático Médio (MSE) por Região", x = "Índice da Região", y = "MSE") + theme_bw()

ggsave("plots_output/w_estimado_vs_true.png", p_w1, width = 8, height = 6)
ggsave("plots_output/w_vies_por_regiao.png", p_w2, width = 10, height = 5)
ggsave("plots_output/w_mse_por_regiao.png", p_w3, width = 10, height = 5)
# Exibir painel combinado
grid.arrange(p_w1, p_w2, p_w3, ncol = 1)
# ==============================================================================
# 6. FINALIZAÇÃO
# ==============================================================================
fim_execucao <- Sys.time()
tempo_total <- fim_execucao - inicio_execucao
cat("\n==================================================\n")
cat(" TEMPO TOTAL: \n")
print(tempo_total)
cat("==================================================\n")