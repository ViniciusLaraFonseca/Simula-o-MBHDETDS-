# Nome: ValoresIniciais_gamma_FIXO_Lambda_it.R
# Objetivo: Simular dados (T=23) com Gamma e Epsilon FIXOS no tempo.

# Carregar dados originais
if(file.exists("_dataCaseStudy.r")) {
  source("_dataCaseStudy.r")
} else {
  stop("Erro: Arquivo '_dataCaseStudy.r' não encontrado.")
}

set.seed(987)
cat("--- Iniciando Simulação (Gamma/Epsilon FIXOS, T=23) ---\n")

# --- PASSO 1: DEFINIR PARÂMETROS E CLUSTERS ---
n_regions <- 75 
n_times   <- 23 
p         <- 3   
K         <- 4   

# Parâmetros Escalares
beta_true  <- c(-1.0, 1.0, 0.5) 
li_w = 0.7
ls_w = 0.99
w_true     <- runif(75, li_w, ls_w)
a0_true    <- 1.0
b0_true    <- 1.0

# --- CORREÇÃO DOS CLUSTERS ---
# Cluster 1 (Melhor) ... Cluster 4 (Pior)
cluster_ids <- 5 - clAI 

# Matriz h CUMULATIVA (Fixa para as regiões)
h_cumulativo <- matrix(0, nrow = n_regions, ncol = K)
for(i in 1:n_regions) {
  k <- cluster_ids[i]
  h_cumulativo[i, 1:k] <- 1
}

# --- PASSO 2: DEFINIR GAMMAS (FIXOS) ---
# Agora é um vetor de tamanho K, não mais uma matriz K x T
gamma_true <- c(0.05, 0.10, 0.10, 0.15) 
# Cluster 1: 0.05
# Cluster 2: 0.05 + 0.10 = 0.15
# Cluster 3: 0.15 + 0.10 = 0.25
# Cluster 4: 0.25 + 0.15 = 0.40

# --- PASSO 3: SIMULAR EPSILON (FIXO POR REGIÃO) ---
cat("--- Calculando Epsilon Fixo ---\n")

epsilon_true <- numeric(n_regions)

for(i in 1:n_regions) {
  # Penalidade é fixa no tempo
  penalidade <- sum(h_cumulativo[i, ] * gamma_true)
  epsilon_true[i] <- 1 - penalidade
}

# Simulando Covariáveis e E
x_true <- array(rnorm(n_regions * n_times * p), dim = c(n_regions, n_times, p))
E_raw  <- matrix(runif(n_regions * n_times, 150, 250), nrow = n_regions)
E_true <- E_raw / mean(E_raw) 

# Calcular g_it (Agora usa epsilon_true[i] fixo)
g_it_true <- array(NA, dim = c(n_regions, n_times))
for(i in 1:n_regions) {
  for(t in 1:n_times) {
    prod_val <- sum(x_true[i, t, ] * beta_true)
    # Epsilon não tem índice t
    g_it_true[i, t] <- E_true[i, t] * epsilon_true[i] * exp(prod_val)
  }
}

# --- PASSO 4: SIMULAR LAMBDA E Y ---
cat("--- Simulando trajetória lambda_it e Y_it ---\n")

lambda_true <- matrix(NA, nrow = n_regions, ncol = n_times)
Y_ini       <- matrix(NA, nrow = n_regions, ncol = n_times)

at_true  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
bt_true  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
att_true <- matrix(NA, nrow = n_regions, ncol = n_times)
btt_true <- matrix(NA, nrow = n_regions, ncol = n_times)

for(i in 1:n_regions) {
  at_true[i, 1] <- a0_true
  bt_true[i, 1] <- b0_true
  
  for(t in 1:n_times) {
    att_true[i, t] <- w_true[i] * at_true[i, t]
    btt_true[i, t] <- w_true[i] * bt_true[i, t]
    
    lambda_true[i, t] <- rgamma(1, shape = att_true[i, t], rate = btt_true[i, t])
    
    mu_it <- lambda_true[i, t] * g_it_true[i, t]
    Y_ini[i, t] <- rpois(1, mu_it)
    
    at_true[i, t+1] <- att_true[i, t] + Y_ini[i, t]
    bt_true[i, t+1] <- btt_true[i, t] + g_it_true[i, t]
  }
}

# --- PASSO 5: PREPARAR OBJETOS PARA O NIMBLE ---

constants_nimble <- list(
  n_regions = n_regions, 
  n_times   = n_times, 
  p         = p, 
  K         = K, 
  h         = h_cumulativo,
  mu_beta   = rep(0, p),
  # Limites para o primeiro gamma (escalares agora)
  a_unif    = 0.0, 
  b_unif    = 0.1,  
  a0        = a0_true, 
  b0        = b0_true, 
  li_w      = li_w,
  ls_w      = ls_w
)

data_nimble <- list(
  Y = Y_ini, 
  E = E_true, 
  x = x_true
)

# Inits agora passam vetores para gamma, não matrizes
inits_nimble_1 <- list(
  beta   = beta_true,
  gamma  = gamma_true,
  lambda = lambda_true,
  w      = w_true
)

gamma_init_2 <- gamma_true * 0.9
inits_nimble_2 <- list(
  beta   = rnorm(p, 0, 0.5),
  gamma  = gamma_init_2,
  lambda = matrix(rgamma(n_regions * n_times, 1, 1), nrow = n_regions),
  w      = runif(n_regions, li_w, ls_w)
)

inits_list_nimble <- list(inits_nimble_1, inits_nimble_2)

cat("--- Simulação Concluída (Gamma Fixo) ---\n")
print("Epsilon por Cluster (Valores Únicos):")
print(tapply(epsilon_true, cluster_ids, mean))