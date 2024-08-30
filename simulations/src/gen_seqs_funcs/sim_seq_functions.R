##################################
# Functions for generating the simulated sequences
###################################

get_n_order_tmat <- function(dat, n){
  #####################
  # Get estimated TPM from a set of sequences
  # dat: matrix whose columns are individual sequeces
  # n: order of Markov chain
  ####################
  aux <- apply(dat, 2, function(col) {
    from <- head(apply(embed(col, n)[, n:1], 1, paste, collapse = ""), -1)
    to <- col[-1:-n]
    rbind(from, to)
  })
  aux <- data.frame(t(matrix(aux, nrow = 2)))
  names(aux) <- c("From", "To")
  TM <- table(aux)
  TM <- TM / rowSums(TM)
  return(TM)
}

################################################################
### The three following functions all simulate a single sequence from a TPM 
### seed: integer, the random generation seed
### P: TPM to generate the sequence from
### num.iters: length of the sequences
### start_state: integer, starting state of the sequences 
###############################################################
##### first order #####
run.mc.sim <- function(seed, P, num.iters = 500, start_state = NULL) {
  # number of possible states
  num.states <- nrow(P)
  
  # stores the states X_t through time
  states <- numeric(num.iters)
  
  # initialize variable for first state
  if(is.null(start_state) == FALSE){
    states[1] <- start_state
  } else if(is.null(start_probs) == FALSE){
    set.seed(seed*num.iters + 1)
    states[1] <- which(rmultinom(1, 1, start_probs) == 1)
  }
  
  for(t in 2:num.iters) {
    set.seed(seed*num.iters + t)
    # probability vector to simulate next state X_{t+1}
    p  <- P[states[t-1], ]
    
    # draw from multinomial and determine state
    states[t] <-  which(rmultinom(1, 1, p) == 1)
  }
  return(states)
}

##### second order #####
run.mc.sim2 <- function(seed, P, num.iters = 500, start_state = NULL) {
  # number of possible states
  num.states <- ncol(P)
  state_dict <- colnames(P)
  
  # stores the states X_t through time
  states <- numeric(num.iters)
  
  states[1:2] <- start_state
  
  for(t in 3:num.iters) {
    set.seed(seed*num.iters + t)
    row_ind <- paste0(state_dict[states[t-2]], state_dict[states[t-1]])
    # probability vector to simulate next state X_{t+1}
    p  <- P[row_ind, ]
    
    # draw from multinomial and determine state
    states[t] <-  which(rmultinom(1, 1, p) == 1)
  }
  return(states)
}

##### fifth order #####
run.mc.sim5 <- function(seed, P, num.iters = 500, start_state = NULL) {
  # number of possible states
  num.states <- ncol(P)
  state_dict <- colnames(P)
  
  # stores the states X_t through time
  states <- numeric(num.iters)
  
  states[1:5] <- start_state
  
  for(t in 6:num.iters) {
    set.seed(seed*num.iters + t)
    row_ind <- paste0(state_dict[states[t-5]], state_dict[states[t-4]], state_dict[states[t-3]], state_dict[states[t-2]], state_dict[states[t-1]])
    #work around if the last 5 generate states are not present in the TPM
    #use the most recent set of 5 previous states that have a row in the TPM
    num = 1
    while(!(row_ind %in% rownames(P))){
      row_ind <- paste0(state_dict[states[t-5 - num]], state_dict[states[t-4 -num]], state_dict[states[t-3-num]], state_dict[states[t-2-num]], state_dict[states[t-1-num]])
      num <- num + 1
    }
    
    # probability vector to simulate next state X_{t+1}
    p  <- P[row_ind, ]
    
    # draw from multinomial and determine state
    states[t] <-  which(rmultinom(1, 1, p) == 1)
  }
  return(states)
}


gen_clust_seq <- function(seed, P_list, order = 1, num.iters = 500, start_state = NULL, clust_prob){
  ################################################################
  ### Simulate a sequences with respect a set of potential cluster
  ### seed: integer, the random generation seed
  ### P: list of TPMs where each TPM is generates sequences for a single cluster
  ### order: input 1,2,5 for first, second, fifth order Markov chain TPMs
  ### num.iters: length of the sequence
  ### start_state: list of integers, each entry is the starting state for the cluster
  ### clust_prob: list of probability of being generated from each cluster (must sum to 1)
  ###############################################################
  num_clust <- length(clust_prob)
  set.seed(seed)
  #select cluster to generate from
  j <- which(rmultinom(1, 1, clust_prob) == 1) 
  
  #simulate the sequences
  if(order == 1){
    sim_seq_vec <- run.mc.sim(seed, P_list[[j]], num.iters, start_state = start_state[[j]])
  } else if (order == 2){
    sim_seq_vec <- run.mc.sim2(seed, P_list[[j]], num.iters, start_state = start_state[[j]])
  } else if (order ==5){
    sim_seq_vec <- run.mc.sim5(seed, P_list[[j]], num.iters, start_state = start_state[[j]])
  }
  
  sim_seq_vec <- ifelse(sim_seq_vec == 1, 'A',
                        ifelse(sim_seq_vec == 2, 'B',
                               ifelse(sim_seq_vec == 3, 'C', 'D')))
  sim_seq <- paste0(sim_seq_vec, collapse = '')
  return(c('seqs' = sim_seq, 'cluster' = as.character(j)))
}


state_change <- function(mat, z_row_ind, z_col_ind){
  #######################################################
  # Helper function for generating the state-based simulation scenario TPMs
  # mat: starting TPM
  # z_row_ind: the rows in the TPM to be set to zero
  # x_col_ind: the columns in the TPM to be set to zero
  ######################################################
  mat[z_row_ind,] <- 0
  nz_row_ind <- seq(1, nrow(mat))[-z_row_ind]
  nz_col_ind <- seq(1, 4)[-z_col_ind]
  for(k in nz_row_ind){
    distr_num <- sum(mat[k, z_col_ind])
    #specifics for certain scenarios with high duration states so all states are reasonably present in the resulting sequences
    if((length(intersect(c(1,4), nz_col_ind)) == 2) & (length(nz_col_ind) == 2)){
      mat[k, nz_col_ind] <- mat[k, nz_col_ind] + c(distr_num*0.25, distr_num*0.75)
    } else if((length(intersect(c(1,4), nz_col_ind)) == 2) & (length(nz_col_ind) == 3)){
      mat[k, nz_col_ind] <- mat[k, nz_col_ind] + c(distr_num*0, distr_num*0.75, distr_num*0.25)
    } else{
      mat[k, nz_col_ind] <- mat[k, nz_col_ind] + rep(distr_num/length(nz_col_ind), length(nz_col_ind))
    }
    mat[k, z_col_ind] <- 0
  }
  return(mat)
}





