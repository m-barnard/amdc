source('src/gen_seqs_funcs/sim_seq_functions.R')
source('src/run_functions/new_mat_funcs.R')


get_sim_TPMs <- function(order){
  ##############################################################
  # Get the TPMs for generating simulated sequences
  # order: input 1,2,5 for first, second, fifth order Markov chain TPMs
  # returns list of the cluster TPMs for all simulation scenarios
  ##############################################################
  
  ### get the estimated TPM from Daynamica data
  if(order == 1){
    df <- read_csv('data/clean_day_seqs.csv', show_col_types = FALSE)
    indices <- list(c(1))
    type <- c('keep')
    trans_res <- final_gen(df$seqs, indices, type)
    mats <- trans_res$mats$`1_keep`
    total_mat <- Reduce('+', mats)
    sums <- rowSums(total_mat)
    avg_TPM <- t(sapply(seq(1, 4), function(x){total_mat[x,]/sums[x]}))
    rownames(avg_TPM) <- c('A', 'B', 'C', 'D')
  } else if((order == 2) | (order == 5)){
    df <- read_csv('data/clean_day_seqs.csv', show_col_types = FALSE)
    seq_data <- sapply(df$seqs, function(x){str_split(x, '')[[1]]})
    colnames(seq_data) <- paste0('s',seq(1, 2229))
    avg_TPM <- get_n_order_tmat(seq_data, order)
    name_switch <- c(H = 'A', O = 'B',W = 'D', `T` = 'C')
    colnames(avg_TPM) <- str_replace_all(colnames(avg_TPM), name_switch)
    rownames(avg_TPM) <- str_replace_all(rownames(avg_TPM), name_switch)
  }
  
  ### STATE LOW ###
  zero_ind1 <- which(str_detect(rownames(avg_TPM), 'C') == T | str_detect(rownames(avg_TPM), 'D') == T)
  state_low1 <- state_change(avg_TPM, zero_ind1, c(3,4))
  zero_ind2 <- which(str_detect(rownames(avg_TPM), 'A') == T | str_detect(rownames(avg_TPM), 'B') == T)
  state_low2 <- state_change(avg_TPM, zero_ind2, c(1,2))
  P_state_low <- list(state_low1, state_low2)
  
  ### STATE MED ###
  zero_ind1 <- which(str_detect(rownames(avg_TPM), 'D') == T)
  state_med1 <- state_change(avg_TPM, zero_ind1, c(4))
  zero_ind2 <- which(str_detect(rownames(avg_TPM), 'C') == T | str_detect(rownames(avg_TPM), 'B') == T)
  state_med2 <- state_change(avg_TPM, zero_ind2, c(2,3))
  zero_ind3 <- which(str_detect(rownames(avg_TPM), 'A') == T | str_detect(rownames(avg_TPM), 'B') == T)
  state_med3 <- state_change(avg_TPM, zero_ind3, c(1,2))
  P_state_med <- list(state_med1, state_med2, state_med3)
  
  ### STATE HIGH ###
  zero_ind1 <- which(str_detect(rownames(avg_TPM), 'D') == T)
  state_high1 <- state_change(avg_TPM, zero_ind1, c(4))
  zero_ind2 <- which(str_detect(rownames(avg_TPM), 'A') == T)
  state_high2 <- state_change(avg_TPM, zero_ind2, c(1))
  zero_ind3 <- which(str_detect(rownames(avg_TPM), 'C') == T)
  state_high3 <- state_change(avg_TPM, zero_ind3, c(3))
  zero_ind4 <- which(str_detect(rownames(avg_TPM), 'B') == T)
  state_high4 <- state_change(avg_TPM, zero_ind4, c(2))
  P_state_high <- list(state_high1, state_high2, state_high3, state_high4)
  
  ### DUR1 LOW ###
  dur1_low1 <- avg_TPM
  dur1_low2 <- avg_TPM
  dur1_low2[paste0(rep('A', order), collapse =''),] <- dur1_low2[paste0(rep('A', order), collapse =''),] + c(-0.02, 0, 0, 0.02)
  dur1_med3 <- avg_TPM
  dur1_low3 <-avg_TPM
  dur1_low3[paste0(rep('A', order), collapse =''),] <- dur1_low3[paste0(rep('A', order), collapse =''),] + c(-0.08, 0, 0, .08)
  P_dur1_low <- list(dur1_low1, dur1_low2, dur1_low3)
  
  ### DUR1 MED ###
  dur1_med1 <- avg_TPM
  dur1_med2 <- avg_TPM
  dur1_med2[paste0(rep('A', order), collapse =''),] <- dur1_med2[paste0(rep('A', order), collapse =''),] + c(-0.015, 0, 0, 0.015)
  dur1_med3 <- avg_TPM
  dur1_med3[paste0(rep('A', order), collapse =''),] <- dur1_med3[paste0(rep('A', order), collapse =''),] + c(-0.05, 0, 0, 0.05)
  P_dur1_med <- list(dur1_med1, dur1_med2, dur1_med3)
  
  ### DUR1 HIGH ###
  dur1_high1 <- avg_TPM
  dur1_high2 <- avg_TPM
  dur1_high2[paste0(rep('A', order), collapse =''),] <- dur1_high2[paste0(rep('A', order), collapse =''),] + c(-0.01, 0, 0, 0.01)
  dur1_high3 <- avg_TPM
  dur1_high3[paste0(rep('A', order), collapse =''),] <- dur1_high3[paste0(rep('A', order), collapse =''),] +  c(-0.02, 0, 0, 0.02)
  P_dur1_high <- list(dur1_high1, dur1_high2, dur1_high3)
  
  ### DUR2 LOW ###
  dur2_low1 <- avg_TPM
  dur2_low1[paste0(rep('A', order), collapse =''),] <- dur2_low1[paste0(rep('A', order), collapse =''),] + c(-0.08, 0, 0, 0.08)
  dur2_low2 <- avg_TPM
  dur2_low2[paste0(rep('D', order), collapse =''),] <- dur2_low2[paste0(rep('D', order), collapse =''),] + c(0.08, 0, 0, -0.08)
  dur2_low3 <- dur2_low2
  dur2_low3[paste0(rep('A', order), collapse =''),] <- dur2_low3[paste0(rep('A', order), collapse =''),] + c(-0.08, 0, 0, 0.08)
  dur2_low3[paste0(rep('D', order), collapse =''),] <- dur2_low3[paste0(rep('D', order), collapse =''),] + c(0.08, 0, 0, -0.08)
  P_dur2_low <- list(dur2_low1, dur2_low2, dur2_low3)
  
  ### DUR2 MED ###
  dur2_med1 <- avg_TPM
  dur2_med1[paste0(rep('A', order), collapse =''),] <- dur2_med1[paste0(rep('A', order), collapse =''),] + c(-0.05, 0, 0, 0.05)
  dur2_med2 <- avg_TPM
  dur2_med2[paste0(rep('D', order), collapse =''),] <- dur2_med2[paste0(rep('D', order), collapse =''),] +  c(0.05, 0, 0, -0.05)
  dur2_med3 <- dur2_med2
  dur2_med3[paste0(rep('A', order), collapse =''),] <- dur2_med3[paste0(rep('A', order), collapse =''),] + c(-0.05, 0, 0, 0.05)
  dur2_med3[paste0(rep('D', order), collapse =''),] <- dur2_med3[paste0(rep('D', order), collapse =''),] + c(0.05, 0, 0, -0.05)
  P_dur2_med <- list(dur2_med1, dur2_med2, dur2_med3)
  
  ### DUR2 HIGH ###
  dur2_high1 <- avg_TPM
  dur2_high1[paste0(rep('A', order), collapse =''),] <- dur2_high1[paste0(rep('A', order), collapse =''),] +  c(-0.02, 0, 0, 0.02)
  dur2_high2 <- avg_TPM
  dur2_high2[paste0(rep('D', order), collapse =''),] <- dur2_high2[paste0(rep('D', order), collapse =''),] +  c(0.02, 0, 0, -0.02)
  dur2_high3 <- dur2_high2
  dur2_high3[paste0(rep('A', order), collapse =''),] <- dur2_high3[paste0(rep('A', order), collapse =''),] + c(-0.02, 0, 0, 0.02)
  dur2_high3[paste0(rep('D', order), collapse =''),] <- dur2_high3[paste0(rep('D', order), collapse =''),] + c(0.02, 0, 0, -0.02)
  P_dur2_high <- list(dur2_high1, dur2_high2, dur2_high3)
  
  ### DUR3 LOW ###
  dur3_low1 <- avg_TPM
  dur3_low1[paste0(rep('A', order), collapse =''),] <- dur3_low1[paste0(rep('A', order), collapse =''),] + c(-0.08, .08/2, 0, .08/2)*2
  dur3_low1[paste0(rep('B', order), collapse =''),] <- dur3_low1[paste0(rep('B', order), collapse =''),] + c(.08/2, -0.08, 0, .08/2)
  dur3_low2 <- avg_TPM
  dur3_low2[paste0(rep('B', order), collapse =''),] <- dur3_low2[paste0(rep('B', order), collapse =''),] +  c(.08/2, -0.08, 0, .08/2)*2
  dur3_low2[paste0(rep('D', order), collapse =''),] <- dur3_low2[paste0(rep('D', order), collapse =''),] +  c(.08/2, .08/2, 0, -0.08)
  dur3_low3 <- avg_TPM
  dur3_low3[paste0(rep('A', order), collapse =''),] <- dur3_low3[paste0(rep('A', order), collapse =''),] + c(-0.08, .08/2, 0, .08/2)
  dur3_low3[paste0(rep('D', order), collapse =''),] <- dur3_low3[paste0(rep('D', order), collapse =''),] + c(.08/2, .08/2, 0, -0.08)*2
  P_dur3_low <- list(dur3_low1, dur3_low2, dur3_low3)
  
  ### DUR3 MED ###
  dur3_med1 <- avg_TPM
  dur3_med1[paste0(rep('A', order), collapse =''),] <- dur3_med1[paste0(rep('A', order), collapse =''),] + c(-0.05, .05/2, 0, .05/2)*2
  dur3_med1[paste0(rep('B', order), collapse =''),] <- dur3_med1[paste0(rep('B', order), collapse =''),] + c(.05/2, -0.05, 0, .05/2)
  dur3_med2 <- avg_TPM
  dur3_med2[paste0(rep('B', order), collapse =''),] <- dur3_med2[paste0(rep('B', order), collapse =''),] + c(.05/2, -0.05, 0, .05/2)*2
  dur3_med2[paste0(rep('D', order), collapse =''),] <- dur3_med2[paste0(rep('D', order), collapse =''),] + c(.05/2, .05/2, 0, -0.05)
  dur3_med3 <- avg_TPM
  dur3_med3[paste0(rep('A', order), collapse =''),] <- dur3_med3[paste0(rep('A', order), collapse =''),] + c(-0.05, .05/2, 0, .05/2)
  dur3_med3[paste0(rep('D', order), collapse =''),] <- dur3_med3[paste0(rep('D', order), collapse =''),] + c(.05/2, .05/2, 0, -0.05)*2
  P_dur3_med <- list(dur3_med1, dur3_med2, dur3_med3)
  
  ### DUR3 HIGH ###
  dur3_high1 <- avg_TPM
  dur3_high1[paste0(rep('A', order), collapse =''),] <- dur3_high1[paste0(rep('A', order), collapse =''),] + c(-0.02, 0.02/2, 0, 0.02/2)*2
  dur3_high1[paste0(rep('B', order), collapse =''),] <- dur3_high1[paste0(rep('B', order), collapse =''),] + c(0.02/2, -0.02, 0, 0.02/2)
  dur3_high2 <- avg_TPM
  dur3_high2[paste0(rep('B', order), collapse =''),] <- dur3_high2[paste0(rep('B', order), collapse =''),] + c(0.02/2, -0.02, 0, 0.02/2)*2
  dur3_high2[paste0(rep('D', order), collapse =''),] <- dur3_high2[paste0(rep('D', order), collapse =''),] + c(0.02/2, 0.02/2, 0,-0.02)
  dur3_high3 <- avg_TPM
  dur3_high3[paste0(rep('A', order), collapse =''),] <- dur3_high3[paste0(rep('A', order), collapse =''),] + c(-0.02, 0.02/2, 0, 0.02/2)
  dur3_high3[paste0(rep('D', order), collapse =''),] <- dur3_high3[paste0(rep('D', order), collapse =''),] + c(0.02/2, 0.02/2, 0,-0.02)*2
  P_dur3_high <- list(dur3_high1, dur3_high2, dur3_high3)
  
  
  return(list('state_low' = P_state_low, 'state_med' = P_state_med, 'state_high' = P_state_high,
              'dur1_low' = P_dur1_low, 'dur1_med' = P_dur1_med, 'dur1_high' = P_dur1_high,
              'dur2_low' = P_dur2_low, 'dur2_med' = P_dur2_med, 'dur2_high' = P_dur2_high,
              'dur3_low' = P_dur3_low, 'dur3_med' = P_dur3_med, 'dur3_high' = P_dur3_high))
}

