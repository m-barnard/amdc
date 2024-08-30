trans.matrix <- function(chain, states)  {
  #########################################
  # Get adjacency matrix for a numeric sequence
  # chain: numeric vector for the sequences 
  # states: number of states in the sequences
  #########################################
  #from https://stats.stackexchange.com/questions/103914/markov-chain-state-transition-probability-in-r
  states <- states
  p <- matrix(nrow = states, ncol = states, 0) # initialise transition matrix
  for (t in 1:(length(chain) - 1)){ p[chain[t], chain[t + 1]] <- p[chain[t], chain[t + 1]] + 1}
  return (p)
}

gen_one_set_vecs <- function(vec){
  ###############################
  # Get set of adjacency matrices for a vector of character sequences
  # vec: vector of character sequences (e.g., c('AABBBBAAA', 'ABBBCCCA'))
  ###############################
  states <- sort(unique(unlist(lapply(vec, function(seq){unique(unlist(str_split(seq, '')))}))))
  vec_list <- lapply(vec, function(x){str_split(x, '')[[1]]})
  vec_list_num <- lapply(vec_list, function(x){for(j in 1:length(states)){x[x == states[j]] <- j}
    return(as.numeric(x))}) #have to turn into numerical for the trans.matrix() function to work
  trans_mats <- lapply(vec_list_num, trans.matrix, states = length(states))
  return(trans_mats)
}
get_correct_vecs <- function(vec, i, type){
  ###############################
  # Select portion of sequence to get adjacency matrices for
  # vec: vector of character sequences (e.g., c('AABBBBAAA', 'ABBBCCCA'))
  # i: indices of the sequence to select or remove (i = c(1) to use the entire sequence)
  # type: 'keep' if i indicates indices to keep or if you want to use the entire sequence
  #       'remove' if i indicates indicates to remove from the sequence
  ###############################
  if(length(i) == 1 & type == 'keep'){
    #this takes the full sequences
    vec_new <- vec
  } else if(length(i) == 2 & type == 'keep'){
    #this takes a portion of the sequences
    vec_new <- substr(vec,i[1], i[2])
  } else if(length(i) == 2 & type == 'remove'){
    #this removes a portion of the sequences, and kepps the rest
    vec1 <- substr(vec, 1, i[1])
    vec2 <- substr(vec, i[2], nchar(vec[1]))
    #want to paste the second vector into the first because that is actually where the sequences join
    vec_new <- paste0(vec2, vec1)
  }
  return(vec_new)
}

final_gen <- function(vec, i, type){
  ###############################
  # Get adjacency matrices for a vector of sequences and any desired subsets of the sequences
  # vec: vector of character sequences (e.g., c('AABBBBAAA', 'ABBBCCCA'))
  # i: list of indices of the sequence to select or remove (i = list(c(1)) to only use the entire sequence)
  # type: list of types, where type options are 'keep' if i indicates indices to keep or if you want to use the entire sequence
  #       'remove' if i indicates indicates to remove from the sequence
  # length of i and type inputs must be the same
  ###############################
  vec_list <- lapply(seq(1, length(i)), function(x){get_correct_vecs(vec, i[x][[1]], type[x])})
  mats <- lapply(vec_list, gen_one_set_vecs)
  vecs <- lapply(mats, function(x){sapply(x, as.vector)})
  names(mats) <- paste0(i, '_', type)
  names(vecs) <- paste0(i, '_', type)
  return(list('vecs'=vecs, 'mats'=mats))
}




