library(tidyverse) 
library(DiscreteWeibull)
library(pracma)
tic()
clear()
cumulative_probability=0.995
FuncDistrSum <- function(n, q, beta){
  n_y <- 0
  max_cum_prob <- 0
  cont <- 0
  while(isTRUE(max_cum_prob <= cumulative_probability) ==TRUE){
    cat("searching n_y:", n_y)
    cont <- cont+1
    y <- seq(0,n_y,1) 
    d <- pmap_dbl(.l= list(x=y, q=list(q), beta=list(beta), 
                           zero = list(TRUE)), .f = ddweibull)
    df_probs <- data.frame(y,d)     
    # the grid considering n distributions
    tib1 <- expand.grid(as.data.frame(
      matrix(y, n_y+1, n))) %>%     
      # naming y1, y2...yn
      setNames(paste0("y",1:n)) %>%       
      # creating id for each row of the dataframe
      mutate(id = 1:n()) %>%      
      # realocate the column  "id"
      dplyr::relocate(any_of("id")) %>% rowwise() %>%       
      # creating the column yk = sum of ys 
      # through the values of the reference row "id"
      mutate(y_k = sum(c_across(-id)))     
    # the grid considering  "n" distributions
    tib2 <- expand.grid(as.data.frame
                        (matrix(y, n_y+1, n))) %>%       
      # naming prob1, prob2...probn
      setNames(paste0("prob",1:n)) %>% mutate(id = 1:n()) %>%       
      # creating  id for each row 
      dplyr::relocate(any_of("id")) %>%       
      # gathering y in groups
      gather(grupo, y, -id) %>% left_join(df_probs,by="y") %>%
      dplyr::select(-y) %>% spread(grupo, d) %>% rowwise() %>%
      mutate(prob_soma = prod(c_across(-id)))         
    # joing dataframes tib1 and tib2 by"id"    
    tib_final <- left_join(tib1, tib2, by = "id")%>%
      arrange(y_k) %>% group_by(y_k) %>%
      summarise(events = n(), prob_soma = sum(prob_soma)) %>%
      mutate(cum_prob_1 = cumsum(prob_soma)) %>% 
      mutate(cut_prob = cumsum(prob_soma)[cont]) %>% 
      mutate(n = n) %>% mutate(q = q) %>% 
      mutate(beta = beta) %>% 
      dplyr::relocate(any_of("beta")) %>%
      dplyr::relocate(any_of("q")) %>%
      dplyr::relocate(any_of("n"))      
    max_cum_prob <- max(tib_final$cut_prob)
    cat(", prob max:", max_cum_prob)
    cat(", y_k:", tib_final$y_k[tib_final$cum_prob_1 == 
                                  max_cum_prob], "\n")
    n_y <- n_y + 1    
  }  
  # showing the results obtained
  results <- tib_final %>% 
    # filter maximum value of interest
    filter(cum_prob_1<=cumulative_probability) %>% 
    filter(cum_prob_1 == max(cum_prob_1))  
  cat("results (1):prob max <= 0.995:",results$cum_prob_1)
  cat(", y_k:", results$y_k, "\n")
  cat("results (2):cutoff probability / 
  nearest value _ prob > 0.995:", 
      tib_final$cum_prob_1[tib_final$cum_prob_1 == 
                             max_cum_prob], ", sum of" ,n ,"distributions")
  cat(", y_k:", tib_final$y_k[tib_final$cum_prob_1 == 
                                max_cum_prob], "\n")  
  return(tib_final)
}
FuncDistrSum(n=3, q=0.5, beta=0.5) #Change according to the desired case
toc()
