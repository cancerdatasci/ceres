#' Calculate area under the curve for ROC and PR curves.
#' 
#' @param predicted_values Vector of continuous-valued prediction scores. 
#' Assumes that positive prediction scores correspond to 1 outcome in the true values.
#' @param true_values Vector of logical or binary values of ground truth
#' @return A data.frame of containing the AUROC, AUPRC, and F1 measure
#' @examples
#' true_values <- sample(0:1, 1000, replace=T)
#' predicted_values <- rnorm(1000)
#' evaluate_prediction_auc(predicted_values, true_values)
#' @export evaluate_prediction_auc
evaluate_prediction_auc <- function(predicted_values, true_values,
                                    min_resolution=1000, plot=T, fdr_threshold=0.05){
  
  
  non_na_indices <- !is.na(predicted_values) & !is.na(true_values)
  predicted_values <- predicted_values[non_na_indices]
  true_values <- true_values[non_na_indices]
  
  thresholds <- seq(min(predicted_values), 
                    max(predicted_values), 
                    length.out = min(sum(!is.na(predicted_values)), min_resolution))
  
  predicted_values <- predicted_values[!is.na(predicted_values) & !is.na(true_values)]
  true_values <- true_values[!is.na(predicted_values) & !is.na(true_values)]
  
  auc_mat <- sapply(thresholds, function(x){
    
    tpr <- sum(predicted_values >= x & true_values) / sum(true_values)
    fpr <- sum(predicted_values >= x & !true_values) / sum(!true_values)
    ppv <- sum(predicted_values >= x & true_values) / sum(predicted_values >= x)
    
    return(c(tpr, fpr, ppv) %>% set_names(c("tpr", "fpr", "ppv")))
    
  }) %>%
    t()
  
  
  id_auroc <- order(1-auc_mat[,"fpr"])
  avg_diff_auroc <- (auc_mat[id_auroc, "tpr"][-length(auc_mat[, "tpr"])] + auc_mat[id_auroc, "tpr"][-1])/2
  auroc <- sum(diff(1-auc_mat[id_auroc,"fpr"])*avg_diff_auroc)
  
  id_auprc <- order(auc_mat[,"tpr"])
  avg_diff_auprc <- (auc_mat[id_auprc, "ppv"][-length(auc_mat[,"ppv"])] + auc_mat[id_auprc, "ppv"][-1])/2
  auprc <- sum(diff(auc_mat[id_auprc, "tpr"])*avg_diff_auprc)
  
  eval_point <- data.frame(pred=predicted_values, tar=true_values) %>%
    dplyr::arrange(pred) %>%
    dplyr::mutate(prob_yes=cumsum(tar) / sum(tar),
                  prob_no=1 - cumsum(abs(1-tar)) / sum(abs(1-tar))) %>%
    dplyr::mutate(prob_diff=(prob_yes - prob_no)^2) %>%
    dplyr::filter(prob_diff == min(prob_diff)) %>%
    magrittr::extract2("pred") %>%
    unique()
  
  precis <- sum(predicted_values >= eval_point & true_values) / sum(predicted_values >= eval_point)
  recall <- sum(predicted_values >= eval_point & true_values) / sum(true_values)
  f1 <- 2*(precis*recall) / (precis + recall)
  recall_eval_point <- auc_mat[id_auroc, "ppv"] %>%
    magrittr::subtract(1-fdr_threshold) %>%
    magrittr::is_less_than(0) %>%
    ifelse(-1, 1) %>%
    diff() %>%
    abs() %>%
    magrittr::is_greater_than(0) %>%
    which() %>%
    max()
  recall_at_fdr <- auc_mat[id_auroc, "tpr"][recall_eval_point]
  
  if(plot){
    auroc_plot <- data.frame(tpr=auc_mat[id_auroc, "tpr"],
                             fpr=auc_mat[id_auroc,"fpr"]) %>%
      ggplot2::ggplot(aes(x=1-fpr, y=tpr)) +
      geom_line() +
      geom_ribbon(aes(ymin=0, ymax=tpr), alpha=0.25, fill="#377eb8") +
      geom_ribbon(aes(ymax=1, ymin=tpr), alpha=0.25, fill="grey50") +
      labs(title="ROC curve")
    
    auprc_plot <- data.frame(tpr=auc_mat[id_auroc, "tpr"],
                             ppv=auc_mat[id_auroc,"ppv"]) %>%
      ggplot2::ggplot(aes(x=tpr, y=ppv)) +
      geom_line() +
      geom_ribbon(aes(ymin=0, ymax=ppv), alpha=0.25, fill="#377eb8") +
      geom_ribbon(aes(ymax=1, ymin=ppv), alpha=0.25, fill="grey50") +
      labs(title="Precision-recall curve", x="Recall", y="Precision") +
      geom_segment(aes(y=1-fdr_threshold, yend=1-fdr_threshold,
                       x=0, xend=recall_at_fdr), linetype=3) +
      geom_segment(aes(x=recall_at_fdr, xend=recall_at_fdr,
                       y=1-fdr_threshold, yend=0), linetype=3)
  }
  
  res <- list(performance_df=data.frame(auc_roc=auroc,
                                        auc_pr=auprc,
                                        f1=f1,
                                        recall_at_fdr=recall_at_fdr),
              performance_plots=list(auroc=auroc_plot,
                                     auprc=auprc_plot))
  
  return(res)
  
}