#' Calculates d score.
#' @param x A numeric vector or matrix.
#' @return If input is a numeric vector, a number is returned. For a matrix, a numeric vector is returned.
#' @export
calcDis <- function(x) {
  if (class(x)== "numeric") return(2*(atan(abs(sd(y)/mean(y))))/pi )
  return(apply(x, 1, function(y) (2*(atan(abs(sd(y)/mean(y))))/pi ) ))
}

#' Normalizes reactivity vector.
#' @description Normalizes raw reactivities using 2-8 \% method.
#' @param raw.estimates A vector of raw reactivities.
#' @return A vector of normalized reactivities.
#' @export
two.eight.normalize <- function(raw.estimates) {
  normalizer <- normalizer(raw.estimates)
  raw.estimates[which(raw.estimates != -999)] <- raw.estimates[which(raw.estimates != -999)] / normalizer
  return(raw.estimates)
}


#' Returns normalizer for reactivity vector.
#' @description Assesses normalization factor for raw reactivities using 2-8 \% method.
#' @param raw.estimates A vector of raw reactivities.
#' @return The normalization factor.
#' @export
normalizer <- function(raw.estimates) {
  raw.estimates[which(raw.estimates == -999)] <- NA
  sorted <- raw.estimates[order(raw.estimates)]
  if (any(is.na(sorted))) {
    normalize.range <- c(round((min(which(is.na(sorted)))-1) * .9), round((min(which(is.na(sorted)))-1) * .98))
  } else {
    normalize.range <- c(round(length(sorted) * .9), round(length(sorted)* .98))
  }
  normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  return(normalizer)
}



#' Identifies subgroupings of replicates for assessing within-group and between-group variation.
#' @param reps_A Number of replicates of group A.
#' @param reps_B Number of replicates of group B.
#' @param batches Logical suggesting if replicates of group A and B were performed in batches and are labelled accordingly. If TRUE, a heterogeneous/homogeneous subset may not have multiple samples from the same batch.
#' @param between_combs Dataframe with each column containing groupings of replicates of groups A and B, which will be used to assess between-group variation.
#' @param within_combs Data.frame with each column containing groupings of replicates of groups A or B, which will be used to assess within-group variation.
#' @return List of two dataframes, containing groupings for within-group and between-group variation.
#' @export
getCombs <- function(reps_A, reps_B, batches = F, between_combs= NULL, within_combs= NULL) {
  if ((is.null(within_combs) & !is.null(between_combs)) | (is.null(between_combs) & !is.null(within_combs))) stop("Homogeneous and heterogeneous sets should either both be null or both be specified with equal number of samples in each set.")
  if (!is.null(within_combs)) if (nrow(between_combs)!=nrow(within_combs)) stop("Heterogeneous and homogeneous sets should have equal number of samples.")
  
  set_membership = if (is.null(within_combs)) max(c(reps_A, reps_B)) else nrow(within_combs)
  while ( (set_membership/2 > min(c(reps_A, reps_B))) | (set_membership %% 2 != 0) ) {
    
    if (set_membership  == 3) break
    
    set_membership = set_membership - 1
    
  }
  
  all_combs = combn(c(paste0("A", 1:reps_A), paste0("B", 1:reps_B)),
                    set_membership)
  

  if (batches) {
    invalid_due_to_batch = apply(all_combs, 2, function(x) {
      curr_reps = mapply(function(y) y[2], strsplit(x, split= ""))
      return(length(unique(curr_reps)) < length(curr_reps))
    })

    all_combs = all_combs[, !invalid_due_to_batch] #Removes sets that have multiple samples from the same batch.
  }

  if (is.null(within_combs)) {
    if (min(c(reps_A, reps_B)) >= set_membership) {
      within_combs = cbind(combn(paste0("A", 1:reps_A), set_membership),
                           combn(paste0("B", 1:reps_B), set_membership))
    } else if (reps_A >= set_membership) {
      within_combs = cbind(combn(paste0("A", 1:reps_A), set_membership))
    } else if (reps_B >= set_membership) {
      within_combs = cbind(combn(paste0("B", 1:reps_B), set_membership))
    }
      
  }


  if (is.null(between_combs)) {
    between_combs_not_in = which(apply(all_combs, 2, paste0, collapse= "") %in% apply(within_combs, 2, paste0, collapse= ""))
    between_combs = data.frame(all_combs[, -between_combs_not_in])
    deg_of_heterogeneity = apply(between_combs, 2, function(x) {
      curr_reps = mapply(function(y) y[1], strsplit(x, split= ""))
      gA = length(which(curr_reps == "A"))
      gB = length(which(curr_reps == "B"))
      return(gA*gB/(gA + gB))
    })

    to_keep = which(deg_of_heterogeneity == max(deg_of_heterogeneity))
    between_combs = data.frame(between_combs[, to_keep])
  }

  return(list(between_combs= between_combs, within_combs= within_combs))
}


#' Performs guided discovery of differentially reactive regions.
#' @param reps_A Number of replicates of group A.
#' @param reps_B Number of replicates of group B.
#' @param rdf Dataframe of reactivities for each sample. Each column must be labelled as A1, A2, ..., B1, B2, ...
#' @param batches Logical suggesting if replicates of group A and B were performed in batches and are labelled accordingly. If TRUE, a heterogeneous/homogeneous subset may not have multiple samples from the same batch.
#' @param between_combs Dataframe with each column containing groupings of replicates of groups A and B, which will be used to assess between-group variation.
#' @param within_combs Data.frame with each column containing groupings of replicates of groups A or B, which will be used to assess within-group variation.
#' @param check_quality Logical, if TRUE, check regions for quality.
#' @param quality Worst allowed quality for a region to be tested.
#' @param evidence Minimum evidence of increase in variation from within-group comparisons to between-group comparisons for a region to be tested.
#' @return p-value for the tested region, estimated using one-sided Wilcoxon signed rank test.
#' @export
dStruct.guided <- function(rdf, reps_A, reps_B, batches = F,
                           within_combs = NULL, between_combs= NULL, check_quality = TRUE,
                           quality = "auto", evidence = 0) {

  if ((quality == "auto") & min(c(reps_A, reps_B)) != 1) quality = 0.5 else if (quality== "auto") quality = 0.2
  if (is.null(between_combs) | is.null(within_combs)) idcombs = getCombs(reps_A, reps_B, batches,
                                                                         between_combs, within_combs)
  if (is.null(between_combs)) between_combs = idcombs$between_combs
  if (is.null(within_combs)) within_combs = idcombs$within_combs

  d_within = dCombs(rdf, within_combs)
  d_between = dCombs(rdf, between_combs)

  if (mean(d_within, na.rm = T) > quality) return(NA)
  if (median(d_between - d_within, na.rm = T) < evidence) return(NA)

  result <- tryCatch({
    wilcox.test(d_within, d_between, alternative = "less", paired= T)$p.value
  }, error= function(e) {
    #Place holder for those transcripts that can't be tested due to error.
    result = NA
  })

  return(result)

}



#' Constructs potential differentially reactive regions.
#' @param d_within Nucleotide-wise d score for within-group variation.
#' @param d_spec Nucleotide-wise d score for between-group variation.
#' @param rdf Dataframe of reactivities for each sample.
#' @param min_length Minimum length of constructed regions.
#' @param check_signal_strength Logical, if TRUE, construction of regions must be based on nucleotides that have a minimum absolute value of reactivity.
#' @param check_nucs Logical, if TRUE, constructed regions must have a minimum number of nucleotides participating in Wilcoxon signed rank test.
#' @param check_quality Logical, if TRUE, check constructed regions for quality.
#' @param quality Worst allowed quality for a region to be tested.
#' @param evidence Minimum evidence of increase in variation from within-group comparisons to between-group comparisons for a region to be tested.
#' @param signal_strength Threshold for minimum signal strength.
#' @export
getRegions <- function(d_within, d_spec, rdf, min_length= 11,
                       check_signal_strength = T, check_nucs = T, check_quality = T,
                       quality = 0.5, evidence = 0, signal_strength = 0.1) {

  if (check_signal_strength) {
    insufficient_signal = apply(rdf, 1, function(x) all(abs(x) < signal_strength, na.rm=T))
    d_within[insufficient_signal] = NA
    d_spec[insufficient_signal] = NA
  }

  #Set min_length to next odd integer.
  min_length = if (min_length %% 2 == 0) min_length +1 else min_length

  to_test = c()

  #How many nucleotides with information must a constructed region have.
  min_nucs = min(min_length, 6)-1

  smooth_evidence = zoo::rollapply(d_spec - d_within, width = min_length,
                                   FUN= function(x) mean(x, na.rm=T))

  where_evidence = smooth_evidence > evidence
  smooth_evidence[where_evidence] = 1
  smooth_evidence[!where_evidence] = 0

  evidence_rle = rle(smooth_evidence)
  consider_regs = which((evidence_rle$values == 1) & (evidence_rle$lengths >= min_length))

  for (i in consider_regs) {
    check_n = check_nucs
    check_q = check_quality

    if (i == 1) {
      start = which(!is.na(d_within - d_spec))[1] + ((min_length-1)/2)
      end = evidence_rle$lengths[1]+((min_length-1)/2)
    } else {
      start = sum(evidence_rle$lengths[1:(i-1)]) + ((min_length-1)/2) + 1
      end = sum(evidence_rle$lengths[1:(i)]) + ((min_length-1)/2)
    }

    #Trim nucleotides from end with very low reactivities, only if regions detected are not too short.
    if (end - start + 1 >= 11) {
      trim_reac = apply(round(rdf[start:end, ], 2), 1,
                        function(x) all(abs(x) <= signal_strength, na.rm= T))
      if (all(trim_reac)) next

      start = start + which(!trim_reac)[1] - 1
      end = end - length(trim_reac) + tail(which(!trim_reac), 1)

      trim_reac = is.na(d_within[start:end]-d_spec[start:end])
      if (all(trim_reac)) next
      start = start + which(!trim_reac)[1] - 1
      end = end - length(trim_reac) + tail(which(!trim_reac), 1)
    }

    if (end - start + 1 < min_length) next


    if (check_n) {
      check_n = sum(!is.na(c(d_within[start:end] - d_spec[start:end]))) > min_nucs

      #check_q can be NA if for example, reactivity were all 0s in a region.
      if (is.na(check_n)) check_n = FALSE
    } else check_n = TRUE

    if (check_n) curr_test = start:end

    if (check_q) {
      check_q = mean(d_within[start:end], na.rm= T) <= quality

      #check_q can be NA if for example, reactivity were all 0s in a region.
      if (!is.na(check_q)) {
        #If the quality is not good for entire region, search for smaller regions that are good quality.
        if (!check_q) {
          wi_d = zoo::rollapply(d_within[start:end], width = min_length,
                                FUN= function(x) mean(x), align= "left")
          if (any(wi_d <= quality, na.rm= T))  {
            which_good = start + which(wi_d <= quality) -1
            curr_test = unique(unlist(purrr::map(which_good, function(x) x:(x+min_length-1))))
            check_q = TRUE
          }
        }
      } else check_q = FALSE
    } else check_q = TRUE

    if (check_n & check_q) to_test = c(to_test, start:end)

  }

  return(to_test)
}


#' Assesses within-group or between-group variation.
#' @param rdf Data.frame of reactivities for each sample.
#' @param combs Data.frame with each column containing groupings of samples.
#' @return Nucleotide-wise d score.
#' @export
dCombs <- function(rdf, combs) {
  d = matrix(, nrow(rdf), ncol(combs))
  for (i in 1:ncol(combs)) {
    curr_comb = as.character(combs[, i])
    curr_dat = rdf[, curr_comb]
    d[, i] = calcDis(curr_dat)
  }

  d = apply(d, 1, mean, na.rm=T)
  return(d)
}

#' Identifies contiguous regions from a list of nucleotide indices.
#' @param x A vector of integers.
#' @param gap Allowed gap to merge regions.
#' @return Dataframe storing start and stop sites of continguous regions.
#' @export
getContigRegions <- function(x, gap = 0) {
  if ((gap %% 1 != 0) | (gap < 0)) stop("parameter \'gap\' supplied to getContigRegions must be positive integer.")
  x = sort(x)
  return(data.frame(Start= c(x[1], x[which(diff(x) > 1 + gap) +1]),
                    Stop = c(x[which(diff(x) > 1 + gap)], tail(x, 1))))
}



#' Performs de novo discovery of differentially reactive regions.
#' @param reps_A Number of replicates of group A.
#' @param reps_B Number of replicates of group B.
#' @param rdf Dataframe of reactivities for each sample.
#' @param min_length Minimum length of constructed regions.
#' @param check_signal_strength Logical, if TRUE, construction of regions must be based on nucleotides that have a minimum absolute value of reactivity.
#' @param check_nucs Logical, if TRUE, constructed regions must have a minimum number of nucleotides participating in Wilcoxon signed rank test.
#' @param check_quality Logical, if TRUE, check constructed regions for quality.
#' @param batches Logical suggesting if replicates of group A and B were performed in batches and are labelled accordingly. If TRUE, a heterogeneous/homogeneous subset may not have multiple samples from the same batch.
#' @param between_combs Dataframe with each column containing groupings of replicates of groups A and B, which will be used to assess between-group variation.
#' @param within_combs Data.frame with each column containing groupings of replicates of groups A or B, which will be used to assess within-group variation.
#' @param quality Worst allowed quality for a region to be tested.
#' @param evidence Minimum evidence of increase in variation from within-group comparisons to between-group comparisons for a region to be tested.
#' @param signal_strength Threshold for minimum signal strength.
#' @param ind_regions Logical, if TRUE, test each region found in the transcript separately.
#' @param gap Integer. Join regions if they are separated by these many nucleotides.
#' @param get_FDR Logical, if FALSE, FDR is not reported.
#' @param proximity_assisted Logical, if TRUE, proximally located regions are tested together.
#' @param proximity Maximum distance between constructed regions for them to be considered proximal.
#' @param proximity_defined_length If performing a "proximity-assisted" test, minimum end-to-end length of a region to be tested.
#' @return Constructs regions, reports p-values and FDR for them.
#' @export
dStruct <- function(rdf, reps_A, reps_B, batches = F, min_length = 11,
                    check_signal_strength = T, check_nucs = T, check_quality = T,
                    quality = "auto", evidence = 0, signal_strength = 0.1,
                    within_combs = NULL, between_combs= NULL, ind_regions = T, gap = 1,
                    get_FDR = T, proximity_assisted = F, proximity = 10,
                    proximity_defined_length = 30) {

  if ((quality == "auto") & min(c(reps_A, reps_B)) != 1) quality = 0.5 else if (quality == "auto") quality = 0.2
  if (is.null(between_combs) | is.null(within_combs)) idcombs = getCombs(reps_A, reps_B, batches,
                                                                         between_combs, within_combs)
  if (is.null(between_combs)) between_combs = idcombs$between_combs
  if (is.null(within_combs)) within_combs = idcombs$within_combs

  d_within = dCombs(rdf, within_combs)
  d_between = dCombs(rdf, between_combs)

  to_test = getRegions(d_within, d_between, rdf, min_length,
                       check_signal_strength, check_nucs, check_quality,
                       quality, evidence, signal_strength)

  if (is.null(to_test) | !length(to_test)) return(NULL)

  contigs_test = getContigRegions(to_test, gap)

  if (!ind_regions & !proximity_assisted) {

    result <- tryCatch({
      wilcox.test(d_within[to_test], d_between[to_test],
                  alternative = "less", paired= T)$p.value
    }, error= function(e) {
      #Place holder for those transcripts that can't be tested due to insufficient data points.
      result = NA
    })

    result = list(regions= contigs_test, pval = result)

  } else if (!proximity_assisted) {

    pvals = c()
    for (i in 1:nrow(contigs_test)) {
      curr_res = tryCatch({
        wilcox.test(d_within[contigs_test$Start[i]:contigs_test$Stop[i]],
                    d_between[contigs_test$Start[i]:contigs_test$Stop[i]],
                    alternative = "less", paired= T)$p.value

      }, error= function(e) {
        #Place holder for those transcripts that can't be tested due to insufficient data points.
        result = NA
      })
      pvals = c(pvals, curr_res)


    }

    contigs_test = cbind(contigs_test, pval = pvals)
    if (get_FDR) contigs_test = cbind(contigs_test, FDR = p.adjust(pvals, "BH"))
    result = contigs_test
  } else {

    proximal_contigs = which(contigs_test$Start[-1]-
                               contigs_test$Stop[-nrow(contigs_test)] < proximity)+1

    if (length(proximal_contigs)) {
      proximally_tied_regs = getContigRegions(proximal_contigs)
      proximally_tied_regs$Start = proximally_tied_regs$Start-1
      proximal_contigs = unlist(apply(proximally_tied_regs, 1,
                                      function(x) x[1]:x[2]))
      non_proximal = setdiff(1:nrow(contigs_test),
                             proximal_contigs)
    } else {
      non_proximal = 1:nrow(contigs_test)
    }

    which_too_short = apply(contigs_test[non_proximal, ], 1,
                            function(x) x[2]- x[1] + 1 < proximity_defined_length)

    pvals = c()
    for (i in non_proximal) {

      curr_nucs = contigs_test$Start[i]:contigs_test$Stop[i]
      if (length(curr_nucs) < proximity_defined_length) {
        curr_res = NA
      } else {
        curr_res = tryCatch({
          wilcox.test(d_within[curr_nucs],
                      d_between[curr_nucs],
                      alternative = "less", paired= T)$p.value

        }, error= function(e) {
          #Place holder for those transcripts that can't be tested due to insufficient data points.
          result = NA
        })
      }

      pvals = c(pvals, curr_res)
    }

    if (length(proximal_contigs)) {
      proximity_defined_regs = list(Start = c(), Stop = c())
      for (i in 1:nrow(proximally_tied_regs)) {
        curr = contigs_test[proximally_tied_regs$Start[i]:proximally_tied_regs$Stop[i], ]
        curr_nucs = unlist(apply(curr, 1,
                                 function(x) x[1]:x[2]))
        curr_length = max(curr_nucs) - min(curr_nucs) + 1

        if (curr_length >= proximity_defined_length) {
          curr_res = tryCatch({
            wilcox.test(d_within[curr_nucs],
                        d_between[curr_nucs],
                        alternative = "less", paired= T)$p.value

          }, error= function(e) {
            #Place holder for those transcripts that can't be tested due to insufficient data points.
            result = NA
          })
        } else {
          curr_res = NA
        }

        # print(paste(curr$Start, collapse = ","))
        proximity_defined_regs$Start = c(proximity_defined_regs$Start,
                                         paste(curr$Start, collapse = ","))
        proximity_defined_regs$Stop = c(proximity_defined_regs$Stop,
                                        paste(curr$Stop, collapse = ","))
        pvals = c(pvals, curr_res)
      }

      contigs_test = rbind(contigs_test[non_proximal, ],
                           as.data.frame(proximity_defined_regs))
    }

    contigs_test = cbind(contigs_test, pval = pvals)
    if (get_FDR) contigs_test = cbind(contigs_test, FDR = p.adjust(pvals, "BH"))
    result = contigs_test
  }

  return(result)

}



#' Performs de novo discovery of differentially reactive regions for transcriptome-wide data.
#' @param reps_A Number of replicates of group A.
#' @param reps_B Number of replicates of group B.
#' @param rl List of dataframes of reactivities for each sample.
#' @param min_length Minimum length of constructed regions.
#' @param check_signal_strength Logical, if TRUE, construction of regions must be based on nucleotides that have a minimum absolute value of reactivity.
#' @param check_nucs Logical, if TRUE, constructed regions must have a minimum number of nucleotides participating in Wilcoxon signed rank test.
#' @param check_quality Logical, if TRUE, check constructed regions for quality.
#' @param batches Logical suggesting if replicates of group A and B were performed in batches and are labelled accordingly. If TRUE, a heterogeneous/homogeneous subset may not have multiple samples from the same batch.
#' @param between_combs Dataframe with each column containing groupings of replicates of groups A and B, which will be used to assess between-group variation.
#' @param within_combs Data.frame with each column containing groupings of replicates of groups A or B, which will be used to assess within-group variation.
#' @param quality Worst allowed quality for a region to be tested.
#' @param evidence Minimum evidence of increase in variation from within-group comparisons to between-group comparisons for a region to be tested.
#' @param signal_strength Threshold for minimum signal strength.
#' @param ind_regions Logical, if TRUE, test each region found in the transcript separately.
#' @param gap Integer. Join regions if they are separated by these many nucleotides.
#' @param processes Number of parallel processes to use.
#' @param method Character specifying either guided or de novo discovery approach.
#' @param proximity_assisted Logical, if TRUE, proximally located regions are tested together.
#' @param proximity Maximum distance between constructed regions for them to be considered proximal.
#' @param proximity_defined_length If performing a "proximity-assisted" test, minimum end-to-end length of a region to be tested.
#' @return Constructs regions, reports p-values and FDR for them.
#' @export
dStructome <- function(rl, reps_A, reps_B, batches= F, min_length = 11,
                       check_signal_strength = T, check_nucs = T, check_quality = T,
                       quality = "auto", evidence = 0, signal_strength = 0.1,
                       within_combs = NULL, between_combs= NULL, ind_regions = T, gap = 1,
                       processes = "auto", method = "denovo",
                       proximity_assisted = F, proximity = 10,
                       proximity_defined_length = 30) {

  if (is.null(names(rl)) | (length(names(rl)) != length(rl))) stop("List \'rl\' supplied to dStructome without transcript names.")

  if (processes == "auto") {
    processes = parallel::detectCores()-1
  } else {
    processes = as.numeric(processes)
    if ((processes %% 1 != 0) | (processes < 0)) stop("parameter \'processes\' supplied to dStructome must be positive integer.")
  }

  if (is.null(between_combs) | is.null(within_combs)) idcombs = getCombs(reps_A, reps_B, batches,
                                                                         between_combs, within_combs)
  if (is.null(between_combs)) between_combs = idcombs$between_combs
  if (is.null(within_combs)) within_combs = idcombs$within_combs

  if (method == "guided") {

    result <- parallel::mcmapply(function(x) {
      dStruct.guided(x, reps_A, reps_B, batches,
                     within_combs, between_combs, check_quality,
                     quality, evidence)
    }, rl, mc.cores=processes)

    res_df = data.frame(t = names(rl), pval = result)
    res_df = subset(res_df, !is.na(pval))
    res_df$FDR = p.adjust(res_df$pval, "BH")

  } else if (method == "denovo") {

    result <- parallel::mcmapply(function(x) {
      dStruct(x, reps_A, reps_B, batches, min_length,
              check_signal_strength, check_nucs, check_quality,
              quality, evidence, signal_strength,
              within_combs, between_combs, ind_regions, gap,
              get_FDR = F, proximity_assisted, proximity,
              proximity_defined_length)
    }, rl, mc.cores=processes, SIMPLIFY = F)

    if (ind_regions) {
      res_df = data.frame(t= NA, Start= NA, Stop = NA, pval= NA)
      for (i in 1:length(result)) {
        if (is.null(result[[i]])) next
        res_df = rbind(res_df, data.frame(t= names(result)[i], result[[i]]))
      }

      res_df = res_df[-1, ]
      res_df$FDR = p.adjust(res_df$pval, "BH")
    } else {
      res_df = data.frame(t= NA, Start= NA, Stop = NA)
      pvals = data.frame(t= NA, pval= NA)
      for (i in 1:length(result)) {
        if (is.null(result[[i]])) next
        pvals = rbind(pvals, data.frame(t =  names(result)[i], pval = result[[i]]$pval), stringsAsFactors= F)
        res_df = rbind(res_df, data.frame(t= names(result)[i], result[[i]]$regions))
      }

      pvals = pvals[-1, ]
      res_df = res_df[-1, ]

      pvals$FDR = p.adjust(pvals$pval, "BH")

      res_df$pval = apply(res_df, 1, function(x) subset(pvals, t == as.character(x[1]))$pval)
      res_df$FDR = apply(res_df, 1, function(x) subset(pvals, t == as.character(x[1]))$FDR)
    }
  }

  return(res_df)

}



#' Plots differentially reactive regions.
#' @param rl List of dataframes of reactivities for each sample.
#' @param diff_regions Dataframe of regions with significance of differentially reactivity.
#' @param outfile The name for pdf file which will be saved.
#' @param fdr FDR threshold for plotted regions.
#' @param ylim Y-axis limits for plots.
#' @return Saves a PDF for all differentially reactive regions. Returns NULL.
#' @export
plot_dStructurome <- function(rl, diff_regions, outfile, fdr = 0.05, ylim = c(-0.05, 3)) {
  diff_regions = subset(diff_regions, FDR < fdr)
  diff_t = unique(diff_regions$t)

  pdf(paste0(outfile, ".pdf"), width=7,height=6)
  for (i in 1:length(diff_t)) {
    curr_t = as.character(diff_t[i])
    curr_regs = subset(diff_regions, t == curr_t)
    curr_df = rl[[curr_t]]

    print(ggplot2::ggplot(data=data.frame(x=0,y=0), ggplot2::aes(x=x, y=y)) +
            ggplot2::annotate("text", x = 4, y = 25,
                     label = curr_t) +
            ggplot2::theme_classic()+
            ggplot2::theme(axis.line=ggplot2::element_blank(),
                  axis.text.x=ggplot2::element_blank(),
                  axis.text.y=ggplot2::element_blank(),
                  axis.ticks=ggplot2::element_blank(),
                  axis.title.x=ggplot2::element_blank(),
                  axis.title.y=ggplot2::element_blank(),
                  legend.position="none",
                  panel.background=ggplot2::element_blank(),
                  panel.border=ggplot2::element_blank(),
                  panel.grid.major=ggplot2::element_blank(),
                  panel.grid.minor=ggplot2::element_blank(),
                  plot.background=ggplot2::element_blank())
    )

    dat = data.frame(curr_df, n = 1:nrow(curr_df))
    dat = reshape2::melt(dat, id.vars = "n")

    print(ggplot2::ggplot(dat, ggplot2::aes(x= n, y= value)) + 
            ggplot2::geom_bar(stat="identity") +ggplot2::facet_grid(variable~.)+
            ggplot2::coord_cartesian(ylim=ylim) +
            ggplot2::geom_rect(ggplot2::aes(NULL, NULL, xmin=Start-0.5, xmax=Stop+0.5),
                      ymin= -Inf, ymax= Inf, data= curr_regs, fill= "red", 
                      color= NA, alpha= 0.3) +
            ggplot2::theme(strip.text = ggplot2::element_text(size = 6), 
                           legend.position = "none") +
            ggplot2::ggtitle(curr_t)+
            ggplot2::xlab("Nucleotide") + ggplot2::ylab("Reactivity"))


    for (r in 1:nrow(curr_regs)) {
      dat = data.frame(curr_df[curr_regs$Start[r]:curr_regs$Stop[r], ],
                       n = curr_regs$Start[r]:curr_regs$Stop[r])
      dat = reshape2::melt(dat, id.vars = "n")

      print(ggplot2::ggplot(dat, ggplot2::aes(x= n, y= value)) + 
              ggplot2::geom_bar(stat="identity") +ggplot2::facet_grid(variable~.)+
              ggplot2::coord_cartesian(ylim=ylim) +
              ggplot2::theme(strip.text = ggplot2::element_text(size = 6), 
                             legend.position = "none") +
              ggplot2::ggtitle(curr_t)+
              ggplot2::xlab("Nucleotide") + ggplot2::ylab("Reactivity"))
    }

  }
  dev.off()
  return()
}

