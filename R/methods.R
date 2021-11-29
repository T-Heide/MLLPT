#' MLE assignment of (low-depth) samples to a tree
#'
#' @param tree A tree object (of class phylo).
#' @param phydata Phylogenetic data used for the tree construction (of class phyDat).
#'
#' Note: The attribute attr(*, "id") containng mutation ids has to be added (e.g., c("chr1:5_A/T", ...)
#' @param sample_data A list of mutation data of samples to add the the tree.
#'
#' Each element must contain a data frame with the following columns:
#'  \itemize{
#'  \item{'id': }{The mutation id as set in the 'phydata' argument.}
#'  \item{'alt' (or 'alt_count'): }{The number of mutated reads.}
#'  \item{'depth' (or 'dp'): }{The coverage of the site.}
#'  \item{'cn_total' (or 'cn'): }{The copy-number of the site.}
#'  \item{'cn_mutated' (or 'mm'): }{The number of mutated alleles (this is assumed to be 1 if missing).}
#'  }
#' @param min_confidence (optional) minimum confidence in the the edge a sample is assigned to (default: 0).
#' @param vaf_bkgr (optional) expected background mutation rate of unmutated sites (default 0.01).
#' @param purity_estimates  (optional) purity of the samples.
#' Must be the same length as the sample_data list (default 1.0, pure samples).
#' @param min_edge_length  (optional) numeric value indicating the length of the edges
#' to the added tips relative to the tree height (default: 0.01).
#' @param return_details  (optional) logical flag indicating if detailed
#' information should be returned (default: false).
#' @param optimize_values  (optional) logical flag indicating if purity, background mutation rate,
#' (and loss fraction) should be optimized (default: true).
#' @param max_vaf_bkgr  (optional) maximum background mutation rate (default: 0.1).
#' @param max_loss_frac  (optional) maximum loss fraction (default: 0).
#' @param loss_frac_init  (optional) initial value of loss fraction (default: 0).
#' @param rescale_tree (optional) flag indicating if the tree should be rescaled (default: true).
#' @param control (optional) control passed to optimizer (default: NULL).
#' @param n_bootstraps (optional) Number of bootstraps to perform (default: 0).
#' @param n_cores (optional) Number of cores to use for bootstrapping (default: 0).
#' @param ...  unused arguments.
#' @return A tree object with the samples in 'sample_data' added to it.
#' @export
#' @import treeman
add_lowpass_sampled =
  function(
    tree,
    phydata,
    sample_data,
    min_confidence = 0,
    vaf_bkgr = 0.01,
    purity_estimates = rep(0.75, length(sample_data)),
    min_edge_length = 0,
    return_details = TRUE,
    optimize_values = TRUE,
    max_vaf_bkgr = 0.1,
    max_loss_frac = 0,
    loss_frac_init = 0,
    rescale_tree = TRUE,
    control = NULL,
    n_bootstraps = 0,
    n_cores = 1,
    ...
  ) {

  if (utils::packageVersion("treeman") < "1.1.5") {
    stop("Please install treeman > 1.1.5 using remotes::install_github(\"DomBennett/treeman\")\n")
  }

  initial_values = mget(ls())

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Helper function
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  # helper function
  get_data = function(data, ids) {

    data %>%
      dplyr::filter(id %in% ids) %>%
      dplyr::filter(!is.na(cn_total)) %>%
      dplyr::filter(cn_total <= 4) %>%
      dplyr::filter(cn_total > 0) %>%
      dplyr::filter(depth > 0) %>%
      dplyr::select(alt_count, depth, cn_mutated, cn_total) %>%
      dplyr::count(alt_count, depth, cn_mutated, cn_total) %>%
      magrittr::set_colnames(c("n", "N", "mm", "cn", "weight"))

  }

  correct_mdata_structure = function(d) {

    # check columns
    if (!"alt_count" %in% colnames(d)) { # missing alt_count column
      if ("alt" %in% colnames(d)) { # use alt if available else stop
        d$alt_count = d$alt
      } else {
        stop("missing 'alt_count' column")
      }
    }

    if (!"depth" %in% colnames(d)) { # missing depth column
      if ("dp" %in% colnames(d)) { # use dp if available else stop
        d$depth = d$dp
      } else {
        stop("missing 'depth' column")
      }
    }

    if (!"cn_total" %in% colnames(d)) { # missing cn_total column
      if ("cn" %in% colnames(d)) { # use cn if available else stop
        d$cn_total = d$cn
      } else if ("copy_number" %in% colnames(d)) {
        d$cn_total = d$copy_number
      } else {
        stop("missing 'cn_total' column")
      }
    }

    if (!"cn_mutated" %in% colnames(d)) { # missing cn_mutated column
      if ("mm" %in% colnames(d)) { # use mm if available else use 1
        d$cn_mutated = d$mm
      } else {
        d$cn_mutated = rep(1, NROW(d)) # defaults to 1
      }
    }

    return(d)
  }

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  al = checkmate::makeAssertCollection()
  checkmate::assertClass(tree, "phylo", add = al)
  checkmate::assertClass(phydata, "phyDat", add = al)
  checkmate::assertSetEqual(names(phydata), tree$tip.label, add = al)
  checkmate::assertList(sample_data, names = "named", any.missing = FALSE, types = "data.frame", min.len = 1, add = al)
  checkmate::assertList(lapply(sample_data, correct_mdata_structure), add = al)
  checkmate::assertNumeric(min_confidence, lower = 0, upper = 1, len = 1, any.missing = FALSE, add = al)
  checkmate::assertNumeric(vaf_bkgr, lower = 0, upper = 1, len = 1, any.missing = FALSE, add = al)
  checkmate::assertNumeric(purity_estimates, lower = 0, upper = 1, len = length(sample_data), any.missing = FALSE, add = al)
  checkmate::assertNumeric(min_edge_length, lower = 0, len = 1, finite = TRUE, any.missing = FALSE, add = al)
  checkmate::assertFlag(return_details, add = al)
  checkmate::assertFlag(optimize_values, add = al)
  checkmate::assertNumeric(max_vaf_bkgr, lower = 0, upper = 0.5, len = 1, any.missing = FALSE, add = al)
  checkmate::assertNumeric(max_loss_frac, lower = 0, upper = 0.5, len = 1, any.missing = FALSE, add = al)
  checkmate::assertTRUE(max_loss_frac >= loss_frac_init, add = al)
  checkmate::assertTRUE(max_vaf_bkgr >= vaf_bkgr, add = al)
  checkmate::assertFlag(rescale_tree, add = al)
  checkmate::assertIntegerish(n_bootstraps, add = al)
  checkmate::assertIntegerish(n_cores, add = al)
  checkmate::reportAssertions(al)

  if (max_vaf_bkgr > 0.1)
    warning(
      "Do you really want to set a background rate of ",
      max_vaf_bkgr,
      "?\n",
      sep = ""
    )

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # General data for analysis
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if (!is.null(names(purity_estimates)))
    purity_estimates = as.numeric(purity_estimates[names(sample_data)])

  # get mutation ids for each node of the tree
  sample_data = lapply(sample_data, correct_mdata_structure)
  initial_values$sample_data = sample_data # more compact version of data

  root_node = phangorn::getRoot(tree)
  tip_rooted_on = min(tree$edge[tree$edge[, 1] == root_node, 2])
  root_to_tip_edge = which(tree$edge[, 1] == root_node & tree$edge[, 2] == tip_rooted_on)

  mids_per_edge = get_mutation_ids_per_edge(tree, phydata)
  stopifnot(is.null(unlist(mids_per_edge[[as.character(root_to_tip_edge)]])))

  lost_ids = unlist(lapply(mids_per_edge, "[[", "lost"))
  mids_per_edge_gain = lapply(lapply(mids_per_edge, "[[", "added"), function(x) x[!x %in% lost_ids])

  if (rescale_tree) {
    tree$edge.length =
      vapply(
        mids_per_edge_gain,
        length,
        integer(1),
        USE.NAMES = FALSE
      ) # rescale the tree?
  }

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Adding of samples to the tree
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  added_tips = logical()
  per_edge_max_ll_data = list()
  per_edge_ll_data = list()
  per_edge_mdata = list()
  tree$node.label = paste0("N", seq_len(tree$Nnode))
  tree_tm = methods::as(tree, "TreeMan") #
  purity_optim = numeric()
  vaf_optim = numeric()
  loss_frac_optim = numeric()
  mll = numeric()
  s_idx = numeric()
  edges_samples = integer()
  pi_samples = numeric()
  bootstrap_results = list()

  # get tree and species ages:
  tree_age_o =  treeman::getAge(tree_tm)
  spns = treeman::getSpnsAge(tree_tm, tree_tm@all, tree_age_o)
  rownames(spns) = spns$spn

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Lookup tables
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  edge_to_node = tree$edge[, 2]
  names(edge_to_node) = seq_along(edge_to_node)
  edge_to_node = edge_to_node[-root_to_tip_edge]

  node_to_edge = as.numeric(names(edge_to_node))
  names(node_to_edge) = edge_to_node
  node_to_edge = node_to_edge[order(as.numeric(names(node_to_edge)))]

  anc_edge = lapply(names(edge_to_node), function(x) {
    # map edge to all edges before it
    next_node = edge_to_node[as.character(x)]
    anc_nodes = phangorn::Ancestors(tree, next_node)
    anc_edges = as.character(node_to_edge[as.character(anc_nodes)])
    stats::na.omit(anc_edges)
  })

  not_anc_edge = lapply(names(edge_to_node), function(x) {
    names(edge_to_node)[!names(edge_to_node) %in% c(anc_edge[[x]], x)]
  })


  # add all samples to the tree:
  for (j in seq_along(sample_data)) {

    cat("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=\n\n")
    cat(
      "Processing sample: ",
      names(sample_data)[j],
      " (",
      j,
      "/",
      length(sample_data),
      ")\n\n",
      sep = ""
    )

    # current data
    sample = names(sample_data)[j]
    sample_purity = purity_estimates[j]
    data = sample_data[[j]]
    vaf_bkgr_sample = vaf_bkgr
    if (max_loss_frac) loss_frac = loss_frac_init else loss_frac = 0

    # skip in some cases
    if (is.na(sample_purity)) next ()
    if (NROW(data) == 0) next ()
    if (sample %in% tree$tip.label) {
      n_dup = sum(c(tree$tip.label, added_tips) == sample)
      sample = paste0(sample, paste0(rep("_", n_dup), collapse = ""))
    }

    # mutation data of each edge
    per_edge_data = lapply(mids_per_edge_gain, get_data, data = data)
    per_edge_mdata[[sample]] = per_edge_data

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Internal functions for likelihoods ####
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    .ll_edge =
      function(d, edge, pi_m, purity, bkg_rate) {
        v = purity
        d = d[[edge]]

        if (pi_m > 0) {
          vaf_m = (d$mm * v) / (v * d$cn + (1 - v) * 2) # expected vaf for each mutation
          vaf_m = vaf_m * (1 - bkg_rate) + bkg_rate
          d1 = stats::dbinom(d$n, d$N, vaf_m)
        } else {
          d1 = 0
        }

        d0 = stats::dbinom(d$n, d$N, bkg_rate)

        val = sum(log(d1 * pi_m + d0 * (1 - pi_m)) * d$weight)
        return(val)
      }

    .ll_path = function(d, edge, pi_m, purity, bkg_rate, loss_frac) {

      edge = as.character(edge)

      #
      ll_mut =
        vapply(
          anc_edge[[edge]],
          .ll_edge,
          numeric(1),
          d = d,
          pi_m = 1 - loss_frac,
          purity = purity,
          bkg_rate = bkg_rate
        )

      ll_bkg =
        vapply(
          not_anc_edge[[edge]],
          .ll_edge,
          numeric(1),
          d = d,
          pi_m = 0,
          purity = purity,
          bkg_rate = bkg_rate
        )

      ll_self =
        .ll_edge(
          edge,
          d = d,
          pi_m = pi_m,
          purity = purity,
          bkg_rate = bkg_rate
        )

      # sum log-lik
      sum(c(unlist(ll_mut), unlist(ll_bkg), unlist(ll_self)))
    }

    .optimize_params = function(d, init, lower, upper, optim) {

      edges = names(edge_to_node)

      results = list(
        per_edge_max_ll = magrittr::set_names(rep(NA, length(edges)), edges),
        per_edge_opt_params = list(),
        params = magrittr::set_names(init, c("pos", "purity", "bkg", "loss")),
        max_ll = -Inf,
        edge_opt = NA
      )

      # initial values
      .f_optim = function(e, p) {
        args = init; args[optim] = p
        do.call(.ll_path, c(list(d, e), as.list(args)))
      }

      for (e in edges) {

        optim[1] = sum(d[[e]]$n) > 0

        opt_new = optim(
          par = init[optim],
          function(p) -1 * .f_optim(e, p),
          lower = lower[optim],
          upper = upper[optim],
          method = "L-BFGS-B",
          control = c(list(ndeps = rep(1e-6, sum(optim))), control)
        )

        if (opt_new$convergence != 0) {
          stop(opt_new$message)
        }

        # store results of this edge
        results$per_edge_max_ll[e] = -opt_new[["value"]]
        results$per_edge_opt_params[[e]] = init
        results$per_edge_opt_params[[e]][optim] = opt_new[["par"]]
        names(results$per_edge_opt_params[[e]]) = c("pos", "purity", "bkg", "loss")

        if (-opt_new$value > results$max_ll) {
          # if generally better then other update results
          results$edge_opt = e
          results$max_ll = -opt_new[["value"]]
          results$params[optim] = opt_new[["par"]]
        }
      }

      return(results)
    }

    .bootstrap = function(d, ...) {

      # permutate data
      i = as.numeric(unlist(lapply(d, function(x) seq_len(NROW(x)))))
      e = rep(names(d), vapply(d, NROW, integer(1)))
      n = as.numeric(unlist(lapply(d, "[", "weight")))

      idx = sample(seq_along(i), sum(n), replace = TRUE, prob = n)
      n_dash = table(idx)[as.character(seq_along(i))]; n_dash[is.na(n_dash)] = 0
      #n_dash = MCMCpack::rdirichlet(1, n) * sum(n)
      for (j in seq_along(i)) d[[e[j]]]$weight[i[j]] = n_dash[j]

      res = tryCatch({
        .optimize_params(d, ...)},
        error = function(e) return(NULL)
      )

      if (is.null(res)) return(NULL)

      cbind(
        data.frame(edge = res$edge_opt),
        t(data.frame(res$params))
      )
    }

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Optimize the samples parameters? ####
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    optim = c(TRUE, optimize_values, optimize_values, max_loss_frac > 0 & optimize_values)
    lower = c(0, 0, .Machine$double.eps, 0)
    upper = c(1, 1 - 1e-3, max_vaf_bkgr, max_loss_frac)
    init = c(0, scales::squish(sample_purity, 0, 1 - 1e-3), vaf_bkgr_sample, loss_frac)
    optimiser_results = .optimize_params(per_edge_data, init, lower, upper, optim)

    if (n_bootstraps) {
      cat("Bootstrapping data", n_bootstraps, "times...\n")
      withr::with_options(list(width = 60), {
        bootstrap_results[[sample]] =
          pbmcapply::pbmclapply(seq_len(n_bootstraps), function(i) {
            .bootstrap(per_edge_data, init, lower, upper, optim)
          }, mc.cores = n_cores) %>% do.call(what = rbind)
      })
    } else {
      bootstrap_results[sample] = list(NULL)
    }

    if (return_details) {
      per_edge_max_ll_data[[sample]] = optimiser_results$per_edge_max_ll
      mll[sample] = optimiser_results$max_ll
      edges_samples[sample] = optimiser_results$edge_opt
      pi_samples[sample] = optimiser_results$params["pos"]
      s_idx[sample] = j
    }

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if (optimize_values) {

      orig_mll = .ll_path(
        per_edge_data,
        edge = optimiser_results$edge_opt,
        pi_m = optimiser_results$params["pos"],
        purity = sample_purity,
        bkg_rate = vaf_bkgr_sample,
        loss_frac = loss_frac
      )

      with(optimiser_results, {
        cat("\n")
        cat("New values:\n")
        cat("\n")
        cat("- Background rate:", vaf_bkgr_sample, "->", params["bkg"], "\n")
        cat("- Purity:", sample_purity, "->", params["purity"], "\n")
        if (max_loss_frac)
          cat("- Loss fraction:", loss_frac, "->", params["loss"], "\n")
        cat("- MLL:", orig_mll, "->", max_ll, "\n")
        cat("\n")
      })


      vaf_optim[sample] = optimiser_results$params["bkg"]
      purity_optim[sample] = optimiser_results$params["purity"]
      loss_frac_optim[sample] = optimiser_results$params["loss"]

    } else {
      vaf_optim[sample] = NA
      purity_optim[sample] = NA
      loss_frac_optim[sample] = NA
    }

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Grid data for summary plots ####
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if (return_details) {

      # sum of log-likelihoods along the tree
      vals_pi_per_edge = seq(0, 1, length = 100)

      per_edge_ll =
        lapply(names(edge_to_node), function(edge) {
          vapply(vals_pi_per_edge, function(pi) {

            params = optimiser_results$per_edge_opt_params[[edge]]

            .ll_path(
                d = per_edge_data,
                edge =  edge,
                pi_m = pi,
                purity = params["purity"],
                bkg_rate = params["bkg"],
                loss_frac = params["loss"]
              )
          }, numeric(1))
        })

      names(per_edge_ll) = names(edge_to_node)
      per_edge_ll_data[[sample]] = per_edge_ll

    }

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    # Add sample to tree
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    pi_opt = optimiser_results$params["pos"]
    edge_opt = optimiser_results$edge_opt

    if (n_bootstraps) {
      max_lik_edge = mean(bootstrap_results[[sample]]$edge == edge_opt)
    } else {
      # normalize edge data and use as confidence value
      log_sum_exp = with(optimiser_results, exp(per_edge_max_ll - max_ll))
      lik_edge = log_sum_exp / sum(log_sum_exp, na.rm = TRUE)
      max_lik_edge = max(lik_edge, na.rm = TRUE)
    }

    if (max_lik_edge > min_confidence) {  # assign sample to tree if high confidence edge found.

      tryCatch({

        # find position to add the sample
        wh_node_idx = edge_to_node[edge_opt]

        sid = ifelse(length(tree$tip.label) >= wh_node_idx,
                     tree$tip.label[wh_node_idx],
                     tree$node.label[wh_node_idx - length(tree$tip.label)])


        # relative position to add to
        range_sid = spns[sid, ]
        rel_pos = 1 - pi_opt
        if (rel_pos <= 0) rel_pos = min_edge_length
        if (rel_pos >= 1) rel_pos = 1 - min_edge_length


        # location in tree
        pos_add = range_sid[["end"]] + (range_sid[["start"]] - range_sid[["end"]]) * rel_pos
        pos_end = pos_add - tree_age_o * min_edge_length

        if (pos_end < 0) pos_end = 0 # don't increase tree age!

        # find correct_modified sid:
        correct_nsid = FALSE
        nsid = sid
        while (!correct_nsid) {

          tree_age_cur = treeman::getAge(tree_tm)
          stopifnot(all.equal(tree_age_cur, tree_age_o))
          age_nsid = treeman::getSpnsAge(tree_tm, nsid, tree_age_cur)

          correct_nsid = (pos_add >= age_nsid[["end"]] |
                            isTRUE(all.equal(age_nsid[["end"]], pos_add))) &
            (age_nsid[["start"]] >= pos_add |
               isTRUE(all.equal(age_nsid[["start"]], pos_add)))

          if (!correct_nsid) {
            nnsid = treeman::getNdPrids(tree_tm, nsid)[1]
            if (nnsid == nsid)
              stop("Stuck at root ...")
            nsid = nnsid
          }
        }


        # add the tip to the tree
        tree_tm =
          treeman::addTip(
            tree_tm,
            sample,
            nsid,
            strt_age = pos_add,
            end_age = pos_end,
            tree_age = tree_age_cur
          )


        # store labels for later change
        cat("=> Added sample (confidence: ", max_lik_edge, ")\n\n", sep = "")
        new_tip = paste0(sample, " (Added confidence: ", signif(max_lik_edge, 6), ")")
        added_tips[[sample]] = new_tip

      })#, error=function(e) {
        #cat("=> Error: Failed to add sample!\n")
        #print(e)
      #})
    }
  }
  cat("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=\n\n")


  # convert the tree back to phylo class
  tree_mod = methods::as(tree_tm, "phylo")
  wh_relabel = tree_mod$tip.label %in% names(added_tips)
  tree_mod$tip.label[wh_relabel] = added_tips[tree_mod$tip.label[wh_relabel]]
  tree_mod$node.label = NULL


  # return extensive result set?
  if (return_details) {

    per_sample_results =
      data.frame(
        idx = s_idx,
        sample = names(edges_samples),
        edge = edges_samples,
        pi = pi_samples[names(edges_samples)],
        background_vaf = vaf_optim[names(edges_samples)],
        purity = purity_optim[names(edges_samples)],
        loss_frac = loss_frac_optim[names(edges_samples)],
        mll = mll[names(edges_samples)]
      )

    results =
      list(
        tree = tree_mod,
        inital_values = initial_values,
        per_sample_results = per_sample_results,
        max_ll_per_edge = per_edge_max_ll_data,
        ll_per_edge = per_edge_ll_data,
        mutation_data = per_edge_mdata,
        mutations_per_edge = mids_per_edge_gain,
        bootstrap_results = bootstrap_results
      )

    return(results)

  }

  return(tree_mod)
}


#' Internal function used to get mutation ids for each edge of the tree.
#'
#' @param tree A tree object (of class phylo).
#' @param phydata Phylogenetic data used for the tree construction (of class phyDat).
#'
#' @return A list of mutations gained and loosed on each edge.
get_mutation_ids_per_edge = function(tree, phydata) {

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Check arguments
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  al = checkmate::makeAssertCollection()
  checkmate::assertClass(tree, "phylo", add = al)
  checkmate::assertClass(phydata, "phyDat", add = al)
  checkmate::reportAssertions(al)

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  # determine ancestral states:
  anc_state = phangorn::ancestral.pars(tree, phydata, "MPR")
  mids_per_site_pattern = split(attr(phydata, "id"), attr(phydata, "index"))

  root_node = phangorn::getRoot(tree)
  wh_rooted_on = min(tree$edge[tree$edge[, 1] == root_node, 2])
  anc_state[[root_node]] = anc_state[[wh_rooted_on]]

  # assign mutations to each
  muts_per_edge =
    lapply(seq_len(NROW(tree$edge)), function(idx_edge) {

      node_from = tree$edge[idx_edge, 1]
      node_to = tree$edge[idx_edge, 2]

      state_mut_from = anc_state[[node_from]][, 2]
      state_mut_to = anc_state[[node_to]][, 2]

      muts_added = state_mut_to > state_mut_from
      muts_lossed = state_mut_to < state_mut_from

      mids_added = unlist(mids_per_site_pattern[which(muts_added)])
      mids_lossed = unlist(mids_per_site_pattern[which(muts_lossed)])

      return(list(added = mids_added, lost = mids_lossed))
    })

  names(muts_per_edge) = seq_len(NROW(muts_per_edge))

  return(muts_per_edge)
}


#' Function to load genotyping data from a file
#'
#' @param file Input file
#' @param use Optional vector of mutation ids to use.
#' @param cn_data Either a Genomic Ranges object containing cn and mm metadata columns or
#' a data.frame with a id, cn and mm column.
#'
#' @return A tibble containing the genotyping data.
#' @export
load_genotyping_file = function(file, use=NULL, cn_data=NULL) {

  stopifnot(file.exists(as.character(file)))
  stopifnot(is.null(use) | is.character(use))

  if (!is.null(cn_data)) {
    if (inherits(cn_data, "GenomicRanges")) {
      stopifnot("cn" %in% names(S4Vectors::elementMetadata(cn_data)))
    } else if (is.data.frame(cn_data)) {
      stopifnot(c("id", "cn") %in% colnames(cn_data))
    } else {
      stop("'cn_data' has to be a data.frame or GenomicRanges object.\n")
    }
  }

  #---------------------------------------------

  .parse_data_chunks =
    function(x, pos) {

      chunk_data =
        x %>%
        dplyr::filter(gsub("chr", "", chr) %in% 1:22) %>%
        dplyr::filter(source == "S") %>%
        dplyr::mutate(id = sprintf("%s:%d_%s/%s", chr, pos, ref, alt)) %>%
        dplyr::mutate(vaf = alt_count / depth) %>%
        dplyr::filter(id %in% use | is.null(use)) %>%
        dplyr::mutate(copy_number = NA, mm = NA)


      if (!is.null(cn_data)) {
        if (inherits(cn_data, "GenomicRanges")) {
          # 1) GRanges

          mdata_gr = with(chunk_data, GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos)))
          ol = GenomicRanges::findOverlaps(mdata_gr, cn_data, type = "any", select = "first")
          chunk_data$copy_number = cn_data$cn[ol]

          # mm col
          if ("mm" %in% names(S4Vectors::elementMetadata(cn_data))) {
            chunk_data$mm = cn_data$mm[ol]
          } else {
            warning(
              "For file ",
              basename(file),
              ": mutation multiplicity data missing! Assuming mm = 1.\n",
              sep = ""
            )
            chunk_data$mm[!is.na(chunk_data$copy_number)] = 1
          }

        } else if (is.data.frame(cn_data)) {
          # 2) data.frame

          # cn col
          mt = match(chunk_data$id, as.character(cn_data$id))
          chunk_data$copy_number = cn_data$cn[mt]

          # mm col
          if ("mm" %in% colnames(cn_data)) {
            chunk_data$mm = cn_data$mm[mt]
          } else {
            warning(
              "For file ",
              basename(file),
              ": mutation multiplicity data missing! Assuming mm = 1.\n",
              sep = ""
            )
            chunk_data$mm[!is.na(chunk_data$copy_number)] = 1
          }
        }
      }  else {
        # 3) No CN data
        warning(
          "For file ",
          basename(file),
          ": CN data missing, assuming CN 1+1.\n",
          sep = ""
        )
        chunk_data$copy_number = 2
        chunk_data$mm = 1
      }

      parsed_data <<- rbind(parsed_data, chunk_data)

      return(!all(x[, "source"] == "GL"))
    }

  #---------------------------------------------

  # parse data in chunks, stop once reading germline mutations
  parsed_data = NULL

  tryCatch({
    readr::read_tsv_chunked(
      as.character(file),
      .parse_data_chunks,
      chunk_size = 10000,
      col_types = "-cicc--iic",
     progress = FALSE
    )
  }, error = function(e) {
    stop(
      "Failed to load file '",
      basename(file),
      "'. Wrong file format?\n",
      sep = ""
    )
  })

  if (!all(is.na(parsed_data$copy_number))) {
    parsed_data =
      parsed_data %>%
      dplyr::filter(!is.na(copy_number))
  }

  return(parsed_data)
}
