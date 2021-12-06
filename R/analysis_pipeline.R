phydat_to_tree = function(phydat, all_trees = FALSE, ...) {
  phydat %>%
    phangorn::pratchet(trace = FALSE, all = all_trees, ...) %>%
    ape::unroot() %>%
    ape::root(out = "GL", resolve.root = TRUE) %>%
    phangorn::acctran(phydat)
}

construct_deep_tree = function(d, vaf_threshold = 0.05, minit = 100) {
  phydat = get_phydat_from_genotyping(d, vaf_threshold = vaf_threshold)
  tree = phydat_to_tree(phydat, all_trees = TRUE, minit = minit, k = minit)
  return(list(tree = tree, phydat = phydat))
}

get_phydat_from_genotyping = function(d, vaf_threshold=0.05) {

  # check input
  checkmate::assertNumber(vaf_threshold, lower = 0, upper = 1)
  checkmate::assertList(d, names = "named")
  for (i in seq_along(d)) {
    checkmate::assertDataFrame(d[[i]])
    checkmate::assertSubset(c("id", "vaf"), colnames(d[[i]]))
    checkmate::assertTRUE(all(d[[1]]$id == d[[i]]$id))
  }


  # convert to phyDat object:
  py_data =
    lapply(d, function(x) x$vaf >= vaf_threshold) %>%
    do.call(what = cbind) %>%
    data.frame() %>%
    dplyr::mutate(GL = 0) %>%
    t() %>%
    as.matrix() %>%
    phangorn::phyDat(type = "USER", levels = c(0, 1))

  attr(py_data, "id") = d[[1]]$id

  return(py_data)
}

parse_data_files =
  function(
    files,
    cutoff_deep = 20,
    vaf_threshold = 0.05,
    minit = 100,
    cna_data = NULL,
    estimate_mm = TRUE,
    purity_estimates = NULL,
    use = NULL,
    ...
  ) {

  checkmate::assertFileExists(files, "r")
  checkmate::assertNumber(cutoff_deep, lower = 0)
  checkmate::assertNumber(vaf_threshold, lower = 0, upper = 1)
  checkmate::assertNumber(minit, lower = 0, finite = TRUE)
  checkmate::assert_list(cna_data, names = "named")
  checkmate::assertFlag(estimate_mm)
  checkmate::assertNumeric(purity_estimates, lower = 0, upper = 1, null.ok = TRUE)

  # load data and find deep and shallow sequencing data
  data = Map(
    MLLPT::load_genotyping_file,
    file = as.list(files),
    use = list(use),
    cn_data = magrittr::set_names(cna_data[names(files)], names(files))
  )

  # check that all ids are equal
  ids = data[[1]]$id
  for (i in seq_along(data)[-1])
    checkmate::assertTRUE(all(data[[i]]$id == ids))

  # split into deep and shallow
  coverage = vapply(data, function(x) mean(x$depth), numeric(1))
  data_deep = data[coverage > cutoff_deep]
  data_lp = data[coverage <= cutoff_deep]

  # identify regions with equal cn in all samples:
  cn_matrix_deep = do.call(cbind, lapply(data_deep, "[[", "copy_number"))
  rownames(cn_matrix_deep) = ids
  all_equal = apply(cn_matrix_deep, 1, function(x) all(x == x[1]))
  wh_use = all_equal & dplyr::between(cn_matrix_deep[, 1], 1, 4)
  for (i in seq_along(data_lp)) data_lp[[i]] = data_lp[[i]][wh_use, ]
  for (i in seq_along(data_deep)) data_deep[[i]] = data_deep[[i]][wh_use, ]
  cn_deep = cn_matrix_deep[wh_use, 1]

  if (!is.null(purity_estimates) & estimate_mm) {
    # ML estimate of multiplicity
    # assuming that all samples have same MM if CN is equal
    x = do.call(cbind, lapply(data_deep, "[[", "alt_count"))
    n = do.call(cbind, lapply(data_deep, "[[", "depth"))
    cn = do.call(cbind, rep(list(cn_deep), length(data_deep)))
    v = do.call(rbind, rep(list(purity_estimates[names(data_deep)]), NROW(x)))

    wh_mask = x / n < vaf_threshold
    x[wh_mask] = n[wh_mask] = 0

    .lik = function(mm) {
      p = (mm * v) / (2 * (1 - v) + v * cn)
      p[mm > cn] = NA
      stats::dbinom(x, n, p, log = 1)
    }

    mm_lik = abind::abind(lapply(seq(1, max(cn)), .lik), along = 3)
    mm_lik_m = apply(mm_lik, c(1, 3), sum, na.rm = TRUE)
    mm_lik_m[mm_lik_m == 0] = NA
    mm = apply(mm_lik_m, 1, function(x) which(x == max(x, na.rm = TRUE))[1])

    multiplicity_data = data.frame(
      n_alt = apply(x, 1, sum),
      depth = apply(n, 1, sum),
      multiplicity = mm,
      copy_number = cn[, 1],
      row.names = data_deep[[1]]$id
    )

    for (i in seq_along(data_lp))
      data_lp[[i]]$mm = mm

  } else {
    multiplicity_data = NULL
  }

  # Only keep LP samples with CN equal to main CN sate
  for (i in seq_along(data_lp))
    data_lp[[i]] = data_lp[[i]] %>%
      dplyr::filter(copy_number == cn_deep & !is.na(mm))

  for (i in seq_along(data_lp))
    stopifnot(all(data_lp[[i]]$copy_number >= data_lp[[i]]$mm))


  # print warning if very few sites are keeped
  frac_equal_cn = mean(wh_use)
  if (frac_equal_cn < 0.1) {
    paste0(
      "Only ",
      round(frac_equal_cn * 100),
      "% of sites with equal CN between 1 and 4!\n"
    ) %>% warning()
  }


  # get tree and phyDat
  deep_tree_results =
    construct_deep_tree(
      data_deep,
      vaf_threshold = vaf_threshold,
      minit = minit
    )

  if (inherits(deep_tree_results$tree, "multiPhylo")) {
    trees_deep = lapply(deep_tree_results$tree, function(x) x)
  } else {
    trees_deep = list(deep_tree_results$tree)
  }

  # add samples to lp tree(s):
  lp_tree_results = list()
  for (i in seq_along(trees_deep)) {
    lp_tree_results[[i]] =
      MLLPT::add_lowpass_sampled(
        tree = trees_deep[[i]],
        phydata = deep_tree_results$phydat,
        sample_data = data_lp,
        min_confidence = 0,
        min_edge_length = 0,
        return_details = TRUE,
        ...
      )
  }

  return(list(
    tree_deep = trees_deep,
    tree_lp = lapply(lp_tree_results, "[[", "tree"),
    tree_deep_results = deep_tree_results,
    tree_lp_results = lp_tree_results,
    cn_matrix = cn_matrix_deep,
    frac_snvs_used = frac_equal_cn,
    multiplicity_data = multiplicity_data,
    mutation_data_deep = data_deep,
    mutation_data_lp = data_lp,
    purity_estimates = purity_estimates
  ))
}


remove_clonal_variants_tree =
  function(tree) {
    root = phangorn::getRoot(ape::root.phylo(tree, "GL"))
    wh_edge = tree$edge[, 1] == root
    tree$edge.length[wh_edge] = 0
    return(tree)
  }

get_id_mutated_before_edge = function(tree, muts_edge) {

  muts_before_edge = list()
  r_node = phangorn::getRoot(tree)

  for (i in seq_len(NROW(tree$edge))) {
    a_node = tree$edge[i, 1]
    path = ape::nodepath(tree, r_node, a_node)
    .get_edge = function(a, b) which(tree$edge[, 1] == a & tree$edge[, 2] == b)
    edges = sapply(seq_along(path)[-1], function(i) .get_edge(path[i - 1], path[i]))
    muts_before_edge[[as.character(i)]] = unique(unlist(muts_edge[as.character(edges)]))
  }

  return(muts_before_edge)
}

calc_post_lik_data = function(data_set) {

  results = data_set$tree_lp_results[[1]]
  mdata = c(data_set$mutation_data_deep, data_set$mutation_data_lp)

  # tree data:
  tree = results$inital_values$tree
  mutated_on_edge = results$mutations_per_edge
  mutated_before_edge = get_id_mutated_before_edge(tree, mutated_on_edge)


  # lookup tables
  edge_to_node = tree$edge[,2]
  names(edge_to_node) = seq_along(edge_to_node)

  node_to_edge = as.numeric(names(edge_to_node))
  names(node_to_edge) = edge_to_node
  node_to_edge = node_to_edge[order(as.numeric(names(node_to_edge)))]

  muts_to_edge =
    mutated_on_edge[sapply(mutated_on_edge, length)>1] %>%
    reshape::melt() %>%
    dplyr::filter(!value %in% value[duplicated(value)]) %>%
    (function(x) magrittr::set_names(x$L1, x$value))


  # check for missin mutation data
  missing =  mean(!unique(unlist(mutated_on_edge)) %in% mdata[[1]]$id)
  if (sum(missing)) {
    cat(paste0("Missing ", signif(missing*100, 2), "% of mutation data.\n"))
  }


  # calculate post-prop for each sample:
  results_all = list()

  for (sample in names(mdata)) {

    # get all variant data:
    mdata_smp = mdata[[sample]] %>%
      dplyr::mutate(edge = muts_to_edge[id]) %>%
      dplyr::mutate(prior_m = NA)

    # get ml estimates for position:
    if (sample %in% rownames(results$per_sample_results)) {
      prior_m = results$per_sample_results[sample,"pi"]
      edge = results$per_sample_results[sample,"edge"]
      purity = results$per_sample_results[sample,"purity"]
      vaf_um = results$per_sample_results[sample,"background_vaf"]
      depth = results$per_sample_results[sample,"background_vaf"]
    } else {
      prior_m = 1
      edge = node_to_edge[as.character(which(tree$tip.label == sample))]
      purity = results$purity_estimates[sample]
      if (is.null(purity)) purity = 1
      vaf_um = 0.01
    }

    # set prior lik:
    mdata_smp$prior_m[mdata_smp$id %in% c(unlist(mutated_before_edge), unlist(mutated_on_edge))] = 0
    mdata_smp$prior_m[mdata_smp$id %in% mutated_before_edge[[edge]]] = 1
    mdata_smp$prior_m[mdata_smp$id %in% mutated_on_edge[[edge]]] = prior_m

    # calculate posterior lik:
    vaf_m = with(mdata_smp, mm * purity / (purity * copy_number + (1-purity) * 2)) # expected vaf for each mutation
    mdata_smp$log_lik_m = dbinom(mdata_smp$alt_count, mdata_smp$depth, prob=vaf_m, log = TRUE)
    mdata_smp$log_lik_um = dbinom(mdata_smp$alt_count, mdata_smp$depth, prob=vaf_um, log = TRUE)
    mdata_smp$pD = with(mdata_smp, exp(log_lik_m) * prior_m + exp(log_lik_um) * (1 - prior_m))
    mdata_smp$post_lik_m = with(mdata_smp, exp(log_lik_m) * prior_m / pD)

    results_all[[sample]] = list(data=mdata_smp, purity=purity, vaf_bkg=vaf_um, edge=edge, sample=sample)
  }

  return(results_all)
}


calc_expected_vaf_density = function(x) {

  # predict denisty
  mdata_smp = x$data %>% dplyr::count(depth, copy_number, prior_m, mm, edge)
  all_data = NULL

  for (i in seq_len(NROW(mdata_smp))) {
    if (is.na(mdata_smp$edge[i])) next()
    for (n in seq(0, mdata_smp$depth[i])) {

      N = mdata_smp$depth[i]
      cn = mdata_smp$copy_number[i]
      p_m = mdata_smp$prior_m[i]
      mm = mdata_smp$mm[i]
      edge = mdata_smp$edge[i]
      v = x$purity

      vaf_m = (mm * v) / (v * cn + (1-v) * 2) # expected vaf for each mutation
      d0 = dbinom(n, N, scales::squish(x$vaf_bkg, c(1e-8, 1-1e-8)))
      d0[d0 == 0] = .Machine$double.xmin
      d1 = dbinom(n, N, scales::squish(vaf_m, c(1e-8, 1-1e-8)))
      d1[d1 == 0] = .Machine$double.xmin

      d = (d1 * p_m + d0 * (1 - p_m)) * mdata_smp$n[i]
      dc = data.frame(alt_count=n, depth=N, d=d, edge=edge, row.names = NULL)
      all_data = rbind(all_data, dc)
    }
  }

  all_data %>%
    split(., .$edge) %>%
    lapply(dplyr::mutate, d=d / sum(d)) %>%
    do.call(what=rbind) %>%
    dplyr::mutate(d = d * tapply(mdata_smp$n, mdata_smp$edge, sum)[as.character(edge)]) %>%
    dplyr::mutate(edge = factor(edge, 1:100, paste0("Edge ", 1:100)))
}


plot_vaf_post_lik_data = function(x, show_expected=TRUE) {

  caption = paste0(
    "Avg. depth: ", signif(mean(x$data$depth), 2),
    ", Purity: ",  signif(x$purity, 3),
    "; Background VAF: ", signif(x$vaf_bkg, 3)
  )

  if (mean(x$data$depth) < 10 & show_expected) {
    expected_density =
      calc_expected_vaf_density(x) %>%
      dplyr::mutate(vaf=alt_count / depth) %>%
      dplyr::group_by(vaf, edge) %>%
      dplyr::summarise(d=sum(d)) %>%
      dplyr::arrange(-vaf) %>%
      dplyr::filter(!is.na(vaf))
  } else {
    expected_density =
      NULL
  }


  plot_of_data =
    x$data %>%
    dplyr::filter(depth > 0) %>%
    dplyr::mutate(edge = factor(edge, 1:100, paste0("Edge ", 1:100))) %>%
    ggplot(aes(x=alt_count / depth)) +
      cowplot::theme_cowplot() +
      facet_wrap(~edge, scales = "free_y") +
      geom_histogram(bins=20, aes(fill=signif(prior_m, 3))) +
      ggtitle(x$sample) +
      xlab("VAF") +
      ylab("Count") +
      labs(fill="p(M = 1)", color="", caption = caption) +
      scale_fill_viridis_c(na.value="gray20", end = 0.8) +
      scale_color_manual(values = "red")

  if (!is.null(expected_density)) {
    plot_of_data =
      plot_of_data +
      geom_freqpoly(
        data = expected_density,
        aes(x = vaf, weight = d, color = "Expected"),
        alpha = 0.8,
        bins = 20
      )
  }

  return(plot_of_data)
}

plot_post_lik_data = function(x) {

  caption = paste0(
    "Avg. depth: ", signif(mean(x$data$depth), 2),
    ", Purity: ",  signif(x$purity, 3),
    "; Background VAF: ", signif(x$vaf_bkg, 3)
  )

  plot_of_data =
    x$data %>%
    dplyr::mutate(edge = factor(edge, 1:100, paste0("Edge ", 1:100))) %>%
    dplyr::filter(!is.na(edge)) %>%
    ggplot(aes(x = post_lik_m)) +
      cowplot::theme_cowplot() +
      facet_wrap( ~ edge, scales = "free_y") +
      geom_histogram(bins = 20) +
      ggtitle(x$sample) +
      xlab("Posterior likelihood p(M=1|D)") +
      ylab("Count") +
      scale_fill_brewer(palette = "Set1", na.value = "gray20") +
      labs(caption = caption)

  return(plot_of_data)
}

