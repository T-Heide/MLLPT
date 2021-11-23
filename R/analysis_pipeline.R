phyDat_to_tree = function(phyDat, all_trees=FALSE, ...) {
  phyDat %>%
    phangorn::pratchet(trace=FALSE, all=all_trees, ...) %>%
    ape::unroot() %>%
    ape::root(out="GL", resolve.root=TRUE) %>%
    phangorn::acctran(phyDat)
}

construct_deep_tree = function(d, vaf_threshold=0.05, minit=100) {
  phydat = get_phydat_from_genotyping(d, vaf_threshold=vaf_threshold)
  tree = phyDat_to_tree(phydat, all_trees=TRUE, minit=minit, k=minit)
  return(list(tree=tree, phydat=phydat))
}

get_phydat_from_genotyping = function(d, vaf_threshold=0.05) {

  # check input
  checkmate::assertNumber(vaf_threshold, lower = 0, upper = 1)
  checkmate::assertList(d, names="named")
  for (i in seq_along(d)) checkmate::assertDataFrame(d[[i]])
  for (i in seq_along(d)) checkmate::assertSubset(c("id","vaf"), colnames(d[[i]]))
  for (i in seq_along(d)) checkmate::assertTRUE(all(d[[1]]$id == d[[i]]$id))


  # convert to phyDat object:
  pyData =
    lapply(d, function(x) x$vaf >= vaf_threshold) %>%
    do.call(what=cbind) %>%
    data.frame() %>%
    dplyr::mutate(GL=0) %>%
    t() %>% as.matrix() %>%
    phangorn::phyDat(type="USER", levels=c(0,1))

  attr(pyData, "id") = d[[1]]$id

  return(pyData)
}

parse_data_files = function(files, cutoff_deep=10, vaf_threshold=0.05, minit = 100, cna_data=NULL, estimate_mm=TRUE, purity_estimates=NULL, ...) {

  checkmate::assertFileExists(files, "r")
  checkmate::assertNumber(cutoff_deep, lower = 0)
  checkmate::assertNumber(vaf_threshold, lower = 0, upper=1)
  checkmate::assertNumber(minit, lower = 0, finite = TRUE)
  checkmate::assert_list(cna_data, names = "named")
  checkmate::assertFlag(estimate_mm)
  checkmate::assertNumeric(purity_estimates, lower = 0, upper=1, null.ok = TRUE)

  # load data and find deep and shallow sequencing data
  data = Map(
    MLLPT::load_genotyping_file,
    file = as.list(files),
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
  all_equal = apply(cn_matrix_deep, 1, function(x) all(x == x[1]))
  wh_use = all_equal & dplyr::between(cn_matrix_deep[,1], 1, 4)
  for (i in seq_along(data_lp)) data_lp[[i]] = data_lp[[i]][wh_use,]
  for (i in seq_along(data_deep)) data_deep[[i]] = data_deep[[i]][wh_use,]
  cn_deep = cn_matrix_deep[wh_use,1]

  if (!is.null(purity_estimates) & estimate_mm) {
    # ML estimate of multiplicity
    # assuming that all samples have same MM if CN is equal
    x = do.call(cbind, lapply(data_deep, "[[", "alt_count"))
    n = do.call(cbind, lapply(data_deep, "[[", "depth"))
    cn = do.call(cbind, rep(list(cn_deep), length(data_deep)))
    v = do.call(rbind, rep(list(purity_estimates[names(data_deep)]), NROW(x)))

    wh_mask = x/n < vaf_threshold
    x[wh_mask] = n[wh_mask] = 0

    .lik = function(mm) {
      p = (mm * v) / (2 * (1 - v) + v * cn)
      p[mm > cn] = NA
      dbinom(x, n, p, log = 1)
    }

    mm_lik = abind::abind(lapply(seq(1, max(cn)), .lik), along = 3)
    mm_lik_m = apply(mm_lik, c(1,3), sum, na.rm=TRUE)
    mm_lik_m[mm_lik_m == 0] = NA
    mm = apply(mm_lik_m, 1, function(x) which(x==max(x, na.rm=TRUE))[1])

    multiplicity_data = data.frame(
      n_alt = apply(x, 1, sum),
      depth = apply(n, 1, sum),
      multiplicity = mm,
      copy_number = cn[, 1]
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
      vaf_threshold=vaf_threshold,
      minit=minit
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
        return_details = TRUE
      )
  }

  return(list(
    tree_deep = trees_deep,
    tree_lp = lapply(lp_tree_results, "[[", "tree"),
    tree_deep_results = deep_tree_results,
    tree_lp_results = lp_tree_results,
    cn_matrix = cn_matrix_deep,
    frac_snvs_used = frac_equal_cn,
    multiplicity_data = multiplicity_data
  ))
}


remove_clonal_variants_tree =
  function(tree) {
    root = phangorn::getRoot(ape::root.phylo(tree, "GL"))
    wh_edge = tree$edge[,1] == root
    tree$edge.length[wh_edge] = 0
    return(tree)
  }
