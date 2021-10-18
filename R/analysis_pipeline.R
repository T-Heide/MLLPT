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
    mutate(GL=0) %>%
    t() %>% as.matrix() %>%
    phangorn::phyDat(type="USER", levels=c(0,1))

  attr(pyData, "id") = d[[1]]$id

  return(pyData)
}

parse_data_files = function(files, cutoff_deep=10, vaf_threshold=0.05, minit = 100, cna_data=cna_data, ...) {

  checkmate::assertFileExists(files, "r")
  checkmate::assertNumber(cutoff_deep, lower = 0)

  # load data and find deep and shallow sequencing data
  data = mapply(
    MLLPT::load_genotyping_file,
    cn_data = magrittr::set_names(cna_data[names(files)], names(files)),
    file = as.list(files),
    SIMPLIFY = FALSE
  )

  # check that all ids are equal
  ids = data[[1]]$id
  for (i in seq_along(data)[-1])
    checkmate::assertTRUE(all(data[[i]]$id == ids))


  # identify regions with equal cn in all samples:
  cn_matrix = do.call(cbind, lapply(data, "[[", "copy_number"))
  all_equal = apply(cn_matrix, 1, function(x) all(x == x[1]))
  wh_use = all_equal & between(cn_matrix[,1], 1, 4)
  for (i in seq_along(data)) data[[i]] = data[[i]][all_equal,]

  frac_equal_cn = mean(wh_use)
  if (frac_equal_cn < 0.1) {
    paste0(
      "Only ",
      round(frac_equal_cn * 100),
      "% of sites with equal CN between 1 and 4!\n"
    ) %>% warning()
  }

  # split into deep and shallow
  coverage = sapply(data, function(x) mean(x$depth))
  data_deep = data[coverage > cutoff_deep]
  data_lp = data[coverage <= cutoff_deep]


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
    tree_lp = sapply(lp_tree_results, "[[", "tree"),
    tree_deep_results = deep_tree_results,
    tree_lp_results = lp_tree_results,
    cn_matrix = cn_matrix,
    frac_snvs_used = frac_equal_cn
  ))
}
