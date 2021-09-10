test_that("lp assignment method works", {

  load("data.rda")

  lp_adding_results =
    expect_type({
      add_lowpass_sampled(
        tree = tree,
        phydata = phy_data,
        sample_data = head(lp_data, 2),
        return_details = TRUE
      )
    }, "list")

  lp_adding_results_plus_loss =
    expect_type({
      add_lowpass_sampled(
        tree = tree,
        phydata = phy_data,
        sample_data = head(lp_data, 2),
        return_details = TRUE,
        max_loss_frac = 0.5
      )
    }, "list")

  lp_adding_results_no_optim =
    expect_type({
      add_lowpass_sampled(
        tree = tree,
        phydata = phy_data,
        sample_data = head(lp_data, 2),
        return_details = TRUE,
        optimize_values = FALSE
      )
    }, "list")


  lp_adding_results_no_details =
    expect_s3_class({
      add_lowpass_sampled(
        tree = tree,
        phydata = phy_data,
        sample_data = head(lp_data, 2),
        return_details = FALSE
      )
    }, "phylo")


  # alternative column names
  lp_data_alt = lp_data
  for (i in seq_along(lp_data_alt)) {
    colnames(lp_data_alt[[i]]) = c("clone","alt","dp","id","ccf","vaf","cn")
  }

  lp_adding_results_alt_names =
    expect_type({
      add_lowpass_sampled(
        tree = tree,
        phydata = phy_data,
        sample_data = head(lp_data_alt, 2),
        return_details = TRUE,
        max_vaf_bkgr=0.05
      )
    }, "list")


  # warning when using very high background rate?
  expect_warning({
    add_lowpass_sampled(
      tree = tree,
      phydata = phy_data,
      sample_data = head(lp_data, 1),
      return_details = TRUE,
      max_vaf_bkgr = 0.5
    )
  }, "Do you really want to set a background rate of")


  # failing when not including essential data?
  for (drop in c("alt","depth","cn")) {

    lp_data_alt = lp_data

    for (i in seq_along(lp_data_alt))
      lp_data_alt[[i]][,drop] = NULL

    expect_error({
      add_lowpass_sampled(
        tree = tree,
        phydata = phy_data,
        sample_data = head(lp_data_alt, 2)
      )
    }, "missing")

  }

  expect_s3_class(plot_lp_loglik(lp_adding_results), "ggplot")
  expect_s3_class(plot_lp_loglik_edge(lp_adding_results), "ggplot")

  purity = purity_vals[lp_adding_results$per_sample_results$sample]
  labeller = function(x) gsub("_L1", "", gsub("EPICC_", "", x))
  expect_s3_class(plot_sample_data(lp_adding_results, labeller, purity), "ggplot")
  expect_s3_class(plot_sample_data(lp_adding_results_plus_loss), "ggplot")

  expect_s3_class(plot_tree(lp_adding_results$tree), "ggplot")
  expect_s3_class(plot_tree(lp_adding_results_no_details), "ggplot")


})


test_that("mutation count data files can be read", {

  # objects for cn annotation:
  cn_df = data.frame(
    id = c("chr1:1_A/T", "chr1:2_C/G", "chr1:3_G/C", "chr1:4_T/A"),
    cn = 1:4,
    mm = 1:4
  )

  cn_gr = GenomicRanges::GRanges("chr1", IRanges::IRanges(1:4, 1:4))
  S4Vectors::elementMetadata(cn_gr) = data.frame(cn=1:4, mm=1:4)

  # try loading with default
  tab_file = system.file("extdata", "test_mutation.table", package="MLLPT")
  res = expect_warning(load_genotyping_file(tab_file), ".*CN data missing.*")
  expect_known_hash(res, "437c8be539")

  # modify result so it is equal to the expected using the above annotations
  res$copy_number = 1:4
  res$mm = 1:4

  # 1) Test using a Granges object as input
  expect_equal(load_genotyping_file(tab_file, cn_data=cn_gr), res)

  # 2) Test using a data frame as input
  expect_equal(load_genotyping_file(tab_file, cn_data=cn_df), res)
  expect_equal(load_genotyping_file(tab_file, cn_data=head(cn_df, 3)), dplyr::filter(res, pos %in% 1:3))
  res2 = expect_warning(load_genotyping_file(tab_file, cn_data=dplyr::select(cn_df, -mm)))
  expect_equal(res2, dplyr::mutate(res, mm=1))

})
