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

