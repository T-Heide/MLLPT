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
        max_loss_frac = 0.1
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
    expect_type({
      add_lowpass_sampled(
        tree = tree,
        phydata = phy_data,
        sample_data = head(lp_data, 2),
        return_details = FALSE
      )
    }, "list")


  expect_s3_class(plot_lp_loglik(lp_adding_results), "ggplot")
  expect_s3_class(plot_lp_loglik_edge(lp_adding_results), "ggplot")

  purity = purity_vals[lp_adding_results$per_sample_results$sample]
  labeller = function(x) gsub("_L1", "", gsub("EPICC_", "", x))
  expect_s3_class(plot_sample_data(lp_adding_results, labeller, purity), "ggplot")
  expect_s3_class(plot_sample_data(lp_adding_results_plus_loss), "ggplot")

  expect_s3_class(plot_tree(lp_adding_results$tree), "ggplot")
  expect_s3_class(plot_tree(lp_adding_results_no_details), "ggplot")


})
