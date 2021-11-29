#' Plot showing log-likelihood of all possible edges.
#'
#' @param d Object returned from 'add_lowpass_sampled'.
#' @param labeller_function (optional) function returning a modified version of sample labels.
#' @param subset (optional) character vector defining a subset of samples to plot.
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
#'
plot_lp_loglik =
  function(d, labeller_function=NULL, subset=NULL) {


    if ("min_confidence" %in% names(d$inital_value)) {
      min_confidence = d$inital_values$min_confidence
    } else {
      min_confidence = 1
    }

    if (d$inital_values$n_bootstraps) {
      d_edge =
        lapply(d$bootstrap_results, function(x) {
          vapply(
            names(d$max_ll_per_edge[[1]]),
            function(y) log((sum(x$edge == y) + 1) / (length(x$edge) + 1)),
            numeric(1)
          )
        })
    } else {
      d_edge = d$max_ll_per_edge
    }

    if (!is.null(subset)) {
      d_edge = d_edge[names(d_edge) %in% subset]
    }

    d_plot =
      lapply(d_edge, function(x) {
        p = exp(x - max(x))
        p / sum(p)
      }) %>%
      do.call(what = cbind) %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("edge", "sample", "p")) %>%
      dplyr::mutate(is_max = c(p == tapply(p, sample, max)[as.character(sample)])) %>%
      dplyr::mutate(is_accept = p > min_confidence) %>%
      dplyr::mutate(edge = factor(edge, unique(sort(edge)), ordered = TRUE)) %>%
      dplyr::mutate(shape = ifelse(is_accept & is_max, 8, ifelse(is_max, 20, NA)))

    if (!is.null(labeller_function)) {
      d_plot$sample = labeller_function(as.character(d_plot$sample))
    }

    plot_p_edge =
      d_plot %>%
      ggplot(aes(x = sample, y = edge, fill = p + 1e-16)) +
      cowplot::theme_cowplot() +
      geom_tile() +
      geom_point(aes(shape = shape), size = 1.2) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      xlab("") +
      ylab("Tree edge") +
      labs(fill = "P(Edge)") +
      scale_fill_continuous(trans = "log10", low = "gray98", high = "red") +
      scale_shape_identity()

    plot_p_edge
  }


#' Plot showing log-likelihood of fraction mutated from most likely edge.
#'
#' @param d Object returned from 'add_lowpass_sampled'.
#' @param labeller_function (optional) function returning a modified version of sample labels.
#' @param subset (optional) character vector defining a subset of samples to plot.
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
#'
plot_lp_loglik_edge =
  function(d, labeller_function = NULL, subset = NULL) {

    per_sample_res = d$per_sample_results

    if (!is.null(subset)) {
      per_sample_res = per_sample_res %>% dplyr::filter(sample %in% subset)
    }

    if (!is.null(labeller_function)) {
      per_sample_res$sample_mod = labeller_function(per_sample_res$sample)
    } else {
      per_sample_res$sample_mod = per_sample_res$sample
    }

    # extract ll data of best edge
    lld =
      lapply(seq_len(NROW(per_sample_res)), function(i) {
        smp = per_sample_res$sample[i]
        smp_m = per_sample_res$sample_mod[i]
        edge =  per_sample_res$edge[i]
        lapply(names(d$ll_per_edge[[smp]]), function(edge_c) {
          cdat = d$ll_per_edge[[smp]][[edge_c]]
          data.frame(
            sample_mod = smp_m,
            edge = edge_c,
            best = (edge_c == edge),
            frac = seq_along(cdat) / length(cdat),
            llik = cdat
          )
        }) %>% do.call(what = rbind)
      }) %>% do.call(what = rbind)

    #
    edge_data_per_smp =
      lld %>%
      split(.$sample) %>%
      lapply(dplyr::mutate, p = exp(llik - max(llik))) %>%
      lapply(dplyr::mutate, p = p / sum(p)) %>%
      do.call(what = rbind)

    per_case_max_pos =
      per_sample_res %>%
      dplyr::mutate(sample_label = paste0("Edge ", edge, " - ", sample_mod)) %>%
      dplyr::select(frac = pi, sample_mod = sample_label)

    anc_list = phangorn::Ancestors(d$inital_values$tree, type = "all")
    for (i in seq_along(anc_list)) anc_list[[i]] = c(rev(anc_list[[i]]), i)
    node_order = unique(unlist(anc_list))
    edge_order = order(match(d$inital_values$tree$edge[, 2], node_order), decreasing = FALSE)

    edge_label =
      paste0(
        seq_along(d$inital_values$tree$edge[, 1]),
        " (", d$inital_values$tree$edge[, 1],
        "-> ", d$inital_values$tree$edge[, 2], ")"
      ) %>% magrittr::set_names(seq_along(.))

    edge_label = edge_label[edge_order]

    per_sample_res =
      per_sample_res %>%
      dplyr::mutate(edge = factor(edge, names(edge_label), edge_label, ordered = TRUE))

    d_plot =
      edge_data_per_smp %>%
      dplyr::mutate(edge = factor(edge, names(edge_label), edge_label, ordered = TRUE)) %>%
      ggplot(aes(x = frac, y = p + 1e-16, color = best)) +
        cowplot::theme_cowplot() +
        geom_line() +
        facet_grid(sample_mod~edge, scales = "free_y") +
        xlab("Position on edge") +
        ylab("Likelihood") +
        geom_vline(data = per_sample_res, aes(xintercept = pi), linetype = 3) +
        scale_color_manual(values = c("black", "red"), breaks = c("FALSE", "TRUE")) +
        guides(color = "none") +
        scale_x_continuous(n.breaks = 2)

    d_plot
  }


#' Plot showing distribution of edge positions during boostraping.
#'
#' @param d Object returned from 'add_lowpass_sampled'.
#' @param labeller_function (optional) function returning a modified version of sample labels.
#' @param subset (optional) character vector defining a subset of samples to plot.
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
#'
plot_lp_position_edge =
  function(d, labeller_function=NULL, subset=NULL) {

    if (d$inital_values$n_bootstraps == 0) {
      warning(
        "Need bootstrap data to use",
        sQuote("plot_lp_position_edge"),
        ".\n",
        " => Please call ",
        sQuote("add_lowpass_sampled"),
        " with ",
        sQuote("n_bootstraps > 0"),
        ".\n"
      )


      pl =
        ggplot(NULL) +
        cowplot::theme_cowplot()

      return(pl)
    }

    data_ml = d$per_sample_result
    data_bootstrap = reshape2::melt(d$bootstrap_results, measure.vars = vector())

    if (!is.null(subset)) {
      data_bootstrap = data_bootstrap %>% dplyr::filter(L1 %in% subset)
      data_ml = data_ml %>% dplyr::filter(sample %in% subset)
    }

    if (!is.null(labeller_function)) {
      data_bootstrap$sample_mod = labeller_function(data_bootstrap$L1)
      data_ml$sample_mod = labeller_function(data_ml$sample)
    } else {
      data_bootstrap$sample_mod = data_bootstrap$L1
      data_ml$sample_mod = data_ml$sample
    }

    # modify edge labels
    edge_lab =
      paste0(
        seq_along(d$inital_values$tree$edge[, 1]),
        " (", d$inital_values$tree$edge[, 1],
        "-> ", d$inital_values$tree$edge[, 2], ")"
      ) %>% magrittr::set_names(seq_along(.))


    anc_list = phangorn::Ancestors(d$inital_values$tree, type = "all")
    for (i in seq_along(anc_list)) anc_list[[i]] = c(rev(anc_list[[i]]), i)
    node_order = unique(unlist(anc_list))
    edge_order = order(match(d$inital_values$tree$edge[, 2], node_order), decreasing = FALSE)
    edge_lab = edge_lab[edge_order]

    data_bootstrap = data_bootstrap %>%
      dplyr::mutate(edge = factor(edge, names(edge_lab), edge_lab, ordered = TRUE))

    data_ml = data_ml %>%
      dplyr::mutate(edge = factor(edge, names(edge_lab), edge_lab, ordered = TRUE))

    d_plot =
      data_bootstrap %>%
      ggplot(aes(x = pos, after_stat(count) / d$inital_values$n_bootstraps)) +
      cowplot::theme_cowplot() +
      geom_histogram(breaks = seq(0, 1, by = 0.05)) +
      facet_grid(sample_mod ~ edge, drop = FALSE, scales = "free_y") +
      xlab("Position on edge") +
      ylab("Density") +
      geom_vline(data = data_ml, aes(xintercept = pi), linetype = 3, color = "red") +
      guides(color = "none") +
      scale_x_continuous(n.breaks = 2, limits = c(0, 1))

    d_plot
  }



#' Plot showing ML estimate of sample parameters.
#'
#' @param d Object returned from 'add_lowpass_sampled'.
#' @param labeller_function (optional) function returning a modified version of sample labels.
#' @param external_purity_estimate (optional) numeric vector of external purity estimates of each sample.
#' @param label_external_purity (optional) label to be used in the legend of the external purity estimates.
#' @param subset (optional) character vector defining a subset of samples to plot.
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
#'
plot_sample_data =
  function(
    d,
    labeller_function = NULL,
    external_purity_estimate = NULL,
    label_external_purity = "Independent estimate",
    subset = NULL
  ) {

    if (is.null(labeller_function)) {
      labeller_function = function(x) x
    }

    if (!d$inital_values$optimize_values) {
      return(NULL)
    }

    per_sample_results = d$per_sample_results

    if (!is.null(external_purity_estimate)) {
      stopifnot(length(external_purity_estimate) == NROW(per_sample_results))
      per_sample_results$external_estimate = external_purity_estimate
    }

    per_sample_results =
      per_sample_results %>%
      dplyr::mutate(initial_purity = d$inital_values$purity[idx]) %>%
      dplyr::filter(sample %in% subset | is.null(subset)) %>%
      dplyr::mutate(sample = labeller_function(sample)) %>%
      dplyr::mutate(sample = factor(sample, rev(sort(unique(sample))), ordered = TRUE))

    if (d$inital_values$n_bootstraps) {
      ci_est = d$bootstrap_results %>%
        lapply(dplyr::select, -edge) %>%
        lapply(apply, 2, stats::quantile, c(0.025, 0.975), names = FALSE) %>%
        reshape2::melt() %>%
        dplyr::mutate(Var1 = factor(Var1, c(1, 2), c("lower", "upper"))) %>%
        reshape2::dcast(L1 ~ Var1 + Var2) %>%
        dplyr::mutate(sample = labeller_function(L1))
    } else {
      ci_est = NULL
    }

    bkgr_vaf_data = data.frame(x = d$inital_values$vaf_bkgr)
    limit_y = min(c(per_sample_results$background_vaf, bkgr_vaf_data$x), na.rm = TRUE) / 10

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    plot_vaf =
      per_sample_results %>%
      ggplot(aes(x = sample, y = background_vaf + 1e-6)) +
      cowplot::theme_cowplot()

    if (!is.null(ci_est)) {
      plot_vaf = plot_vaf +
        geom_errorbar(
          data = ci_est,
          aes(x = sample, ymin = lower_bkg + 1e-6, ymax = upper_bkg + 1e-6),
          inherit.aes = FALSE,
          alpha = 0.8,
          width = 0.2
        )
    }

    plot_vaf = plot_vaf +
      geom_point() +
      geom_hline(
        data = bkgr_vaf_data,
        aes(
          yintercept = x,
          linetype = "Initial value")
        ) +
      coord_flip() +
      scale_y_log10(
        limits = c(limit_y + 1e-6, d$inital_values$max_vaf_bkgr + 1e-6),
        n.breaks = 3
      ) +
      scale_linetype_manual(values = 4) +
      theme(legend.position = "bottom") +
      labs(linetype = "") +
      ylab("Error rate") +
      xlab("") +
      cowplot::background_grid()


    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    plot_purity =
      per_sample_results %>%
      ggplot(aes(x = sample, y = purity)) +
      cowplot::theme_cowplot() +
      geom_segment(aes(xend = sample, yend = initial_purity), linetype = 3)

    if (!is.null(ci_est)) {
      plot_purity = plot_purity +
        geom_errorbar(
          data = ci_est,
          aes(x = sample, ymin = lower_purity, ymax = upper_purity),
          inherit.aes = FALSE,
          alpha = 0.8,
          width = 0.2
        )
    }

    plot_purity = plot_purity +
      geom_point(size = 2) +
      geom_point(
        aes(
          shape = "Initial value",
          color = "Initial value",
          y = initial_purity
        ), size = 2) +
      scale_y_continuous(n.breaks = 4, limits = c(0, 1)) +
      coord_flip() +
      labs(color = "", shape = "") +
      ylab("Purity") +
      xlab("") +
      theme(legend.position = "bottom") +
      scale_shape_manual(
        values = c(4, 3),
        breaks = c("Initial value", label_external_purity)
      ) +
      scale_color_manual(
        values = c("red", "red"),
        breaks = c("Initial value", label_external_purity)
      ) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      cowplot::background_grid()

    if (!is.null(external_purity_estimate)) {

      plot_purity = plot_purity +
        geom_segment(
          aes(
            xend = sample,
            yend = external_estimate
          ), linetype = 3
        ) +
        geom_point(
          aes(
            shape = label_external_purity,
            color = label_external_purity,
            y = external_estimate
          ),  size = 2
        )

    }


    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


    plot_loss_frac =
      per_sample_results %>%
      ggplot(aes(x = sample, y = loss_frac)) +
      cowplot::theme_cowplot()

    if (!is.null(ci_est)) {
      plot_loss_frac = plot_loss_frac +
        geom_errorbar(
          data = ci_est,
          aes(x = sample, ymin = lower_loss, ymax = upper_loss),
          inherit.aes = FALSE,
          alpha = 0.8,
          width = 0.2
        )
    }

    plot_loss_frac = plot_loss_frac +
      cowplot::background_grid() +
      geom_point(size = 2) +
      coord_flip() +
      labs(color = "", shape = "") +
      ylab("Loss rate") +
      xlab("") +
      theme(
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      scale_y_continuous(
        n.breaks = 3,
        limits = c(0, NA)
      )

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    plot_dp =
      d$mutation_data %>%
      lapply(do.call, what = rbind) %>%
      lapply(function(x) stats::weighted.mean(x$N, x$weight)) %>%
      reshape2::melt() %>%
      dplyr::mutate(L1 = labeller_function(L1)) %>%
      dplyr::mutate(sample = factor(L1, levels(per_sample_results$sample), ordered = TRUE)) %>%
      dplyr::filter(!is.na(sample)) %>%
      ggplot(aes(x = sample, y = value)) +
      cowplot::theme_cowplot() +
      geom_point() +
      ylim(0, NA) +
      coord_flip() +
      ylab("Coverage") +
      xlab("") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      cowplot::background_grid() +
      scale_y_continuous(n.breaks = 3, limits = c(0, NA))


    if (d$inital_values$max_loss_frac > 0) {
      pg =
        cowplot::plot_grid(
          cowplot::as_grob(plot_vaf),
          cowplot::as_grob(plot_loss_frac),
          cowplot::as_grob(plot_purity),
          cowplot::as_grob(plot_dp),
          align = "h",
          axis = "b",
          rel_widths = c(0.75, 0.45, 0.45, 0.45),
          nrow = 1
        )
    } else {
      pg =
        cowplot::plot_grid(
          cowplot::as_grob(plot_vaf),
          cowplot::as_grob(plot_purity),
          cowplot::as_grob(plot_dp),
          align = "h",
          axis = "b",
          rel_widths = c(0.8, 0.5, 0.5),
          nrow = 1
        )
    }

    return(pg)
  }


#' A internal function used for the plotting of trees
#'
#' @param tree object of class phylo to plot.
#' @param HI (optional) HI of the tree (unused).
#' @param color_by (optional) lookup table with color labels (NULL).
#' @param linewidth (optional) numeric value defining linewidth of tree edges (1).
#' @param pointsize (optional)  numeric value defining point size of tree tips  (1).
#' @param textsize (optional)  numeric value defining linewidth of tree edges (2.5).
#' @param labeller_function (optional) HI of the tree (unused).
#' @param keep_conf_label (optional) Flat indicating if `confidence label' should be keeped (default: FALSE).
#' @param subset (optional) character vector defining a subset of samples to plot.
#'
#' @return ggplot object
#' @export
#' @import ggplot2
#'
#' @examples plot_tree(ape::rtree(5))
plot_tree =
  function(
    tree,
    HI = NULL,
    color_by = NULL,
    linewidth = 1,
    pointsize = 1,
    textsize = 2.5,
    labeller_function = NULL,
    keep_conf_label = FALSE,
    subset = NULL
  ) {


  if ("tree" %in% names(tree)) tree = tree$tree
  checkmate::assertClass(tree, "phylo", null.ok = FALSE)
  checkmate::assertNumber(HI, null.ok = TRUE)
  if (is.factor(color_by)) checkmate::assertFactor(color_by, names = "named", null.ok = TRUE)
    else checkmate::assertCharacter(color_by, names = "named", null.ok = TRUE)
  checkmate::assertNumber(linewidth, null.ok = FALSE, na.ok = FALSE, lower = 0, finite = TRUE)
  checkmate::assertNumber(pointsize, null.ok = FALSE, na.ok = FALSE, lower = 0, finite = TRUE)
  checkmate::assertNumber(textsize, null.ok = FALSE, na.ok = FALSE, lower = 0, finite = TRUE)
  checkmate::assertFunction(labeller_function, null.ok = TRUE)

  if (is.null(labeller_function)) {
    labeller_function = function(x) x
  }

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plot tree
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  added_node = grepl("Added", tree$tip.label)
  tree$tip.label = unlist(tree$tip.label)

  if (!is.null(subset)) {
    wh_keep = gsub(" [(]Added.*", "", tree$tip.label) %in% subset | !added_node
    tree = ape::keep.tip(tree, tree$tip.label[wh_keep])
  }

  #tree$tip.label = gsub(" [(]Added.*[)]", " ", tree$tip.label)
  shape_lookup = ifelse(added_node, "A", "U") %>% magrittr::set_names(tree$tip.label)
  rep_label = ifelse(keep_conf_label, " [(]Added", " [(]Added.*[*]")
  rep_with = ifelse(keep_conf_label, ":", "")

  new_label =
    gsub(")", " *", labeller_function(tree$tip.label)) %>%
    gsub(pattern = rep_label, replacement = rep_with) %>%
    magrittr::set_names(tree$tip.label)

  new_label_base =
    gsub(" [(]Added.*", "", tree$tip.label) %>%
    magrittr::set_names(tree$tip.label)

  if (is.null(color_by)) {
    ids = unique(new_label_base[tree$tip.label])
    color_by = magrittr::set_names(rep("gray10", length(ids)), ids)
    color_scale = scale_color_identity()
  } else {
    color_scale = scale_color_brewer(palette = "Set1", drop = FALSE)
  }

  plt =
    tree %>%
    ggtree::ggtree(
      layout = "rectangular",
      color = "gray10",
      size = linewidth
    ) +
    ggtree::geom_tiplab(
      aes(
        label = new_label[label],
        color = color_by[new_label_base[label]]
      ),
      size = textsize,
      hjust = -0.25
    ) +
    ggtree::geom_tippoint(
      aes(
        shape = shape_lookup[label],
        color = color_by[new_label_base[label]]
      ),
      size = pointsize
    ) +
    ggtree::geom_treescale(
      y = -1,
      x = 0,
      fontsize = textsize,
      linesize = linewidth,
      offset = 0.9
    ) +
    scale_shape_manual(
      breaks = c("U", "A"),
      values = c(U = 16, A = 8)
    ) +
    labs(color = "") +
    guides(shape = "none") +
    guides(color = guide_legend(nrow = 1)) +
    theme(
      title = element_text(size = 8, color = "gray20"),
      legend.position = "bottom",
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA)
    ) +
    color_scale

  if (!is.null(HI)) {
    plt = plt + labs(caption = paste0("HI = ", round(HI, 2)))
  }

  plt = plt + xlim(0, 1.6 * max(plt$data$x))

  return(plt)
}


#' Internal function that removes the reference tip from a tree
#'
#' @param tree object of class phylo.
#' @param out character vector with the tip name.
#'
#' @return modified tree
#'
remove_root_tip = function(tree, out) {
  data = tree %>% tidytree::as_tibble()

  wh_o = which(data$label == out)
  wh_ep = which(data$parent == data$parent[wh_o] & data$node != data$node[wh_o] & data$parent != data$node)

  data$branch.length[wh_ep] = data$branch.length[wh_ep] + data$branch.length[wh_o]

  data = data[-wh_o, ]

  data$parent = match(data$parent, data$node)
  data$node = match(data$node, data$node)

  tidytree::as.phylo(data)
}


#' Scale edge length of added LP samples for plotting
#'
#' @param x A phylo object.
#' @param frac_height The edge length relative to the tree height.
#'
#' @return A phylo object.
#' @export
#'
set_lp_tiplength = function(x, frac_height=0.01) {

  checkmate::assertClass(x, "phylo")
  checkmate::assertNumber(frac_height, lower = 0, upper = 1)

  wh_added = grep("Added", x$tip.label)
  idx_added = x$edge[, 2] %in% wh_added
  x$edge.length[idx_added] = max(ape::node.depth.edgelength(x)) * frac_height

  return(x)
}



#' Plot showing tree and CNA data site-by-site.
#'
#' @param x Object returned from 'add_lowpass_sampled' or a object of type 'phylo'.
#' @param cn_data CNA data as a named list of genomic ranges objects with a 'cn' metadata column.
#' @param title (optional) title for the plot.
#' @param dropped_samples (optional) vector of dropped samples.
#' @param plot_tree_function (optional) function to plot the tree object.
#' @param labeller_function (optional) function returning a modified version of sample labels.
#' @param ... (optional) arguments passed to the 'plot_tree_function' function.
#' @return A grob.
#' @export
#' @import ggplot2
#'
plot_tree_and_cn_data =
  function(
    x,
    cn_data,
    title = "",
    dropped_samples = NULL,
    plot_tree_function = NULL,
    labeller_function = function(x) x,
    ...
  ) {

  if (try("tree" %in% names(x))) x = x$tree
  if (is.null(plot_tree_function)) plot_tree_function = MLLPT::plot_tree

  checkmate::assertClass(x, "phylo", null.ok = FALSE)
  checkmate::assertList(cn_data, names = "named")
  checkmate::assertString(title)
  checkmate::assertCharacter(dropped_samples, null.ok = TRUE)
  checkmate::assert_function(plot_tree_function, null.ok = TRUE)
  checkmate::assert_function(labeller_function)


  cols_cn = c(
    "0" = "#142132",
    "1" = "#256D7C",
    "2" = "#9CD098",
    "3" = "#F36749",
    "4" = "#CB290C"
  )

  alt_ids = c(
    "chromosome" = "seqnames",
    "start.pos" = "start",
    "end.pos" = "end",
    "CNt" = "cn",
    "CN" = "cn"
  )


  # tree data:
  if (!is.null(dropped_samples))
    x = ape::drop.tip(x, dropped_samples)

  if (is.list(cn_data))
    cn_data =
      lapply(cn_data, as.data.frame) %>%
        reshape2::melt(measure.vars = c()) %>%
        dplyr::mutate(sample = L1)

  for (i in seq_along(alt_ids))
    if (!alt_ids[i] %in% colnames(cn_data))
      if (names(alt_ids)[i] %in% colnames(cn_data))
        cn_data[, alt_ids[i]] = cn_data[, names(alt_ids)[i]]


  # plot of tree
  tree_plot =
    plot_tree_function(x) +
    ggtitle("LP tree") +
    theme(title = element_text(face = "bold"))

  # cn plot
  new_labs =
    magrittr::set_names(
      with(tree_plot$data, label[isTip]),
      x$tip.label[x$tip.label != "GL"]
    )

  lab_order = with(tree_plot$data, label[isTip][order(y[isTip])])
  o_labs = gsub(" [(].*", "", names(new_labs)[match(lab_order, new_labs)])

  cna_plot =
    cn_data %>%
    dplyr::filter(sample %in% o_labs) %>%
    dplyr::filter(!seqnames %in% c("chrX", "chrY")) %>%
    dplyr::mutate(seqnames = factor(seqnames, paste0("chr", 1:22))) %>%
    dplyr::mutate(cn = factor(scales::squish(cn, c(0, 4)))) %>%
    dplyr::mutate(sample = match(sample, o_labs)) %>%
    ggplot(aes(
      x = (start + end) / 2,
      width = end - start,
      y = as.numeric(sample),
      height = 0.75,
      fill = factor(cn)
    )) +
      geom_tile() +
      facet_grid(
        .~seqnames,
        space = "free",
        scales = "free",
        switch = "x"
      ) +
      scale_y_continuous(
        limits = c(0.5, length(o_labs) + 0.5),
        breaks = seq_along(o_labs),
        labels = labeller_function(o_labs)
      ) +
      theme(
        strip.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 0),
        strip.background.x = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 10),
        title = element_text(size = 8, color = "gray20"),
        legend.position = "bottom"
      ) +
      scale_fill_manual(
        values = cols_cn,
        breaks = names(cols_cn),
        drop = FALSE
      ) +
      xlab("") +
      ylab("") +
      labs(fill = "CN") +
      ggtitle("CN data")


  # merge both plots
  merged_plot =
    cowplot::plot_grid(
      ggplot(NULL) +
        theme(axis.line = element_blank()) +
        ggtitle(title),
      cowplot::plot_grid(
        tree_plot +
          ylim(0.5, length(o_labs) + 0.5) + ggtitle(""),
        cna_plot +
          ggtitle(""),
        nrow = 1,
        align = "h",
        axis = "bt",
        labels = c("A", "B"),
        vjust = c(2, 2),
        rel_widths = c(1, 5)
      ),
      hjust = c(1, 1),
      rel_heights = c(0.1, 1),
      ncol = 1
    )

  merged_plot

}
