colors_overload <- union(ggthemes::tableau_color_pal('Tableau 20')(20), RColorBrewer::brewer.pal(12, 'Set3'))
colors_overload <- c(colors_overload, 'black')

fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
}



do_scatter <- function (umap_use, meta_data, label_name, facet_var, no_guides = TRUE, 
    do_labels = TRUE, nice_names, palette_use = colors_overload, 
    pt_size = 4, point_size = 0.5, pt_shape = ".", base_size = 20, 
    do_points = TRUE, do_density = FALSE, h = 3, w = 4, 
                        alpha_fore=1, alpha_back=.3, color_back='lightgrey', 
                       nrow = 1, do_raster = FALSE) 
{
    if (do_raster) {
        geom_point_fxn <- function(...) geom_point_rast(..., width = w, height = h)
    } else {
        geom_point_fxn <- geom_point
    }
    
    plt_df <- data.frame(umap_use)[, 1:2]
    colnames(plt_df) <- c('X1', 'X2')
    plt_df <- plt_df %>% 
        cbind(meta_data) %>% 
        dplyr::sample_frac(1L)
    plt_df$given_name <- plt_df[[label_name]]
    if (!missing(nice_names)) {
        plt_df %<>% dplyr::inner_join(nice_names, by = "given_name") %>% 
            subset(nice_name != "" & !is.na(nice_name))
        plt_df[[label_name]] <- plt_df$nice_name
    }
    
    plt <- plt_df %>% ggplot(aes_string("X1", "X2", col = label_name, 
        fill = label_name)) + 
#         theme_tufte(base_size = base_size) + 
#         theme(panel.background = element_rect(fill = NA, color = "black")) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, 
            alpha = 1, shape = 16, size = 4)), alpha = FALSE) + 
        scale_color_manual(values = palette_use) + scale_fill_manual(values = palette_use) + 
        theme(plot.title = element_text(hjust = 0.5)) + labs(x = "UMAP 1", 
        y = "UMAP 2")
    if (do_points) {
        ## this facets while keeping non-facet points in the background
        if (!missing(facet_var)) {
            if (!is(facet_var, 'quosure')) {
                stop('facet_var must be a quosure. e.g. quo(\'donor\')')
            }            

            plt <- plt + geom_point_fxn(
                data = dplyr::select(plt_df, -!!facet_var), 
                shape = pt_shape, size = point_size,
                color = color_back, fill = color_back, alpha = alpha_back
            ) +
                facet_wrap(vars(!!facet_var), nrow = nrow)
        }
        plt <- plt + geom_point_fxn(shape = pt_shape, size = point_size, alpha = alpha_fore)
    }
    if (do_density) 
        plt <- plt + geom_density_2d()
    if (no_guides) 
        plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)
    if (do_labels) {
        plt <- plt + 
#             geom_text_repel(
#                 data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
#                 label.size = NA, aes_string(label = label_name), 
#                 color = "black", 
#                 size = pt_size, alpha = 1, segment.size = 0
#             ) + 
            geom_label(
                data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
                label.size = NA, aes_string(label = label_name, color = label_name), 
#                 color = "black", 
                fill = 'white', 
                size = pt_size, alpha = .6, segment.size = 0
            ) + 
            geom_text(
                data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
                label.size = NA, aes_string(label = label_name, color = label_name), 
#                 color = "black", 
                size = pt_size, alpha = 1, segment.size = 0
            ) + 
            guides(col = FALSE, fill = FALSE)
    }
    return(plt)
}


setupVals <- function(data_mat, feature, qlo, qhi) {
    .x <- data_mat[feature, , drop = FALSE] %>% as("dgTMatrix")
    cutoffs <- quantileSparse(.x, c(qlo, qhi))
    cutoffs[2] <- max(cutoffs[2], min(.x@x))
    if (qlo == 0 & qhi == 1) {
        return(.x)
    } 
    
    if (qlo > 0) {
        .x@x[.x@x < cutoffs[1]] <- cutoffs[1]
#         message(sprintf("For %s, lo = %.3f", feature, ifelse(length(.x@x) == ncol(.x), cutoffs[1], NA)))
    }
    if (qhi < 1) {
        .x@x[.x@x > cutoffs[2]] <- cutoffs[2]
#         message(sprintf("For %s, hi = %.3f", feature, cutoffs[2]))
        
    }
    return(.x)
}


quantileSparse <- function(.x, qlist) {
    ratio_zero <- 1 - (length(.x@x) / ncol(.x))
    q_nz <- which(qlist > ratio_zero)
    q_adj <- (qlist[q_nz] - ratio_zero) / (1 - ratio_zero)
    res <- rep(0, length(qlist))
    res[q_nz] <- quantile(.x@x, q_adj)
    res
}

## TODO: test is feature is present
## TODO: allow for different cutoffs, for each marker
## TODO: somehow draw canvas first, then do plotting? 
library(patchwork)
library(ggthemes)

plotFeatures <- function(data_mat, dim_df, features, nrow = 1, 
                         qlo = 0.05, qhi = 1, order_by_expression = FALSE, 
                         pt_shape = 16, pt_size = .5, no_guide = FALSE,
                         .xlim = c(NA, NA), .ylim = c(NA, NA), color_high = muted("blue")) {
    plt_df <- data.frame(dim_df[, 1:2])
    colnames(plt_df) <- c("X1", "X2")


    plt_list <- lapply(features, function(feature) {
        .x <- setupVals(data_mat, feature, qlo, qhi)
        plt_df$value <- 0
        plt_df[.x@j + 1, "value"] <- .x@x
        if (order_by_expression) {
            plt_df %<>% dplyr::arrange(value)             
        } else {
            plt_df %<>% dplyr::sample_frac(1L)
        }

        plt <- plt_df %>% 
            ggplot(aes(X1, X2, color = value)) + 
#             geom_point_rast(dpi = 300, width = 6, height = 4, size = .5, shape = pt_shape) + 
            geom_point(shape = ".") + 
            scale_color_gradient2(na.value = "lightgrey", mid = "lightgrey", midpoint = 0, high = color_high) + 
#             theme_tufte(base_size = 14, base_family = "Helvetica") + 
#             theme(panel.background = element_rect(), plot.title = element_text(hjust = .5)) +
            theme(plot.title = element_text(hjust = .5)) +
            labs(x = "UMAP 1", y = "UMAP 2", title = feature) + 
            NULL
        if (no_guide) {
            plt <- plt + 
            guides(color = FALSE) 
        }
        
        if (sum(is.na(.xlim)) < 2) 
            plt <- plt + xlim(.xlim)
        if (sum(is.na(.ylim)) < 2) 
            plt <- plt + ylim(.ylim)
        plt

    })
    if (length(plt_list) > 1) {
        Reduce(`+`, plt_list) + patchwork::plot_layout(nrow = nrow)
    } else {
        plt_list[[1]]
    }
}


cbind2_fill <- function(A, B) {
    rownames_all <- union(rownames(A), rownames(B))
    
    add_to_A <- setdiff(rownames_all, rownames(A))
    add_to_B <- setdiff(rownames_all, rownames(B))

    A2 <- Matrix::rsparsematrix(length(add_to_A), ncol(A), 0)
    B2 <- Matrix::rsparsematrix(length(add_to_B), ncol(B), 0)

    A3 <- Matrix::rbind2(A, A2)
    rownames(A3) <- c(rownames(A), add_to_A)
    B3 <- Matrix::rbind2(B, B2)
    rownames(B3) <- c(rownames(B), add_to_B)

    res <- Matrix::cbind2(A3[rownames_all, ], B3[rownames_all, ])
    return(res)    
}





get_markers <- function(counts, meta_data, cluster_colname) {
    meta_data$cluster <- meta_data[[cluster_colname]]
    pb <- presto::collapse_counts(
        counts, meta_data, 
        c('library', 'cluster'), 
        min_cells_per_group = 10, 
        get_norm = TRUE
    )  
    # pb$meta_data$cluster <- paste0('C', pb$meta_data$cluster)

    
    system.time({
        suppressWarnings({
            presto_res <- presto.presto(
                y ~ 1 + (1|cluster) + (1|cluster:library) + (1|library) + offset(logUMI), 
                pb$meta_data, 
                pb$counts_mat,
                size_varname = 'logUMI', 
                effects_cov = c('cluster'),
                ncore = 20, 
                min_sigma = .05, 
                family = 'poisson',
                nsim = 1000
            ) 
        })
    })

    contrasts_mat <- make_contrast.presto(presto_res, 'cluster')
    effects_marginal <- contrasts.presto(presto_res, contrasts_mat, one_tailed = TRUE) %>% 
        dplyr::mutate(cluster = contrast) %>% 
        dplyr::mutate(
            ## convert stats to log2 for interpretability 
            logFC = sign(beta) * log2(exp(abs(beta))),
            SD = log2(exp(sigma)),
            zscore = logFC / SD
        ) %>% 
        dplyr::select(cluster, feature, logFC, SD, zscore, pvalue) %>% 
        arrange(pvalue)
    effects_marginal$fdr <- p.adjust(effects_marginal$pvalue, method = 'BH')
    return(effects_marginal)
}


tt <- function(effects_marginal, fdr_max = .2, logFC_min = .2) {
    # cluster_max = max(as.integer(gsub('.*?(\\d+)', '\\1', unique(effects_marginal$cluster))))    
    effects_marginal %>% 
        subset(fdr < fdr_max & logFC > logFC_min) %>% 
        # dplyr::mutate(cluster = factor(cluster, paste0('C', 0:cluster_max))) %>% 
        split(.$cluster) %>% 
        map(function(.SD) {
            .SD %>% 
                # dplyr::arrange(-zscore) %>% 
                dplyr::arrange(-logFC) %>% 
                head(10) %>% 
                dplyr::select(feature, cluster) %>% 
                tibble::rowid_to_column('rank')
        }) %>% 
        bind_rows() %>% 
        tidyr::spread(cluster, feature) %>% 
        dplyr::select(-rank) %>% 
        t() %>% 
        identity()
}
