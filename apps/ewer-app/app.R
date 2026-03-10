library(shiny)
library(ggplot2)
library(dplyr)

# ── Color palette (consistent with MA255 course style) ──
custom_blue   <- "#1A5276"
custom_red    <- "#C0392B"
custom_green  <- "#27AE60"
custom_orange <- "#E67E22"
custom_gray   <- "#95A5A6"
light_blue    <- "#D6EAF8"
light_red     <- "#FADBD8"
custom_purple <- "#8E44AD"
custom_teal   <- "#16A085"

# ── Helper: pairwise adjusted p-values for any method ──
# Returns a vector of adjusted p-values for all pairwise comparisons
get_pairwise_pvals <- function(dat, method) {
  if (method %in% c("none", "bonferroni", "holm", "BH")) {
    pw <- pairwise.t.test(dat$y, dat$group, p.adjust.method = method)
    return(pw$p.value[!is.na(pw$p.value)])
  } else if (method == "tukey") {
    mod <- aov(y ~ group, data = dat)
    tukey_res <- TukeyHSD(mod)$group
    return(tukey_res[, "p adj"])
  } else if (method == "scheffe") {
    g <- nlevels(dat$group)
    N <- nrow(dat)
    df_err <- N - g
    # Get unadjusted p-values from pooled-variance t-tests
    pw <- pairwise.t.test(dat$y, dat$group, p.adjust.method = "none", pool.sd = TRUE)
    raw_p <- pw$p.value[!is.na(pw$p.value)]
    # Recover |t| from two-sided p, then compute Scheffe-adjusted p
    t_vals <- qt(pmax(raw_p, 1e-15) / 2, df_err, lower.tail = FALSE)
    scheffe_p <- pf(t_vals^2 / (g - 1), g - 1, df_err, lower.tail = FALSE)
    return(scheffe_p)
  }
}

# ═══════════════════════════════════════════════════════════
# UI
# ═══════════════════════════════════════════════════════════
ui <- fluidPage(

  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; background-color: #f8f9fa; }
    .navbar { background-color: #1A5276 !important; border: none; }
    .navbar .navbar-brand { color: white !important; font-weight: bold; font-size: 18px; }
    .navbar .navbar-nav > li > a { color: #D6EAF8 !important; font-size: 14px; }
    .navbar .navbar-nav > .active > a { color: white !important; background-color: #2980B9 !important; }
    .well { background-color: white; border: 1px solid #ddd; border-radius: 8px; }
    h4 { color: #1A5276; font-weight: bold; margin-bottom: 15px; }
    .btn-primary { background-color: #1A5276; border-color: #1A5276; }
    .btn-primary:hover { background-color: #2980B9; border-color: #2980B9; }
    .info-box { background-color: #D6EAF8; border-left: 4px solid #1A5276;
                padding: 12px 15px; margin: 10px 0; border-radius: 4px; font-size: 14px; }
    .warning-box { background-color: #FADBD8; border-left: 4px solid #C0392B;
                   padding: 12px 15px; margin: 10px 0; border-radius: 4px; font-size: 14px; }
    .stat-card { background: white; border-radius: 8px; padding: 15px; text-align: center;
                 box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 10px; }
    .stat-value { font-size: 28px; font-weight: bold; }
    .stat-label { font-size: 12px; color: #666; margin-top: 5px; }
    .formula-box { background-color: #fafafa; border: 1px solid #ddd; border-radius: 6px;
                   padding: 15px; text-align: center; margin: 10px 0; font-size: 16px; }
  ")))
  ,

  navbarPage(
    title = "MA255: Experiment-Wise Error Rate",
    id = "main_nav",

    # ── Tab 1: The EWER Formula ──
    tabPanel("1. The Formula",
      fluidRow(
        column(4,
          wellPanel(
            h4("Settings"),
            sliderInput("formula_alpha",
                        HTML("Individual significance level (&alpha;)"),
                        min = 0.01, max = 0.20, value = 0.05, step = 0.01),
            sliderInput("formula_groups",
                        "Number of groups (g)",
                        min = 2, max = 15, value = 5, step = 1),
            div(class = "info-box",
              HTML("<b>Number of pairwise comparisons:</b>"),
              br(),
              HTML("<i>k</i> = <i>g</i>(<i>g</i>&minus;1)/2 = "),
              textOutput("n_comparisons_text", inline = TRUE)
            ),
            div(class = "formula-box",
              HTML("EWER = 1 &minus; (1 &minus; &alpha;)<sup><i>k</i></sup>")
            ),
            div(class = "warning-box",
              HTML("<b>Key insight:</b> Even with &alpha; = 0.05 for each test,
                   the probability of at least one false rejection grows rapidly
                   with the number of comparisons.")
            )
          )
        ),
        column(8,
          plotOutput("ewer_curve", height = "420px"),
          br(),
          fluidRow(
            column(4,
              div(class = "stat-card",
                div(class = "stat-value", style = paste0("color:", custom_blue),
                    textOutput("ewer_value")),
                div(class = "stat-label", "Experiment-Wise Error Rate")
              )
            ),
            column(4,
              div(class = "stat-card",
                div(class = "stat-value", style = paste0("color:", custom_red),
                    textOutput("n_comparisons_val")),
                div(class = "stat-label", "Pairwise Comparisons (k)")
              )
            ),
            column(4,
              div(class = "stat-card",
                div(class = "stat-value", style = paste0("color:", custom_orange),
                    textOutput("expected_false")),
                div(class = "stat-label", "Expected False Rejections")
              )
            )
          )
        )
      )
    ),

    # ── Tab 2: Simulation ──
    tabPanel("2. Simulate It!",
      fluidRow(
        column(4,
          wellPanel(
            h4("Simulation Settings"),
            sliderInput("sim_groups", "Number of groups (g)",
                        min = 2, max = 10, value = 5, step = 1),
            sliderInput("sim_n", "Observations per group (n)",
                        min = 5, max = 50, value = 15, step = 5),
            sliderInput("sim_alpha",
                        HTML("Significance level (&alpha;)"),
                        min = 0.01, max = 0.20, value = 0.05, step = 0.01),
            numericInput("sim_reps", "Number of simulated experiments",
                         value = 100, min = 100, max = 1000, step = 100),
            actionButton("run_sim", "Run Simulation",
                         class = "btn-primary btn-block",
                         icon = icon("play")),
            br(),
            div(class = "info-box",
              HTML("<b>How it works:</b> Each simulated experiment generates
                   data from identical populations (all group means are equal, so
                   every null hypothesis is true). We then run all pairwise
                   <i>t</i>-tests and check if <i>any</i> produce a false rejection.")
            )
          )
        ),
        column(8,
          fluidRow(
            column(6,
              div(class = "stat-card",
                div(class = "stat-value", style = paste0("color:", custom_red),
                    textOutput("sim_ewer_obs")),
                div(class = "stat-label", "Observed EWER (from simulation)")
              )
            ),
            column(6,
              div(class = "stat-card",
                div(class = "stat-value", style = paste0("color:", custom_blue),
                    textOutput("sim_ewer_theory")),
                div(class = "stat-label", "Theoretical EWER (from formula)")
              )
            )
          ),
          plotOutput("sim_barplot", height = "300px"),
          plotOutput("sim_histogram", height = "250px")
        )
      )
    ),

    # ── Tab 3: Adjustments ──
    tabPanel("3. Adjustments Help",
      fluidRow(
        column(4,
          wellPanel(
            h4("Adjustment Comparison"),
            sliderInput("adj_groups", "Number of groups (g)",
                        min = 3, max = 10, value = 5, step = 1),
            sliderInput("adj_n", "Observations per group (n)",
                        min = 5, max = 50, value = 15, step = 5),
            sliderInput("adj_alpha",
                        HTML("Significance level (&alpha;)"),
                        min = 0.01, max = 0.20, value = 0.05, step = 0.01),
            numericInput("adj_reps", "Number of simulated experiments",
                         value = 100, min = 100, max = 1000, step = 100),
            actionButton("run_adj", "Run Simulation",
                         class = "btn-primary btn-block",
                         icon = icon("play")),
            br(),
            div(class = "info-box",
              HTML("<b>Methods compared:</b>
                   <ul style='margin-bottom:0; padding-left:20px;'>
                   <li><b>None</b>: Raw <i>p</i>-values (no correction)</li>
                   <li><b>Bonferroni</b>: Multiply each <i>p</i> by <i>k</i></li>
                   <li><b>Holm</b>: Sequential step-down (less conservative)</li>
                   <li><b>BH (FDR)</b>: Controls false discovery rate</li>
                   <li><b>Tukey HSD</b>: Based on studentized range distribution</li>
                   <li><b>Scheff&eacute;</b>: Most conservative; valid for any contrast</li>
                   </ul>")
            )
          )
        ),
        column(8,
          plotOutput("adj_barplot", height = "350px"),
        )
      )
    ),

    # ── Tab 4: Interactive Example ──
    tabPanel("4. Try It Yourself",
      fluidRow(
        column(4,
          wellPanel(
            h4("Your Experiment"),
            HTML("<p style='font-size:13px;'>Imagine you ran an experiment with several groups and found
                  a significant overall <i>F</i>-test. Now you want to check which
                  pairs are different. Explore what happens!</p>"),
            sliderInput("try_groups", "Number of groups (g)",
                        min = 2, max = 12, value = 4, step = 1),
            radioButtons("try_method", "P-value adjustment method:",
                         choices = c("None" = "none",
                                     "Bonferroni" = "bonferroni",
                                     "Holm" = "holm",
                                     "BH (FDR)" = "BH",
                                     "Tukey HSD" = "tukey",
                                     "Scheffé" = "scheffe"),
                         selected = "none"),
            actionButton("try_run", "Generate New Experiment",
                         class = "btn-primary btn-block",
                         icon = icon("flask")),
            br(),
            div(class = "info-box",
              HTML("<b>Note:</b> Data are generated under H<sub>0</sub> &mdash;
                   all groups have the <i>same</i> true mean. Any 'significant'
                   result is a <b>Type I error</b>!")
            )
          )
        ),
        column(8,
          plotOutput("try_boxplot", height = "250px"),
          plotOutput("try_pvalues", height = "320px"),
          fluidRow(
            column(6,
              div(class = "stat-card",
                div(class = "stat-value", style = paste0("color:", custom_red),
                    textOutput("try_n_sig")),
                div(class = "stat-label", "False Rejections (Type I Errors)")
              )
            ),
            column(6,
              div(class = "stat-card",
                div(class = "stat-value", style = paste0("color:", custom_blue),
                    textOutput("try_n_total")),
                div(class = "stat-label", "Total Comparisons Made")
              )
            )
          )
        )
      )
    )
  )
)

# ═══════════════════════════════════════════════════════════
# SERVER
# ═══════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # ── Tab 1: The Formula ──

  output$n_comparisons_text <- renderText({
    k <- input$formula_groups * (input$formula_groups - 1) / 2
    as.character(k)
  })

  output$ewer_value <- renderText({
    k <- input$formula_groups * (input$formula_groups - 1) / 2
    ewer <- 1 - (1 - input$formula_alpha)^k
    sprintf("%.1f%%", ewer * 100)
  })

  output$n_comparisons_val <- renderText({
    k <- input$formula_groups * (input$formula_groups - 1) / 2
    as.character(k)
  })

  output$expected_false <- renderText({
    k <- input$formula_groups * (input$formula_groups - 1) / 2
    sprintf("%.1f", k * input$formula_alpha)
  })

  output$ewer_curve <- renderPlot({
    g_current <- input$formula_groups
    alpha <- input$formula_alpha
    k_current <- g_current * (g_current - 1) / 2
    ewer_current <- 1 - (1 - alpha)^k_current

    # Build curve over range of groups
    g_seq <- 2:20
    k_seq <- g_seq * (g_seq - 1) / 2
    ewer_seq <- 1 - (1 - alpha)^k_seq

    df_curve <- data.frame(groups = g_seq, k = k_seq, ewer = ewer_seq)
    df_point <- data.frame(groups = g_current, ewer = ewer_current)

    ggplot(df_curve, aes(x = groups, y = ewer)) +
      geom_hline(yintercept = alpha, linetype = "dashed", color = custom_green, linewidth = 0.8) +
      annotate("text", x = 19, y = alpha + 0.03,
               label = paste0("Individual α = ", alpha),
               color = custom_green, size = 4, hjust = 1) +
      geom_line(linewidth = 1.2, color = custom_blue) +
      geom_point(data = df_point, aes(x = groups, y = ewer),
                 color = custom_red, size = 5) +
      geom_segment(data = df_point,
                   aes(x = groups, xend = groups, y = 0, yend = ewer),
                   color = custom_red, linetype = "dotted", linewidth = 0.7) +
      geom_segment(data = df_point,
                   aes(x = 2, xend = groups, y = ewer, yend = ewer),
                   color = custom_red, linetype = "dotted", linewidth = 0.7) +
      annotate("text", x = g_current + 0.5, y = ewer_current,
               label = sprintf("%.1f%%", ewer_current * 100),
               color = custom_red, fontface = "bold", size = 5, hjust = 0) +
      scale_x_continuous(breaks = seq(2, 20, 2)) +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      labs(title = "Experiment-Wise Error Rate vs. Number of Groups",
           subtitle = bquote("EWER = 1 - (1 -" ~ alpha ~ ")" ^ {k} ~
                             ", where k = g(g-1)/2 pairwise comparisons"),
           x = "Number of Groups (g)",
           y = "Experiment-Wise Error Rate") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold", color = custom_blue),
            plot.subtitle = element_text(color = "#555"),
            panel.grid.minor = element_blank())
  })

  # ── Tab 2: Simulation ──

  sim_results <- eventReactive(input$run_sim, {
    g <- input$sim_groups
    n <- input$sim_n
    alpha <- input$sim_alpha
    n_reps <- input$sim_reps

    withProgress(message = "Running simulation...", value = 0, {
      false_rejections <- numeric(n_reps)
      any_rejection <- logical(n_reps)

      for (i in 1:n_reps) {
        # Generate data under H0: all means equal
        dat <- data.frame(
          y = rnorm(g * n, mean = 0, sd = 1),
          group = factor(rep(1:g, each = n))
        )
        # All pairwise t-tests (no adjustment)
        pw <- pairwise.t.test(dat$y, dat$group, p.adjust.method = "none")
        p_vals <- pw$p.value[!is.na(pw$p.value)]

        false_rejections[i] <- sum(p_vals < alpha)
        any_rejection[i] <- any(p_vals < alpha)

        if (i %% 100 == 0) incProgress(100 / n_reps)
      }

      list(false_rejections = false_rejections,
           any_rejection = any_rejection,
           obs_ewer = mean(any_rejection),
           theory_ewer = 1 - (1 - alpha)^(g*(g-1)/2))
    })
  })

  output$sim_ewer_obs <- renderText({
    req(sim_results())
    sprintf("%.1f%%", sim_results()$obs_ewer * 100)
  })

  output$sim_ewer_theory <- renderText({
    req(sim_results())
    sprintf("%.1f%%", sim_results()$theory_ewer * 100)
  })

  output$sim_barplot <- renderPlot({
    req(sim_results())
    res <- sim_results()

    df <- data.frame(
      Category = c("At Least One\nFalse Rejection", "No False\nRejections"),
      Proportion = c(res$obs_ewer, 1 - res$obs_ewer)
    )
    df$Category <- factor(df$Category, levels = df$Category)

    ggplot(df, aes(x = Category, y = Proportion, fill = Category)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.1f%%", Proportion * 100)),
                vjust = -0.5, fontface = "bold", size = 5) +
      geom_hline(yintercept = input$sim_alpha, linetype = "dashed",
                 color = custom_green, linewidth = 0.8) +
      annotate("text", x = 2.4, y = input$sim_alpha + 0.02,
               label = paste0("α = ", input$sim_alpha), color = custom_green, size = 4) +
      scale_fill_manual(values = c(custom_red, custom_blue)) +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      labs(title = "Simulation Results: How Often Did We Get a False Positive?",
           subtitle = sprintf("%d experiments simulated under H₀ (all means equal)",
                              input$sim_reps),
           y = "Proportion of Experiments") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", color = custom_blue),
            axis.title.x = element_blank(),
            panel.grid.minor = element_blank())
  })

  output$sim_histogram <- renderPlot({
    req(sim_results())
    res <- sim_results()

    df <- data.frame(n_false = res$false_rejections)

    ggplot(df, aes(x = n_false)) +
      geom_histogram(binwidth = 1, fill = custom_orange, color = "white", alpha = 0.85) +
      scale_x_continuous(breaks = 0:max(c(5, max(df$n_false)))) +
      labs(title = "Distribution of False Rejections per Experiment",
           x = "Number of False Rejections (out of k pairwise tests)",
           y = "Count") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", color = custom_blue),
            panel.grid.minor = element_blank())
  })

  # ── Tab 3: Adjustments ──

  adj_results <- eventReactive(input$run_adj, {
    g <- input$adj_groups
    n <- input$adj_n
    alpha <- input$adj_alpha
    n_reps <- input$adj_reps

    methods <- c("none", "bonferroni", "holm", "BH", "tukey", "scheffe")
    method_labels <- c("None", "Bonferroni", "Holm", "BH (FDR)", "Tukey HSD", paste0("Scheff", "\u00e9"))

    n_methods <- length(methods)

    withProgress(message = "Simulating adjustments...", value = 0, {
      any_reject_mat <- matrix(FALSE, nrow = n_reps, ncol = n_methods)

      for (i in 1:n_reps) {
        dat <- data.frame(
          y = rnorm(g * n, mean = 0, sd = 1),
          group = factor(rep(1:g, each = n))
        )

        for (j in seq_along(methods)) {
          p_vals <- get_pairwise_pvals(dat, methods[j])
          any_reject_mat[i, j] <- any(p_vals < alpha)
        }

        if (i %% 50 == 0) incProgress(50 / n_reps)
      }

      ewer_by_method <- colMeans(any_reject_mat)

      data.frame(
        method = factor(method_labels, levels = method_labels),
        ewer = ewer_by_method
      )
    })
  })

  output$adj_barplot <- renderPlot({
    req(adj_results())
    df <- adj_results()
    alpha <- input$adj_alpha

    bar_colors <- c(custom_red, custom_blue, custom_orange, custom_green, custom_purple, custom_teal)

    ggplot(df, aes(x = method, y = ewer, fill = method)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.1f%%", ewer * 100)),
                vjust = -0.5, fontface = "bold", size = 4.5) +
      geom_hline(yintercept = alpha, linetype = "dashed",
                 color = "black", linewidth = 0.8) +
      annotate("text", x = 6.4, y = alpha + 0.015,
               label = paste0("Target α = ", alpha), size = 4) +
      scale_fill_manual(values = bar_colors) +
      scale_y_continuous(labels = scales::percent_format(),
                         limits = c(0, max(df$ewer) * 1.2 + 0.05)) +
      labs(title = "Observed EWER by Adjustment Method",
           subtitle = sprintf("%d groups, %d obs/group, %d simulations under H₀",
                              input$adj_groups, input$adj_n, input$adj_reps),
           y = "Experiment-Wise Error Rate") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", color = custom_blue),
            axis.title.x = element_blank(),
            panel.grid.minor = element_blank())
  })

  # ── Tab 4: Try It Yourself ──

  try_data <- eventReactive(input$try_run, {
    g <- input$try_groups
    n <- 15

    dat <- data.frame(
      y = rnorm(g * n, mean = 5, sd = 1.5),
      group = factor(paste("Group", LETTERS[1:g]),
                     levels = paste("Group", LETTERS[1:g]))
    )

    method <- input$try_method

    if (method %in% c("tukey", "scheffe")) {
      # For Tukey and Scheffe, use helper to get p-values
      # and extract pair labels from TukeyHSD (consistent naming)
      mod <- aov(y ~ group, data = dat)
      tukey_res <- TukeyHSD(mod)$group
      pair_names <- rownames(tukey_res)
      # Replace "-" with " vs. " for display
      pair_names <- gsub("-", " vs. ", pair_names)

      p_vals <- get_pairwise_pvals(dat, method)

      pairs_df <- data.frame(
        pair = pair_names,
        p_value = as.numeric(p_vals),
        significant = as.numeric(p_vals) < 0.05
      )
    } else {
      pw <- pairwise.t.test(dat$y, dat$group, p.adjust.method = method)

      # Extract p-values into tidy format
      p_mat <- pw$p.value
      pairs_list <- list()
      for (i in 1:nrow(p_mat)) {
        for (j in 1:ncol(p_mat)) {
          if (!is.na(p_mat[i, j])) {
            pairs_list[[length(pairs_list) + 1]] <- data.frame(
              pair = paste(colnames(p_mat)[j], "vs.", rownames(p_mat)[i]),
              p_value = p_mat[i, j],
              significant = p_mat[i, j] < 0.05
            )
          }
        }
      }
      pairs_df <- bind_rows(pairs_list)
    }

    pairs_df$pair <- factor(pairs_df$pair, levels = rev(pairs_df$pair))

    list(dat = dat, pairs = pairs_df, method = method)
  })

  output$try_boxplot <- renderPlot({
    req(try_data())
    dat <- try_data()$dat

    ggplot(dat, aes(x = group, y = y, fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = 21) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
      scale_fill_brewer(palette = "Set3") +
      labs(title = "Simulated Data (All Groups Have the Same True Mean!)",
           x = "", y = "Response") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", color = custom_blue, size = 13),
            axis.text.x = element_text(angle = 30, hjust = 1))
  })

  output$try_pvalues <- renderPlot({
    req(try_data())
    pairs <- try_data()$pairs
    method_label <- switch(try_data()$method,
                           "none" = "No Adjustment",
                           "bonferroni" = "Bonferroni",
                           "holm" = "Holm",
                           "BH" = "BH (FDR)",
                           "tukey" = "Tukey HSD",
                           "scheffe" = paste0("Scheff", "\u00e9"))

    ggplot(pairs, aes(x = p_value, y = pair, fill = significant)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.8) +
      annotate("text", x = 0.07, y = nrow(pairs) * 0.95,
               label = "α = 0.05", size = 3.5, hjust = 0) +
      scale_fill_manual(values = c("TRUE" = custom_red, "FALSE" = custom_gray),
                        labels = c("TRUE" = "Significant (Type I Error!)",
                                   "FALSE" = "Not Significant"),
                        name = "") +
      scale_x_continuous(limits = c(0, 1.05)) +
      labs(title = paste0("Pairwise P-Values (", method_label, ")"),
           x = "P-Value", y = "") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", color = custom_blue),
            legend.position = "bottom",
            panel.grid.minor = element_blank())
  })

  output$try_n_sig <- renderText({
    req(try_data())
    as.character(sum(try_data()$pairs$significant))
  })

  output$try_n_total <- renderText({
    req(try_data())
    as.character(nrow(try_data()$pairs))
  })
}

shinyApp(ui, server)

