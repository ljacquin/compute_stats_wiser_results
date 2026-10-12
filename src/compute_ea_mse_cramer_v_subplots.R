# load libraries
library(mvtnorm)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(MASS)
library(data.table)
library(stringr)
library(whitening)
library(lme4)
library(tidyr)
# library(rstudioapi)
# setwd(dirname(getActiveDocumentContext()$path))
setwd("compute_stats_wiser_results/src/")
source("functions.R")

# get results for all species
results_df_rice <- as.data.frame(
    fread('../results/wiser_blue_ea_gain/Rice_sim_results_cramer_v_h2_0.5.csv')
)
results_df_Maize <- as.data.frame(
    fread('../results/wiser_blue_ea_gain/Maize_sim_results_cramer_v_h2_0.5.csv')
)
results_df_Apple <- as.data.frame(
    fread('../results/wiser_blue_ea_gain/Apple_sim_results_cramer_v_h2_0.5.csv')
)
results_df_Pine <- as.data.frame(
    fread('../results/wiser_blue_ea_gain/Pine_sim_results_cramer_v_h2_0.5.csv')
)

# add species names
results_df_rice$Species <- "Rice"
results_df_Maize$Species <- "Maize"
results_df_Apple$Species <- "Apple"
results_df_Pine$Species <- "Pine"

# combine results for all species
results_df <- rbind(
    results_df_rice,
    results_df_Maize,
    results_df_Apple,
    results_df_Pine
)

# define species order
results_df$Species <- factor(
    results_df$Species,
    levels = c(
        "Rice",
        "Maize",
        "Apple",
        "Pine"
    )
)

# EA versus Cramer's V
df_rho <- results_df |>
    select(
        Species,
        Observed_Cramer_V,
        WISER_Correlation,
        BLUE_Correlation
    ) |>
    pivot_longer(
        cols = c(
            WISER_Correlation,
            BLUE_Correlation
        ),
        names_to = "Method",
        values_to = "Correlation"
    ) |>
    mutate(
        Method = recode(
            Method,
            WISER_Correlation = "WISER",
            BLUE_Correlation = "BLUE"
        )
    ) |>
    drop_na()

ggplot_rho <- ggplot(
    df_rho,
    aes(
        x = Observed_Cramer_V,
        y = Correlation,
        color = Method
    )
) +
    geom_point(
        alpha = 0.6,
        size = 2
    ) +
    geom_smooth(
        method = "loess",
        se = TRUE
    ) +
    facet_wrap(
        ~Species,
        ncol = 2,
        scales = "free"
    ) +
    scale_color_brewer(
        palette = "Set2"
    ) +
    theme_minimal(
        base_family = "Calibri"
    ) +
    labs(
        x = "Cramer's V",
        y = "Estimation accuracy (EA)",
        color = NULL
    ) +
    theme(
        plot.background = element_rect(
            fill = "white",
            color = NA
        ),
        panel.background = element_rect(
            fill = "white",
            color = NA
        ),
        plot.title = element_text(
            size = 14,
            face = "bold",
            hjust = 0.5
        ),
        strip.text = element_text(
            size = 14,
            face = "bold"
        ),
        axis.line = element_line(
            colour = "black",
            linewidth = 1
        ),
        axis.ticks = element_line(
            linewidth = 1,
            color = "black"
        ),
        axis.text = element_text(
            color = "black"
        ),
        axis.ticks.length = unit(
            0.2,
            "cm"
        ),
        legend.position = "bottom"
    ) +
    font(
        "xylab",
        size = 15
    ) +
    font(
        "xy.text",
        size = 13
    ) +
    font(
        "legend.text",
        size = 13
    )
ggplot_rho

ggsave(
    filename = "../results/wiser_blue_ea_gain/sim_rho_cramer_v_all_species_h2_0.5.png",
    plot = ggplot_rho,
    width = 12,
    height = 9,
    dpi = 300
)

# MSE versus Cramer's V
df_mse <- results_df |>
    select(
        Species,
        Observed_Cramer_V,
        WISER_MSE,
        BLUE_MSE
    ) |>
    pivot_longer(
        cols = c(
            WISER_MSE,
            BLUE_MSE
        ),
        names_to = "Method",
        values_to = "MSE"
    ) |>
    mutate(
        Method = recode(
            Method,
            WISER_MSE = "WISER",
            BLUE_MSE = "BLUE"
        )
    ) |>
    drop_na()

ggplot_mse <- ggplot(
    df_mse,
    aes(
        x = Observed_Cramer_V,
        y = MSE,
        color = Method
    )
) +
    geom_point(
        alpha = 0.6,
        size = 2
    ) +
    geom_smooth(
        method = "loess",
        se = TRUE
    ) +
    facet_wrap(
        ~Species,
        ncol = 2,
        scales = "free"
    ) +
    scale_color_brewer(
        palette = "Set2"
    ) +
    theme_minimal(
        base_family = "Calibri"
    ) +
    labs(
        x = "Cramer's V",
        y = "Mean squared error (MSE)",
        color = NULL
    ) +
    theme(
        plot.background = element_rect(
            fill = "white",
            color = NA
        ),
        panel.background = element_rect(
            fill = "white",
            color = NA
        ),
        plot.title = element_text(
            size = 14,
            face = "bold",
            hjust = 0.5
        ),
        strip.text = element_text(
            size = 14,
            face = "bold"
        ),
        axis.line = element_line(
            colour = "black",
            linewidth = 1
        ),
        axis.ticks = element_line(
            linewidth = 1,
            color = "black"
        ),
        axis.text = element_text(
            color = "black"
        ),
        axis.ticks.length = unit(
            0.2,
            "cm"
        ),
        legend.position = "bottom"
    ) +
    font(
        "xylab",
        size = 15
    ) +
    font(
        "xy.text",
        size = 13
    ) +
    font(
        "legend.text",
        size = 13
    )
ggplot_mse

ggsave(
    filename = "../results/wiser_blue_ea_gain/sim_mse_cramer_v_all_species_h2_0.5.png",
    plot = ggplot_mse,
    width = 12,
    height = 9,
    dpi = 300
)

# EA gain versus Cramer's V
df_ea_gain <- results_df |>
    select(
        Species,
        Observed_Cramer_V,
        EA_gain
    ) |>
    drop_na()

ggplot_ea_gain <- ggplot(
    df_ea_gain,
    aes(
        x = Observed_Cramer_V,
        y = EA_gain
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "grey40"
    ) +
    geom_point(
        alpha = 0.6,
        size = 2
    ) +
    geom_smooth(
        method = "loess",
        se = TRUE
    ) +
    facet_wrap(
    ~Species,
    ncol = 2,
    scales = "free"
    ) +
    theme_minimal(
        base_family = "Calibri"
    ) +
    labs(
        x = "Cramer's V",
        y = "EA gain (WISER - BLUE)"
    ) +
    theme(
        plot.background = element_rect(
            fill = "white",
            color = NA
        ),
        panel.background = element_rect(
            fill = "white",
            color = NA
        ),
        plot.title = element_text(
            size = 14,
            face = "bold",
            hjust = 0.5
        ),
        strip.text = element_text(
            size = 14,
            face = "bold"
        ),
        axis.line = element_line(
            colour = "black",
            linewidth = 1
        ),
        axis.ticks = element_line(
            linewidth = 1,
            color = "black"
        ),
        axis.text = element_text(
            color = "black"
        ),
        axis.ticks.length = unit(
            0.2,
            "cm"
        )
    ) +
    font(
        "xylab",
        size = 15
    ) +
    font(
        "xy.text",
        size = 13
    )
ggplot_ea_gain

ggsave(
    filename = "../results/wiser_blue_ea_gain/sim_ea_gain_cramer_v_all_species_h2_0.5.png",
    plot = ggplot_ea_gain,
    width = 12,
    height = 9,
    dpi = 300
)