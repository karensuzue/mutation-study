# This script plots fitness over generations for the best and worst const mutation runs

library(fs)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

root_const <- "/mnt/d/MINE/sse-staticenv-constmut/"
best_U <- "3.1623e\\-04"
worst_U <- "1.0000e\\+02"

# Function to parse U (genome-wide mutation rate), rep, kind from filename
# Expects: history_<startU>_<rep>_(fitness|mutation).csv
# startU refers to the starting rate in runs with evolving mutation
get_meta <- function(file) {
    nm <- basename(file)
    m <- regexec("^history_([^_]+)_([0-9]+)_(fitness|mutation)\\.csv$", nm)
    parts <- regmatches(nm, m)[[1]]
    if (length(parts) == 0) return(NULL)
    return(list(U=as.numeric(parts[2]), rep=as.integer(parts[3]), kind=parts[4]))
}

read_fit_file <- function(file) {
    meta <- get_meta(file)
    if (is.null(meta)) {
        stop(paste("Filename didn't match expected pattern:", basename(file)))
    }
    df <- read_csv(file, show_col_types=FALSE) %>%
          select(-Fittest_ID) %>% 
          mutate(Start_U = meta$U)
    return(df)
}

# Grab all history files we need
files_fit_best <- dir_ls(
    root_const, 
    recurse=TRUE, 
    regexp=paste0("history_", best_U, "_.*_fitness\\.csv$")
)
files_fit_best

files_fit_worst <- dir_ls(
    root_const, 
    recurse=TRUE, 
    regexp=paste0("history_", worst_U, "_.*_fitness\\.csv$")
)
files_fit_worst

# -------------------------
# Best and worst summary
# -------------------------
all_best_df <- list()
for (f in files_fit_best) {
    all_best_df[[length(all_best_df)+1]] <- read_fit_file(f)
}
all_best_df <- bind_rows(all_best_df)

all_worst_df <- list()
for (f in files_fit_worst) {
    all_worst_df[[length(all_worst_df)+1]] <- read_fit_file(f)
}
all_worst_df <- bind_rows(all_worst_df)

avg_best_fit_df <- all_best_df %>%
              group_by(Start_U, Generation) %>%
              summarize(
                Avg_Mean_F=mean(Mean_F, na.rm=TRUE),
                SD_Mean_F=sd(Mean_F, na.rm=TRUE),
                .groups="drop"
              )

avg_worst_fit_df <- all_worst_df %>%
              group_by(Start_U, Generation) %>%
              summarize(
                Avg_Mean_F=mean(Mean_F, na.rm=TRUE),
                SD_Mean_F=sd(Mean_F, na.rm=TRUE),
                .groups="drop"
              )

avg_fit_df <- bind_rows(avg_best_fit_df, avg_worst_fit_df)

median_best_fit_df <- all_best_df %>%
                 group_by(Start_U, Generation) %>%
                 summarize(
                    Median_F=median(Mean_F, na.rm=TRUE),
                    Q1_F=quantile(Mean_F, 0.25, na.rm=TRUE),
                    Q3_F=quantile(Mean_F, 0.75, na.rm=TRUE),
                    .groups="drop"
                )
median_worst_fit_df <- all_worst_df %>%
                 group_by(Start_U, Generation) %>%
                 summarize(
                    Median_F=median(Mean_F, na.rm=TRUE),
                    Q1_F=quantile(Mean_F, 0.25, na.rm=TRUE),
                    Q3_F=quantile(Mean_F, 0.75, na.rm=TRUE),
                    .groups="drop"
                )

median_fit_df <- bind_rows(median_best_fit_df, median_worst_fit_df)

# -------------------------
# Safe values for log plots
# -------------------------
eps <- 1e-12

# FOR SSE - use error instead of fitness
avg_err_df <- avg_fit_df %>%
            mutate(
                Avg_Error = -Avg_Mean_F,
                SD_Error = SD_Mean_F
            )


avg_err_log <- avg_err_df %>%
               mutate(
                y=pmax(Avg_Error, eps),
                ymin=pmax(Avg_Error - SD_Error, eps),
                ymax=pmax(Avg_Error + SD_Error, eps)
               )

median_err_df <- median_fit_df %>%
    mutate(
        Median_Error=-Median_F,
        Q1_Error=-Q3_F,
        Q3_Error=-Q1_F
    )

median_err_log <- median_err_df %>%
    # mutate(
    #     y = pmax(Median_Error, eps),
    #     ymin = pmax(Q1_Error, eps),
    #     ymax = pmax(Q3_Error, eps)
    # )
    mutate(
        y = pmax(Median_Error, eps),
        ymin = ifelse(Q1_Error > 0, Q1_Error, NA_real_),
        ymax = ifelse(Q1_Error > 0 & Q3_Error > 0, Q3_Error, NA_real_)
    )

# -------------------------
# Linear and log
# -------------------------
p_fit_linear <- ggplot(avg_fit_df, aes(x=Generation, y=Avg_Mean_F, color=factor(Start_U), group=Start_U)) +
                geom_ribbon(
                    aes(ymin=Avg_Mean_F - SD_Mean_F,
                        ymax=Avg_Mean_F + SD_Mean_F,
                        group=Start_U,
                        fill=factor(Start_U)),
                    # inherit.aes=TRUE,
                    alpha=0.15,
                    color=NA
                ) +
                geom_line(linewidth=0.8) +
                labs(
                    x="Generation",
                    y="Mean fitness (across replicates)",
                    color="Start U",
                    fill="Start U"
                ) +
                theme_bw()

p_err_logy <- ggplot(
    avg_err_log, aes(x=Generation, y=y, color=factor(Start_U), group=Start_U)) +
    geom_ribbon(
        aes(ymin=ymin, ymax=ymax, fill=factor(Start_U)),
        alpha=0.15,
        color=NA
    ) +
    geom_line(linewidth=0.8) + 
    scale_y_log10() +
    labs(
        x="Generation",
        y="Mean error (-fitness)",
        color="Start U",
        fill="Start U"
    ) + 
    theme_bw()


p_err_logy_median <- ggplot(
    median_err_log,
    aes(x=Generation, y=y, color=factor(Start_U), group=Start_U)) +
    geom_ribbon(
        aes(
            ymin=ymin,
            ymax=ymax,
            fill=factor(Start_U)
        ),
        alpha=0.15,
        color=NA
    ) +
    geom_line(linewidth=0.8) +
    scale_y_log10() +
    labs(
        x="Generation",
        y="Median error (-fitness)",
        color="Start U",
        fill="Start U"
    ) +
    theme_bw()


ggsave("fig_generation_constmut_linear.pdf", p_fit_linear, width=7, height=5)
ggsave("fig_generation_constmut_logy.pdf", p_err_logy, width=7, height=5)
ggsave("fig_generation_constmut_logy_median.pdf", p_err_logy_median, width=7, height=5)