# This script plots evolving mutation rate and fitness over generations

library(fs)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

root_evolve <- "./Data/evolve-diag/"
root_const <- "./Data/const-diag2/"
best_const_U <- 10

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

read_mut_file <- function(file) {
    meta <- get_meta(file)
    if (is.null(meta)) {
        stop(paste("Filename didn't match expected pattern:", basename(file)))
    }
    df <- read_csv(file, show_col_types=FALSE) %>%
          select(-Best_mu, -Mean_mu, -Highest_ID) %>%
          mutate(Start_U = meta$U)
    return(df)
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
files_mut <- dir_ls(root_evolve, recurse=TRUE, regexp="history_.*_mutation\\.csv$")
files_fit <- dir_ls(root_evolve, recurse=TRUE, regexp="history_.*_fitness\\.csv$")
files_best_U_fit <- dir_ls(root_const, recurse=TRUE, regexp="history_1.0000e\\+01_.*_fitness\\.csv$")

# --- DATAFRAMES FROM EVOLVING RUNS ---
all_mut_df <- list()
for (f in files_mut) {
    all_mut_df[[length(all_mut_df)+1]] <- read_mut_file(f)
}
all_mut_df <- bind_rows(all_mut_df)
# all_mut_df

all_fit_df <- list()
for (f in files_fit) {
    all_fit_df[[length(all_fit_df)+1]] <- read_fit_file(f)
}
all_fit_df <- bind_rows(all_fit_df)
# all_fit_df

avg_mut_df <- all_mut_df %>%
              group_by(Start_U, Generation) %>%
              summarize(
                Avg_Mean_U=mean(Mean_U, na.rm=TRUE),
                SD_Mean_U=sd(Mean_U, na.rm=TRUE),
                .groups="drop"
              )
# avg_mut_df

avg_fit_df <- all_fit_df %>%
              group_by(Start_U, Generation) %>%
              summarize(
                Avg_Mean_F=mean(Mean_F, na.rm=TRUE),
                SD_Mean_F=sd(Mean_F, na.rm=TRUE),
                .groups="drop"
              )
# avg_fit_df

median_mut_df <- all_mut_df %>%
                 group_by(Start_U, Generation) %>%
                 summarize(
                    Median_Mean_U = median(Mean_U, na.rm=TRUE),
                    Q1_Mean_U = quantile(Mean_U, 0.25, na.rm=TRUE),
                    Q3_Mean_U = quantile(Mean_U, 0.75, na.rm=TRUE),
                    .groups="drop"
                )

# --- DATAFRAMES FROM CONSTANT RUN WITH BEST U ---
all_best_U_fit_df <- list()
for (f in files_best_U_fit) {
    all_best_U_fit_df[[length(all_best_U_fit_df)+1]] <- read_fit_file(f)
}
all_best_U_fit_df <- bind_rows(all_best_U_fit_df)
# all_best_U_fit_df

avg_best_U_fit_df <- all_best_U_fit_df %>%
                     group_by(Generation) %>%
                     summarize(
                        Avg_Mean_F=mean(Mean_F, na.rm=TRUE),
                        SD_Mean_F=sd(Mean_F, na.rm=TRUE),
                        .groups="drop"
                     )
# avg_best_U_fit_df


# --- SAFE FOR LOG PLOTS ---
eps <- 1e-12
avg_mut_log <- avg_mut_df %>%
               mutate(
                y=pmax(Avg_Mean_U, eps),
                ymin=pmax(Avg_Mean_U - SD_Mean_U, eps),
                ymax=pmax(Avg_Mean_U + SD_Mean_U, eps)
               )

avg_fit_log <- avg_fit_df %>%
               mutate(
                y=pmax(Avg_Mean_F, eps),
                ymin=pmax(Avg_Mean_F - SD_Mean_F, eps),
                ymax=pmax(Avg_Mean_F + SD_Mean_F, eps)
               )

avg_best_U_fit_log <- avg_best_U_fit_df %>%
                      mutate(
                      y=pmax(Avg_Mean_F, eps),
                      ymin=pmax(Avg_Mean_F - SD_Mean_F, eps),
                      ymax=pmax(Avg_Mean_F + SD_Mean_F, eps)
                      )
# --- PLOTS ---
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
                # For best constant U
                geom_ribbon(
                    data=avg_best_U_fit_df,
                    aes(x=Generation,
                        ymin=Avg_Mean_F - SD_Mean_F,
                        ymax=Avg_Mean_F + SD_Mean_F),
                    inherit.aes=FALSE,
                    fill="grey",
                    alpha=0.35
                ) +
                geom_line(
                    data=avg_best_U_fit_df,
                    aes(x=Generation, y=Avg_Mean_F),
                    inherit.aes=FALSE,
                    linewidth=0.8,
                    linetype="dashed",
                    color="black"
                ) +
                labs(
                    x="Generation",
                    y="Mean fitness (across replicates)",
                    color="Start U"
                ) +
                theme_bw()

p_mut_linear <- ggplot(avg_mut_df, aes(x=Generation, y=Avg_Mean_U, color=factor(Start_U), group=Start_U)) +
                geom_ribbon(
                    aes(ymin=Avg_Mean_U - SD_Mean_U,
                        ymax=Avg_Mean_U + SD_Mean_U,
                        group=Start_U,
                        fill=factor(Start_U)),
                    # inherit.aes=TRUE,
                    alpha=0.15,
                    color=NA
                ) +
                geom_line(linewidth=0.8) +
                geom_hline(yintercept=best_const_U, linetype="dashed", linewidth=0.8, color="black") +
                labs(
                    x="Generation",
                    y="Mean U (across replicates)",
                    color="Start U"
                ) +
                theme_bw()


p_fit_logy <- ggplot(avg_fit_log, aes(x=Generation, y=y, color=factor(Start_U), group=Start_U)) +
                geom_ribbon(
                    aes(ymin=ymin,
                        ymax=ymax,
                        group=Start_U,
                        fill=factor(Start_U)),
                    # inherit.aes=TRUE,
                    alpha=0.15,
                    color=NA
                ) +
                geom_line(linewidth=0.8) +
                geom_ribbon(
                    data=avg_best_U_fit_log,
                    aes(x=Generation,
                        ymin = ymin,
                        ymax = ymax),
                    inherit.aes=FALSE,
                    fill="grey",
                    alpha=0.35
                ) +
                geom_line(
                    data=avg_best_U_fit_log,
                    aes(x=Generation, y=y),
                    inherit.aes=FALSE,
                    linewidth=0.8,
                    linetype="dashed",
                    color="black"
                ) +
                scale_y_log10() +
                labs(
                    x="Generation",
                    y="Mean fitness (across replicates)",
                    color="Start U"
                ) +
                theme_bw()

p_mut_logy <- ggplot(avg_mut_log, aes(x=Generation, y=y, color=factor(Start_U), group=Start_U)) +
                geom_ribbon(
                    aes(ymin=ymin,
                        ymax=ymax,
                        group=Start_U,
                        fill=factor(Start_U)),
                    # inherit.aes=TRUE,
                    alpha=0.15,
                    color=NA
                ) +
                geom_line(linewidth=0.8) +
                geom_hline(yintercept=best_const_U, linetype="dashed", linewidth=0.8, color="black") +
                scale_y_log10() +
                labs(
                    x="Generation",
                    y="Mean U (across replicates)",
                    color="Start U"
                ) +
                theme_bw()

p_mut_median_logy <- ggplot(median_mut_df, aes(x=Generation, y=Median_Mean_U, color=factor(Start_U), group=Start_U)) +
                       geom_ribbon(
                            aes(ymin=Q1_Mean_U, 
                                ymax=Q3_Mean_U,
                                group=Start_U,
                                fill=factor(Start_U)), 
                            alpha=0.15,
                            color=NA
                       ) +
                       geom_line(linewidth=0.8) +
                       geom_hline(yintercept=best_const_U, linetype="dashed", linewidth=0.8, color="black") +
                       scale_y_log10() +
                       labs(
                            x="Generation",
                            y="Median U (across replicates)",
                            color="Start_U"
                       )
ggsave("fig2_fit_linear_diag.pdf", p_fit_linear, width=7, height=5)
ggsave("fig2_fit_logy_diag.pdf", p_fit_logy, width=7, height=5)
ggsave("fig2_mut_linear_diag.pdf", p_mut_linear, width=7, height=5)
ggsave("fig2_mut_logy_diag.pdf", p_mut_logy, width=7, height=5)
ggsave("fig2_mut_median_logy_diag.pdf", p_mut_median_logy, width=7, height=5)