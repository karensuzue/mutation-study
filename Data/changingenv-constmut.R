# This script produces:
# 1) A heatmap of final mean fitness, 
#    showing the interaction effects of environment-change step (X)
#    and genomic mutation rate (Y), kept constant.
# 2) A timeseries 

library(tidyverse)
library(patchwork)   # for multi-panel layouts
library(scales)      # for scientific notation on axes
library(fs)

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------
DATA_DIR <- "/mnt/d/MINE/sse-changingenv-constmut/" 
OUT_DIR <- "."

# Values to plot
U_RATES    <- c("2.5119e-04", "3.1623e-04", "6.6667e-01", "1.0000e+00", "1.0000e+02")
CHANGE_STEPS <- c(2000, 1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01)

# These appear in filenames
CHANGE_STEPS_FILE <- c(2000, 1000, 500, 200, 100, 50, 20, 10, 5, 2, 1) 
CHANGE_PER_FILE <- c(1, 2, 5, 10, 20, 50, 100)

# Toggle best or mean fitness
WHICH_FITNESS <- "Best" # "Best", "Mean"
FITNESS_NAME <- paste0(WHICH_FITNESS, "_F")

# Reading every generation is heavy
# Set to 1 to read everything, or e.g. 1000 to keep every 1000th generation.
SUBSAMPLE_STEP <- 1000

# Filename pattern:
# history_U{tag}_change{change_per_step}_per{change_env_step}_{rep}_fitness.csv
# history_U{tag}_change{change_per_step}_per{change_env_step}_{rep}_mutation.csv (ignored since U is kept constant)
REGEX_PATTERN <- "^history_U([^_]+)_change([0-9]+)_per([0-9]+)_([0-9]+)_fitness\\.csv$"

# ---------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------
parse_filename <- function(file) {
    nm <- basename(file)
    m <- regexec(REGEX_PATTERN, nm)
    parts <- regmatches(nm, m)[[1]]
    if (length(parts) == 0) return(NULL)
    parsed <- list(
        U = parts[2],
        change_per = as.integer(parts[3]),
        change_env_step = as.integer(parts[4]),
        rep = as.integer(parts[5])
        # type = parts[6]
    )
    return(parsed)
}

# Function to grab the last row of a file
read_last_row <- function(file) {
    df <- read_csv(file, show_col_types=FALSE)
    if (nrow(df) == 0) return(NULL)
    return(df[nrow(df), , drop=FALSE])
}

read_file <- function(file) {
    meta <- parse_filename(file)
    df <- read_csv(file, show_col_types=FALSE) %>%
            filter(Generation %% SUBSAMPLE_STEP == 0) %>%
            select(-Fittest_ID) %>%
            mutate(
                U = meta$U,
                effective_step = meta$change_env_step / meta$change_per
            )
    return(df)
}

# Function to compile all final generation fitness results
compile_final <- function(files) { # 'files' is a vector of filenames
    rows <- list()
    # One row corresponds to one file
    for (f in files) {
        # Get metadata
        meta <- parse_filename(f)
        if (is.null(meta)) next

        # Grab last row
        last_row <- read_last_row(f)
        if (is.null(last_row)) next

        # Get the "Mean/Best_*" value from that row
        val <- last_row[[FITNESS_NAME]]

        # Append a new row to 'rows'
        rows[[length(rows)+1]] <- data.frame(U = meta$U, 
                                            change_per = meta$change_per, 
                                            change_env_step = meta$change_env_step,
                                            rep = meta$rep, 
                                            fitness = as.numeric(val))
    }
    rows_df <- bind_rows(rows) # rows is a list of rows, so we bind
    return(rows_df) 
}

# ---------------------------------------------------------------------
# THEME
# ---------------------------------------------------------------------
custom_theme <- function() {
    theme_bw(base_size = 11) +
    theme(
        strip.background  = element_rect(fill = "grey92", colour = NA),
        strip.text        = element_text(size = 10, face = "bold"),
        panel.grid.minor  = element_blank(),
        legend.position   = "right",
        legend.key.size   = unit(0.5, "cm")
    )
}

# Colour palette: one colour per change_env_step (11 levels)
step_palette <- colorRampPalette(c("#1a6ca8", "#1d9e75", "#f0a500", "#c0392b"))(
  length(CHANGE_STEPS)
)
names(step_palette) <- as.character(sort(CHANGE_STEPS, decreasing = FALSE))

# ---------------------------------------------------------------------
# LOCATE & SUMMARIZE
# ---------------------------------------------------------------------
files <- dir_ls(DATA_DIR, recurse=TRUE)
files <- files[str_detect(basename(files), REGEX_PATTERN)]

all_df <- list()
for (f in files) {
    all_df[[length(all_df) + 1]] <- read_file(f)
}
all_df <- bind_rows(all_df)
all_df$U_fac <- factor(all_df$U, levels = U_RATES)
all_df$step_fac <- factor(all_df$effective_step, 
                                levels = sort(unique(CHANGE_STEPS), 
                                                decreasing = FALSE)) # slow to fast change                               

# Include 'error' column, collapse change_per x change_env_step pairs into a single value
final_rows <- compile_final(files) %>% 
                mutate(
                    error = -fitness,
                    effective_step = change_env_step / change_per
                )
final_summary <- final_rows %>%
                    group_by(U, effective_step) %>%
                    summarize(
                        avg_fit = mean(fitness),
                        sd_fit = sd(fitness),
                        sem_fit = sd_fit / sqrt(n()),

                        avg_err = mean(error),
                        sd_err = sd(error),
                        sem_err = sd_err / sqrt(n()),
                        .groups = "drop"
                    )

final_summary$U_fac <- factor(final_summary$U, levels = U_RATES)
final_summary$step_fac <- factor(final_summary$effective_step, 
                                levels = sort(unique(CHANGE_STEPS), 
                                                decreasing = FALSE)) # slow to fast change

ts_summary <- all_df %>%
            group_by(U_fac, step_fac, Generation) %>%
            summarize(
                med_fit = median(.data[[FITNESS_NAME]]),
                q25_fit = quantile(.data[[FITNESS_NAME]], 0.25),
                q75_fit = quantile(.data[[FITNESS_NAME]], 0.75),

                med_err = -med_fit,
                q25_err = -q75_fit,
                q75_err = -q25_fit,

                .groups = "drop"
            )

eps <- 1e-8
ts_summary_log <- ts_summary %>%
  mutate(
    med_err = pmax(med_err, eps),
    q25_err = pmax(q25_err, eps),
    q75_err = pmax(q75_err, eps)
  )

# ---------------------------------------------------------------------
# PLOT 1: HEATMAP - mean final fitness, U x change_env_step
# ---------------------------------------------------------------------
plot_heatmap <- ggplot(final_summary,
                        aes(x = step_fac, y = U_fac, fill = avg_err)) +
                geom_tile(color = "white", linewidth = 0.5) +
                geom_text(aes(label = sprintf("%.2e\n±%.2e", avg_err, sd_err)),
                            size = 2.0, color = "white") +
                scale_fill_gradientn( 
                    colours  = c("#1a6ca8", "#1d9e75", "#f0a500"),
                    trans    = "log10", # error values span many orders of magnitude
                    labels   = label_scientific(),
                    name     = paste0(WHICH_FITNESS, " error (SSE)\n(final gen)")
                ) +
                scale_x_discrete(guide = guide_axis(angle = 45)) +
                labs(
                    title = paste0(WHICH_FITNESS, " final error across 20 replicates"),
                    x = "Rate of environment change (fast to slow)",
                    y = "Genome-wide mutation rate"
                ) +
                custom_theme() +
                theme(axis.text.x = element_text(size = 9))

# ---------------------------------------------------------------------
# PLOT 2: TIME SERIES - fitness over generations, faceted by U
# ---------------------------------------------------------------------
plot_ts <- ggplot(ts_summary,
                    aes(x = Generation, color = step_fac, fill = step_fac)) +
            geom_ribbon(aes(ymin = q25_err, ymax = q75_err), alpha = 0.15) +
            geom_line(aes(y = med_err), linewidth = 0.6) +
            facet_wrap(~ U_fac, ncol=2, scales = "free_y") +
            scale_colour_manual(values = step_palette, name = "change_env_step") +
            scale_fill_manual(values = step_palette,   name = "change_env_step") +
            scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
            labs(
                title    = paste0(WHICH_FITNESS, " error over time"),
                subtitle = "Median ± IQR across 20 replicates, faceted by U rate",
                x        = "Generation",
                y        = paste0(WHICH_FITNESS, " Error")
            ) +
            custom_theme()

plot_ts_log <- ggplot(ts_summary_log,
                    aes(x = Generation, color = step_fac, fill = step_fac)) +
            geom_ribbon(aes(ymin = q25_err, ymax = q75_err), alpha = 0.15) +
            geom_line(aes(y = med_err), linewidth = 0.6) +
            facet_wrap(~ U_fac, ncol=2, scales = "free_y") +
            scale_colour_manual(values = step_palette, name = "change_env_step") +
            scale_fill_manual(values = step_palette,   name = "change_env_step") +
            scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
            scale_y_log10() + 
            labs(
                title    = paste0(WHICH_FITNESS, " error over time"),
                subtitle = "Median ± IQR across 20 replicates, faceted by U rate",
                x        = "Generation",
                y        = paste0(WHICH_FITNESS, " Error")
            ) +
            custom_theme()

# ---------------------------------------------------------------------
# SAVE
# ---------------------------------------------------------------------
ggsave(file.path(OUT_DIR, "changingenv_constmut_heatmap.pdf"), plot_heatmap, width = 9, height = 5)
ggsave(file.path(OUT_DIR, "changingenv_constmut_ts.pdf"), plot_ts, width = 9, height = 5)
ggsave(file.path(OUT_DIR, "changingenv_constmut_tslog.pdf"), plot_ts_log, width = 9, height = 5)