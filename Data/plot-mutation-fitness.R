# For each genome-wide mutation rate U, this script collects the final generation
# fitness and mutation rate from each replicate and plots log fitness as a
# function of log mutation rate (with replicate aggregation).

library(fs)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)


root_const <- "/mnt/d/MINE/sse-staticenv-constmut/"
root_evolve <- "./diag-exploit-vc/evolve-diag-final/"

# Include results from evolving mutation runs?
include_evolve <- FALSE

# Function to grab the last row of a file
read_last_row <- function(file) {
    df <- read_csv(file, show_col_types=FALSE)
    if (nrow(df) == 0) return(NULL)
    return(df[nrow(df), , drop=FALSE])
}

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

# Function to compile all final generation results
# Takes in a vector of filenames
# Outputs a single data frame where each row corresponds to one file (one rep + kind)
# i.e., rep | kind | value | U 
compile_all <- function(files) {
    rows <- list()
    for (f in files) {
        # Get metadata
        meta <- get_meta(f)
        if (is.null(meta)) next
        # Grab last row
        last_row <- read_last_row(f)
        if (is.null(last_row)) next
        # Get the "Mean_*" value from that row
        val <- if (meta$kind == "fitness") last_row$Mean_F else last_row$Mean_U
        # Append a new row to 'rows'
        rows[[length(rows)+1]] <- data.frame(rep=meta$rep, kind=meta$kind, value=as.numeric(val), U=meta$U)
    }
    # rows_df <- do.call(rbind, rows)
    rows_df <- bind_rows(rows) # rows is a list of rows, so we bind
    return(rows_df)
}

# Grab all history files we need
files_const <- dir_ls(root_const, recurse=TRUE, regexp="history_.*_(fitness|mutation)\\.csv$")
# files_const

if (include_evolve) {
   files_evolve <- dir_ls(root_evolve, recurse=TRUE, regexp="history_.*_(fitness|mutation)\\.csv$")
   # files_evolve
}

# -------------------------
# Constant mutation summary
# -------------------------
const_long <- compile_all(files_const)
const_wide <- const_long %>% 
              pivot_wider(names_from=kind, values_from=value) %>%
              filter(is.finite(fitness), is.finite(mutation)) %>%
              mutate(error = -fitness)

const_summary <- const_wide %>%
                 group_by(U) %>%
                 summarize(
                    avg_fit=mean(fitness),
                    sd_fit=sd(fitness),
                    sem_fit=sd_fit / sqrt(n()), # n() is # of replicates, SEM gets smaller as you have more reps

                    avg_err=mean(error),
                    sd_err=sd(error),
                    sem_err=sd_err / sqrt(n()),
                    .groups="drop"
                 ) %>%
                 arrange(U) # sort rows by increasing U
# const_summary
write_csv(const_summary, "const_summary.csv")


best_U_avg <- const_summary %>%
              slice_max(order_by=avg_fit, n=1, with_ties=TRUE)

best_U_avg

# These are numbers!!
best_U_x_log <- log10(best_U_avg$U[1])
best_U_x_lin <- best_U_avg$U[1]


worst_U_avg <- const_summary %>%
               slice_min(order_by=avg_fit, n=1, with_ties=TRUE)
worst_U_avg
worst_U_x_log <- log10(worst_U_avg$U[1])
worst_U_x_lin <- worst_U_avg$U[1]


# -------------------------
# Evolving mutation summary
# -------------------------
if (include_evolve) {
   evo_long <- compile_all(files_evolve)
   evo_wide <- evo_long %>%
               pivot_wider(names_from=kind, values_from=value) %>%
               rename(endU=mutation, startU=U) %>%
               filter(is.finite(fitness), is.finite(endU)) %>%
               mutate(error = -fitness)

   evo_summary <- evo_wide %>%
                  group_by(startU) %>%
                  summarize(
                     avg_fit=mean(fitness),
                     sd_fit=sd(fitness),
                     sem_fit=sd_fit / sqrt(n()),

                     avg_err=mean(error),
                     sd_err=sd(error),
                     sem_err=sd_err / sqrt(n()),

                     avg_endU=mean(endU),
                     sd_endU=sd(endU),
                     sem_endU=sd_endU / sqrt(n()),
                     .groups = "drop"
                  ) %>%
                  arrange(startU)
   # evo_summary

   write_csv(evo_summary, "evo_summary.csv")
}

# -------------------------
# Log-log plot
# -------------------------
p <- ggplot() +
     geom_point(
        data=const_wide,
        aes(x=log10(mutation), y=log10(pmax(error, 1e-12))),
        alpha=0.25, size=1
     ) +
     geom_line(
        data=const_summary,
        aes(x=log10(U), y=log10(pmax(avg_err, 1e-12))),
        linewidth=1
     ) +
     # best constant U vertical line
     geom_vline(
         xintercept = best_U_x_log,
         linetype="dashed",
         linewidth=0.8
     ) +
     # worst constant U vertical line
     geom_vline(
         xintercept = worst_U_x_log,
         linetype="dotted",
         linewidth=0.8
     ) +
     annotate(
         "text",
         x=best_U_x_log, # horizontal pos of text
         y=min(log10(pmax(const_wide$error, 1e-12)), na.rm=TRUE), # vertical pos
         label=paste0("best=", signif(best_U_x_lin, 3)), # show raw/non-log U
         vjust=1.5,
         size=3
     ) +
     annotate(
         "text",
         x=worst_U_x_log, # horizontal pos of text
         y=min(log10(pmax(const_wide$error, 1e-12)), na.rm=TRUE), # vertical pos
         label=paste0("worst=", signif(worst_U_x_lin, 3)), # show raw/non-log U
         vjust=1.5,
         size=3
     ) +
     labs(
        x="log10(Mutation Rate)",
        y="log10(Error)"
     ) +
     theme_bw()

if (include_evolve) {
   p <- p + 
      # evolving-mutation run vertical error bars (fitness sem)
      geom_errorbar(
         data=evo_summary,
         aes(
            x=log10(avg_endU), 
            ymin=log10(pmax(avg_err - sem_err, 1e-12)),
            ymax=log10(pmax(avg_err + sem_err, 1e-12)),
            color=factor(startU)
         )
      ) +
      # evolving-mutation run horizontal error bars (mutation sem)
      geom_errorbarh(
        data=evo_summary,
        aes(
            y=log10(pmax(avg_err, 1e-12)),
            xmin=log10(pmax(avg_endU - sem_endU, 1e-12)),
            xmax=log10(pmax(avg_endU + sem_endU, 1e-12)),
            color=factor(startU)
         )
      ) +
      geom_point(
        data=evo_summary,
        aes(x=log10(avg_endU), y=log10(pmax(avg_err, 1e-12)), color=factor(startU)),
        size=3
      )
}

ggsave(
   if (include_evolve) "fig_loglog_evolvemut.pdf" else "fig_loglog_constmut.pdf",
   p,
   width = 7,
   height = 5
)

# -------------------------
# Linear plot
# -------------------------
p2 <- ggplot() +
     geom_point(
        data=const_wide,
        aes(x=mutation, y=fitness),
        alpha=0.25, size=1
     ) +
     geom_line(
        data=const_summary,
        aes(x=(U), y=(avg_fit)),
        linewidth=1
     ) +
     # best constant U vertical line
     geom_vline(
         xintercept = best_U_x_lin,
         linetype="dashed",
         linewidth=0.8
     ) +
     # worst constant U vertical line
     geom_vline(
         xintercept = worst_U_x_lin,
         linetype="dotted",
         linewidth=0.8
     ) +
      annotate(
         "text",
         x=best_U_x_lin, # horizontal pos of text
         y=min(const_wide$fitness, na.rm=TRUE), # vertical pos
         label=paste0("best=", signif(best_U_x_lin, 3)), # show raw/non-log U
         vjust=1.5,
         size=3
     ) +
      annotate(
         "text",
         x=worst_U_x_lin, # horizontal pos of text
         y=min(const_wide$fitness, na.rm=TRUE), # vertical pos
         label=paste0("worst=", signif(worst_U_x_lin, 3)), # show raw/non-log U
         vjust=1.5,
         size=3
     ) +

     labs(
        x="(Mutation Rate)",
        y="(Fitness)"
     ) +
     theme_bw()

if (include_evolve) {
   p2 <- p2 +
      # evolving-mutation run vertical error bars (fitness sem)
      geom_errorbar(
         data=evo_summary,
         aes(
               x=(avg_endU), 
               ymin=(avg_fit - sem_fit),
               ymax=(avg_fit + sem_fit),
               color=factor(startU)
         )
      ) +
      # evolving-mutation run horizontal error bars (mutation sem)
      geom_errorbarh(
         data=evo_summary,
         aes(
               y=(avg_fit),
               xmin=(avg_endU - sem_endU),
               xmax=(avg_endU + sem_endU),
               color=factor(startU)
         )
      ) +
      geom_point(
         data=evo_summary,
         aes(x=(avg_endU), y=(avg_fit), color=factor(startU)),
         size=3
      )
}

ggsave(
   if (include_evolve) "fig_linear_evolvemut.pdf" else "fig_linear_constmut.pdf",
   p2,
   width = 7,
   height = 5
)