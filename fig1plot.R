# For each genome-wide mutation rate U, this script collects the final generation
# fitness and mutation rate from each replicate and plots log fitness as a
# function of log mutation rate (with replicate aggregation).

library(fs)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

root_const <- "./fig1/"
root_evolve <- "./fig1-evolve/"

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
files_evolve <- dir_ls(root_evolve, recurse=TRUE, regexp="history_.*_(fitness|mutation)\\.csv$")
# files_evolve


# Constant mutation stats summary
const_long <- compile_all(files_const)
const_wide <- const_long %>% 
              pivot_wider(names_from=kind, values_from=value) %>%
              filter(is.finite(fitness), is.finite(mutation))
const_summary <- const_wide %>%
                 group_by(U) %>%
                 summarize(
                    avg_fit=mean(fitness),
                    sd_fit=sd(fitness),
                    sem_fit=sd(fitness) / sqrt(n()) # n() is # of replicates, SEM gets smaller as you have more reps
                 ) %>%
                 arrange(U) # sort rows by increasing U
const_summary

# Evolving mutation stats summary
evo_long <- compile_all(files_evolve)
evo_wide <- evo_long %>%
            pivot_wider(names_from=kind, values_from=value) %>%
            rename(endU=mutation, startU=U) %>%
            filter(is.finite(fitness), is.finite(endU))
evo_summary <- evo_wide %>%
               group_by(startU) %>%
               summarize(
                    avg_fit=mean(fitness),
                    sd_fit=sd(fitness),
                    sem_fit=sd(fitness) / sqrt(n()),
                    avg_endU=mean(endU),
                    sd_endU=sd(endU),
                    sem_endU=sd(endU) / sqrt(n())
               ) %>%
               arrange(startU)
evo_summary

p <- ggplot() +
     geom_point(
        data=const_wide,
        aes(x=log10(mutation), y=log10(fitness)),
        alpha=0.25, size=1
     ) +
     geom_line(
        data=const_summary,
        aes(x=log10(U), y=log10(avg_fit)),
        linewidth=1
     ) +
     # evolving-mutation run vertical error bars (fitness sem)
     geom_errorbar(
        data=evo_summary,
        aes(
            x=log10(avg_endU), 
            ymin=log10(avg_fit - sem_fit),
            ymax=log10(avg_fit + sem_fit)
        )
     ) +
     # evolving-mutation run horizontal error bars (mutation sem)
     geom_errorbarh(
        data=evo_summary,
        aes(
            y=log10(avg_fit),
            xmin=log10(avg_endU - sem_endU),
            xmax=log10(avg_endU + sem_endU)
        )
     ) +
     geom_point(
        data=evo_summary,
        aes(x=log10(avg_endU), y=log10(avg_fit)),
        size=3
     ) +
     # labels for evolve points
     geom_text(
        data=evo_summary,
        aes(x=log10(avg_endU), y=log10(avg_fit),
            label=paste0("start=", format(startU, scientific=TRUE))),
            size=3, vjust=-1
     ) +
     labs(
        x="log10(Mutation Rate)",
        y="log10(Fitness)"
     ) +
     theme_bw()

    
    ggsave("fig1_loglog.pdf", p, width=7, height=5)


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
     # evolving-mutation run vertical error bars (fitness sem)
     geom_errorbar(
        data=evo_summary,
        aes(
            x=(avg_endU), 
            ymin=(avg_fit - sem_fit),
            ymax=(avg_fit + sem_fit)
        )
     ) +
     # evolving-mutation run horizontal error bars (mutation sem)
     geom_errorbarh(
        data=evo_summary,
        aes(
            y=(avg_fit),
            xmin=(avg_endU - sem_endU),
            xmax=(avg_endU + sem_endU)
        )
     ) +
     geom_point(
        data=evo_summary,
        aes(x=(avg_endU), y=(avg_fit)),
        size=3
     ) +
     # labels for evolve points
     geom_text(
        data=evo_summary,
        aes(x=(avg_endU), y=(avg_fit),
            label=paste0("start=", format(startU, scientific=TRUE))),
            size=3, vjust=-1
     ) +
     labs(
        x="(Mutation Rate)",
        y="(Fitness)"
     ) +
     theme_bw()

    
    ggsave("fig1_linear.pdf", p2, width=7, height=5)