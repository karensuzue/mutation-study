# For each genome-wide mutation rate U, this script collects the final generation
# fitness and mutation rate from each replicate and plots log fitness as a
# function of log mutation rate (with replicate aggregation).

library(fs)
library(readr)
library(ggplot2)

root <- "./fig1/"

# Grab the last row of a file
read_tail <- function(file) {
    df <- read_csv(file, show_col_types = FALSE)
    tail(df, n = 1)
}

# Grab all history files we need
files <- dir_ls(root, recurse = TRUE, regexp = "history_.*_(fitness|mutation)\\.csv$")
files[1]

# Parse U, rep, kind from filename
get_meta <- function(file) {
    nm <- basename(file)
    m <- regexec("^history_([^_]+)_([0-9]+)_(fitness|mutation)\\.csv$", nm)
    parts <- regmatches(nm, m)[[1]]
    if (length(parts) == 0) return(NULL)
    list(U=as.numeric(parts[2]), rep=as.integer(parts[3]), kind=parts[4])
}

meta<- get_meta(files[1])
data <- read.csv(files[1])
last <- data[nrow(data), ]
last
val <- if (meta$kind == "fitness") last$Mean_F else last$Mean_U 
val

rows <- list()
for (f in files) {
    meta <- get_meta(f)
    if (is.null(meta)) next

    data <- read.csv(f)
    if (nrow(data) == 0) next 

    # get last generation row
    last <- data[nrow(data), ] 
    # get the "Mean_*" value from that row
    val <- if (meta$kind == "fitness") last$Mean_F else last$Mean_U
    # append a new row to 'rows'
    rows[[length(rows)+1]] <- data.frame(U=meta$U, rep=meta$rep, kind=meta$kind, value=val)
}

# rows is a list of small dfs, so we combine them into a big df
rows_df <- do.call(rbind, rows)

# get rid of the 'kind' column, create 'fitness' and 'mutationrate' columns
wide_rows_df <- reshape(rows_df, idvar=c("U", "rep"), timevar="kind", direction="wide")

# rename 'value.fitness' and 'value.mutation' columns to 'fitness' and 'mutation' respectively
names(wide_rows_df) <- sub("^value\\.", "", names(wide_rows_df))

plot_df <- subset(wide_rows_df, is.finite(fitness) & is.finite(mutation))
plot_df

# Across replicates
mean_fitness <- aggregate(fitness ~ U, data=plot_df, FUN=mean)
sd_fitness <- aggregate(fitness ~ U, data=plot_df, FUN=sd)

# Find the U with the largest mean fitness
best_idx <- which.max(mean_fitness$fitness)
best_U <- mean_fitness$U[best_idx]
best_mean_fit <- mean_fitness$fitness[best_idx]

cat("Best mean fitness:", best_mean_fit, "at U =", format(best_U, scientific=TRUE), "\n")

# Linear scale plot
pdf("fig1_linear.pdf", width=7, height=5)
plot(plot_df$mutation, plot_df$fitness,
     xlab="Mutation rate (linear)", ylab="Fitness (linear)",
     pch=16, col=rgb(0,0,0,0.25))

lines(mean_fitness$U, mean_fitness$fitness, lwd=2)
points(mean_fitness$U, mean_fitness$fitness, pch=16, cex=1.1)

dev.off()

# Log-x scale plot
pdf("fig1_logx.pdf", width=7, height=5)

plot(log10(plot_df$mutation), plot_df$fitness,
     xlab="log(Mutation rate)", ylab="Fitness (linear)",
     pch=16, col=rgb(0,0,0,0.25))

lines(log10(mean_fitness$U), mean_fitness$fitness, lwd=2)
points(log10(mean_fitness$U), mean_fitness$fitness, pch=16, cex=1.1)

dev.off()


# Log-log plot
if (all(plot_df$fitness > 0)) {
  pdf("fig1_loglog.pdf", width=7, height=5)

  plot(log10(plot_df$mutation), log10(plot_df$fitness),
       xlab="log(Mutation rate)", ylab="log(Fitness)",
       pch=16, col=rgb(0,0,0,0.25))

  lines(log10(mean_fitness$U), log10(mean_fitness$fitness), lwd=2)
  points(log10(mean_fitness$U), log10(mean_fitness$fitness), pch=16, cex=1.1)

  dev.off()
}