# -------------------------------------------------------------------------
# Figure S.10
#
# Manuscript reference : Section S.7
# Description    : Relation between metric curvature and Alexandrov curvature
# -------------------------------------------------------------------------

pdf("FigureS10.pdf", width = 8, height = 8)   # open file device
old_par <- par(no.readonly = TRUE)   # save current settings

f <- function(x) {
  x
}

# Generate x values
x_values <- seq(-1.5, 1.5, by = 0.01)  # Define range from -5 to 5

# Compute f(x)
y_values <- sapply(x_values, f)

par(mfrow=c(1,1),mar = c(5,5,3,2))
# Plot the function
plot(x_values, y_values, type = "p", col = "coral3", lwd = 2,
     xlab = expression(kappa), ylab = expression(rho~"\'"),
     ylim = c(-1.5, 1.5), xaxt = "n", yaxt = "n",cex.lab=2)

# Add axis labels
axis(1, at = seq(-1, 1, by = 1),cex.axis=2)
axis(2, at = c(-1, 0, 1), labels = c(expression(rho~"\'"[-1]), "0", expression(rho~"\'"[1])), cex.axis=2)

# Add grid
grid()

par(old_par)                          # restore original settings
dev.off()