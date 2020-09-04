library(exactICLforLCA)
library(BayesLCA)
library(ggplot2)

# Read data
data(Alzheimer)

## set CONSTANTS
alpha_var <- 1
beta_var <- 1
delta_var <- 1
# G <- as.integer(readline(prompt="Enter the number of groups: "))
G <- 7

# LCA fit
fit <- blca.em(Alzheimer, G)
Z <- unMAP(MAP(Zscore(Alzheimer, fit)))

Y <- as.matrix(Alzheimer)

res <- ICLFit(Z, Y, G, alpha_var, beta_var, delta_var)
ICLsummary()


# extra - BEGIN

z_new <- res$Z_max
z_new_sum <- apply(z_new, 2, sum)

apply(z_new, 2, sum)/nrow(z_new)

theta_hat_new <- t(z_new) %*% as.matrix(Alzheimer)
theta_hat_new <- theta_hat_new / z_new_sum
theta_hat_new

r <- barplot(theta_hat_new, main="'exactICLforLCA' model fit",
             xlab="SYMPTOMS",col=c('#004c6d','#7fa9c9'),
             legend = rownames(theta_hat_new), beside=TRUE)

legend("topleft",
       legend = c("Group 1"),
       col = c('#004c6d'),
       pch = c(15,25),
       bty = "n",
       pt.cex = 2,
       cex = 1.2,
       text.col = "black",
       horiz = F ,
       inset = c(0.05, 0.05))

legend("topleft",
       legend = c("Group 2"),
       col = c('#7fa9c9'),
       pch = c(15,25),
       bty = "n",
       pt.cex = 2,
       cex = 1.2,
       text.col = "black",
       horiz = F ,
       inset = c(0.05, 0.1))

# extra - EnD

#extra - BEGIN

fit <- blca.em(Alzheimer, 2)
Z <- unMAP(MAP(Zscore(Alzheimer, fit)))

Zsum <- apply(Z, 2, sum)

apply(Z, 2, sum)/nrow(Z)

theta_hat <- t(Z) %*% as.matrix(Alzheimer)
theta_hat <- theta_hat / Zsum
theta_hat

r <- barplot(theta_hat, main="Expectation–Maximization algorithm (blca.em) fit",
             xlab="SYMPTOMS",col=c('#004c6d','#7fa9c9'),
             legend = rownames(theta_hat), beside=TRUE)
usr <- par("usr")
par(usr=c(usr[0:2], 0, 1))
axis(2,at=seq(0,1,1))


legend("topleft",
       legend = c("Group 1"),
       col = c('#004c6d'),
       pch = c(15,25),
       bty = "n",
       pt.cex = 2,
       cex = 1.2,
       text.col = "black",
       horiz = F ,
       inset = c(0.05, 0.05))

legend("topleft",
       legend = c("Group 2"),
       col = c('#7fa9c9'),
       pch = c(15,25),
       bty = "n",
       pt.cex = 2,
       cex = 1.2,
       text.col = "black",
       horiz = F ,
       inset = c(0.05, 0.1))


## extra - END


