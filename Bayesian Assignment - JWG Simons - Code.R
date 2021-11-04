############################################################################
###################### Assignment Bayesian statistics ######################
############################################################################

##############################################################
###################### Data preparation ######################
##############################################################

# If necessary, install package "Stat2Data" from which to import the dataset. 
if(!require(Ecdat)) install.packages("Ecdat")
library(Ecdat)
# Load the data into R. 
data(Computers)
# Store the variables of interest in the dataset "Computers" into a generic data object.
dat <- Computers[, c("price", "speed", "hd", "ram", "ads", "premium")]
# Assign column names. 
colnames(dat) <- c("price","speed", "hd", "ram", "ads", "premium")
# Convert 'premium' variable from character to numeric. 
dat$premium <- as.numeric(as.factor(dat$premium)) - 1 

# Grand-mean centering of independent variables. 
dat$speed <- dat$speed - mean(dat$speed)
dat$hd <- dat$hd - mean(dat$hd)
dat$ram <- dat$ram - mean(dat$ram)
dat$ads <- dat$ads - mean(dat$ads)
dat$premium <- dat$premium - mean(dat$premium)
# Summary and standard deviations of grand-mean centered data. 
summary(dat)
apply(dat, 2, sd)

# MLE analysis second model. 
output.m2 <- lm(price ~ speed + hd + ram + ads + premium, data = dat)
# Request MLE output second model. 
summary(output.m2)
#Prerequisites to MH step in second model.
coef.m2 <- as.vector(output.m2$coefficients) # Store MLE coefficients for MH proposal density.
sd.m2 <- as.vector(sqrt(diag(vcov(output.m2)))) # Store MLE standard deviations for MH proposal density.

#########################################################################################################
###################### Estimation - Gibbs Sampler & Metropolis Hastings Algorithm  ######################
#########################################################################################################

# Gibbs sampler algorithm for estimating first model.  
gibbs.sampler.m1 <- function(){ # Define a function for executing a gibbs sampler. 
  
  # Store sample size length. 
  n <- nrow(dat)
  # Define matrix for storing conditional posterior parameter draws for two chains. 
  iter <- matrix(data = NA, nrow = 12000, ncol = 10)
  # Assign column names for respective parameters of the two chains. 
  colnames(iter) <- c("b0.c1", "b0.c2", "b1.c1", "b1.c2", "b2.c1", "b2.c2", "b3.c1", "b3.c2", "s2.c1", "s2.c2")
  # Specify initial values for the two chains. 
  inits <- c(1, -1, 0.6, -2, -1, 1, 0.5, -0.8, 1, 5)
  # Store initial values in storage matrices.
  iter[1, ] <- inits
  # Specify and store uninformative prio mean for intercepts and slopes. 
  mu.prior <- 0
  # Specify and store uninformative prior variance for slopes.
  tau2.prior <- 1000
  # Specify and store uninformative shape and scale parameter for variance.
  alpha.prior <- 0.001
  beta.prior <- 0.001
  
  # Set seed for reproducibility.
  set.seed(3791)
  # Estimation for the first proposed model. 
  for (j in 1 : 2){ # for 2 chains; 
    
    for (i in 1 : 11999){ # for 1 to 11999 instances;
      
      ## Conditional posterior of intercept (b0). 
      # Determine mean.
      sum.b0 <- sum(dat[, 1] - (iter[i, j + 2] * dat[, 2]) - (iter[i, j + 4] * dat[, 3]) - (iter[i, j + 6] * dat[, 4])) 
      mean.01 <- ((sum.b0 / as.numeric(iter[i, j + 8])) + (mu.prior / tau2.prior)) / ((n / as.numeric(iter[i, j + 8])) + (1 / tau2.prior))
      # Determine variance.
      tau2.01 <- 1 / ((n / as.numeric(iter[i, j + 8])) + (1 / tau2.prior))
      # Draw from conditional posterior.
      iter[i + 1, j] <- rnorm(n = 1, mean = mean.01, sd = tau2.01)
      
      ## Conditional posterior of predictor "speed' (b1). 
      # Determine mean.
      sum.b1 <- sum((dat[ , 2] * (dat[ , 1] - iter[i + 1, j] - (iter[i, j + 4] * dat[ , 3]) - (iter[i, j + 6] * dat[ , 4])))) 
      mean.11 <- ((sum.b1 / as.numeric(iter[i, j + 8])) + (mu.prior / tau2.prior))  / ((sum((dat[ , 2] ^ 2) / as.numeric(iter[i, j + 8]))) + (1 / tau2.prior))
      # Determine variance.
      tau2.11 <- 1 / ((sum(dat[ , 2] ^ 2) / as.numeric(iter[i, j + 8])) + (1 / tau2.prior))
      # Draw from conditional posterior.
      iter[i + 1, j + 2] <- rnorm(n = 1, mean = mean.11, sd = tau2.11) 
      
      ## Conditional posterior of predictor "hd" (b2).
      # Determine mean for the first chain.
      sum.b2 <- sum(dat[ , 3] * (dat[ , 1] - iter[i + 1, j] - (iter[i + 1, j + 2] * dat[ , 2]) - (iter[i, j + 6] * dat[ , 4]))) 
      mean.21 <- ((sum.b2 / as.numeric(iter[i, j + 8])) + (mu.prior / tau2.prior)) / ((sum((dat[ , 3] ^ 2) / as.numeric(iter[i, j + 8]))) + (1 / tau2.prior))
      # Determine variance for the first chain. 
      tau2.21 <- 1 / ((sum(dat[ , 3] ^ 2) / as.numeric(iter[i, j + 8])) + (1 / tau2.prior))
      # Draw from conditional posterior for the first chain.
      iter[i + 1, j + 4] <- rnorm(n = 1, mean = mean.21, sd = tau2.21)
      
      ## Conditional posterior of predictor "ram" (b3).
      # Determine mean for the first chain.
      sum.b3 <- sum(dat[ , 4] * (dat[ , 1] - iter[i + 1, j] - (iter[i + 1, j + 2] * dat[ , 2]) - (iter[i + 1, j + 4] * dat[ , 3]))) 
      mean.31 <- ((sum.b3 / as.numeric(iter[i, j + 8])) + (mu.prior / tau2.prior)) / ((sum((dat[ , 4] ^ 2) / as.numeric(iter[i, j + 8]))) + (1 / tau2.prior))
      # Determine variance for the first chain. 
      tau2.31 <- 1 / ((sum(dat[ , 4] ^ 2) / as.numeric(iter[i, j + 8])) + (1 / tau2.prior))
      # Draw from conditional posterior for the first chain.
      iter[i + 1, j + 6] <- rnorm(n = 1, mean = mean.31, sd = tau2.31)
      
      ## Conditional posterior of variance (s2).
      # Define alpha. 
      alpha <- (n / 2) + alpha.prior
      # Residuals. 
      res <- dat[ , 1] - (iter[i + 1, j] + (iter[i + 1, j + 2] * dat[ , 2]) + (iter[i + 1, j + 4] * dat[ , 3]) + (iter[i + 1, j + 6] * dat[ , 4])) 
      # Sum of squared residuals. 
      S <- sum((res ^ 2))
      # Define Beta based on sum of squared residuals. 
      beta <- (S / 2) + beta.prior
      # Draw frrom conditional posterior. 
      iter[i + 1, j + 8] <- 1 / (rgamma(n = 1, shape = alpha, rate = beta))
    }
  }
  
  # Define burn-in. 
  T <- 2000
  # Remove burn-in from iteration matrix.
  iter.burned <-  iter[-(1:T), ]
  # Subset iteration matrix into seperate matrices for respective chains, excluding burn-in. 
  iter.chain1 <- subset(iter.burned, TRUE, c("b0.c1", "b1.c1", "b2.c1", "b3.c1", "s2.c1"))
  iter.chain2 <- subset(iter.burned, TRUE, c("b0.c2", "b1.c2", "b2.c2", "b3.c2", "s2.c2"))
  # Return iteration matrix. 
  return(list(iter = iter, iter.burned = iter.burned, iter.chain1 = iter.chain1, iter.chain2 = iter.chain2)) 
  
}

# Metropolis-Hastings algorithm for estimating second model.  
metropolis.hastings.m2 <- function(){ # Define a function for gibbs sampler / MH sampler. 
  
  # Store sample size length. 
  n <- nrow(dat)
  # Define matrix for storing conditional posterior parameter draws for two chains. 
  iter <- matrix(data = NA, nrow = 12000, ncol = 14)
  # Assign column names for respective parameters of the two chains. 
  colnames(iter) <- c("b0.c1", "b0.c2", "b1.c1", "b1.c2", "b2.c1", "b2.c2", "b3.c1", "b3.c2", "b4.c1", "b4.c2", "b5.c1", "b5.c2", "s2.c1", "s2.c2")
  # Specify initial values for the two chains. 
  inits <- c(1, -1, 0.6, -2, -1, 1, 0.5, -0.8, 1, 0.2, 0.4, -0.2, 1, 5)
  # Store initial values in storage matrices.
  iter[1, ] <- inits
  # Specify and store uninformative prior mean for intercept and slopes.    
  mu.prior <- 0
  # Specify and store uninformative prior variance for intercept and slopes.
  tau2.prior <- 1000
  # Specify and store uninformative shape and scale parameter for variance.
  alpha.prior <- 0.001
  beta.prior <- 0.001
  # Specify and store uninformative m.0, nu.0, s2.0 for specification of prior t-distribution. 
  m.0 <- 0
  nu.0 <- 1
  s2.0 <- 0.001
  
  #Define a function for computing the conditional posterior of the slope of 'ads' in an MH step. 
  cond.ads <- function(B){
    # First term of density. 
    denst.1 <- -(B ^ 2) * (sum(dat[, 5] ^ 2) / (2 * as.numeric(iter[i, j + 12]))) 
    # # Second term of density. 
    denst.2 <- B * (sum(dat[ , 5] * (dat[ , 1] - iter[i + 1, j] - iter[i + 1, j + 2] * dat[ , 2] - iter[i + 1, j + 4] * dat[ , 3] - iter[i + 1, j + 6] * dat[ , 4] - iter[i, j + 10] * dat[ , 6]) / as.numeric(iter[i, j + 12]))) 
    # Prior t-distribution.
    tdist <- log((1 + ((prop - m.0) ^ 2 / (nu.0 * s2.0))) ^ (-(nu.0 + 1) / 2)) 
    # Conditional distribution. 
    cond <- (denst.1 + denst.2) + tdist 
    return(cond)
  }
  
  #Define a function for computing the conditional posterior of the slope of 'premium' in an MH step. 
  cond.prem <- function(B){
    # First term of density. 
    denst.1 <- -(B ^ 2) * (sum(dat[, 6] ^ 2) / (2 * as.numeric(iter[i, j + 12]))) 
    # # Second term of density. 
    denst.2 <- B * (sum(dat[ , 6] * (dat[ , 1] - iter[i + 1, j] - iter[i + 1, j + 2] * dat[ , 2] - iter[i + 1, j + 4] * dat[ , 3] - iter[i + 1, j + 6] * dat[ , 4] - iter[i + 1, j + 8] * dat[ , 5]) / as.numeric(iter[i, j + 12]))) 
    # Prior t-distribution.
    tdist <- log((1 + ((prop - m.0) ^ 2 / (nu.0 * s2.0))) ^ (-(nu.0 + 1) / 2)) 
    # Conditional distribution. 
    cond <- (denst.1 + denst.2) + tdist 
    return(cond)
  }
  
  # Set seed for reproducibility.
  set.seed(123)
  # Estimation for the first proposed model. 
  for (j in 1 : 2){ # for 2 chains; 
    
    for (i in 1 : 11999){ # for 1 to 9999 instances;
      
      ## Conditional posterior of intercept (b0). 
      # Determine mean.
      sum.b0 <- sum(dat[, 1] - (iter[i, j + 2] * dat[, 2]) - (iter[i, j + 4] * dat[, 3]) - (iter[i, j + 6] * dat[, 4]) - (iter[i, j + 8] * dat[, 5]) - (iter[i, j + 10] * dat[, 6])) 
      mean.01 <- ((sum.b0 / as.numeric(iter[i, j + 12])) + (mu.prior / tau2.prior)) / ((n / as.numeric(iter[i, j + 12])) + (1 / tau2.prior))
      # Determine variance.
      tau2.01 <- 1 / ((n / as.numeric(iter[i, j + 12])) + (1 / tau2.prior))
      # Draw from conditional posterior.
      iter[i + 1, j] <- rnorm(n = 1, mean = mean.01, sd = tau2.01)
      
      ## Conditional posterior of predictor "speed' (b1). 
      # Determine mean.
      sum.b1 <- sum((dat[ , 2] * (dat[ , 1] - iter[i + 1, j] - (iter[i, j + 4] * dat[ , 3]) - (iter[i, j + 6] * dat[ , 4]) - (iter[i, j + 8] * dat[ , 5]) - (iter[i, j + 10] * dat[ , 6])))) 
      mean.11 <- ((sum.b1 / as.numeric(iter[i, j + 12])) + (mu.prior / tau2.prior))  / ((sum((dat[ , 2] ^ 2) / as.numeric(iter[i, j + 12]))) + (1 / tau2.prior))
      # Determine variance.
      tau2.11 <- 1 / ((sum(dat[ , 2] ^ 2) / as.numeric(iter[i, j + 12])) + (1 / tau2.prior))
      # Draw from conditional posterior.
      iter[i + 1, j + 2] <- rnorm(n = 1, mean = mean.11, sd = tau2.11) 
      
      ## Conditional posterior of predictor "hd" (b2).
      # Determine mean for the first chain.
      sum.b2 <- sum(dat[ , 3] * (dat[ , 1] - iter[i + 1, j] - (iter[i + 1, j + 2] * dat[ , 2]) - (iter[i, j + 6] * dat[ , 4]) - (iter[i, j + 8] * dat[ , 5]) - (iter[i, j + 10] * dat[ , 6]))) 
      mean.21 <- ((sum.b2 / as.numeric(iter[i, j + 12])) + (mu.prior / tau2.prior)) / ((sum((dat[ , 3] ^ 2) / as.numeric(iter[i, j + 12]))) + (1 / tau2.prior))
      # Determine variance for the first chain. 
      tau2.21 <- 1 / ((sum(dat[ , 3] ^ 2) / as.numeric(iter[i, j + 12])) + (1 / tau2.prior))
      # Draw from conditional posterior for the first chain.
      iter[i + 1, j + 4] <- rnorm(n = 1, mean = mean.21, sd = tau2.21)
      
      ## Conditional posterior of predictor "ram" (b3).
      # Determine mean for the first chain.
      sum.b3 <- sum(dat[ , 4] * (dat[ , 1] - iter[i + 1, j] - (iter[i + 1, j + 2] * dat[ , 2]) - (iter[i + 1, j + 4] * dat[ , 3]) - (iter[i, j + 8] * dat[ , 5]) - (iter[i, j + 10] * dat[ , 6]))) 
      mean.31 <- ((sum.b3 / as.numeric(iter[i, j + 12])) + (mu.prior / tau2.prior)) / ((sum((dat[ , 4] ^ 2) / as.numeric(iter[i, j + 12]))) + (1 / tau2.prior))
      # Determine variance for the first chain. 
      tau2.31 <- 1 / ((sum(dat[ , 4] ^ 2) / as.numeric(iter[i, j + 12])) + (1 / tau2.prior))
      # Draw from conditional posterior for the first chain.
      iter[i + 1, j + 6] <- rnorm(n = 1, mean = mean.31, sd = tau2.31)
      
      ## Metropolis Hastings Step for conditional posterior of predictor "ads" (b4).
      # Sample candidate value from proposal distribution and random value from uniform distribution on interval 0 to 1. 
      prop <- rnorm(n = 1, mean = coef.m2[5], sd = 2.5 * sd.m2[5])  
      u <- runif(n = 1, min = 0, max = 1) 
      # Determine acceptance ratio.
      r <- cond.ads(prop) - cond.ads(as.numeric(iter[i, j + 8])) 
      if(r > log(u)){ # If acceptance ratio larger than log of u; 
        iter[i + 1, j + 8] <- prop # store proposal;
      }else { # otherwise; 
        iter[i + 1, j + 8] <- iter[i, j + 8] # store previous value.
      }
      
      ## Metropolis Hastings Step for conditional posterior of predictor "premium" (b5).
      # Sample candidate value from proposal distribution and random value from uniform distribution on interval 0 to 1. 
      prop <- rnorm(n = 1, mean = coef.m2[6], sd = 2.5 * sd.m2[6])  
      u <- runif(n = 1, min = 0, max = 1) 
      # Determine acceptance ratio.
      r <- cond.prem(prop) - cond.prem(as.numeric(iter[i, j + 10])) 
      if(r > log(u)){ # If acceptance ratio larger than log of u; 
        iter[i + 1, j + 10] <- prop # store proposal;
      }else { # otherwise; 
        iter[i + 1, j + 10] <- iter[i, j + 10] # store previous value.
      }
      
      ## Conditional posterior of variance (s2).
      # Define alpha. 
      alpha <- (n / 2) + alpha.prior
      # Residuals. 
      res <- dat[ , 1] - (iter[i + 1, j] + (iter[i + 1, j + 2] * dat[ , 2]) + (iter[i + 1, j + 4] * dat[ , 3]) + (iter[i + 1, j + 6] * dat[ , 4]) + (iter[i + 1, j + 8] * dat[ , 5]) + (iter[i + 1, j + 10] * dat[ , 6])) 
      # Sum of squared residuals. 
      S <- sum((res ^ 2))
      # Define Beta based on sum of squared residuals. 
      beta <- (S / 2) + beta.prior
      # Draw frrom conditional posterior. 
      iter[i + 1, j + 12] <- 1 / (rgamma(n = 1, shape = alpha, rate = beta))
    }
  }
  
  # Define burn-in. 
  T <- 2000
  # Remove burn-in from iteration matrix.
  iter.burned <-  iter[-(1:T), ]
  # Subset iteration matrix into seperate matrices for respective chains, excluding burn-in. 
  iter.chain1 <- subset(iter.burned, TRUE, c("b0.c1", "b1.c1", "b2.c1", "b3.c1", "b4.c1", "b5.c1", "s2.c1"))
  iter.chain2 <- subset(iter.burned, TRUE, c("b0.c2", "b1.c2", "b2.c2", "b3.c2", "b4.c2", "b5.c2", "s2.c2"))
  # Return iteration matrix. 
  return(list(iter = iter, iter.burned = iter.burned, iter.chain1 = iter.chain1, iter.chain2 = iter.chain2)) 
  
}

# Store raw iteration matrix, burn-in corrected iteration matrix, and seperate chain iteration matrices for first model. 
model1 <- gibbs.sampler.m1()
# Store raw iteration matrix, burn-in corrected iteration matrix, and seperate chain iteration matrices for second model. 
model2 <- metropolis.hastings.m2()

########################################################################
###################### Assessing Model Convergence #####################
########################################################################

# If necessary, install package "gridExtra" for combining caterpillar plots.  
if(!require(ggplot2)) install.packages("ggplot2")
# Require package "gridExtra". 
library("ggplot2")

### Caterpillar plots. 
# Define a function that overlays the iterations of both chains. 
caterpillar <- function(data, param.c1, param.c2, main){ # Parameter input from respective chains, and associated main title; 
  ggplot(data = as.data.frame(data)) +  # Input is iteration matrix with burn-in removed; 
    geom_line(mapping = aes(y = data[ , param.c1], x = c(1 : length(data[ , param.c1])), colour = "red", alpha = 0.2)) + # insert lines of first chain in red;
    geom_line(mapping = aes(y = data[ , param.c2], x = c(1 : length(data[ , param.c2])), colour = "blue", alpha = 0.2)) + # insert lines of second chain in blue; 
    ylab("") + # empty y-axis title;
    xlab("Iteration") + # name x-axis; 
    theme(legend.position = "none") + # remove legend;
    ggtitle("", main) # main title. 
}

# If necessary, install package "gridExtra" for combining caterpillar plots.  
if(!require(gridExtra)) install.packages("gridExtra")
# Require package "gridExtra". 
library("gridExtra")
# Arrange respective caterpillar plots into one plot for first model. 
grid.arrange(caterpillar(data = model1$iter.burned, param.c1 = 1, param.c2 = 2, main = "B0 - Intercept"), caterpillar(data = model1$iter.burned, param.c1 = 3, param.c2 = 4, main = "B1 - Speed"),
             caterpillar(data = model1$iter.burned, param.c1 = 5, param.c2 = 6, main = "B2 - Ram"), caterpillar(data = model1$iter.burned, param.c1 = 7, param.c2 = 8, main = "B3 - Hd"),
             caterpillar(data = model1$iter.burned, param.c1 = 9, param.c2 = 10, main = "S2 - Residual Variance"),
             ncol = 2, nrow = 4)
# Arrange respective caterpillar plots into one plot for second model. 
grid.arrange(caterpillar(data = model2$iter.burned, param.c1 = 1, param.c2 = 2, main = "B0 - Intercept"), caterpillar(data = model2$iter.burned, param.c1 = 3, param.c2 = 4, main = "B1 - Speed"),
             caterpillar(data = model2$iter.burned, param.c1 = 5, param.c2 = 6, main = "B2 - Ram"), caterpillar(data = model2$iter.burned, param.c1 = 7, param.c2 = 8, main = "B3 - Hd"),
             caterpillar(data = model2$iter.burned, param.c1 = 9, param.c2 = 10, main = "B4 - Ads"), caterpillar(data = model2$iter.burned, param.c1 = 11, param.c2 = 12, main = "B5 - Premium"), 
             caterpillar(data = model2$iter.burned, param.c1 = 13, param.c2 = 14, main = "S2 - Residual Variance"),
             ncol = 2, nrow = 4)

### Auto-correlation Plots and MH Acceptance ratio. 

## Auto-correlation Plots. 
autocorr <- function(dat, lag, par){ # Function with iteration data, autocorrelation lag, and parameter of interest as input. 
  ac.matrix <- matrix(data = NA, nrow = lag, ncol = par) # Define autocorrelation storage matrix;
  for (j in 1: par){ # For j number of parameters;
    t <- 0 # set time to 0 and; 
    ac.matrix[1 , j] <- 1 # assign correlation value 1 for lag 0;
    for (i in 1 : (lag - 1)){ # for i number of lags; 
      ac.matrix[i + 1, j] <- cor(dat[-((nrow(dat) - t) : nrow(dat)), j], dat[-(1 : (t + 1)), j]) # calculate autocorrelation between time t and t - 1;
      t <- t + 1 # Set next iteration in time; 
    }
  }
  return(ac.matrix) # return autocorrelation matrix.
}

# Store autocorrelations first model in dataframe. 
post.autocorr.m1 <- as.data.frame(autocorr(dat = model1$iter.burned, lag = 60, par = 10)) 
# Store autocorrelations second model in dataframe. 
post.autocorr.m2 <- as.data.frame(autocorr(dat = model2$iter.burned, lag = 60, par = 14)) 
# Column names first model. 
colnames(post.autocorr.m1) <- c("b0.c1", "b0.c2", "b1.c1", "b1.c2", "b2.c1", "b2.c2", "b3.c1", "b3.c2", "s2.c1", "s2.c2")
# Column names second model. 
colnames(post.autocorr.m2) <- c("b0.c1", "b0.c2", "b1.c1", "b1.c2", "b2.c1", "b2.c2", "b3.c1", "b3.c2", "b4.c1", "b4.c2", "b5.c1", "b5.c2", "s2.c1", "s2.c2")
# Define function for auto-correlation plotting. 
autocorr.plot <- function(post.autocorr, param.c1, param.c2, main){ # Function with as its input parameters of respective chains, and main title. 
  ggplot(data = post.autocorr) +  # Input data is posterior parameter autocorrelation; 
    geom_col(aes(x = 1 : nrow(post.autocorr), y = post.autocorr[, param.c1]), colour = "black", fill = "red", alpha = 0.5) + # Lags on x-axis, autocorrelation chain 1 on y-axis,
    geom_col(aes(x = 1 : nrow(post.autocorr), y = post.autocorr[, param.c2]), colour = "black", fill = "blue", alpha = 0.5) + # Lags on x-axis, autocorrelation chain 2 on y-axis,
    ylab("Correlation") + # Define y-axis title;
    xlab("Lag") + # Define x-axis title;
    ggtitle("", main) # Main title. 
}

# Arrange respective autocorrelation plots into one plot for first model. 
grid.arrange(autocorr.plot(post.autocorr = post.autocorr.m1, param.c1 = 1, param.c2 = 2, main = "B0 - Intercept"), autocorr.plot(post.autocorr = post.autocorr.m1, param.c1 = 3, param.c2 = 4, main = "B1 - Speed"),
             autocorr.plot(post.autocorr = post.autocorr.m1, param.c1 = 5, param.c2 = 6, main = "B2 - Ram"), autocorr.plot(post.autocorr = post.autocorr.m1, param.c1 = 7, param.c2 = 8, main = "B3 - Hd"),
             autocorr.plot(post.autocorr = post.autocorr.m1, param.c1 = 9, param.c2 = 10, main = "S2 - Residual Variance"),
             ncol = 2, nrow = 3)
# Arrange respective autocorrelation plots into one plot for second model. 
grid.arrange(autocorr.plot(post.autocorr = post.autocorr.m2, param.c1 = 1, param.c2 = 2, main = "B0 - Intercept"), autocorr.plot(post.autocorr = post.autocorr.m2, param.c1 = 3, param.c2 = 4, main = "B1 - Speed"),
             autocorr.plot(post.autocorr = post.autocorr.m2, param.c1 = 5, param.c2 = 6, main = "B2 - Ram"), autocorr.plot(post.autocorr = post.autocorr.m2, param.c1 = 7, param.c2 = 8, main = "B3 - Hd"),
             autocorr.plot(post.autocorr = post.autocorr.m2, param.c1 = 9, param.c2 = 10, main = "B4 - Ads"), autocorr.plot(post.autocorr = post.autocorr.m2, param.c1 = 11, param.c2 = 12, main = "B5 - Premium"),
             autocorr.plot(post.autocorr = post.autocorr.m2, param.c1 = 13, param.c2 = 14, main = "S2 - Residual Variance"),
             ncol = 2, nrow = 4)

### Markov chain error. 
# Obtain MC error by taking the standard deviation of parameter posterior and dividing by the square root of iterations. 
mc.err.m1 <- apply(rbind(model1$iter.chain1, model1$iter.chain2), 2, sd) / sqrt(20000) # Markov chain error in first model. 
mc.err.m2 <- apply(rbind(model2$iter.chain1, model2$iter.chain2), 2, sd) / sqrt(20000) # Markov chain error in second model.  
# Evaluate whether MC error is smaller than 5% of sample standard deviation in both chains. 
mc.err.m1 < (apply(rbind(model1$iter.chain1, model1$iter.chain2), 2, sd) / 100) * 5 # MC error smaller than 5% of sample SD.
mc.err.m2 < (apply(rbind(model2$iter.chain1, model2$iter.chain2), 2, sd) / 100) * 5 # MC error smaller than 5% of sample SD.

## Acceptance ratio of MH step in second model. 
length(unique(rbind(model2$iter.chain1, model2$iter.chain2)[, 5])) / length(rbind(model2$iter.chain1, model2$iter.chain2)[ , 5]) # Acceptance ratio of ~ 42%.
length(unique(rbind(model2$iter.chain1, model2$iter.chain2)[, 6])) / length(rbind(model2$iter.chain1, model2$iter.chain2)[ , 6]) # Acceptance ratio of ~ 44%.

#######################################################################################
###################### Posterior Predictive Checks ####################################
#######################################################################################

# Combine chains.
post.m1 <- as.data.frame(rbind(model1$iter.chain1, model1$iter.chain2)) # Combine chains for the first model. 
post.m2 <- as.data.frame(rbind(model2$iter.chain1, model2$iter.chain2)) # Combine chains for the second model. 
# Assign names to column of combined posterior matrix. 
colnames(post.m1) <- c("b0", "b1", "b2", "b3", "S2") # Column names for posterior matrix of the first model. 
colnames(post.m2) <- c("b0", "b1", "b2", "b3", "b4", "b5", "S2") # Column names for posterior matrix of the second model. 

ppc <- function(post, model, t){ # Function for computing the posterior predictive p-value for t samples.
  
  n <- nrow(dat) # Store number of rows in data. 
  
  hom.dat <- vector(length = t) # Define vector of length t for storing observed sample statistic homogeneity assumption.
  hom.dat.rep <- vector(length = t) # Define vector of length t for storing replicated sample statistic homogeneity assumption.
  
  norm.dat <- vector(length = t) # Define vector of length t for storing observed sample statistic normality assumption.
  norm.dat.rep <- vector(length = t) # Define vector of length t for storing replicated sample statistic normality assumption.
  
  for (i in 1 : t){ # For t samples; 
    
    if (model == 1){
      # Calculate residuals observed data first model. 
      dat.resid <- dat[ , 1] - (post[i, 1] + post[i, 2] * dat[ , 2] + post[i, 3] * dat[ , 3] + post[i, 4] * dat[ , 4])  
      # Replicate data set under first model. 
      dat.rep <- rnorm(n = n, mean = post[i, 1] + post[i, 2] * dat[ , 2] + post[i, 3] * dat[ , 3] + post[i, 4] * dat[ , 4], sd = sqrt(post[i, 5]))
      # Calculate residuals replicated data. 
      dat.rep.resid <-  dat.rep - (post[i, 1] + post[i, 2] * dat[ , 2] + post[i, 3] * dat[ , 3] + post[i, 4] * dat[ , 4]) 
    } else if (model == 2) {
      # Calculate residuals observed data second model. 
      dat.resid <- dat[ , 1] - (post[i, 1] + post[i, 2] * dat[ , 2] + post[i, 3] * dat[ , 3] + post[i, 4] * dat[ , 4] + post[i, 5] * dat[ , 5] + post[i, 6] * dat[ , 6])  
      # Replicate data set second model. 
      dat.rep <- rnorm(n = n, mean = post[i, 1] + post[i, 2] * dat[ , 2] + post[i, 3] * dat[ , 3] + post[i, 4] * dat[ , 4] + post[i, 5] * dat[ , 5] + post[i, 6] * dat[ , 6], sd = sqrt(post[i, 7]))
      # Calculate residuals replicated data. 
      dat.rep.resid <-  dat.rep - (post[i, 1] + post[i, 2] * dat[ , 2] + post[i, 3] * dat[ , 3] + post[i, 4] * dat[ , 4] + post[i, 5] * dat[ , 5] + post[i, 6] * dat[ , 6])     
    } 
    
    ## Inspect homogeneity of residuals assumption.
    # Test statistic observed data. Ratio of variance in two seperate halves of the observed sample.  
    sample.obs <- sample.int(n = length(dat.resid), size = floor(0.50 * length(dat.resid)), replace = F) # Halving indicator. 
    var.dat.1 <- var(dat.resid[sample.obs]) # Split data in one half;
    var.dat.2 <- var(dat.resid[-sample.obs]) # Split data in other half;
    hom.dat[i] <- var.dat.1 / var.dat.2 # Compute ratio of half variances.
    # Test statistic replicated data. Ratio of variance in two seperate halves of the replicated sample.  
    sample.rep <- sample.int(n = length(dat.rep.resid), size = floor(0.50 * length(dat.rep.resid)), replace = F)
    var.dat.rep.1 <- var(dat.rep.resid[sample.rep]) # Split data in one half;
    var.dat.rep.2 <- var(dat.rep.resid[-sample.rep]) # Split data in other half;
    hom.dat.rep[i] <- var.dat.rep.1 / var.dat.rep.2 # Compute ratio of half variances.
    
    ## Inspect normality of residuals assumption.
    # Use skewness as a statistic to determine the PPC for normality of residuals.
    norm.dat[i] <- (sum((dat.resid - mean(dat.resid)) ** 3) / n) / (((sum((dat.resid - mean(dat.resid)) ** 2)) / n) ** (3 / 2)) # Skewness statistic observed data. 
    norm.dat.rep[i] <-   (sum((dat.rep.resid - mean(dat.rep.resid)) ** 3) / n) / (((sum((dat.rep.resid - mean(dat.rep.resid)) ** 2)) / n) ** (3 / 2)) # skewness statistic replicated data. 
    
  }
  
  ppc.res.hom <- sum(((hom.dat.rep >  hom.dat) == TRUE), na.rm = T) / t # Calculate PPC homogeneity assumption. 
  ppc.res.norm <- sum(((abs(norm.dat.rep) >  abs(norm.dat)) == TRUE), na.rm = T) / t # Calculate PPC normality assumption. 
  
  return(list(ppc.res.hom = ppc.res.hom, ppc.res.norm = ppc.res.norm)) # Return linearity, homogeneity, and normality PPC values. 
  
}

# Request and print regression assumption PPC values. 
(ppc.m1 <- ppc(post = post.m1, model = 1, t = nrow(post.m1))) # PPC values homogeneity and normality assumptions first model. 
(ppc.m2 <- ppc(post = post.m2, model = 2, t = nrow(post.m2))) # PPC values homogeneity and normality assumptions second model. 

#########################################################################################
################# Model comparison and selection with the DIC  #######################
#########################################################################################

### DIC. 
DIC <- function(post, model){
  
  if (model == 1){
    ## Calculate Dhat. 
    n <- nrow(dat) # Store n. 
    # Calculate residual evaluated at posterior mean. 
    res.Dhat <- dat[ , 1] - (mean(post[ , 1]) + mean(post[ , 2]) * dat[ , 2] + mean(post[ , 3]) * dat[ , 3] + mean(post[ , 4]) * dat[ , 4])  
    loglik.Dhat <- (-(n / 2) * (log(2 * pi))) - (n * log(sqrt(mean(post[ , 5])))) - ((1 / (2 * mean(post[ , 5]))) * sum(res.Dhat ** 2))
    Dhat <- -2 * loglik.Dhat
    
    ## Calculate Dbar.
    loglik.Dbar <- vector(length = nrow(post)) 
    for (i in 1 : nrow(post)){ 
      res.Dbar <- dat[ , 1] - (post[i, 1] + post[i , 2] * dat[ , 2] + post[i , 3] * dat[ , 3] + post[i , 4] * dat[ , 4])  
      loglik.Dbar[i] <- (-(n / 2) * (log(2 * pi))) - (n * log(sqrt(post[i , 5]))) - ((1 / (2 * post[i , 5])) * sum(res.Dbar ** 2))
    }
    Dbar <- 1 / nrow(post) * sum(-2 * loglik.Dbar)
    pD <- 2 * (Dbar - Dhat) 
    DIC <- Dhat + pD
    
    cat("Mean deviance:", Dhat, "\n")
    cat("Penalty:", pD, "\n")
    cat("Mean deviance:", DIC, "\n")
    
  } else if (model == 2){
    ## Calculate Dhat. 
    n <- nrow(dat) # Store n. 
    # Calculate residual evaluated at posterior mean. 
    res.Dhat <- dat$price - (mean(post[, "b0"]) + mean(post[ , "b1"]) * dat$speed + mean(post[ , "b2"]) * dat$hd + mean(post[ , "b3"]) * dat$ram + mean(post[ , "b4"]) * dat$ads + mean(post[ , "b5"]) * dat$premium)  
    loglik.Dhat <- (-(n / 2) * (log(2 * pi))) - (n * log(sqrt(mean(post[ , "S2"])))) - ((1 / (2 * mean(post[ , "S2"]))) * sum(res.Dhat ** 2))
    Dhat <- -2 * loglik.Dhat
    
    ## Calculate Dbar.
    loglik.Dbar <- vector(length = nrow(post)) 
    for (i in 1 : nrow(post)){ 
      res.Dbar <- dat$price - (post[i, "b0"] + post[i , "b1"] * dat$speed + post[i , "b2"] * dat$hd + post[i , "b3"] * dat$ram + post[i , "b4"] * dat$ads + post[i , "b5"] * dat$premium)  
      loglik.Dbar[i] <- (-(n / 2) * (log(2 * pi))) - (n * log(sqrt(post[i , "S2"]))) - ((1 / (2 * post[i , "S2"])) * sum(res.Dbar ** 2))
    }
    Dbar <- 1 / nrow(post) * sum(-2 * loglik.Dbar)
    pD <- 2 * (Dbar - Dhat) 
    DIC <- Dhat + pD
    
    cat("Mean deviance:", Dhat, "\n")
    cat("Penalty:", pD, "\n")
    cat("Mean deviance:", DIC, "\n")
  }
  
}

(DIC.model1 <- DIC(post = post.m1, model = 1)) # DIC for first model.
(DIC.model2 <- DIC(post = post.m2, model = 2)) # DIC for second model. 

###################################################################################
################# Hypothesis evaluation with the Bayes factor #####################
###################################################################################

# Standardized MLE analysis. 
standard.output <- lm(scale(price, center = F) ~ scale(speed, center = F) + scale(hd, center = F) + scale(ram, center = F) 
                      + scale(ads, center = F) + scale(premium, center = F), data = dat)
# Store standardized coefficients and . 
standard.coef <- as.vector(standard.output$coefficients)
# Specify fraction b. 
b <- (5 / nrow(dat))
# Store vector with mean regression coefficients. 
mu <- standard.coef[1:6]
# Store a positive-definite symmetric matrix specifying the covariance matrix of the regression coefficients.
Sigma.raw <- vcov(standard.output)
# Remove intercept from covariance matrix. 
Sigma <- Sigma.raw[1:6, 1:6]
# If necessary, install package "MASS" for simulating from multivariate normal distribution.  
if(!require(MASS)) install.packages("MASS")
# Require package "gridExtra". 
library("MASS")
# Simulate a posterior from a multivariate normal distribution. 
set.seed(3791)
mvr.post <- mvrnorm(n = 10000, mu = mu, Sigma = Sigma) 
colnames(mvr.post) <- c("b0", "b1", "b2", "b3", "b4", "b5") 
# Center prior distribution on boundary of hypothesis. 
prior.mu <- c(rep(0, 6))
# Base prior variance on minimal training sample.  
prior.Sigma <- Sigma / b 
# Simulate a prior from a multivariate normal distribution. 
set.seed(3791)
mvr.prior <- mvrnorm(n = 10000, mu = prior.mu, Sigma = prior.Sigma)
colnames(mvr.prior) <- c("b0", "b1", "b2", "b3", "b4", "b5") 

## Hypothesis 1: The effect of each coefficient is positive.  
# Determine the fit based on the posterior for which b1 > 0, b2 > 0, b3 > 0, b4 > 0, and b5 > 0.
f_i1 <- mean(mvr.post[, "b1"] > 0 & mvr.post[, "b2"] > 0 & mvr.post[, "b3"] > 0 & mvr.post[, "b4"] > 0 & mvr.post[, "b5"] > 0)
# Determine the complexity proportion in the prior for which b1 > 0, b2 > 0, b3 > 0, b4 > 0, and b5 > 0.
c_i1 <- mean(mvr.prior[, "b1"] > 0 & mvr.prior[, "b2"] > 0 & mvr.prior[, "b3"] > 0 & mvr.prior[, "b4"] > 0 & mvr.prior[, "b5"] > 0) 
# Computing the Bayes factor.
(bf1u <- f_i1 / c_i1)

## Hypothesis 2: The standardized effect of each of the coefficients is larger than 0.1. 
# Determine the fit based on the posterior for which b1, b2, b3 > b4, b5.
f_i2 <- mean(mvr.post[, "b1"] > 0.1 & abs(mvr.post[, "b2"]) > 0.1 & mvr.post[, "b3"] > 0.1 & mvr.post[, "b4"] > 0.1 & abs(mvr.post[, "b5"]) > 0.1)
# Determine the complexity proportion in the prior for which b1, b2, b3 > b4, b5.
c_i2 <- mean(mvr.prior[, "b1"] > 0.1 & abs(mvr.prior[, "b2"]) > 0.1 & mvr.prior[, "b3"] > 0.1 & mvr.prior[, "b4"] > 0.1 & abs(mvr.prior[, "b5"]) > 0.1) 
# Computing the Bayes factor.
(bf2u <- f_i2 / c_i2)

## Hypothesis 3: Objective price determinants have a stronger effect than normative price determinants.  
f_i3 <- mean(abs(mvr.post[, "b1"]) > abs(mvr.post[, "b4"]) & abs(mvr.post[, "b1"]) > abs(mvr.post[, "b5"]) & 
               abs(mvr.post[, "b2"]) > abs(mvr.post[, "b4"]) & abs(mvr.post[, "b2"]) > abs(mvr.post[, "b5"]) &
               abs(mvr.post[, "b3"]) > abs(mvr.post[, "b4"]) & abs(mvr.post[, "b3"]) > abs(mvr.post[, "b5"]))
# Determine the complexity proportion in the prior for which b1, b2, b3 > b4, b5.
c_i3 <- mean(abs(mvr.prior[, "b1"]) > abs(mvr.prior[, "b4"]) & abs(mvr.prior[, "b1"]) > abs(mvr.prior[, "b5"]) & 
               abs(mvr.prior[, "b2"]) > abs(mvr.prior[, "b4"]) & abs(mvr.prior[, "b2"]) > abs(mvr.prior[, "b5"]) &
               abs(mvr.prior[, "b3"]) > abs(mvr.prior[, "b4"]) & abs(mvr.prior[, "b3"]) > abs(mvr.prior[, "b5"])) 
# Computing the Bayes factor.
(bf3u <- f_i3 / c_i3)

###########################################################################################
###################### Parameter Estimates & Credible Intervals ###########################
###########################################################################################

# Request means of posterior of each parameter for second model. 
(round(param.ests <- apply(post.m2, 2, mean), 2))
# Request 95% credible interval of posterior of each parameter for second model. 
round((param.cred <- sapply(post.m2, quantile, probs = c(0.025, 0.975))), 2)

post.hist <- function(post, param, main){ # Function for plotting parameter posteriors.
  ggplot(data = post, mapping = aes(x = post[ , param])) + # Input posterior data, with iterative parameters on x-axis.
    geom_histogram(colour = "black", fill = "white", bins = 30) + # Define histogram. 
    geom_vline(xintercept = mean(post[ , param]), linetype = "dashed") + # Insert mean as vertical line, indicated by dashed line. 
    geom_vline(xintercept = quantile(post[,param], probs = 0.025), linetype = "twodash") + # Insert 0.025 quantile, indicated by twodashed line. 
    geom_vline(xintercept = quantile(post[,param], probs = 0.975), linetype = "twodash") + # Insert 0.975 quantile, indicated by twodashed line. 
    scale_x_continuous(breaks = c(as.vector(quantile(post[ , param], probs = 0.025)), mean(post[ , param]), 
                                  as.vector(quantile(post[ , param], probs = 0.975)))) + # Add values associated with vertical lines. 
    
    xlab("") + # Remove x-axis. 
    ylab("") + # Remove y-axis. 
    ggtitle("", main) # Add main title. 
}

# If necessary, install package "gridExtra" for combining catterplots.  
if(!require(gridExtra)) install.packages("gridExtra")
# Require package "gridExtra". 
library("gridExtra")
# Plotting posterior distribution parameters second model 
grid.arrange(arrangeGrob(post.hist(post = post.m2, param = "b0", main = "B0 - Intercept"), post.hist(post = post.m2, param =  "b1", main = "B1 - Speed"), post.hist(post = post.m2, param = "b2", main = "B2 - Hd"),
                         post.hist(post = post.m2, param = "b3", main = "B3 - Ram"), post.hist(post = post.m2, param  = "b4", main = "B4 - Ads"), post.hist(post = post.m2, param  = "b5", main = "b5 - Premium"),
                         post.hist(post = post.m2, param = "S2", main = "S2 - Residual Variance"),
                         ncol = 2, nrow = 4))


