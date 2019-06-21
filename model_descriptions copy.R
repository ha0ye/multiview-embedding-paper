## model 1: 3sp coupled logistic map (lotka-volterra)
make_model_1 <- function()
{
    run_model_1 <- function(r, alpha, run_time = 1000)
    {
        n <- length(r)
        data <- matrix(0, nrow = run_time, ncol = n)
        data[1,] <- rep.int(0.2, times = n)
        
        for(i in 2:run_time)
        {
            data[i,] <- r * data[i-1,] * (1 - alpha %*% data[i-1,])
        }
        
        return(data)
    }
    
    alpha <- matrix(c(1, 0.2, 0.2, 
                      0.2, 1, -0.2, 
                      0.2, -0.2, 1), nrow = 3, byrow = TRUE)
    r <- c(3.6, 3, 3)

    # generate data
    data <- run_model_1(r, alpha, run_time = 3000)
    save(data, file = "model_1.Rdata")
    return()
}

## model 2: Hastings-Powell 3sp model (process error not done properly, so set to 0)
make_model_2 <- function()
{
    run_model_2 <- function(a1, b1, a2, b2, d1, d2, run_time = 1000)
    {
        f <- function(v)
        {
            x <- v[1]
            y <- v[2]
            z <- v[3]
            f1x_y <- a1*x/(1+b1*x) * y
            f2y_z <- a2*y/(1+b2*y) * z
            return(c(x*(1-x) - f1x_y, 
                     f1x_y - f2y_z - d1*y, 
                     f2y_z - d2*z))
        }
        
        rk_step <- function(v, dt)
        {
            k1 <- dt * f(v)
            k2 <- dt * f(v + k1/2)
            k3 <- dt * f(v + k2/2)
            k4 <- dt * f(v + k3)
            return(v + (k1 + 2*k2 + 2*k3 + k4)/6)
        }
        
        n <- 3
        data <- matrix(0, nrow = run_time, ncol = n)
        data[1,] <- c(0.8, 0.2, 8)
        
        sample_freq <- 800
        dt <- 0.01
        
        for(i in 2:run_time)
        {
            v <- data[i-1,]
            for(k in 1:sample_freq)
            {
                v <- rk_step(v, dt)
            }
            data[i,] <- v
        }
        return(data)
    }
    
    a1 <- 2.5
    b1 <- 3.2
    a2 <- 0.1
    b2 <- 2.0
    d1 <- 0.2
    d2 <- 0.015
    
    # generate data
    data <- run_model_2(a1, b1, a2, b2, d1, d2, run_time = 3000)
    save(data, file = "model_2.Rdata")
    return()
}

## model 3: Dennis LPA beetle model
make_model_3 <- function()
{
    run_model_3 <- function(params, run_time = 1000)
    {
        b <- params$b
        c_el <- params$c_el
        c_ea <- params$c_ea
        c_pa <- params$c_pa
        u_a <- params$u_a
        u_l <- params$u_l
        
        data <- matrix(0, nrow = run_time, ncol = 3)
        data[1,] <- c(250, 5, 100) # initial conditions from Dennis 2001 (L, P, A)
        
        for(i in 2:run_time)
        {
            L <- data[i-1,1]
            P <- data[i-1,2]
            A <- data[i-1,3]
            
            data[i,1] <- b * A * exp(-c_el * L - c_ea * A)
            data[i,2] <- L * (1 - u_l)
            data[i,3] <- P * exp(-c_pa * A) + A * (1 - u_a)
        }
        
        return(data)
    }
    
    # maximum likelihood (ML) estimates with manipulated u_a and c_pa for chaotic behavior
    params <- data.frame(b = 10.67, u_l = 0.1955, u_a = 0.96, 
                         c_el = 0.01647, c_ea = 0.01313, c_pa = 0.35)

    # generate data
    data <- run_model_3(params, run_time = 3000)
    save(data, file = "model_3.Rdata")
    return()
}

## Huisman-Weissing resource competition model
make_HW_model <- function()
{
    run_HW_model <- function(run_time = 1000)
    {
        # step function
        step <- function(x, dt)
        {
            f <- function(x)
            {
                N <- head(x, num_species)
                R <- tail(x, num_resources)

                mu <- r * apply(R / (K + R), 2, min)
                N_prime <- N * (mu - m) # vector calculation (each of N, mu, m is 1:num_species)
                R_prime <- D * (S - R) - C %*% (mu * N)
                return(c(N_prime, R_prime))
            }
            
            k1 <- f(x)
            k2 <- f(x + dt/2 * k1)
            k3 <- f(x + dt/2 * k2)
            k4 <- f(x + dt * k3)
            return(x + dt / 6 * (k1 + 2*k2 + 2*k3 + k4))
        }
        
        # initial values & coefficients (Huisman & Weissing 2001, Am. Nat. Fig. 1)
        S <- c(10, 10, 10)
        N0 <- c(0.1, 0.1, 0.1, 0.1, 0.1)
        num_resources <- length(S)
        num_species <- length(N0)
        R <- matrix(S, byrow = TRUE, nrow = run_time, ncol = num_resources)
        N <- matrix(N0, byrow = TRUE, nrow = run_time, ncol = num_species)
        K <- matrix(c(0.20, 0.05, 0.50, 0.05, 0.50,
                      0.15, 0.06, 0.05, 0.50, 0.30,
                      0.15, 0.50, 0.30, 0.06, 0.05), nrow = num_resources, byrow = TRUE)
        C <- matrix(c(0.20, 0.10, 0.10, 0.10, 0.10,
                      0.10, 0.20, 0.10, 0.10, 0.20,
                      0.10, 0.10, 0.20, 0.20, 0.10), nrow = num_resources, byrow = TRUE)

        r <- rep.int(1, num_species) # maximum growth rate
        m <- rep.int(0.25, num_species) # mortality
        D <- 0.25
        dt <- 0.01
        sample_freq <- 500
        
        # loop
        for(t in 1:(run_time-1))
        {
            x <- c(N = N[t,], R = R[t,])
            for(k in 1:sample_freq)
            {
                x <- step(x, dt)
            }
            N[t+1,] <- head(x, num_species)
            R[t+1,] <- tail(x, num_resources)
        }

        temp <- data.frame(cbind(N, R))
        names(temp) <- c(paste("N", 1:NCOL(N), sep = ""), 
                         paste("R", 1:NCOL(R), sep = ""))
        return(temp)
    }
    
    data <- run_HW_model(run_time = 3000)
    data <- data[, grep("^N[0-9]+", names(data))]
    save(data, file = "model_hw.Rdata")
}

## helper functions
# add additive observation error (normal)
add_observation_error <- function(data, obs_err_sd = 0.1)
{
    return(data + rnorm(length(data), sd = obs_err_sd))
}

normalize <- function(data)
{
    column_means <- colMeans(data, na.rm = TRUE)
    column_sds <- apply(data, 2, function(x) {sd(x, na.rm = TRUE)})
    num_rows <- NROW(data)
    return((data - matrix(rep(column_means, each = num_rows), nrow = num_rows)) / 
               matrix(rep(column_sds, each = num_rows), nrow = num_rows))
}