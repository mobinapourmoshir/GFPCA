library(R6)
library(Matrix)

# Define the MFPC class
MFPCClass <- R6Class("MFPC",
                     public = list(
                       mvfd_obj = NULL,
                       centerfns = TRUE,
                       initialize = function(mvfd_obj, centerfns = TRUE) {
                         # Constructor method
                         self$mvfd_obj <- mvfd_obj
                         self$centerfns <- centerfns
                       },
                       
                       MFPC = function() {
                         
                         ############################### get.pen function ###############################
                         get.pen <- function(td, alpha=0) {
                           m = length(td);
                           h = td[2:m] - td[1:(m-1)]; 
                           Q = matrix(0, m, m-1);
                           R = matrix(0, m-1, m-1);
                           for(k in 2:(m-1))
                           {
                             Q[k-1,k] = 1/h[k-1];
                             Q[k,k] = -1/h[k-1] - 1/h[k];     
                             Q[k+1,k] = 1/h[k]
                           }
                           for(j in 2:(m-2))
                           {
                             R[j,j] = 1/3 * (h[j-1] + h[j]);
                             R[j,j+1] = 1/6 * h[j];
                             R[j+1,j] = 1/6 * h[j]
                           }
                           R[m-1,m-1] = 1/3 * (h[m-2] + h[m-1]);
                           s <- solve(R[2:(m-1), 2:(m-1)]) %*% t(Q[1:m, 2:(m-1)]);
                           OMEGA = Q[1:m, 2:(m-1)] %*% s;
                           EIG.O <- eigen(OMEGA); GAMMA=EIG.O$vectors; LAMBDA=diag(EIG.O$values);
                           S.alpha <- GAMMA%*%diag((1/(1+alpha*diag(LAMBDA))))%*%t(GAMMA);
                           return(list(OMEGA=OMEGA, S.alpha=S.alpha))
                         }
                         
                         ############################### optimal number of mfrow for par charts ############################### 
                         
                         get_optimal_factors <- function(n) {
                           factors <- NULL
                           
                           for (i in 1:sqrt(n)) {
                             if (n %% i == 0) {
                               factors <- rbind(factors, c(i, n / i))
                             }
                           }
                           
                           # Find the pair of factors closest in value
                           index <- which.min(abs(factors[,1] - factors[,2]))
                           return(as.integer(factors[index, ]))
                         }
                         
                         
                         ############################### norm of any vector ############################### 
                         norm_vec = function(x) sqrt(sum(x^2))
                         
                         ############################### Considering some values for alpha ############################### 
                         get.Alphas = function(n=101,a=2,s=-75) {return(a^seq(s,s+n))}
                         
                         
                         ############################### Optimal number of alpha using GCV ############################### 
                         opt_alpha = function(X, u, S_alpha) {
                           alphas = get.Alphas()
                           n = nrow(X)
                           m = ncol(X)
                           GCV = c()
                           
                           for (S in S_alpha) {
                             GCV_alpha = (1/m) * (norm_vec((diag(m)-S)%*%t(X)%*%u)^2/(norm_vec(u)^2 %*% (1-(1/m)*sum(diag(S)))^2))
                             GCV = c(GCV, GCV_alpha)}
                           
                           opt.alpha = alphas[which.min(GCV)]
                           opt_s.alpha = S_alpha[[which.min(GCV)]]
                           return(list(GCV = GCV, opt.alpha = opt.alpha, opt_s.alpha = opt_s.alpha))}
                         
                         
                         
                         
                         ############# Multivariate FPCA #############
                         n_var = length(self$mvfd_obj)
                         n = nrow(self$mvfd_obj[[1]])
                         n_cols = as.vector(as.data.frame(sapply(self$mvfd_obj, dim))[2,])
                         
                         ############# Pre-processing: Centralizing the data #############
                         if (self$centerfns == TRUE) {
                           mu = lapply(self$mvfd_obj, function(x){apply(x,2,mean)}) #mean on grid-points (columns)
                           self$mvfd_obj = lapply(self$mvfd_obj, function(x){scale(x , scale = FALSE)}) # centralized/ de-mean data 
                         } 
                         
                         
                         ############# Visualization: matplot for matplot_mvfd_obj #############
                         #dev.off()
                         par(mfrow=c(as.vector(get_optimal_factors(n_var))))
                         # Chart 1: matplot_mvfd_obj
                         matplot_mvfd_obj = for (i in 1:n_var) {
                           matplot(t(self$mvfd_obj[[i]]),type = 'l', main = i)} 
                         
                         
                         ############# initial S_alpha #############
                         
                         GridPoints = list()
                         for (i in 1:n_var) {
                           # We should be able to extract grid points from column names! 
                           
                           # If column names are characters:
                           if(any(is.na(as.numeric(colnames(self$mvfd_obj[[i]])))) ==TRUE){
                             cycle = seq(1:length(colnames(self$mvfd_obj[[i]])))}
                           # otherwise:  
                           else{cycle = as.numeric(colnames(self$mvfd_obj[[i]]))}
                           GridPoints[[i]] = cycle
                           #S_alpha [[i]] = get.pen(td = cycle, alpha = 1e-5)$S.alpha
                         }
                         
                         ############## S_alpha for all alphas 
                         alphas = get.Alphas()
                         S_alpha = list()
                         index = 0
                         
                         for (alpha in alphas ) {
                           index = index +1
                           S = list()
                           for (i in 1:n_var) {
                             S[[i]] = get.pen(td = GridPoints[[i]], alpha = alpha)$S.alpha}
                           
                           S_alpha[[index]] = as.matrix(bdiag(S)) 
                         }
                         
                         
                         ############# Power Algorithm #############  
                         
                         #############  Step1: initializing v 
                         v = list()
                         for (i in 1:n_var) {
                           v[[i]] = svd(self$mvfd_obj[[i]])$v[,1]}
                         
                         ###### data side by side:
                         comb_data = do.call(cbind, self$mvfd_obj)
                         
                         #####  Step2. Updating u and v (Repeat until convergence)
                         iter = 0 
                         u_prev = numeric(n)
                         comb_S_alpha = S_alpha[[60]] #initial S_alpha with initial alpha = 1.525879e-05
                         
                         while (1){
                           iter = iter+1
                           
                           #### v side by side 
                           Combv= t(matrix(c(unlist(v)),nrow = 1))
                           
                           
                           u = comb_data %*% Combv 
                           v_new = comb_S_alpha %*% t(comb_data) %*% u
                           v_new = v_new/norm_vec(v_new)
                           
                           
                           # optimal alpha
                           opt = opt_alpha(X = comb_data, u = u, S_alpha)
                           comb_S_alpha = opt$opt_s.alpha
                           alpha = opt$opt.alpha
                           
                           
                           if(norm_vec(u)-norm_vec(u_prev)<= 1e-15){break}
                           u_prev = u
                         }
                         v = split(v_new, rep(seq_along(n_cols), n_cols))
                         
                         # Chart 2: matplot_PC1
                         matplot_PC = for (i in 1:n_var) {
                           matplot(v[[i]],type = 'l', main = i)} 
                         
                         return(S_alpha)
                       }
                     )
)

# Instantiate an object of the MFPC class
library(fda)
df  = gait
data = list(t(df[,,1]),t(df[,,2]))
mfpc_object <- MFPCClass$new(mvfd_obj = data, centerfns = FALSE)

# Call the MFPC method on the object
mfpc_result <- mfpc_object$MFPC()
