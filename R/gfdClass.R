####################### Define a Set of Multivariate Grid Functional Data objects ####################### 


mvgfd <- R6::R6Class("mvgfd",
                     
                     ############################# Public list #############################
                     
                     #' @description
                     #' Constructor for `mvgfd` objects (Multivariate Grid Functional Data objects)
                     public = list(
                       initialize = function(argval = NULL, mvgfd_obj, centerfns = TRUE, num_pcs = 1) {
                         
                         #' @param argval A list of numeric vectors of argument values at which the `mfd` object is to be evaluated
                         private$.argval <- if (is.list(argval)) argval else NULL
                         
                         
                         #' @param centerfns logical. If TRUE  the input data undergoes centralization or demeaning prior to being processed by the function.
                         #' Demean the data (centerfns)
                         private$.centerfns <- centerfns
                         
                         
                         #' @param mvfd_obj  List of matrices or arrays (Multivariate Grid Functional Data )
                         # mvgfd object: either list of matrices or arrays 
                         private$.mvfd_obj <- mvgfd_obj
                         
                         
                         #' @param num_pcs  The number of PCs. The default is one (The first principal component only). But the user is able to see higher PCs
                         private$.num_pcs <- num_pcs
                         
                         
                       },
                       
                       
                       ########### Check for a valid data 
                       check_data = function() {
                         
                         # The functional data should be a list of matrices or arrays or just one matrix
                         if (!(is.list(private$.mvfd_obj) || is.matrix(private$.mvfd_obj))) {
                           stop("Input data must be a list of matrices and arrays or just one matrix!")
                         }
                         
                         # Check if the list contains only matrices or arrays
                         if (is.list(private$.mvfd_obj)){
                           valid_entries <- all(sapply(private$.mvfd_obj, function(x) is.matrix(x) || is.array(x)))
                           if (!valid_entries) {
                             stop("The list contains entries other than matrices or arrays.")
                           }
                           
                         }
                         

                         # Check if the list contains arrays
                         #Also, checks for univariate case! 
                         if(is.matrix(private$.mvfd_obj)){all_matrices <- 'Univariate case!'
                         } else {all_matrices <- all(sapply(private$.mvfd_obj, is.matrix))}
                         
                         
                         
                         # Check for NA or non-numeric entries in matrices or arrays
                         contains_na_non_numeric <- any(sapply(private$.mvfd_obj, function(x) any(is.na(x) | !is.numeric(x))))
                         
                         # Display report on the data 
                         cat("Report on the data:\n")
                         cat("All matrices: ", all_matrices, "\n")
                         cat("Contains NA or non-numeric entries: ", contains_na_non_numeric, "\n")
                         
                         return(list(all_matrices = all_matrices, contains_na_non_numeric = contains_na_non_numeric))
                       }
                       
                     ),
                     
                     ############################# Private list #############################
                     
                     private = list(
                       .argval = NULL,
                       .mvfd_obj = NULL, 
                       .centerfns = NULL,
                       .num_pcs = NULL, 
                       all_matrices = NULL, 
                       contains_na_non_numeric = NULL
                     ),
                     
                     ############################# Active list #############################
                     active = list(
                       
                       argval = function(value) {
                         if (missing(value)) {
                           private$.argval
                         } else {
                           stop("`$argval` is read only", call. = FALSE)}
                       }, 
                       
                       mvfd_obj = function(value) {
                         if (missing(value)) {
                           private$.mvfd_obj
                         } else {
                           stop("`$mvfd_obj` is read only", call. = FALSE)}
                       }, 
                       
                       
                       centerfns = function(value) {
                         if (missing(value)) {
                           private$.centerfns
                         } else {
                           stop("`$centerfns` is read only", call. = FALSE)}
                       }, 
                       
                       num_pcs = function(value) {
                         if (missing(value)) {
                           private$.num_pcs
                         } else {
                           stop("`$centerfns` is read only", call. = FALSE)}
                       }

                     )
)


############### TEST ###############

library(fda)
df  = gait
datagait = list(t(df[,,1]),t(df[,,2]))
data_checker_new = mvgfd$new(mvgfd_obj = datagait,centerfns =FALSE )
# Check the data
#result2 <- data_checker_new$check_data()

