# R tools for manipulating a multidimensional linear coordinate system by means of elementary rotations of its basis vectors

mgRotateA <- function(L, alpha)
{
    return(cos(alpha) * L[,1] - sin(alpha) * L[,2])
}

mgRotateB <- function(L, alpha)
{
    return(sin(alpha) * L[,1] + cos(alpha) * L[,2])
}

mgRotate <- function(L, alpha)
{
    return(cbind(mgRotateA(L, alpha), mgRotateB(L, alpha)))
}

# Rotational optimization of an orthonormal basis vector system (columns of L) using elementary rotations.
# This is a very simple trust-region algorithm with fixed stepwidth $alpha, doing a fixed number of $iterations.
# Rotational directions bringing L[$zero, $j_target] to zero best are chosen first.
# $k is the dimensionality of the descent: the $k best-scoring directions are applied in each iteration.
# Set $visualization to show a line plot after each iteration. 
# Returns the final basis vectors L and the rotation matrix R

# Example:
# 
# # inital factorization
# P <- prcomp(X, center = FALSE)
# 
# # do the first optimization using the first 15 principal components and give graphical feedback
# O <- mgElemRotOpt(P$rotation[,1:15], c(1:90, 1095:1650), pi/1000, 200, 1, visualization=TRUE, k=3)
#
# # further refine the last result, testing different parameters for $alpha and $k
# O = mgElemRotOpt(O$L[,1:15], c(1:90, 1095:1650), pi/2000, 600, 1, visualization=TRUE, k=1)
#
# # get the coordinates from original data by projecting it onto the new basis
# S = X %*% O$L
#
# # compare reconstruction of observation 1 to original data
# plot(X[1,], type = "l")
# lines((S %*% t(O$L))[1,], type = "l", col="red")

mgElemRotOpt <- function(L, i_zeroes, alpha, iterations = 1, j_target = 1, k = 1, visualization = FALSE)
{
    n_comp = ncol(L)
    j_test = (1:n_comp)[-j_target]

    # remember the rotations
    R = diag(n_comp)

    for(it in 1:iterations)
    {
        P = NULL

        for(i in 1:length(j_test)) for(sgn in c(-1,1))
        {   
            # score of target criterium in current state
            p_current = sum(L[i_zeroes, j_target]^2)

            # pair to test
            j_rot = c(j_target, j_test[i])

            p_test = sum(mgRotateA(L[i_zeroes, j_rot], sgn * alpha)^2)

            # remember performance and the pair and sign
            P <- rbind(P, c(p_current - p_test, j_rot, sgn))
        }

        o = order(P[,1], decreasing = TRUE)

        for(i in 1:k)
        {
            j_rot = P[o[i],2:3]
            sgn = P[o[i],4]

            L[,j_rot] = mgRotate(L[,j_rot], sgn * alpha)
            R[,j_rot] = mgRotate(R[,j_rot], sgn * alpha)
        }

        if(visualization)
            plot(L[,j_target], type = "l")

    } 
 
  return(list(L = L, R = R))
}
