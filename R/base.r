# Collection of commonly used functions.

# stemString() returns stem id of string. Decapitalizes and removes all non a-z and non 0-9.
# Ex: MC-7 -> mc7
stemString = function(str) {
	require(stringr)

	out = str_replace_all(str, "[^[:alnum:]]", "")  # removes all non[a-z.0-9]
	out = tolower(out)
	return(out)
}

# Calculates the indices of matching and ordered elements. Handles input of different length by truncated the
# longest element.
matchingIds = function(names1, names2) {
	common_length = min(length(names1), length(names2))
	ids = which(names1[1:common_length] == names2[1:common_length])
	return(ids)
}

# Counts number of occurances of a character in a string
count = function(str, letter) {
	if (!is.character(letter)) {
		stop("Letter must be a character")
	}
	if (nchar(letter) != 1) {
		stop("Letter must be single character")
	}
	enforced_string = as.character(str)
	string_vector = strsplit(enforced_string, NULL)[[1]]
	return(sum(string_vector == letter))
}

# Distance measure based on relative ratios. Tolerance value is added to zero inputs
relRatio = function(x, y, tol=0.00001) {

	x[x == 0] = tol
	y[y == 0] = tol
	# if (abs(x) < tol) {
	# 	x = x + tol
	# }
	# if (abs(y) < tol) {
	# 	y = y + tol
	# }
	ratio = x/y
	ratio[ratio < 1] = 1/ratio[ratio < 1]
	return(ratio)
}

# Generator of binning functions.
# Convert vector of numerical values to closest bin according to some distance measure.
# USAGE:
# conc2str = conc2strGen()  # generate binning function with default parameters
# conc2str = conc2strGen(my_bins, my_distance_fun)
# assigned_bin_strings = conc2str(values)
conc2strGen = function(bin_def=c(0.1, 1.0, 10.0, 100.0), distance=relRatio)
{
	Vectorize(
		function(num, bins=bin_def) {
			# Find closes bin numerical value
			bin_num = bins[which.min(distance(num, bins))]

			return(format(bin_num, nsmall=1, trim=TRUE))
		}, "num"
	)
}



# Calculates cosine similarity of the columns of provided matrix.
# If another matrix is provided, Y, then the cosine similarity between the columns of the two matrices are calculated.
# This is similar to the cor() {stats} interface.
cosineSimil = function(X, Y=NULL, na.value=0.0) {
	if (is.null(Y)) {
		Y = X
	}

	if (!is.null(na.value)) {
		Y[is.na(Y)] = na.value
		X[is.na(X)] = na.value
	}

	mat = matrix(NA, nrow=ncol(X), ncol=ncol(Y))
	for (i in 1:ncol(X)) {
		for (j in 1:ncol(Y)) {
			mat[i, j] = sum(X[,i] * Y[,j]) / (sqrt(sum(X[,i]^2)) * sqrt(sum(Y[,j]^2)))  # cosine similarity
		}
	}
	return(mat)
}


# Calculates similarity between signed adjacency matrices. Can be either symmetric or non-symmetric corresponding to
# a bipartite graph.
networkSimilarity = function(cor_list, name, corFun=cor, ...) {
	net_sim = matrix(NA, ncol=length(cor_list), nrow=length(cor_list))
	for (i in 1:length(cor_list)) {
		for (j in 1:length(cor_list)) {
			cmat1 = cor_list[[cell_lines[[i]]]][[name]]
			cmat2 = cor_list[[cell_lines[[j]]]][[name]]

			cmat1[is.na(cmat1)] = 0.0
			cmat2[is.na(cmat2)] = 0.0
			if (nrow(cmat1) == ncol(cmat1) & nrow(cmat2) == ncol(cmat2)) {
				# symmetrical only count once
				net_sim[i, j] = corFun(cmat1[lower.tri(cmat1)], cmat2[lower.tri(cmat2)], ...)
			} else {
				# assymetrical, all entries assumed to be unique
				net_sim[i, j] = corFun(as.vector(cmat1), as.vector(cmat2), ...)
			}
		}
	}
	return(net_sim)
}


# Plots
regPlot = function(x, y, ...) {
	linreg = lm(y~x)
	plot(x, y, cex=0.5, col=rgb(0, 0, 0, 0.5), xlim=c(-1.0, 1.0), ylim=c(-1.0, 1.0),
		main=paste0(
			"cor=", format(cor(x, y, use="pairwise.complete.obs"), digits=3, nsmall=3),
			" , slope=", format(linreg$coefficients[2], digits=3, nsmall=3)
			), 
			...)
	# plot(x, y, cex=0.5, col=rgb(0, 0, 0, 0.5), ...)
	abline(linreg, col="red")
}

