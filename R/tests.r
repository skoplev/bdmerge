
# Creates a random dlist for testing purposes.
# n samples, m features, n_exp_ids experimental ids labeled "id_#"
makeRandDlist = function(n, m, n_exp_ids) {
	dlist = list()
	dlist$data = matrix(rnorm(n*m), ncol=n, nrow=m)
	dlist$meta_col = data.frame(id=paste0("id_", sample(x=n_exp_ids, size=n, replace=TRUE)))
	dlist$meta_row = data.frame(probe=paste0("pr_"), 1:m)

	return(dlist)
}


# A collection of test functions. Returns true if test passed and false if not.
testMergeDataListsByCol = function() {

	# Testing parameters
	n1 = 5  # number of samples of data set 1
	m1 = 10  # number of features of data set 2

	n2 = 30
	m2 = 20

	n_exp_ids = 10  # number of string based experimental ids

	# init two random data lists with matching experimental ids
	dlists = list()
	dlists[[1]] = makeRandDlist(n1, m1, n_exp_ids)
	dlists[[2]] = makeRandDlist(n2, m2, n_exp_ids)

	checkDataList(dlists[[1]])
	checkDataList(dlists[[2]])

	unique(dlists[[1]]$meta_col)
	unique(dlists[[2]]$meta_col)

	dcomb = mergeDataListsByCol(
		dlists=dlists, 
		match_conditions=list("id", "id"),
		matrix_selections=list("data", "data"),
		method="cyclical"
		)

	dlists[[1]]$meta_col
	dlists[[2]]$meta_col

	dcomb[[1]]$meta_col
	dcomb[[2]]$meta_col

	dcomb[[1]]$data
	dcomb[[2]]$data
}

