# Connects to the LINCS Clous API, returns JSON-like R object: data frames of comma separated data
lincsApi = function(query, table, max_requests=10) {
	require(stringr)
	require(RJSONIO)

	url = paste0('http://api.lincscloud.org/a2/', table, '?q=', query, '&user_key=b740ab796a4adf3f0f15eb7b9ae32d2b')
	# url = paste0('http://api.lincscloud.org/a2/', table, '?q=', query, '&user_key=lincsdemo')
	url = str_replace_all(url, '\"', "%22")
	print(url)

	count = 0
	result = NULL

	while(is.null(result) && count < max_requests) {
		count = count + 1

		try({
			result = fromJSON(url)
		})
	}

	if (is.null(result)) {
		warning("LINCS request not found, query:", query)
	}

	return(result)
}


# Retrieves perturbation information based on provided pert_id. Uses the LINCS REST API
lincsPertId = function(id) {
	query = paste0('{"pert_id":"', id, '"}')
	return(lincsApi(query, "pertinfo"))
}

lincsGeneId = function(id) {
	query = paste0('{"pr_id":"', id, '"}')
	return(lincsApi(query, "geneinfo"))
}
