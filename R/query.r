# Returns column metadata from target that matches experimental conditions of query.
# Assumes that the provided meta_col data frames has a $exp_id -- a experimental condition id
sampleMatch = function(meta_col_query, meta_col_target) {
	if (!("exp_id" %in% colnames(meta_col_query))) {
		stop("Query meta data does not have an exp_id column")
	}

	if (!("exp_id" %in% colnames(meta_col_target))) {
		stop("Target meta data does not have an exp_id column")
	}

	query_exp = unique(meta_col_query$exp_id)  # query experiments

	target_match = which(meta_col_target$exp_id %in% query_exp)  #inclusion criteria

	return(meta_col_target[target_match,])
}