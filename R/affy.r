
# affy_list: list of Affy expression sets
# returns expression matrix with intersecting
combineAffyExprs = function(affy_list) {
	# require(aff	y)
	# Combine the platform data
	# Get probe names for each platform
	probe_names = lapply(affy_list, function(d) {return(rownames(exprs(d)))})

	# Cross-platform probeset. The genes in common accross platform.
	common_probe_names = Reduce(intersect, probe_names)  # Reduce sequentially applies the intersect operation

	# combine normalized expression matrices into single matrix by using only intersecting probe names.
	# do.call applies cbind operation to list returned by lapply
	combined_exprs = do.call(cbind, 
		lapply(affy_list, function(d) {
				exprs(d)[common_probe_names,]  # slices expression matrix based on common_probe_names
			}
		)
	)

	return(combined_exprs)
}