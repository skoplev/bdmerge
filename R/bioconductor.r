# Converts
# Requires Bioconductor base
eset2dlist = function(eset) {
	out = list()

	out$meta_col = pData(eset)  # Get sample meta data
	out$meta_col$platform = annotation(eset)  # by convention
	out$meta_row = data.frame(id=featureNames(eset))
	out$data = exprs(eset)

	return(out)
}