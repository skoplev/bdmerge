
# Arguments:
#	d: list of data and metadata in $meta_row and $meta_col
#	sample_def: column names of d$meta_col that defines the experimental id's

# Assumes that the metadata for the first entry of each sample is representative for the others.
# TODO: clean up implementation. data.table functionality can be used for evaluating the functions
# more efficiently...
sampleApply = function(d, name="data", sample_def, fun, ...) {
	# Check user input
	checkDataList(d)
	if (!all(sample_def %in% colnames(d$meta_col))) {
		stop("Invalid sample definition: ", paste(sample_def, collapse=", "))
	}

	# construct experimental ids
	out = list()
	out$meta_row = d$meta_row  # the same feature metadata

	if (class(d$meta_col)[1] == "data.table") {
		# data.table
		if (length(sample_def) > 1) {
			exp_ids = apply(d$meta_col[, sample_def, with=FALSE], 1, paste, collapse=":")
		} else {
			exp_ids = d$meta_col[,sample_def, with=FALSE]
		}
	} else {
		# data.frame or matrix
		if (length(sample_def) > 1) {
			exp_ids = apply(d$meta_col[, sample_def], 1, paste, collapse=":")
		} else {
			exp_ids = d$meta_col[,sample_def]
		}
	}

	unique_exp = unique(exp_ids)
	# print(unique_exp)

	# Allocate sample data matrix
	out$meta_col = d$meta_col[match(unique_exp, exp_ids),]

	# Allocate transformed matrices
	out[[name]] = matrix(NA, nrow=nrow(d[[name]]), ncol=length(unique_exp))

	data_is_table = class(d[[name]])[1] == "data.table"

	# Loop over each condition
	for (i in 1:length(unique_exp)) {
		# Get matrix of data for 

		if (data_is_table) {
			sample_data = d[[name]][, exp_ids == unique_exp[i], with=FALSE]
			sample_data = as.matrix(sample_data)  # for single measurements
		} else {
			sample_data = d[[name]][, exp_ids == unique_exp[i]]
			sample_data = as.matrix(sample_data)
		}

		# print(dim(sample_data))
		# print(unique_exp[i])

		if (is.null(dim(sample_data))) {
			warning("Sample data not found for experimental id: ", unique_exp[i])
			next
		}

		sample_transform = apply(sample_data, 1, fun, ...)

		out[[name]][,i] = sample_transform
	}

	return(out)
}

# Transform data list table with some function
# additional arguments are the arguments to the transformation function.
transformDlistCol = function(d, table, W) {
	d$meta_col
	d$meta_row
	d[[table]]

	out = list()
	out$meta_col = d$meta_col
	out$meta_row = data.frame(id=rownames(W))

	out[[table]] = W %*% abs(d[[table]])
	return(out)
}

# Performs ttest on data structure with d$data, d$meta_col, d$meta_row
# $meta_col must specify fields:
#	exp_id, an string identifying the experimental id. e.g. "MC60:erlotinib"
#	neg_ctrl, a boolean modifier of exp_id specifying whether the sample is a negative control or a primary experiment.

# ctrl_code:
#	plate:  find on plate by looking up meta_col$pert_id="DMSO" on meta_col$det_plate
#	neg_ctrl:  1-to-1 correspondance of exp_ids with neg_ctrl boolean value.

# Supports both lists of data.frames, data.tables, or matrices
ttestExpCtrl = function(d, ctrl_code) {
	checkDataList(d)

	# Output data structures.
	exp_unique = unique(d$meta_col$exp_id)
	exp_unique = exp_unique[exp_unique != ""]

	out = list()
	# Copy metadata
	out$meta_row = d$meta_row
	out$meta_col = d$meta_col[match(exp_unique, d$meta_col$exp_id),]  # specific to the experimental ids

	# Initialize output matrices
	out$mean_diff = matrix(NA, nrow=nrow(d$data), ncol=length(exp_unique))  # difference between sample and control
	out$mean_sample = matrix(NA, nrow=nrow(d$data), ncol=length(exp_unique))  # mean sample (postive experiment)
	out$mean_ctrl = matrix(NA, nrow=nrow(d$data), ncol=length(exp_unique))  # mean control (negative experiment)
	# out$mean_sample = matrix(NA, nrow=nrow(d$data), ncol=length(exp_unique))  # sample mean
	out$tstat = matrix(NA, nrow=nrow(d$data), ncol=length(exp_unique))
	out$pval = matrix(NA, nrow=nrow(d$data), ncol=length(exp_unique))

	# data.table indexing
	# setkey(d$meta_col, "exp_id")  # for binary searches

	# if ((class(d$meta_col))[1] == "data.table") {
	# 	setkey(d$meta_col, "exp_id", "neg_ctrl")  # two column binary search
	# }

	data_is_table = (class(d$data) == "data.table")[1]  # main data.table?

	for (i in 1:length(exp_unique)) {
		if (i %% 10 == 0) {
			cat("ttest experiment ", i, "out of ", length(exp_unique), "\n")
		}

		# Get sample ids of samples and associated negative controls.
		# Two different encoding methods for finding negative controls.
		# "neg_ctrl" uses a boolean column along with exp_id.
		# "plate" uses pert_id="DMSO" along with det_plate for the experiments.
		# Note that the exp_id field is used differently either refering to the experimental condition (same id for repeats)
		# or the experimental coordinates (different ids for repeats).
		if (ctrl_code == "exp_type") {
			# Get sample and control data of the experiment. Linear search. Can be extended to binary search via data.table
			sample_cols = which(d$meta_col$exp_id == exp_unique[i] & d$meta_col$exp_type == "sample")  # data.frame
			# sample_cols = d$meta_col[list(exp_unique[i], FALSE), which=TRUE]  # data.table

			ctrl_cols = which(d$meta_col$exp_id == exp_unique[i] & d$meta_col$exp_type == "control")  # data.fame
			# ctrl_cols = d$meta_col[list(exp_unique[i], TRUE), which=TRUE]
		} else if (ctrl_code == "plate") {
			sample_cols = which(d$meta_col$exp_id == exp_unique[i])

			# ctrl_cols 
			plate = unique(d$meta_col$det_plate[sample_cols])
			if (length(plate) > 1) {
				warning("Multiple plates included in experimental id definition.")
			}

			ctrl_cols = which((d$meta_col$det_plate %in% plate) & d$meta_col$pert_id == "DMSO")
		}

		if (length(sample_cols) == 0) {
			warning("Undefiend sample for experiment: ", exp_unique[i])
		}

		if (length(ctrl_cols) == 0) {
			warning("Undefined negative control for experiment: ", exp_unique[i])
		}

		if (data_is_table) {
			# data.table has minor syntax differences
			sample_data = d$data[,sample_cols, with=FALSE]
			control_data = d$data[,ctrl_cols, with=FALSE]
		} else {
			# assume data.frame or matrix 
			sample_data = d$data[,sample_cols]
			sample_data = as.matrix(sample_data)  # force matrix format in case there is only one observation

			control_data = d$data[,ctrl_cols]
			control_data = as.matrix(control_data)
		}

		# Vectorized two-sample t-tests
		sample_mean = apply(sample_data, 1, mean, na.rm=TRUE)
		ctrl_mean = apply(control_data, 1, mean, na.rm=TRUE)

		sample_var = apply(sample_data, 1, var, na.rm=TRUE)
		ctrl_var = apply(control_data, 1, var, na.rm=TRUE)

		# Calculate t statistic
		tstat = (sample_mean - ctrl_mean) / sqrt(sample_var/ncol(sample_data) + ctrl_var/ncol(control_data))

		# Calculate the degrees of freedom:
		# Welch-Satterthwaite approximation of the effective degrees of freedom
		deg_free = (sample_var/ncol(sample_data) + ctrl_var/ncol(control_data))^2 / ((sample_var/ncol(sample_data))^2/(ncol(sample_data) - 1) + (ctrl_var/ncol(control_data))^2/(ncol(control_data) - 1))

		# Calculate pvalues
		pval = 2 * pt(abs(tstat), deg_free, lower=FALSE)  # assume symmetric density

		# Mean difference
		mean_diff = sample_mean - ctrl_mean

		# Save statistics used in t-test
		out$mean_diff[,i] = mean_diff
		out$mean_sample[,i] = sample_mean
		out$mean_ctrl[,i] = ctrl_mean
		out$tstat[,i] = tstat
		out$pval[,i] = pval
	}
	return(out)
}

# Note that the metadata is based on the first matching experimental conditions. Hence,
# if $exp_id contains multiple metadata entries only the first carries over to the output.
chardirExpCtrl = function(d, ctrl_code) {

	# Test environemnt
	# ctrl_code="plate"

	# require(GeoDE)
	checkDataList(d)  # test if valid, throws errors if not

	exp_unique = unique(d$meta_col$exp_id)
	exp_unique = exp_unique[exp_unique != ""]

	out = list()
	# Copy metadata
	out$meta_row = d$meta_row
	out$meta_col = d$meta_col[match(exp_unique, d$meta_col$exp_id),]  # specific to the experimental ids

	# Initialize output matrices
	out$data = matrix(NA, nrow=nrow(d$data), ncol=length(exp_unique))  # difference between sample and control
	out$ctrl = matrix(NA, nrow=nrow(d$data), ncol=length(exp_unique))

	data_is_table = (class(d$data) == "data.table")[1]  # main data.table?

	for (i in 1:length(exp_unique)) {
		if (i %% 10 == 0) {
			cat("chardir experiment ", i, "out of ", length(exp_unique), "\n")
		}

		# Find data column indices for sample and control.
		if (ctrl_code == "exp_type") {
			# Get sample and control data of the experiment. Linear search. Can be extended to binary search via data.table
			sample_cols = which(d$meta_col$exp_id == exp_unique[i] & d$meta_col$exp_type == "sample")  # data.frame
			# sample_cols = d$meta_col[list(exp_unique[i], FALSE), which=TRUE]  # data.table

			ctrl_cols = which(d$meta_col$exp_id == exp_unique[i] & d$meta_col$exp_type == "control")  # data.fame
			# ctrl_cols = d$meta_col[list(exp_unique[i], TRUE), which=TRUE]
		} else if (ctrl_code == "plate") {
			sample_cols = which(d$meta_col$exp_id == exp_unique[i])

			# ctrl_cols 
			plate = unique(d$meta_col$det_plate[sample_cols])
			if (length(plate) > 1) {
				warning("Multiple plates included in experimental id definition.")
			}

			ctrl_cols = which((d$meta_col$det_plate %in% plate) & d$meta_col$pert_id == "DMSO")
		}

		# Duplicate observations if singletons
		if (length(sample_cols) == 1) {
			sample_cols = rep(sample_cols, 2)
		}

		if (length(ctrl_cols) == 1) {
			ctrl_cols = rep(ctrl_cols, 2)
		}

		# Get sample and control data of the experiment
		if (data_is_table) {
			# data.table
			sample_data = d$data[,sample_cols, with=FALSE]
			control_data = d$data[,ctrl_cols, with=FALSE]
		} else {
			# assume data.frame
			sample_data = d$data[,sample_cols]
			sample_data = as.matrix(sample_data)  # force matrix format in case there is only one observation

			control_data = d$data[,ctrl_cols]
			control_data = as.matrix(control_data)
		}

		# Store mean control
		out$ctrl[,i] = apply(control_data, 1, mean, na.rm=TRUE)

		# Make characteristic direction data structure
		# d$meta_row is converted to data frame 
		dchar = data.frame(names=as.data.frame(d$meta_row)[,1])

		dchar = cbind(dchar, as.data.frame(control_data), as.data.frame(sample_data))


		# Truncate missing data (rows)
		complete_rows = apply(!is.na(dchar), 1, all)
		dchar = dchar[complete_rows,]

		# Construct control sample code for the characteristic direction method
		sample_ctrl_class = as.factor(c(rep(1, ncol(control_data)), rep(2, ncol(sample_data))))

		# chdir = chdirAnalysis(dchar, sample_ctrl_class, 
		# 	2.0,   # penalty parameter
		# 	CalculateSig=TRUE, 
		# 	nnull=10)

		# Try charactaristic direction method, save if sucessfull.
		tryCatch(
			{
				chdir = chdirAnalysis(dchar, sample_ctrl_class)  # does not calculate significance levels
			}, error=function(e) {
				warning("chdir not calculated for ", i, "th condition. Possibly singleton control and sample. Ignored error message: ", e)
			}, finally={
				# Test if rownames aggree
				if (!all(as.data.frame(d$meta_row)[complete_rows,1] == rownames(chdir$chdirprops[[1]][[1]]))) {
					stop("Characteristic direction probe name mismatch")
				}

				# Store characteristic direction
				out$data[complete_rows,i] = chdir$chdirprops[[1]][[1]]
			}
		)

	}
	return(out)
}

# correlation merge.
# methodA is a user supplied sequence of methods used to generate the data in sourceA
# transformA is a function that reduces the data or transforms it using some method
# method=c("sample", "bootstrap", "average")
corMerge = function(
	sourceA, descA, transformA=NULL,
	sourceB, descB, transformB=NULL,
	col_selector_field=NULL,  # do a seperate analysis for each factor, returns list for each, including an "all" analysis
	corFun=stats::cor,  # correlation function
	match_def, match_method="sample") 
{
	data_nameA = strsplit(basename(sourceA), "\\.")[[1]][1]
	data_nameB = strsplit(basename(sourceB), "\\.")[[1]][1]

	# Read data
	dA = readData(dirname(sourceA), data_nameA)
	dB = readData(dirname(sourceB), data_nameB)

	# Rename, $data
	dA = list(data=dA[[data_nameA]], meta_col=dA$meta_col, meta_row=dA$meta_row)
	dB = list(data=dB[[data_nameB]], meta_col=dB$meta_col, meta_row=dB$meta_row)

	# Transform according to provided callback function
	if (!is.null(transformA)) {
		dA = transformA(dA)
	}

	if (!is.null(transformB)) {
		dB = transformB(dB)
	}

	# if (!is.null(col_selector_field)) {
	# 	sep_values = unique(c(dA$meta_col[[col_selector_field]], dB$meta_col[[col_selector_field]]))

	# 	for (i in 1:length(sep_values)) {

	# 	}
	# }

	if (match_method == "average") {
		# Apply mean to data list
		dA = sampleApply(d=dA, name="data", sample_def=match_def, fun=mean, na.rm=TRUE)
		dB = sampleApply(d=dB, name="data", sample_def=match_def, fun=mean, na.rm=TRUE)
		dcomb = mergeDataListsByCol(list(dA, dB),
			# list(
			# 	# make sure the data name is the same
			# 	list(data=dA[[data_nameA]], meta_row=dA$meta_row, meta_col=dA$meta_col),
			# 	list(data=dB[[data_nameB]], meta_row=dB$meta_row, meta_col=dB$meta_col)), 
			match_condition=match_def)  # 1-to-1 match
	} else if (match_method == "bootstrap") {
		dcomb = mergeDataListsByCol(list(dA, dB),
			# list(
			# 	list(data=dA[[data_nameA]], meta_row=dA$meta_row, meta_col=dA$meta_col),
			# 	list(data=dB[[data_nameB]], meta_row=dB$meta_row, meta_col=dB$meta_col)), 
			match_condition=match_def, bootstrap=TRUE)
	} else if (match_method == "sample") {
		dcomb = mergeDataListsByCol(list(dA, dB),
			# list(
			# 	list(data=dA[[data_nameA]], meta_row=dA$meta_row, meta_col=dA$meta_col),
			# 	list(data=dB[[data_nameB]], meta_row=dB$meta_row, meta_col=dB$meta_col)), 
			match_condition=match_def, bootstrap=FALSE)
	} else {
		stop("invalid match_method", match_method)
	}

	# Remove entries that were not matched
	allmissing = allMissingValuesCol(dcomb, "data")
	exclude = do.call('|', allmissing)

	dcomb = colSubsetMatchedCollection(dcomb, !exclude)

	# Replace missing values with zeros
	dcomb[[1]]$data[is.na(dcomb[[1]]$data)] = 0.0
	dcomb[[2]]$data[is.na(dcomb[[2]]$data)] = 0.0


	# Calculate correlations
	cor = list()  # output structure
	# sample cor
	cor$A = corFun(dcomb[[1]]$data)
	cor$B = corFun(dcomb[[2]]$data)

	# feature cor
	cor$AT = corFun(t(dcomb[[1]]$data))
	cor$BT = corFun(t(dcomb[[2]]$data))

	cor$ATBT = corFun(t(dcomb[[1]]$data), t(dcomb[[2]]$data))

	# Metadata
	cor$meta_rowA = dcomb[[1]]$meta_row
	cor$meta_colA = dcomb[[1]]$meta_col

	cor$meta_rowB = dcomb[[2]]$meta_row
	cor$meta_colB = dcomb[[2]]$meta_col

	# Invocation record
	cor$invoke = list()
	# cor$invoke$id = 
	cor$invoke$sourceA = sourceA
	cor$invoke$sourceB = sourceB
	cor$invoke$transformA = deparse(transformA)
	cor$invoke$transformB = deparse(transformB)
	cor$invoke$descA = descA
	cor$invoke$descB = descB
	cor$invoke$date = as.character(Sys.Date())
	cor$invoke$match_method = match_method
	return(cor)
}


# Calculates hold-out z-scores across each row. Returns matrix of zscores of the
# same dimensionality as the input matrix.
zscoreHoldout = function(X, separator=NULL) {
	stopifnot(!is.null(separator) & length(separator) == ncol(X))


	zscore = matrix(NA, nrow=nrow(X), ncol=ncol(X))

	if (is.null(separator)) {
		for (i in 1:nrow(X)) {
			for (j in 1:ncol(X)) {
				complement = !((1:ncol(X)) == j)

				zscore[i, j] = (X[i, j] - mean(X[i, complement], na.rm=TRUE)) / sd(X[i, complement], na.rm=TRUE)
			}
		}
	} else {
		for (i in 1:nrow(X)) {
			for (j in 1:ncol(X)) {
				complement = separator[j] == separator & !((1:ncol(X)) == j)
				zscore[i, j] = (X[i, j] - mean(X[i, complement], na.rm=TRUE)) / sd(X[i, complement], na.rm=TRUE)
			}
		}
	}

	return(zscore)
}

