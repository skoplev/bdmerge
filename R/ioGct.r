# Input-output functions handling .gct files.
# -------------------------------------------------------

# reads a gct file and returns data frame 
readGct = function(file_name) {
	gct_header = read.table(file_name, skip=1, nrows=1, sep="\t")  # read .gct file header
	meta_nrow = gct_header[1, 4]  # rows to skip
	meta_ncol = gct_header[1, 3]  # columns to skip

	if (is.null(meta_nrow)) {
		meta_nrow = 0
	}

	if (is.null(meta_ncol)) {
		meta_ncol = 0
	}

	# Reads a .gct file into a data frame. Flat format of both main data, row meta data, and column meta data.
	# main_mat = read.table(file_name, skip=2, sep="\t", header=FALSE)
	main_mat = fread(file_name, skip=2, header=FALSE)

	out = parseAnnotatedMatrix(as.data.frame(main_mat), meta_nrow, meta_ncol)  # seperate main and meta data. out.x, out.meta_row, out.meta_col
	return(out)
}

# Reads
readGctMeta = function(file_name) {
	gct_header = read.table(file_name, skip=1, nrows=1, sep="\t")  # read .gct file header
	meta_nrow = gct_header[1, 4]
	meta_ncol = gct_header[1, 3] 

	meta_data = fread(file_name, skip=2, nrows=meta_nrow+2)  # column meta data and first data entry

	out = parseAnnotatedMatrix(as.data.frame(meta_data), meta_nrow, meta_ncol)

	return(out)
}


# Returns list containing main data in $data and meta data in $meta_row and $meta_col
# meta_nrow, meta_ncol is from .gct file and encodes where the 1-indexed positions where the
# the main data matrix starts.
parseAnnotatedMatrix = function(mat, meta_nrow, meta_ncol) {
	# Input check
	if (!is.integer(meta_nrow) & meta_nrow != 0) {
		stop("meta_nrow is not an integer")
	}

	if (!is.integer(meta_ncol) & meta_ncol != 0) {
		stop("meta_ncol is not an integer")
	}

	# Main data matrix
	x = mat[(meta_nrow+2):nrow(mat), (meta_ncol+2):ncol(mat)]  # indicies are zero indexed and does not account for row and column names (hence +2)
	# Convert to numeric matrix
	x = as.matrix(x)
	class(x) = "numeric"

	# Column meta data, stored in data frame with named columns for each type of metadata. Columns are samples.
	if (meta_nrow == 0) {
		meta_col = mat[1, (meta_ncol+2):ncol(mat), drop=FALSE]  # only the second row
	} else {
		# non-zero.
		meta_col = mat[1:(meta_nrow+1), (meta_ncol+2):ncol(mat)]
	}

	# meta_col = d[2:meta_nrow, (meta_ncol+2):ncol(d)]
	meta_col = t(meta_col)
	meta_col = as.data.frame(meta_col)

	if (meta_nrow == 0) {
		colnames(meta_col) = "id"
	} else {
		colnames(meta_col) = c("id", as.vector(mat[2:(meta_nrow+1), 1]))
	}

	# Row meta data, stored in data frame with named columns for each type of metadata. Rows are peptides
	meta_row = mat[(meta_nrow+2):nrow(mat), 1:(meta_ncol+1)]
	meta_row = as.data.frame(meta_row)

	colnames(meta_row) = as.vector(t(mat[1, 1:(meta_ncol+1)]))  # transpose to column vector

	# out = list(data=data.matrix(x), meta_col=meta_col, meta_row=meta_row)
	out = list(data=x, meta_col=meta_col, meta_row=meta_row)

	return(out)
}
