# Connectivity map (CMAP) functions.

# Returns .cel file name of experiment given a metadata row as provided by the CMAP database (from the .xls file).
getCmapFileName = function(meta_row) {
	# library(stringr)
	require(stringr)

	# Get CMAP pertubation id
	pert_scan_id = as.character(meta_row$perturbation_scan_id)
	pert_scan_id = str_replace_all(pert_scan_id, "'", "")  # remove "'" for correspondende with .cel file names

	file_name = paste0(pert_scan_id, ".CEL")

	return(file_name)
}

# Returns vector of control file names given a CMAP metadata row.
# meta_row has to be a data frame slice
getCmapCtrlFileNames = function(meta_row) {
	# Base filename
	pert_scan_id = as.character(meta_row$perturbation_scan_id)
	pert_scan_id = str_replace_all(pert_scan_id, "'", "")  # remove "'" for correspondende with .cel file names

	# Control file specifications. Can contain multiple files per entry. Dot separated.
	pert_ctrl_scan_id = as.character(meta_row$vehicle_scan_id) 
	pert_ctrl_scan_id = str_replace_all(pert_ctrl_scan_id, "'", "")  # remove "'" for correspondende with .cel file names

	ctrl_file_names = vector()

	num_ctrl_files = count(meta_row$vehicle_scan_id4, ".")  # number of
	# if (cmap_match_sampl$array3[i] == "HT_HG-U133A_EA" || ) {
	if (num_ctrl_files > 1) {
		# Controls have 5-6 .cell files, encoded by dot separation in vehicle_scan_id.
		# Find file names of control exeriments.
		file_extensions = unlist(strsplit(pert_ctrl_scan_id, "[.]"))
		file_extensions = file_extensions[file_extensions != ""]  # remove empty

		for (f_ext in file_extensions) {
			file_name_base = strsplit(pert_scan_id, "[.]")[[1]][1]  # basis filename of perturbation and control experiments

			file_name = paste0(file_name_base, ".", f_ext, ".CEL")

			ctrl_file_names = c(ctrl_file_names, file_name)  # append
		}

	} else {
		# Single file control experiment
		file_name = paste0(pert_ctrl_scan_id, ".CEL")
		ctrl_file_names = c(ctrl_file_names, file_name)
	}

	return(ctrl_file_names)
}
