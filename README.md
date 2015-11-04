# bdmerge
This R library contains a variety of functions that assists in handling multiple types of data. There is an accompanying [website](http://amp.pharm.mssm.edu/bdmerge) where data and a selection of methods are available using a web-based interface.

# Data structures
Data is assumed to consist of sets of matching matrices (features x samples). Each data structure also have data frames describing the meta data of both the features and samples. The field names of the meta data can be freely chosen.

Internally, the library represents the data structure as a list of matrices with two data frames named "meta_col" and "meta_row". Functions are provided which write and read these structures as [HDF5](https://www.hdfgroup.org/HDF5) files.

# Methods
Once data is transformed into the standardized structure described above statistical methods can be formulated which takes as input one data set and returns another.
