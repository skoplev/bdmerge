# Returns the file path of execution.
# Only works if R script is invoked at command line.
# To get absolute paths use:
#	paste0(getwd(), "/", getScriptPath())
# Or simplu:
#   setwd(getScriptPath())
getScriptPath = function() {
    cmd.args = commandArgs()
    m = regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir = dirname(regmatches(cmd.args, m))
    if (length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if (length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(script.dir)
}
