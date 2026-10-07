# ggplot2 is imported whole. The plotting code uses dozens of ggplot2 functions,
# and the explicit importFrom() list it replaced went stale (no geom_col), so
# figure generation failed in every batch run while still working in any
# interactive session that had attached ggplot2 itself.
#' @import data.table
#' @import ggplot2
#' @importFrom SMFnorm normalize_methylation_data find_shared_sites load_data
"_PACKAGE"

