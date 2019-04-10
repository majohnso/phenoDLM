#' MTCI Data
#'
#' A dataset containing MTCI values for 2003 - 2008 on a grid over southern India. Each column corresponds to a grid cell in the region, each row corresponds to an 8-day composite measurement. Data were obtained from the NERC Earth Observation Data Centre (\url{www.neodc.rl.ac.uk/}).
#'
#' @format A data frame with 230 rows and 2163 columns.
"mtci_data"

#' Grid Cell Data
#'
#' A dataset containing latitude, longitude, landcover and cell index complementing the data in \code{mtci_data}. 
#'
#' @format A data frame with 2163 rows and 4 columns:
#' \describe{
#'   \item{lon}{Longitude of the grid cell}
#'   \item{lat}{Latitude of the grid cell}
#'   \item{landcover}{Majority landcover classification for the grid cell, derived from the GLC 2000 dataset (\url{http://forobs.jrc.ec.europa.eu/products/glc2000/products.php})}
#'   \item{loc}{Index for the grid cell, matches column names of \code{mtci_data}}
#' }
"location_attributes"