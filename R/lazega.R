#' Lazega lawyers network data
#'
#' @source This network dataset comes from a network study of corporate law partnership 
#' that was carried out in a Northeastern US corporate law firm 
#' in New England from 1988 to 1991. It represents collaborative relations 
#' among the 36 attorneys (partners and associates) of this firm.
#' Nodal attributes include: Age, Gender, Office, Practice, School, and Years.
#' 
#' @format An oject of class \code{network}.
#' 
#' @references
#' Lazega, E. (2001), "The Collegial Phenomenon: 
#' The Social Mechanisms of Cooperation Among Peers in a 
#' Corporate Law Partnership," Oxford University Press. 
#' 
#' @examples
#' \dontrun{
#'  par(mfrow = c(1, 2), oma = rep(0, 4))
#'  CC <- hcl.colors(3, "Teal")
#'  set.seed(22)
#'  plot(lazega,
#'       vertex.col = CC[lazega %v% "Office"], 
#'      vertex.cex = 2)
#'  legend("topright",
#'         pch    = 21,
#'         pt.bg  = CC,
#'         legend = c("Boston", "Hartford", "Providence"),
#'         title  = "OFFICE")
#' }
"lazega"


