# plotting helpers and useful functions

library(PaletteWoodsHole)
library(ggplot2)
library(tidyverse)

# plotting themes and labels
cols.site <- tulips
labels.site <-c("Fl:Dry Tortugas", "Fl:Grecian", "Fl:Biscayne", "VI:Flat Cay", "VI:Black Point", "VI:Brewer's Bay")
breaks.site <-c("FLA_015", "FLA_049", "FLA_073", "Flat Cay", "Black Point", "Brewer's Bay")
breaks.transect <- c("FLA_015 1","FLA_015 2","FLA_015 3", "FLA_049 1","FLA_049 2","FLA_049 3", "FLA_073 1","Flat Cay 1","Flat Cay 2","Flat Cay 3","Black Point 1","Black Point 2","Black Point 3","Brewer's Bay 1","Brewer's Bay 2","Brewer's Bay 3")
cols.study <- wefa_sun
labels.study <- c("Apprill et al 2021",  "Becker et al 2020", "Becker et al unpub",  "Neave et al 2017", "Weber et al 2020")
theme_set(theme_light()+ theme(legend.text = element_text(size=8),
                               plot.title=element_text(size=9, hjust=0.5, face = "bold")))


# functions
# A funtion to return pairwise bc distances from the phyloseq distance function
# x is a ps
distance <- phyloseq::distance
makePairwiseDistance <- function(x){
  d <- distance(x, method="bray")
  return(d %>% as.matrix() %>% as.data.frame.table(responseName="bray_distance"))
}

# generates an equation ready to insert into a distance-decay plot
# **makes the slope negative because I reverse my y axes in my plots for distance-decay**
# **not suitable for general use**
lm_eqn <- function(lmresult){
  a <- lmresult$coefficients[1]
  a <- ifelse(sign(a) >= 0,
              paste0(" + ", format(a, digits = 4)),
              paste0(" - ", format(-a, digits = 4)))
  pvalue <- ifelse(lmresult$coefficients[2,4]<2e-16,
                   paste0("italic(p) < 0.001"),
                   paste0("italic(p)==",format(lmresult$coefficients[2,4],digits=3)))
  eq1 <- substitute( paste( italic(y) == b, italic(x)),
                     list(#a = a, got rid of intercept as it has no interpretation in this context
                          b = format(-lmresult$coefficients[2], digits = 4))) #note I made coef negative because I made the y axis reversed!
  eq2 <- substitute( paste( italic(R)^2 == r2 ),
                     list(r2 = format(lmresult$r.squared, digits = 3)))
  eq3 <- pvalue
  list(eq1, eq2, eq3)
}

# function for cleaning corncob DA differentialTest output
# retrieves only the mu (abundance) scores for each significant model
# custom removes the extraneous letters in the variable
cleanDA <- function(da){
  df <- tibble()
  n <- da$significant_models[[1]]$np.mu
  for(i in 1:length(da$significant_models)){
    df <- rbind(df, da$significant_models[[i]]$coefficients[1:n,] %>% as.data.frame() %>% rownames_to_column(var="site"))
  }
  df$site <- df$site %>% gsub(x = ., pattern = "mu.site|mu.\\(", replace="") %>% gsub(x=., pattern=")","") %>% factor()
  tax <- da$data %>% tax_table() %>% as("matrix") %>% data.frame() %>% rownames_to_column(var="ASV") %>%
    filter(ASV %in% da$significant_taxa) %>%  slice(rep(1:n(), each = n))
  df <- bind_cols(df, tax) %>% as_tibble() %>% rename(StdE = "Std. Error", t= "t value", p = "Pr(>|t|)")
  return(df)
}

# functions for calculating distance between points on a globe
# taken from https://eurekastatistics.com/calculating-a-distance-matrix-for-geographic-points-using-r/
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.

  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.

  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }

  n.geopoints <- nrow(df.geopoints)

  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints

  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})

  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")

  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name

  return(mat.distances)
}

# borrowed these functions from https://rdrr.io/github/mworkentine/mattsUtils/src/R/microbiome_helpers.R
#' Add taxonomy label
#'
#' add a column to the taxonomy table of a phyloseq object that lists the
#' lowest rank taxonomy assigned to that OTU along with a prefix indicating
#' the taxonomic rank.
#'
#' Example: g:Pseudomonas
#'
#' @param physeq a valid phyloseq object that contains a taxonomy table
#' @param num_species the number of species to retain if more than one are identified
#' @return a phyloseq object with an additional column on the taxonomy
#'         table called "Taxonomy"
#' @export

split_species = function(string, n = 2) {
  splits = str_split(string, "/", n + 1)
  res = map_if(splits, ~length(.x) > 2, ~.x[1:n]) %>%
    map_chr(str_c, collapse = "/")
  return(res)
}

add_taxonomy_column = function(physeq, num_species = 2) {
  tax_df = as.data.frame(tax_table(physeq)@.Data) %>%
    rownames_to_column("OTU") %>%
    mutate(Species = split_species(Species, n = num_species)) %>%
    mutate(Taxonomy =
             case_when(
               is.na(Class)  ~ str_c("p:", Phylum),
               is.na(Order)  ~ str_c("c:", Class),
               is.na(Family)  ~ str_c("o:", Order),
               is.na(Genus)   ~ str_c("f:", Family),
               is.na(Species) ~ str_c("g:", Genus),
               TRUE ~ str_c(Genus, " ", Species)
             )
    )

  tax = as.matrix(tax_df[, -1])
  rownames(tax) = tax_df$OTU
  tax_table(physeq) = tax_table(tax)

  return(physeq)
}
