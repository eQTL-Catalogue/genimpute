suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("readr"))


option_list <- list(optparse::make_option(c("-p", "--filename"), type="character", default=NULL,help="Variant pairs with pearson Rs", metavar = "type"))

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

is_diploid <- function(vector){
  vector %like% "\\|"
}

save_file <- function(df, file_name, has_col_names){
  write.table(df, file = file_name, append = FALSE, quote = F, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = F,
              col.names = has_col_names, qmethod = c("escape", "double"),
              fileEncoding = "")
}
filename = opt$filename

#Import genotypes
gt_matrix = readr::read_tsv(filename, comment = "##")
gt_mat = gt_matrix[,-c(1:9)]

#Identify deploid genotypes
diploid_tibble = dplyr::mutate(gt_mat, across(everything(), is_diploid))
diploid_count = colSums(as.matrix(diploid_tibble))
selection = diploid_count == max(diploid_count) | diploid_count == 0
selected_idvs = gt_mat[,!selection]

#Build final table
problematic_individuals = dplyr::bind_cols(gt_matrix[,c(1:9)], selected_idvs)
problematic_individuals_gt_matrix = problematic_individuals[,-c(1:9)]
if (dim(problematic_individuals_gt_matrix)[2] == 0) {
  save_file(data.frame(), "problematic_individuals.txt", FALSE )
} else{
  save_file(problematic_individuals, "problematic_individuals.txt", TRUE)
}
if (dim(problematic_individuals_gt_matrix)[2] == 0) {
  save_file(data.frame(), "problematic_variants.txt", FALSE )
} else {
  #find problematic variants
  has_diploid = dplyr::mutate(problematic_individuals_gt_matrix, across(everything(), is_diploid))
  has_any_diploid_per_var = has_diploid %>% mutate(has_any_diploid_per_var = reduce(., `|`))
  problematic_indivs_labled_prob_vars = dplyr::bind_cols(problematic_individuals[,c(1:9)], has_any_diploid_per_var)
  problem_vars = filter(problematic_indivs_labled_prob_vars, has_any_diploid_per_var == T)  %>% select(ID)
  save_file(problem_vars, "problematic_variants.txt", FALSE )
}
