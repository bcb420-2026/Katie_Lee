# Exercises

These are intended to be done **after** completing the worked examples.

## Exercise 1 — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119732

Using **GSE119732**, confirm whether the ID column contains Ensembl IDs with version suffixes.

1. Extract the first 20 IDs.
2. Count how many contain a `.`.
3. Create a new column with versions stripped.
4. Map the identifiers to HGNC symbols.

``` r
library(knitr)
library(GEOquery)
library(readr)
library(biomaRt)
library(dplyr)
```


``` r
fetch_geo_supp <- function(gse, destdir = "data") {
  dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
  GEOquery::getGEOSuppFiles(GEO = gse, baseDir = destdir, makeDirectory = TRUE)
  invisible(file.path(destdir, gse))
}

# 1. Extract first 20 IDs
gse <- "GSE119732"
geo_path <- file.path("data", gse)

fetch_geo_supp(gse = gse)
file <- list.files(file.path("data", gse), full.names = TRUE)
df <- readr::read_tsv(file, show_col_types = FALSE)

df <- df[1:20, ]
head(df)
```

```
## # A tibble: 6 × 30
##   gene_id         A1    A2    A3    A4    B1    B2    B3    B4    B5    C1    C2
##   <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
## 1 ENSG0000022…     0     0     0     0     0     0     0     0     0     1     0
## 2 ENSG0000022…    79   119    84    50    80    72    50    22    54    26   110
## 3 ENSG0000027…    17    10    22    19    19    22    16     5    18     9    20
## 4 ENSG0000024…     0     0     0     0     0     0     0     0     0     2     4
## 5 ENSG0000023…     0     0     0     0     0     0     0     0     0     0     0
## 6 ENSG0000026…     0     0     0     0     0     0     0     0     0     0     0
## # ℹ 18 more variables: C3 <dbl>, C4 <dbl>, C5 <dbl>, D1 <dbl>, D2 <dbl>,
## #   D3 <dbl>, D4 <dbl>, D5 <dbl>, E1 <dbl>, E2 <dbl>, E3 <dbl>, E4 <dbl>,
## #   E5 <dbl>, F1 <dbl>, F2 <dbl>, F3 <dbl>, F4 <dbl>, F5 <dbl>
```


``` r
# 2. Count how many contain a `.`
sum(grepl("\\.[0-9]+$", df$gene_id))
```

```
## [1] 20
```

``` r
# 3. Create a new column with versions stripped.
strip_ensembl_version <- function(x) sub("\\..*$", "", x)

id_col <- names(df)[1]
ids <- df[[1]] |> as.character()
if (sum(grepl("\\.[0-9]+$", ids)) >= 1) {
  ids <- strip_ensembl_version(ids)
}
df$versions_stripped <- ids
head(df)
```

```
## # A tibble: 6 × 31
##   gene_id         A1    A2    A3    A4    B1    B2    B3    B4    B5    C1    C2
##   <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
## 1 ENSG0000022…     0     0     0     0     0     0     0     0     0     1     0
## 2 ENSG0000022…    79   119    84    50    80    72    50    22    54    26   110
## 3 ENSG0000027…    17    10    22    19    19    22    16     5    18     9    20
## 4 ENSG0000024…     0     0     0     0     0     0     0     0     0     2     4
## 5 ENSG0000023…     0     0     0     0     0     0     0     0     0     0     0
## 6 ENSG0000026…     0     0     0     0     0     0     0     0     0     0     0
## # ℹ 19 more variables: C3 <dbl>, C4 <dbl>, C5 <dbl>, D1 <dbl>, D2 <dbl>,
## #   D3 <dbl>, D4 <dbl>, D5 <dbl>, E1 <dbl>, E2 <dbl>, E3 <dbl>, E4 <dbl>,
## #   E5 <dbl>, F1 <dbl>, F2 <dbl>, F3 <dbl>, F4 <dbl>, F5 <dbl>,
## #   versions_stripped <chr>
```

``` r
# 4. Map the identifiers to HGNC symbols.
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ids,
    mart = ensembl
  )

map$ensembl_gene_id <- as.character(map$ensembl_gene_id)

df_mapped <- df %>%
  left_join(map, by=c("versions_stripped" = "ensembl_gene_id")) %>%
  dplyr::select(gene_id, hgnc_symbol, everything())
```

``` r
df_mapped
```

```
## # A tibble: 20 × 32
##    gene_id     hgnc_symbol    A1    A2    A3    A4    B1    B2    B3    B4    B5
##    <chr>       <chr>       <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
##  1 ENSG000002… "DDX11L1"       0     0     0     0     0     0     0     0     0
##  2 ENSG000002… "WASH7P"       79   119    84    50    80    72    50    22    54
##  3 ENSG000002… "MIR6859-1"    17    10    22    19    19    22    16     5    18
##  4 ENSG000002… "MIR1302-2…     0     0     0     0     0     0     0     0     0
##  5 ENSG000002… "FAM138A"       0     0     0     0     0     0     0     0     0
##  6 ENSG000002… "OR4G4P"        0     0     0     0     0     0     0     0     0
##  7 ENSG000002… "OR4G11P"       0     0     0     0     0     0     0     0     0
##  8 ENSG000001… "OR4F5"         0     0     0     0     0     0     0     0     0
##  9 ENSG000002…  <NA>           0     0     0     1     1     0     1     0     1
## 10 ENSG000002… ""              0     0     0     0     0     0     0     0     0
## 11 ENSG000002… "CICP27"        0     0     0     0     0     0     0     0     0
## 12 ENSG000002… ""              7    10     9     6    29    30    16     8    32
## 13 ENSG000002… ""             17    11    14     4    27    37    21     9    34
## 14 ENSG000002… ""              0     0     1     0     0     0     0     0     0
## 15 ENSG000002… ""              6    19    12     8    60    45    24    11    30
## 16 ENSG000002… "RNU6-1100…     0     0     0     0     0     0     0     0     0
## 17 ENSG000002… ""              0     0     0     0     0     0     0     0     0
## 18 ENSG000002… "DDX11L17"      8     5     9     6     0     0     0     0     0
## 19 ENSG000002… "WASH9P"      282   298   334   205   156   163   102    42   114
## 20 ENSG000002… "MIR6859-2"     0     0     0     0     0     0     0     0     0
## # ℹ 21 more variables: C1 <dbl>, C2 <dbl>, C3 <dbl>, C4 <dbl>, C5 <dbl>,
## #   D1 <dbl>, D2 <dbl>, D3 <dbl>, D4 <dbl>, D5 <dbl>, E1 <dbl>, E2 <dbl>,
## #   E3 <dbl>, E4 <dbl>, E5 <dbl>, F1 <dbl>, F2 <dbl>, F3 <dbl>, F4 <dbl>,
## #   F5 <dbl>, versions_stripped <chr>
```
## Exercise 2 — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122380

Using **GSE122380**, confirm whether the ID column contains Ensembl IDs with version suffixes.

1. Extract the first 20 IDs.
2. Create a new column with versions stripped.
3. Map the identifiers to HGNC symbols.
4. What is different about this file? 

This file is empty!



