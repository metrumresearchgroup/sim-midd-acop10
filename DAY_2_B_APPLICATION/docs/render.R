library(rmarkdown)

files <- list.files(pattern="\\.R$", full.names=TRUE)
files <- files[!(files %in% c("Exercise_2_EC.html"))]

opt <- list(css = "../common/styles.css")

for(f in files) {
  try(
    rmarkdown::render(
      input=f, 
      output_dir= "docs",
      output_format="html_document",
      output_options = opt, 
      knit_meta=list(author = "MetrumRG")
    )
  )
}

