knitr_output_options = list(
  mathjax = NULL,
  self_contained = FALSE,
  lib_dir = paste0("../../reports/", outdir, "/_libs")
  #lib_dir = "_libs"
)

render("snakemake/scripts/test_render.Rmd",
       output_file = paste0(x, 
                            "_", 
                            template, 
                            ".html"),
       output_dir = paste0("reports/", outdir),
       output_options = knitr_output_options)
