# Make paper and push up to google doc
rmarkdown::render("analysis/paper/manuscript.Rmd", output_format = "bookdown::word_document2")

# upload to google drive for feedback
googledrive::drive_put("analysis/paper/manuscript.docx", 'SD_dynamics_paper')
