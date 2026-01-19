# Define a writable library path in your home directory
lib_path <- Sys.getenv("R_LIBS_USER")

# Create the directory if it doesn't exist
if (!dir.exists(lib_path)) {
  dir.create(lib_path, recursive = TRUE)
}

# Set the library path
.libPaths(lib_path)

# Install renv
if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv", repos='http://cran.rstudio.com/', lib=lib_path)
}
renv::init()

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos='http://cran.rstudio.com/', lib=lib_path)
}

# Install Bioconductor and CRAN packages
BiocManager::install(c("AnnotationDbi", "UniProt.ws", "org.Hs.eg.db", "org.Mm.eg.db"), lib=lib_path)
install.packages(c("data.table", "dplyr"), repos='http://cran.rstudio.com/', lib=lib_path)

# Print library paths to confirm installation
print(.libPaths())

