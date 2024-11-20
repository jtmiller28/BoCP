library(qpdf)
# Native data
pdf_files <- list.files("./outputs/outlier-pdfs/native-data/100km/", full.names = TRUE)
pdf_combine(input = pdf_files, output = "./outputs/outlier-pdfs/native-data/100km/native-combined-100km.pdf")

pdf_files <- list.files("./outputs/outlier-pdfs/native-data/250km/", full.names = TRUE)
pdf_combine(input = pdf_files, output = "./outputs/outlier-pdfs/native-data/250km/native-combined-250km.pdf")

pdf_files <- list.files("./outputs/outlier-pdfs/native-data/500km/", full.names = TRUE)
pdf_combine(input = pdf_files, output = "./outputs/outlier-pdfs/native-data/500km/native-combined-500km.pdf")

# All data
pdf_files <- list.files("./outputs/outlier-pdfs/all-data/100km/", full.names = TRUE)
pdf_combine(input = pdf_files, output = "./outputs/outlier-pdfs/all-data/100km/all-combined-100km.pdf")

pdf_files <- list.files("./outputs/outlier-pdfs/all-data/250km/", full.names = TRUE)
pdf_combine(input = pdf_files, output = "./outputs/outlier-pdfs/all-data/250km/all-combined-250km.pdf")

pdf_files <- list.files("./outputs/outlier-pdfs/all-data/500km/", full.names = TRUE)
pdf_combine(input = pdf_files, output = "./outputs/outlier-pdfs/all-data/500km/all-combined-500km.pdf")

# Introduced + Native Data
pdf_files <- list.files("./outputs/outlier-pdfs/native-and-introduced-data/100km/", full.names = TRUE)
pdf_combine(input = pdf_files, output = "./outputs/outlier-pdfs/native-and-introduced-data/100km/native-and-introduced-combined-100km.pdf")

pdf_files <- list.files("./outputs/outlier-pdfs/native-and-introduced-data/250km/", full.names = TRUE)
pdf_combine(input = pdf_files, output = "./outputs/outlier-pdfs/native-and-introduced-data/250km/native-and-introduced-combined-250km.pdf")

pdf_files <- list.files("./outputs/outlier-pdfs/native-and-introduced-data/500km/", full.names = TRUE)
pdf_combine(input = pdf_files, output = "./outputs/outlier-pdfs/native-and-introduced-data/500km/native-and-introduced-combined-500km.pdf")