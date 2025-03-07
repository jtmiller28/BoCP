### Check max dw time
library(data.table)

one <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl.csv")
two <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl2.csv")
three <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl3.csv")
four <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl4.csv")
five <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl5.csv")
six <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl6.csv")
seven <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl7.csv")
eight <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl8.csv")
nine <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl9.csv")
ten <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl10.csv")
eleven <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl11.csv")
twelve <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl12.csv")
thirteen <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl13.csv")
fourteen <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl14.csv")
fifthteen <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl15.csv")
sixteen <- fread("./data/gbif-dws/raw-dw-info-tbls/raw-info-tbl16.csv")

summary_dt <- rbind(one, two, three, four, five, six, seven, eight, nine, 
                    ten, eleven, twelve, thirteen, fourteen, fifthteen, sixteen)

# Find the names we need
missing_taxa <- setdiff(unique(name_alignment$acceptedNameParent), unique(summary_dt$acceptedParentName))
missing_taxa <- data.frame(acceptedNameParent = missing_taxa)
fwrite(missing_taxa, "/home/millerjared/blue_guralnick/millerjared/BoCP/data/gbif-dws/raw-dw-info-tbls/missing-taxa-v2-07-22-24.csv")
