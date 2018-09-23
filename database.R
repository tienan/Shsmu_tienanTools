Rmysql
#????????

#    2.1 dbConnect??dbDisconnect  #???ݿ?��?Ӻ???
#    2.2 dbListFields ,dbListTables,dbGetInfo,dbListResults,summary,dbGetException,dbExistsTable    #?鿴???ݿ????߲?????Ϣ

#??ѯ????

#    2.3 dbGetQuery    #??ѯ????
#    2.4 dbReadTable   #??ȡ????????

#???º???

#   2.5 dbWriteTable   #???????ݿ??????߽?????д????Ӧ?ı?
#    2.6 dbRemoveTable  # ɾ?????ݿ??еı?

#dbSendQuery????

#    2.7 dbSendQuery ,dbClearResult  #??query???????ݿ?????
#    2.8 dbColumnInfo??dbGetRowsAffected??dbGetRowCount??dbHasCompleted     #?鿴???ݿ?????ִ?н???
#    2.9 dbFetch??fetch   #??dbSendQuery?????Ľ?????ȡ??��
#    2.10 dbNextResult??dbMoreResults   #һ??һ????ȡ????

#????????

#   2.11 dbCommit??dbBegin??dbRollback 

install.packages("RMySQL")
library(RMySQL)

####For example
conn <- dbConnect(RMySQL::MySQL(),dbname="health",username="root",password="feng1234",host="172.29.1.30")


dbListTables(conn)
dbWriteTable(conn, "mtcars", mtcars)
dbWriteTable(conn, "mtcars", mtcars,append=T) 
dbListTables(conn)
dbListFields(conn, "mtcars")
dbReadTable(conn, "mtcars")

# You can fetch all results:
res <- dbSendQuery(conn, "SELECT * FROM mtcars WHERE cyl = 4")
dbFetch(res)
dbClearResult(res)

# Or a chunk at a time
res <- dbSendQuery(con, "SELECT * FROM mtcars WHERE cyl = 4")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}
# Clear the result
dbClearResult(res)

#??ȡ????
dbGetQuery(con, "SELECT * FROM mtcars")


# Disconnect from the database
dbDisconnect(con)


#sql ?????ٴ??ͻ???????

sql = "select l.sampleID,l.cancer,l.patientID,g.gene_id,g.normalized_count,c.gender,c.race,c.pathologic_T,c.pathologic_N,c.pathologic_M,c.pathologic_stage,c.tobacco_smoking_history,c.days_to_last_followup,c.days_to_death from clinic_luad as c,genes2id_link as l,LUAD_genes as g where c.bcr_patient_barcode=l.patientID and g.sample_id=l.sampleID"

sql = " select  patient_id,chromosome, m.methylation_id ,gene_id,  std(m.value)  from LUAD_methylation as m group by m.methylation_id order by std(m.value)"

sql = "select * from LUAD_genes"

sql = "select * from clinic_luad as c, genes2id_link as l where c.bcr_patient_barcode=l.patientID"

clinic = dbGetQuery(conn , sql)

#tmp is the file where 


res = dbGetQuery(conn , sql)




write.table(res,file="LUAD_genes_exp.txt",quote = F,sep="\t",row.names = F,fileEncoding = "utf-8")








