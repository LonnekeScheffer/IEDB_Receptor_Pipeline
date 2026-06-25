This code prepares the OGRDB germline reference set for usage in the IEDB receptor pipeline. 

Running [make_ogrdb_reference.sh](make_ogrdb_reference.sh) requires one argument: the path to the IgBLAST folder `ncbi-igblast-1.22.0`. 
The script will create a new folder `igblast_db_ogrdb` in the working directory. 

A complementary IMGT germline database can be found in the [ Immcantation suite docker v4.6.0](hub.docker.com/layers/immcantation/suite/4.6.0)