default:
	@echo createDB createTables

createDB:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena postgres  -f createDatabase-gh50.sql

psql:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh50


createTables: createElementsTable  createGeneAssociationsTable createTFBsTable  createTissuesTable

createElementsTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh50 -f createElementsTable.sql

createGeneAssociationsTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh50 -f createGeneAssociationsTable.sql

createTFBsTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh50 -f createTFBsTable.sql

createTissuesTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh50 -f createTissuesTable.sql

fillElements:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh50 -f fillElements.sql

# fill:
# 	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh50 -f fillTable-rbp.sql
# 
# index:
# 	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh50 -f createIndex.sql
# 
# 
# 
# 
# 
# 
# 
# create.all: create.elements create.associations create.tfbs create.tissues
# 
# fill.all: fill.elements fill.associations fill.tfbs fill.tissues
# 
# create.elements:
# 	(export PGPASSWORD=trena; psql -h khaleesi -U trena gh411 -f createElementsTable.sql)
# 
# create.associations:
# 	(export PGPASSWORD=trena; psql -h khaleesi -U trena gh411 -f createGeneAssociationsTable.sql)
# create.tfbs:
# 	(export PGPASSWORD=trena; psql -h khaleesi -U trena gh411 -f createTFBsTable.sql)
# create.tissues:
# 	(export PGPASSWORD=trena; psql -h khaleesi -U trena gh411 -f createTissuesTable.sql)
# 
# fill.elements:
# 	(export PGPASSWORD=trena; psql -h khaleesi -U trena gh411 -f fillElements.sql)
# fill.associations:
# 	(export PGPASSWORD=trena; psql -h khaleesi -U trena gh411 -f fillGeneAssociations.sql)
# fill.tfbs:
# 	(export PGPASSWORD=trena; psql -h khaleesi -U trena gh411 -f fillTFBSs.sql)
# 
# fill.tissues:
# 	(export PGPASSWORD=trena; psql -h khaleesi -U trena gh411 -f fillTissues.sql)
# 
# psql.session:
# 	(export PGPASSWORD=trena; psql -h khaleesi -U trena gh411)
# 

