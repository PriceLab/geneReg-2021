default:
	@echo createDB createTables

createDB:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena postgres  -f createDatabase-gh59.sql

psql:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh59


createTables: createElementsTable  createGeneAssociationsTable createTFBsTable  createTissuesTable

createElementsTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh59 -f createElementsTable.sql

createGeneAssociationsTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh59 -f createGeneAssociationsTable.sql

createTFBsTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh59 -f createTFBsTable.sql

createTissuesTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh59 -f createTissuesTable.sql

fillElements:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh59 -f fillElements.sql
