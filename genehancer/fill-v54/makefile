default:
	@echo createDB createTables

createDB:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena postgres  -f createDatabase-gh54.sql

psql:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh54


createTables: createElementsTable  createGeneAssociationsTable createTFBsTable  createTissuesTable

createElementsTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh54 -f createElementsTable.sql

createGeneAssociationsTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh54 -f createGeneAssociationsTable.sql

createTFBsTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh54 -f createTFBsTable.sql

createTissuesTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh54 -f createTissuesTable.sql

fillElements:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena gh54 -f fillElements.sql
