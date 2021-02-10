default:
	@echo targets: createDB  createTables cut fill psql

createDB:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena  postgres -f createDatabase-genereg2021.sql

createTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f createTable-rbp.sql

fill:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f fillTable-rbp.sql

index:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f createIndex.sql

psql:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021

