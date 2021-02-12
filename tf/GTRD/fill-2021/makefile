default:
	@echo targets: createDB  createTable fill index psql

createTable:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f createTable.sql

fill:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f fillTable.sql

index:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f createIndex.sql

psql:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021

