default:
	@echo targets:  createTables  fill psql

createTables:
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f createTable-remapAnno.sql
	 PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f createTable-remapRegions.sql

fill:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f fillTable-rbp.sql

index:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f createIndexRegions.sql
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f createIndexAnno.sql


psql:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021

