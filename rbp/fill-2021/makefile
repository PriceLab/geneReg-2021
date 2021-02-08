default:
	@echo targets: createDB  createTables fill

createDB:
	 PGPASSWORD=trena psql -h khaleesi -U trena  postgres -f createDatabase-genereg2021.sql

createTable:
	 PGPASSWORD=trena psql -h khaleesi -U trena genereg2021 -f createTable-rbp.sql

cut:
	awk -d\t 'NR > 1{print $$1, $$2, $$3, $$4, $$5, $$7, $$8, $$9, $$10, $$11, $$20, $$21}' ../incoming/tblAnno.tsv | tr " " "\t" > tblAnno-ready.tsv

fill:
	PGPASSWORD=trena psql -h khaleesi -U trena genereg2021 -f fillTable-rbp.sql

psql:
	PGPASSWORD=trena psql -h khaleesi -U trena genereg2021

