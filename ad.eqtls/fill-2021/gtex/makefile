default:
	@echo targets: fill-test

fill-test:
	PGPASSWORD=`cat ~/trena.passwd` psql -h khaleesi -U trena genereg2021 -f testLoadBrainCortex.sql
