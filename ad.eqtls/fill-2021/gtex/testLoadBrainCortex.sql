\connect genereg2021;;
\copy eqtls from 'src/unitTests/brain-cortex-ready.csv' delimiter E',' CSV  NULL as 'NULL' HEADER;
select count(*) from eqtls;

