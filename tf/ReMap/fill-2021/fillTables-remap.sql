\connect genereg2021;;
\copy rbp from 'tblAnno-ready.tsv' delimiter E'\t' CSV  NULL as 'NULL';
select count(*) from rbp;

