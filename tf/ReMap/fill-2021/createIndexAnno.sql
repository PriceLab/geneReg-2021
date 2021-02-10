\connect genereg2021;
drop index if exists remapAnno_index;
create index remapAnno_index on remapAnno (id, tf, celltype);

