\connect genereg2021;
drop table rbp;
create table rbp (chrom varchar,
                  start int,
                  endpos int,
                  width int,
		  strand char(1),
		  gene varchar,
		  method varchar,
		  celltype varchar,
		  accession varchar,
		  score numeric,
		  target varchar,
		  targetfeature varchar
                  );
grant all on table rbp to trena;
