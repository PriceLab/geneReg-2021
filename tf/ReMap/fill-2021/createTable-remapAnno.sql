\connect genereg2021;
drop table remapAnno;
create table remapAnno (id varchar,
                        tf varchar,
                        celltype varchar,
                        treatment varchar,
                        experiment varchar,
			url varchar
                        );
grant all on table remapAnno to trena;


