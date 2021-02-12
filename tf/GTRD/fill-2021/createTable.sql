\connect genereg2021;
drop table grtd;
create table grtd (chrom varchar,
                   start int,
                   endpos int,
                   width int,
                   tf varchar
                   );
grant all on table grtd to trena;
