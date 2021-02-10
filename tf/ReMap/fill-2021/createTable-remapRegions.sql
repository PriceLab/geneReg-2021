\connect genereg2021;
drop table remapRegions;
create table remapRegions (chrom varchar,
                           start int,
                           endpos int,
                           name varchar,
                           score numeric,
			   peakstart int,
			   peakend int
                           );
grant all on table remapRegions to trena;

