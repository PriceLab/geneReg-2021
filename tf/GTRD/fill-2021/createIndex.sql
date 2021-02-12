\connect genereg2021;
drop index if exists grtd_index;
create index grtd_index on grtd (chrom, start, endpos, tf)


