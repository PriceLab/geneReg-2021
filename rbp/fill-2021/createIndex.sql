\connect genereg2021;
drop index if exists rbp_index;
create index rbp_index on rbp (chrom, start, endpos, gene, target);

