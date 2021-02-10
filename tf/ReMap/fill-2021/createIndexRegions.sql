\connect genereg2021;
drop index if exists remapRegions_index;
create index remapRegions_index on remapRegions (chrom, start, endpos);

