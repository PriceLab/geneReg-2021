\connect gh59;
drop table associations;

create table associations(GHid varchar,
                          symbol varchar,
			  eqtl_score numeric,
			  erna_score numeric,
			  chic_score numeric,
			  expression_score numeric,
			  distance_score numeric,
			  tss_proximity numeric,
			  combined_score numeric,
                          is_elite boolean
                          );

grant all on table "associations" to trena;

