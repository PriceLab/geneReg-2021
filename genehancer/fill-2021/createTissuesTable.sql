\connect gh50;
drop table tissues;

create table tissues(GHid varchar,
                    source varchar,
		    tissue varchar,
		    category varchar
                    );

grant all on table "tissues" to trena;

