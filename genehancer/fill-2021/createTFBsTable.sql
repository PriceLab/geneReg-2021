\connect gh50;
drop table tfbs;

create table tfbs(GHid varchar,
                  tf varchar,
		  tissues varchar
		   );

grant all on table "tfbs" to trena;
