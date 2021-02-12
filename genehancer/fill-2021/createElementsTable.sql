\connect gh50;
drop table elements;

create table elements(chr varchar,
                      element_start int,
                      element_end int,
                      GHid varchar,
                      is_elite boolean,
                      type varchar
                      );

grant all on table "elements" to trena;

