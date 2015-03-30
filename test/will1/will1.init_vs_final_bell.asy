import "field.asy" as field;

section_t America;
  America.lon_min = -180;
  America.lon_max = 0;
  America.lat_min = -90;
  America.lat_max = 90;

picture pic1, pic2;
size(pic1, 23cm, 0);
size(pic2, 23cm, 0);

draw_field(pic1, 'height', range = new real[] {-10, 1010}, time=0, section=America);
draw_field(pic2, 'height', range = new real[] {-10, 1010}, time=12, section=America);

add(pic1.fit(),(0,0),E);
add(pic2.fit(),(0,0),5W);
