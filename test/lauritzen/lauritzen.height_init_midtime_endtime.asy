import "field.asy" as field;

picture pic1, pic2;
picture pic3, pic4;
size(pic1, 23cm, 0);
size(pic2, 23cm, 0);
size(pic3, 23cm, 0);
size(pic4, 23cm, 0);

draw_field(pic1, 'height', range = new real[] {-0.01, 1.01}, time=0, section=Pacific);
draw_field(pic2, 'height', range = new real[] {-0.01, 1.01}, time=5, section=Pacific);
draw_field(pic3, 'height', range = new real[] {-0.01, 1.01}, time=10, section=Pacific);

picture pic12, pic123;
add(pic12, pic1.fit(),(0,0),E);
add(pic12, pic2.fit(),(0,0),5W);
add(pic123, pic12.fit(),(0,0),E);
add(pic123, pic3.fit(),(0,0),5W);
add(pic123.fit(),(0,0),E);
