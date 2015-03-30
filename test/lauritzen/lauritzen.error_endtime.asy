import "field.asy" as field;
import graph;
import palette;

size(23cm, 0);

draw_field('height', range = new real[] {-0.04, 0.04}, time=20, section=Pacific);

pen[] pal=BWRainbow();
palette("error at t=T", bounds(-0.05, 0.05), (-250,-260), (250,-250), axis=Bottom, pal, 
        PaletteTicks(pTick=fontsize(20)));
