import "field.asy" as field;
import graph;
import palette;

section_t Atlantic;
  Atlantic.lon_min = -109;
  Atlantic.lon_max = 69;
  Atlantic.lat_min = -88;
  Atlantic.lat_max = 88;

pen[] Palette=BWRainbow();
pen[] pal = Gradient(100 ... new pen[] {white, heavymagenta, blue, cyan, green, yellow, red, black});
pen[] pal1= pal[25:675];
pen[] pall = {palegrey, mediumyellow, mediumgreen, lightblue, heavyred, black};

string midticks(real x) {
  if ((((int)(x*2)) % 2) == 0) {
      return format("%d",(int)x);
  } else
      return "";
}

void draw_frame(int time, string fmt="png") {
    picture pic, pic1, pic2;
    size(pic1, 23cm, 0);
    size(pic2, 23cm, 0);
    draw_field(pic1, 'vort', range = new real[] {-0.1, 0.1}, time=time, Atlantic, pale=pal1);
    palette(pic1, "relative vorticity", bounds(-0.1, 0.1), (-250,-260), (250,-252), axis=Bottom, pal1,
        PaletteTicks(pTick=fontsize(40)));
    draw_field(pic2, 'duallevel', range = new real[] {7, 10}, time=time, Atlantic, pale=pall);
    palette(pic2, "grid level", bounds(7-0.5, 10+0.5), (-250,-260), (250,-252), axis=Bottom, pall,
        PaletteTicks(ticklabel=midticks, Step=0.5, pTick=fontsize(40)));
    add(pic,pic1.fit(),(0,0),E);
    add(pic,pic2.fit(),(0,0),5W);
    shipout('frame' + '_' + format("%03d", time), pic, format=fmt);
}

real[][] load_field(int t, int k) {
    file in=input("fort.1" + format("%03d",t) + format("%01d",k));
    real[][] field = in.line().dimension(0,0);
    return transpose(field);
}

/* first arrival and maximum wave height */ 
real TIME_INVALID = 0.9e16;

for (int t = 0; t <= 12; t+=1) {
    draw_frame(t, fmt="png");
}

