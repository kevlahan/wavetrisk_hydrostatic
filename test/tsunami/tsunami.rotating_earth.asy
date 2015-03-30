import "field.asy" as field;
import graph;
import palette;

section_t Indian_Ocean;
  Indian_Ocean.lon_min = -9;
  Indian_Ocean.lon_max = 169;
  Indian_Ocean.lat_min = -88;
  Indian_Ocean.lat_max = 88;

real Hdim=3.344175893265152e+03;

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
    draw_field(pic1, 'height', range = new real[] {-0.7/Hdim, 0.7/Hdim}, time=time, Indian_Ocean, pale=pal1);
    palette(pic1, "height [m]", bounds(-0.69, 0.69), (-250,-260), (250,-252), axis=Bottom, pal1,
        PaletteTicks(pTick=fontsize(40)));
    draw_field(pic2, 'level', range = new real[] {7, 12}, time=time, Indian_Ocean, pale=pall);
    palette(pic2, "grid level", bounds(7-0.5, 12+0.5), (-250,-260), (250,-252), axis=Bottom, pall,
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
int S_ARRIVAL = 1;
int S_WAVE_H = 2;

real[][] arrival = load_field(1,S_ARRIVAL);
real[][] wave_h = load_field(1,S_WAVE_H);
for (int t = 2; t <= 10*2; t+=1) {
    write(t);
    real[][] field = load_field(t,S_WAVE_H);
    wave_h = sequence(new real[](int i) {return max(wave_h[i], field[i]);}, wave_h.length);
    real[][] field = load_field(t,S_ARRIVAL);
    arrival = sequence(new real[](int i) {return min(arrival[i], field[i]);}, arrival.length);
}

picture pic1;
write(max(arrival));
write(min(arrival));
size(pic1, 25cm, 0);
bounds range = image(pic1, arrival/3600.0, Range(0, 8), (-180,-90), (180,90), Palette);
clip(pic1, (25,-55)--(25,35)--(135,35)--(135,-55)--cycle);
palette(pic1, "arrival time [h]", range, (25,-60), (135,-57), axis=Bottom, Palette, 
        PaletteTicks(pTick=fontsize(20)));
shipout("first_arrival_time", pic1, format="png");

picture pic2;
size(pic2, 25cm, 0);
bounds range = image(pic2, wave_h, Range(0,0.7), (-180,-90), (180,90), pal1);
clip(pic2, (25,-55)--(25,35)--(135,35)--(135,-55)--cycle);
palette(pic2, "wave amplitude [m]", range, (25,-60), (135,-57), axis=Bottom, pal1, 
        PaletteTicks(pTick=fontsize(20)));
shipout("max_wave_height", pic2, format="png");

picture pic3;
bounds range = image(pic3, wave_h, Range(0,0.1), (-180,-90), (180,90), Palette);
image(pic3, wave_h, Range(0,0.1), (-180,-90), (180,90), Palette);
palette(pic3, "wave amplitude [m]", range, (-178,-105), (178,-92), axis=Bottom, Palette, 
        PaletteTicks(pTick=fontsize(14)));
shipout("max_wave_height2", pic3, format="pdf");

/*
Indian_Ocean.lon_min -= 0.5 * 116;
Indian_Ocean.lon_max -= 0.4 * 116;
*/
for (int t = 0; t <= 10*12; t+=1) {
    if (t > 96) {
        Indian_Ocean.lon_min = -9 - 0.5*(t-96);
        Indian_Ocean.lon_max = 169 - 0.5*(t-96);
    }
    draw_frame(t, fmt="png");
}
write(Indian_Ocean.lon_min);
write(Indian_Ocean.lon_max);

/*
draw_chi(60, fmt="pdf");
draw_chi(61, fmt="pdf");
draw_chi(62, fmt="pdf");
*/
