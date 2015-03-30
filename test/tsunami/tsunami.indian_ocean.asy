import "field.asy" as field;
import graph;
import palette;

section_t Indian_Ocean;
  Indian_Ocean.lon_min = 20;
  Indian_Ocean.lon_max = 135;
  Indian_Ocean.lat_min = -25;
  Indian_Ocean.lat_max = 25;

void draw_frame(string variable, int time, string fmt="png") {
    picture pic;
    size(pic, 50cm, 0);
    draw_field(pic, variable, time=time, level=0, Indian_Ocean);
    shipout(variable + '_' + format("%02d", time), pic, format=fmt);
}

real[][] load_field(int t, int k) {
    file in=input("fort.1" + format("%03d",t) + format("%01d",k));
    return in.line().dimension(0,0); 
}

/* first arrival and maximum wave height */
real TIME_INVALID = 0.9e16;
int S_ARRIVAL = 1;
int S_WAVE_H = 2;
real[][] arrival = load_field(1,S_ARRIVAL);
real[][] wave_h = load_field(1,S_WAVE_H);
for (int t = 2; t <= 15; t+=1) {
    real[][] field = load_field(t,S_WAVE_H);
    wave_h = sequence(new real[](int i) {return max(wave_h[i], field[i]);}, wave_h.length);
    real[][] field = load_field(t,S_ARRIVAL);
    arrival = sequence(new real[](int i) {return min(arrival[i], field[i]);}, arrival.length);
}

write("ERROR: This file is currently out-dated. Please use tsunami.rotating_earth.asy");

for (int t = 0; t <= 100; t+=1) {
    draw_frame('height', t, fmt="png");
    draw_frame('level', t, fmt="png");
}
draw_frame('chi', 0, fmt="pdf");
