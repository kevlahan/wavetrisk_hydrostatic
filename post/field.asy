import palette;
import three;
import graph;
import "shared.asy" as shared;

int REAL_START = 18;
int PRIMAL = 1;
int DUAL = 2;

int S_DVORT = 9;
int S_DCHI = 10;
int S_DLEVEL = 11;
int S_TOPO = 18;
int S_CHI = 19;
int S_HEIGHT_PERT = 20;
int S_KE = 21;
int S_MASK = 22;
int S_LEVEL = 23;

real earth_radius = 1.0;
real radius=250;

pen[] pal;

int[] corners = {0, 6, 3}; /* unused, primal, dual */

struct section_t {
    real lon_min, lon_max;
    real lat_min, lat_max;
};

section_t GMT;
  GMT.lon_min = -90;
  GMT.lon_max = 90;
  GMT.lat_min = -90;
  GMT.lat_max = 90;

section_t Eurasia;
  Eurasia.lon_min = 0;
  Eurasia.lon_max = 180;
  Eurasia.lat_min = -90;
  Eurasia.lat_max = 90;

section_t America;
  America.lon_min = -180;
  America.lon_max = 0;
  America.lat_min = -90;
  America.lat_max = 90;

section_t Pacific;
  Pacific.lon_min = 90;
  Pacific.lon_max = -90;
  Pacific.lat_min = -90;
  Pacific.lat_max = 90;

pair lonlat(triple c) {
    return (atan2(c.y, c.x), asin(c.z))*180/pi;
}

transform3 rotate_section(section_t section) {
    real lon_c = (section.lon_max + section.lon_min)/2;
    if (section.lon_min > section.lon_max) {
       if (lon_c > 0) lon_c -= 180;
       else           lon_c += 180;
    }
    pair centre = (lon_c, (section.lat_max+section.lat_min)/2);
    write(centre);
    real lonangle = -centre.x;
    real latangle =  centre.y;
    triple yaxis = (0,1,0);
    triple zaxis = (0,0,1);
    return rotate(latangle, yaxis)*rotate(lonangle, zaxis);
}

void draw_hexagon(picture pic, triple[] p, real sval, bool land) {
   guide hex;
   for(int i=0; i < p.length; ++i) {
       pair q = (p[i].y, p[i].z);
       if (p[i].x < 0) return;
       hex = hex--q;
   }
   hex = hex--cycle;
   if (land) {
       fill(pic, hex, gray);
   } else {
       int pidx = round((pal.length-1)*sval);
       fill(pic, hex, pal[pidx]);
   }
}

void draw_errors(picture pic=currentpicture, triple[] errors, section_t section=Eurasia) {
        triple a_node = errors[0];
        real orig_radius = norm(a_node);
        orig_radius = sqrt(a_node.x*a_node.x + a_node.y*a_node.y + a_node.z*a_node.z);
        triple p;
	transform3 rot_transf = rotate_section(section);
        for (int k = 0; k< errors.length; ++k) {
                p = errors[k]/orig_radius;
                p = p*radius;
                p = rot_transf*p;
                if (p.x > 0) dot(pic, (p.y, p.z));
                else write("dot outside picture");
        }
}

void draw_field(picture pic=currentpicture, string variable, real[] range,
                int time=0, section_t section=Eurasia, pen[] pale=BWRainbow(), int max_level=15) {
    int grid, var;
    pal = pale;
    int s = 0;
    file in;
    real[] data;
    real[] mina;
    real[] maxa;
    if (variable == 'level') {
        grid = PRIMAL;
        var = S_LEVEL;
    } else if (variable == 'duallevel') {
        grid = DUAL;
        var = S_DLEVEL;
    } else if (variable == 'height') {
        grid = PRIMAL;
        var = S_HEIGHT_PERT;
    } else if (variable == 'vort') {
        grid = DUAL;
        var = S_DVORT;
    } else if (variable == 'chi') {
        grid = PRIMAL;
        var = S_CHI;
    } else if (variable == 'topo') {
        grid = PRIMAL;
        var = S_TOPO;
    } else if (variable == 'velocity') {
        grid = PRIMAL;
        var = S_KE;
    } else {
        write('ERROR: draw_field: invalid variable name');
        return;
    }

    real minv = range[0];
    real maxv = range[1];
    /* 
    string filename1 = "fort." + (string) (grid*100000 + time*100 + 0);
    file inmm=input(filename1).line();
    mina = inmm;
    maxa = inmm;
    while (mina.length > 0) {
        if (mina.length > var) {
            minv = min(minv, mina[var]);
            maxv = max(maxv, maxa[var]);
        }
        mina = inmm;
        maxa = inmm;
    }
    close(inmm);
    */
    real actualmin = 1e16;
    real actualmax = -1e16;

    triple[] p;
    transform3 rot_transf = rotate_section(section);
    for (int l = 1; l <= max_level; ++l) {
        string filename = "fort." + (string) (grid*100000 + time*100 + l);
        in=input(filename, check=false).line();
        if (error(in)) continue;
        write(filename);

        data = in;
        triple a_node = data[s:s+3];
        real orig_radius = norm(a_node);
        orig_radius = 1.0001*sqrt(a_node.x*a_node.x + a_node.y*a_node.y + a_node.z*a_node.z);
        while (data.length > 0) {
            bool inside = false;
            for (int k = 0; k<corners[grid]; ++k) {
                p[k] = data[s+3*k:s+3*k+3]/orig_radius;
                if (!inside) {
                    pair ll = lonlat(p[k]);
                    inside = section.lat_min < ll.y && ll.y < section.lat_max;
                    if (section.lon_min < section.lon_max)
                        inside = (section.lon_min < ll.x && ll.x < section.lon_max) && inside;
                    else
                        inside = (section.lon_min < ll.x || ll.x < section.lon_max) && inside;
                }
                p[k] = p[k]*radius;
                p[k] = rot_transf*p[k];
            }
            if (inside) {
                real value = data[var];
                if (variable == 'velocity') value = sqrt(2*value); /* compute velocity magnitude from kinetic energy */
                if (var == S_LEVEL && l < 8) value = l;
		actualmin = min(value, actualmin);
		actualmax = max(value, actualmax);
                value = min(maxv, max(minv,value));
		if (grid == PRIMAL)
                	draw_hexagon(pic, p, (value-minv)/(maxv-minv), data[S_CHI]>0.5);
		if (grid == DUAL)
                	draw_hexagon(pic, p, (value-minv)/(maxv-minv), data[S_DCHI]>0.5);
            }
            data = in;
        }
        close(in);
        write((actualmin, actualmax));
    }
    draw(pic, Circle((0,0), radius));
}
