import palette;
import three;
import graph;
import "shared.asy" as shared;

int REAL_START = 18;
int PRIMAL = 1;
int DUAL = 2;

int S_DIVU = 1;
int S_DDIVU = 3;
int S_DVORT = 4;
int S_PHI = 5;
int S_VORT = 9;
int S_TOPO = 18;
int S_CHI = 19;
int S_HEIGHT_PERT = 20;
int S_DHEIGHT = 21;
int S_MASK = 22;
int S_LEVEL = 23;

real earth_radius = 1.0;
real radius=250;
real Hdim=3.344175893265152e+03;

pen[] pal;

int[] corners = {0, 6, 3}; /* unused, primal, dual */

struct section_t {
    real lon_min, lon_max;
    real lat_min, lat_max;
};

/*section_t Eurasia = {0, 180, -90, 90};*/
section_t Eurasia;
  Eurasia.lon_min = 0;
  Eurasia.lon_max = 180;
  Eurasia.lat_min = -90;
  Eurasia.lat_max = 90;

pair lonlat(triple c) {
    return (atan2(c.y, c.x), asin(c.z))*180/pi;
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
        pair centre = (section.lon_max+section.lon_min, section.lat_max+section.lat_min)/2;
        real lonangle = -centre.x;
        real latangle =  centre.y;
        triple yaxis = (0,1,0);
        triple zaxis = (0,0,1);
        for (int k = 0; k< errors.length; ++k) {
                p = errors[k]/orig_radius;
                p = p*radius;
                p = rotate(latangle, yaxis)*rotate(lonangle, zaxis)*p;
                dot(pic, (p.y, p.z));
        }
        if (p.x > 0) write("dot outside picture");
}

void draw_field(picture pic=currentpicture, string variable, 
                int time=0, int level=0, section_t section=Eurasia, pen[] pale=BWRainbow()) {
    int grid, var, min_level, max_level;
    pal = pale;
    int s = 0;
    file in;
    real[] data;
    real[] mina;
    real[] maxa;
    triple[] p;
    pair centre = (section.lon_max+section.lon_min, section.lat_max+section.lat_min)/2;
    if (variable == 'level') {
        grid = PRIMAL;
        var = S_LEVEL;
    } else if (variable == 'height') {
        grid = PRIMAL;
        var = S_HEIGHT_PERT;
    } else if (variable == 'trend') {
        grid = PRIMAL;
        var = S_DHEIGHT;
    } else if (variable == 'chi') {
        grid = PRIMAL;
        var = S_CHI;
    } else if (variable == 'topo') {
        grid = PRIMAL;
        var = S_TOPO;
    } else if (variable == 'vort') {
        grid = DUAL;
        var = S_VORT;
    } else {
        write('ERROR: draw_field: invalid variable name');
        return;
    }
    if (level==0) {
        min_level = 7;
        max_level = 12;
    } else {
        min_level = level;
        max_level = level;
    }
    write(centre);
    real lonangle = -centre.x;
    real latangle =  centre.y;
    triple yaxis = (0,1,0);
    triple zaxis = (0,0,1);
    write(var);
    real minv, maxv;
    string filename;
    if (variable == 'vort') {
        minv = -1.0e-1;
        maxv =  1.0e-1;
    }

    if (variable == 'level') {
        minv = min_level;
        maxv = max_level;
    }

    for (int l = min_level; l <= 8; ++l) {
        filename = "fort." + (string) (grid*100000 + time*100 + l);
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
                    inside = section.lon_min < ll.x && ll.x < section.lon_max && 
                             section.lat_min < ll.y && ll.y < section.lat_max;
                }
                p[k] = p[k]*radius;
                p[k] = rotate(latangle, yaxis)*rotate(lonangle, zaxis)*p[k];
            }
            if (inside) {
                real value;
                if (var >= 9) value = min(maxv, max(minv,data[var]));
                if (var == S_LEVEL && l < 8) value = l;

                draw_hexagon(pic, p, (value-minv)/(maxv-minv), 0==1);
            }
            data = in;
        }
        close(in);
        write((minv, maxv));
    }
    draw(pic, Circle((0,0), radius));
}
