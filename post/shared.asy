pair operator cast(real[] x) {
  return (x[0],x[1]);
}

triple operator cast(real[] x) {
  return (x[0],x[1],x[2]);
}

real[] operator cast(triple x) {
  return new real[] {x.x,x.y,x.z};
}

int[][] operator -(int[][] a, int s) {
  for (int j=0; j<a.length; ++j) {
    for (int i=0; i<a[j].length; ++i) {
      a[j][i] = a[j][i] - s;
    }
  }
  return a;
}

