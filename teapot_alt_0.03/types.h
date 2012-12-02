struct Cell
{
  double w;
  double theta;
  double psi;
};

struct Event
{
  int type;                  // type=1  -  exit
  			     // type=2  -  save dt
  double param;
  double T;
};
