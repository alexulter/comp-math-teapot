#include <math.h>

struct field {
	double omega;
	double psi;
	double teta;
	double J11; 
	double J12; 
	double J21;
	double J_omega;
	double J_teta;
	double T_omega;
	double T_teta;
	double F_omega;
	double F_teta;
	double temperature;
	double oxygen;
	int life_exist;
	int t_resistance;
	int o_resistance;
	double health;
	int gender;

	};

struct field **space_array;
struct field **array_b;
struct field **array_b1;

extern unsigned int numxcells;
extern unsigned int numycells;
extern void (*scheme)();
extern void (*schemefinish)();
extern void (*bounds)();
extern double T;
extern double dt;
extern double dt0;
extern struct Event *events;
extern double numevents;
extern double nextevent;
extern double h_x;
extern double h_z;
extern double K;
extern double beta;
extern long step;
extern const double pi;
extern double lamb1;
extern double lamb2;
extern double Lamb1;
extern double Lamb2;

int progonka (double** koef_array , double* psi , double* f , const int n , const int n_0);
int cycle_progonka (double** koef_array , double* psi , double* f , const int n , const int n_0);
int init_systems (void);
int calc_omega_and_teta (struct field **space_array, struct field **array_b);
int cp_array(struct field **space_array, struct field **array_b);
int cp_array_2(struct field **array_b, struct field **array_b1);
int edge_condition_o (struct field **array_b, long int info);
int edge_condition_t (struct field **array_b, long int info);
int temperature_change (struct field **space_array, struct field **array_b);
int oxygen_change (struct field **space_array, struct field **array_b);
struct field **create_space_array (long int info);
struct field **create_array_b (void);
struct field **create_array_b1 (void);
int status (struct field **space_array, struct field **array_b, long int info);
int status_2 (struct field **space_array, struct field **array_b);
int body_of_program (struct field **space_array, struct field **array_b, long int info);
int red_color_t (double t);
int green_color_t (double t);
int blue_color_t (double t);
int red_color_o (double o);
int green_color_o (double o);
int blue_color_o (double o);
long int reading_from_file(void);
int health_color(double h);
int ec_smpl_t(struct field **array_b, long int info);
int ec_smpl_o(struct field **array_b, long int info);




#define vert_size 100
#define hor_size 100
#define frand() ((double) rand() / (RAND_MAX + 1.0))
#define TEMP_MAX 50
#define TEMP_MIN -50
#define OXY_MAX 100
#define OXY_MIN 0

#define ec_function_t1() ((TEMP_MAX-TEMP_MIN)*frand()+TEMP_MIN)
#define ec_function_o1() ((OXY_MAX)*frand())
//#define ec_function_t_sin(arg) ((TEMP_MAX-TEMP_MIN)*(sin(#arg)+1)/2+TEMP_MIN)


#define SIZE_OF_CELL 4
#define SPEED 10000000

