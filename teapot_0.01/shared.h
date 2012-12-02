#include <math.h>

struct field {
	double temperature;
	double rotor_v;
	double v_potential;
	};

struct field **space_array;
struct field **array_b;
int cp_array(struct field **space_array, struct field **array_b);
int edge_condition_o (struct field **array_b, long int info);
int edge_condition_t (struct field **array_b, long int info);
int temperature_change (struct field **space_array, struct field **array_b);
int oxygen_change (struct field **space_array, struct field **array_b);
struct field **create_space_array (long int info);
struct field **create_array_b (void);
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

