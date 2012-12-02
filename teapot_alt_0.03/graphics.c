#include <gtk/gtk.h>
#include <glade/glade.h>
#include "shared.h"
#include "init.h"
#include "calc.h"
#include "types.h"
#include "defs.h"

static GdkPixmap *pixmap = NULL;
//static GdkPixmap *pixmap_t = NULL;
static gint show_life = 0;
static gint abc = 0;
static gint screen_mode = 1;
static gint LRUN = 0;
// draw one black rectangle
	static void 
draw_brush (GtkWidget *widget, gint x, gint y)
{
	/*GdkGC *gc;
	GdkColor color;
	color.red = 0;
	color.blue = 0;
	color.green = 0;
	gc = gdk_gc_new(widget->window);
	gdk_gc_set_rgb_fg_color (gc, &color);*/
	GdkRectangle update_rect;
	update_rect.x = x - SIZE_OF_CELL / 2;
	update_rect.y = y - SIZE_OF_CELL / 2;
	update_rect.width = SIZE_OF_CELL;
	update_rect.height = SIZE_OF_CELL;
	gdk_draw_rectangle 	(pixmap,
			widget->style->black_gc,
			TRUE,
			update_rect.x, update_rect.y,
			update_rect.width, update_rect.height);
	gdk_window_invalidate_rect(widget->window, &update_rect, FALSE);
}

//takes cordinates, draws rectangle, makes new life in space_array
	static gboolean 
da_button_press_event	(GtkWidget *widget,
		GdkEventButton *event)
{
	int i, j;
	if (event->button == 1 && pixmap != NULL)
	{
		i = (event->x / SIZE_OF_CELL) + 2;
		j = (event->y / SIZE_OF_CELL) + 2;
		if (((j + 1) <= vert_size - 1) &&
				((i + 1) <= hor_size - 1))
		{
			space_array[i][j].life_exist = 1;
			space_array[i][j].health = 100;
			i = (i-1.5)*SIZE_OF_CELL;
			j = (j-1.5)*SIZE_OF_CELL;
			draw_brush(widget, i, j);
		}
	}

	return TRUE;
}

// redraw part of screen
	static gboolean 
da_expose_event	(GtkWidget *widget,
		GdkEventExpose *event)
{	
	GdkGC *gc;
	GdkColor color;
	color.red=65535;
	color.blue=65535;
	color.green=65535;
	gc = gdk_gc_new(widget->window);
	gdk_gc_set_rgb_fg_color(gc, &color);
	gdk_draw_drawable(widget->window,
			widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
			pixmap,
			event->area.x, event->area.y,
			event->area.x, event->area.y,
			event->area.width, event->area.height);
	return FALSE;
}

// clears pixmap
	static gboolean
da_clear (GtkWidget *da)
{
	gint w=0;
	gint h=0;
	if (pixmap) 
		g_object_unref(pixmap);
	gdk_drawable_get_size(da->window,&w,&h);
	pixmap = gdk_pixmap_new(da->window,w,h,-1);
	       gdk_draw_rectangle (pixmap,
			       da->style->white_gc,
			       TRUE,
			       0, 0,
			       da->allocation.width,
			       da->allocation.height);
	return TRUE;
}


//exits
	static gboolean
to_exit (void)

{
	Finish();
	exit(-2);
	return TRUE;
}

// starts life in continious loop
	static gboolean
life_run (void)
{
	LRUN = 1;
	return TRUE;
}

//stops loop
	static gboolean
life_stop (void)
{
	LRUN = 0;
	return TRUE;
}


// helpfunction
static gboolean 
life (GtkWidget *da, struct field **space_array)
{    
	//long int info;
	//info = reading_from_file();
	LRUN = 3;
	//array_b = create_array_b();
	//abc = 1;
	//abc = screen_mode;
	//body_of_program(space_array, array_b,info);
	return TRUE;
}

//shows temperature
static gboolean 
temp (GtkWidget *da)
{

	//array_b = create_array_b();
	//body_of_program(space_array, array_b);
	screen_mode = 2;
	if (show_life != screen_mode) show_life = screen_mode;
	else show_life = 0;
	abc = screen_mode;
	return TRUE;
}


//shows oxygen
static gboolean 
oxy (GtkWidget *da)
{

	//array_b = create_array_b();
	//body_of_program(space_array, array_b);
	screen_mode = 3;
	if (show_life != screen_mode) show_life = screen_mode;
	else show_life = 0;
	abc = screen_mode;
	return TRUE;
}

  
//shows life
static gboolean 
to_life (GtkWidget *da)
{

	//array_b = create_array_b();
	//body_of_program(space_array, array_b);
	screen_mode = 1;
	if (show_life != screen_mode) show_life = screen_mode ;
	else show_life = 0;
	abc = screen_mode;
	return TRUE;
}


/*shows field dependindg on what is needed to show*/
static void  
show_field( GtkWidget *da)
{	
	gint i,j;
	GdkGC *gc; 
	GdkColor color;
	gc = gdk_gc_new (da->window);	
	if (abc==1)
	{	
		for (j=2;j<=hor_size-2;j++)
		{

			for (i=2;i<=vert_size-2;i++)
			{

				//printf("%i",space_array[5][5].life_exist);

				if (space_array[i][j].life_exist == 1)
				{     
					color.red=health_color(space_array[i][j].health);
					color.blue=health_color(space_array[i][j].health);
					color.green=health_color(space_array[i][j].health);
					gdk_gc_set_rgb_fg_color(gc,&color);        
					//printf("%i \n",space_array[5][5]);
					gdk_draw_rectangle (pixmap,
							gc,
							TRUE,
							(i-2)*SIZE_OF_CELL,
							(j-2)*SIZE_OF_CELL,
							SIZE_OF_CELL,
							SIZE_OF_CELL);
				}

				else 
				{	
					color.red=65535;
					color.blue=65535;
					color.green=65535;
					gdk_gc_set_rgb_fg_color(gc,&color);

					gdk_draw_rectangle (pixmap,
							gc,
							TRUE,
							(i-2)*SIZE_OF_CELL,
							(j-2)*SIZE_OF_CELL,
							SIZE_OF_CELL,
							SIZE_OF_CELL);
				}



			}
		}
		/*gdk_draw_drawable(da->window,
		gc,//da->style->fg_gc[GTK_WIDGET_STATE (da)],
		pixmap, 0, 0, 0, 0,
		da->allocation.width, da->allocation.height);*/
	}		
	if(abc == 2)
	{
		for (j=2;j<=hor_size-2;j++)
		{

			for (i=2;i<=vert_size-2;i++)
			{

				//printf("%i",space_array[5][5].life_exist);

				if ((space_array[i][j].life_exist == 1)&&(show_life != 0))
				{
				color.red = health_color(space_array[i][j].health);
				color.blue = health_color(space_array[i][j].health);
				color.green = health_color(space_array[i][j].health);
				gdk_gc_set_rgb_fg_color(gc,&color);
				//printf("%i \n",space_array[5][5]);
				gdk_draw_rectangle (pixmap,
				gc,
				TRUE,
				(i-2)*SIZE_OF_CELL,
				(j-2)*SIZE_OF_CELL,
				SIZE_OF_CELL,
				SIZE_OF_CELL);
				}
				 else	
				{	
				color.red = red_color_t(space_array[i][j].
						omega + TEMP_MAX);
				color.blue = blue_color_t(space_array[i][j].
						omega + TEMP_MAX);
				color.green = green_color_t(space_array[i][j].
						omega+ TEMP_MAX);
				gdk_gc_set_rgb_fg_color(gc,&color);
				gdk_draw_rectangle (pixmap,
						gc,
						TRUE,
						(i-2)*SIZE_OF_CELL,
						(j-2)*SIZE_OF_CELL,
						SIZE_OF_CELL,
						SIZE_OF_CELL);
				}

			}
		}
	}
	if(abc == 3)
	{
		for (j=2;j<=hor_size-2;j++)
		{

			for (i=2;i<=vert_size-2;i++)
			{

				//printf("%i",space_array[5][5].life_exist);
//			Сносить тут !!!
				if ((space_array[i][j].life_exist == 1)&&(show_life != 0))
				{	
				color.red =  health_color(space_array[i][j].health);
				color.blue =  health_color(space_array[i][j].health);
				color.green =  health_color(space_array[i][j].health);
				gdk_gc_set_rgb_fg_color(gc,&color);
				//printf("%i \n",space_array[5][5]);
				gdk_draw_rectangle (pixmap,
				gc,
				TRUE,
				(i-2)*SIZE_OF_CELL,
				(j-2)*SIZE_OF_CELL,
				SIZE_OF_CELL,
				SIZE_OF_CELL);
				}
				else
				{
			//		тут ьыл КИСЛОРОД
				color.red = red_color_o(space_array[i][j].
						psi);
				color.blue = blue_color_o(space_array[i][j].
						psi);
				color.green = green_color_o(space_array[i][j].
						psi);
				gdk_gc_set_rgb_fg_color(gc,&color);
				gdk_draw_rectangle (pixmap,
						gc,
						TRUE,
						(i-2)*SIZE_OF_CELL,
						(j-2)*SIZE_OF_CELL,
						SIZE_OF_CELL,
						SIZE_OF_CELL);
				}

			}
		}
	}	
}			

static gboolean
new(struct field **space_array )
{
	LRUN = 2;

	/*info = reading_from_file();
	space_array =create_space_array(info);//(space_array)
	

	screen_mode = 1;
        abc = screen_mode;*/
	return TRUE;
}



int main (int argc,char **argv )
{

	//	gpointer gp;
	//	gp = (void*)space_array;
	GtkWidget *window;	
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *button_quit;
	GtkWidget *button_run;
	GtkWidget *button_new;
	GtkWidget *button_t;
	GtkWidget *button_o;
	GtkWidget *da;
	GtkWidget *button_step;
	GtkWidget *button_stop;
	GtkWidget *button_l;
	gtk_init(&argc, &argv);
	glade_init();
	long int info;
	int l,k;

	info = reading_from_file();
	space_array = create_space_array(info);
	//
	Init();
	
	
	

	//
	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	g_signal_connect (G_OBJECT (window), "destroy",
			G_CALLBACK (to_exit), NULL);
	gtk_container_set_border_width (GTK_CONTAINER (window), 10);
	hbox = gtk_hbox_new(FALSE,5);
	vbox = gtk_vbox_new(FALSE,0);


	da= gtk_drawing_area_new();
	gtk_widget_set_size_request (da, 
				(hor_size-2)*SIZE_OF_CELL, 
				(vert_size-2)*SIZE_OF_CELL);


	button_quit = gtk_button_new_with_label ("Quit");
	button_step = gtk_button_new_with_label ("Step");
	button_stop = gtk_button_new_with_label ("Stop");
	button_run = gtk_button_new_with_label ("Run");
	button_new = gtk_button_new_with_label ("New");
	button_l = gtk_button_new_with_label ("Show life");
	button_t = gtk_button_new_with_label ("Show temperature");
	button_o = gtk_button_new_with_label ("Show oxygen");


	g_signal_connect (G_OBJECT(button_quit),
			"clicked",
			G_CALLBACK(to_exit),
			NULL);
	g_signal_connect (G_OBJECT(button_step), 
			"clicked", 
			G_CALLBACK(life), 
			space_array);
	g_signal_connect (G_OBJECT(da),
			"button_press_event",
			G_CALLBACK(da_button_press_event),
			NULL);
	g_signal_connect (G_OBJECT(da),
			"configure_event",
			G_CALLBACK(da_clear), 
			NULL);
	g_signal_connect (G_OBJECT(da),
			"expose_event",
			G_CALLBACK(da_expose_event),
			NULL);
	g_signal_connect (G_OBJECT(button_run), 
			"clicked", 
			G_CALLBACK(life_run), 
			space_array);
	g_signal_connect (G_OBJECT(button_stop), 
			"clicked", 
			G_CALLBACK(life_stop), 
			space_array);
	g_signal_connect (G_OBJECT(button_t), 
			"clicked", 
			G_CALLBACK(temp), 
			NULL);
	g_signal_connect (G_OBJECT(button_l), 
			"clicked", 
			G_CALLBACK(to_life), 
			space_array);
	g_signal_connect (G_OBJECT(button_o), 
			"clicked", 
			G_CALLBACK(oxy), 
			space_array);
	g_signal_connect (G_OBJECT(button_new), 
			"clicked", 
			G_CALLBACK(new), 
			da);

	gtk_widget_set_events (da,GDK_EXPOSURE_MASK
			| GDK_LEAVE_NOTIFY_MASK
			| GDK_BUTTON_PRESS_MASK);




	gtk_box_pack_start(GTK_BOX (vbox), button_new,FALSE,FALSE,0);
	gtk_box_pack_start(GTK_BOX (vbox), button_step,FALSE,FALSE,0);
	gtk_box_pack_start(GTK_BOX (vbox), button_run,FALSE,FALSE,0);
	gtk_box_pack_start(GTK_BOX (vbox), button_stop,FALSE,FALSE,0);
	gtk_box_pack_start(GTK_BOX (vbox), button_l,FALSE,FALSE,0);
	gtk_box_pack_start(GTK_BOX (vbox), button_t,FALSE,FALSE,0);
	gtk_box_pack_start(GTK_BOX (vbox), button_o,FALSE,FALSE,0);
	gtk_box_pack_start(GTK_BOX (vbox), button_quit,FALSE,FALSE,0);
	gtk_box_pack_end(GTK_BOX(hbox),vbox,FALSE,FALSE,0);
	gtk_box_pack_end(GTK_BOX(hbox),da,FALSE,FALSE,0);
	gtk_container_add (GTK_CONTAINER (window), hbox);


	gtk_widget_show_all(window);



	gboolean br;
	int i=0,j=0;
	gint w=0,h=0;
	GdkRectangle rect;
	GdkGC *gc;
	GdkColor color;
	color.red = 65535;
	color.blue = 65535;
	color.green = 65535;
	gc = gdk_gc_new (da->window);
	gdk_gc_set_rgb_fg_color (gc, &color);
	rect.x = 0;
	rect.y = 0;
	gdk_drawable_get_size(da->window,&w,&h);

	rect.width = w;// da->allocation.width;
	rect.height = h;// da->allocation.height;

	array_b = create_array_b();
	array_b1 = create_array_b1();
	LRUN = 0;
	while (TRUE)
	{
		while (gtk_events_pending ())
			gtk_main_iteration ();

		if (LRUN == 2)
		{
		space_array =create_space_array(info);//(space_array)
		screen_mode = 1;
        	abc = screen_mode;
		LRUN = 0;
		};

		if ((LRUN == 1)||(LRUN == 3))
		{
			for (k=0; k<SPEED;k++);	
			abc = screen_mode;
			body_of_program(space_array, array_b,info);
			
			if (LRUN == 3)	LRUN = 0;
		}

		if (abc!=0)
		{
			
			
			show_field (da);
			
			gdk_draw_drawable(da->window,
					gc,
					pixmap, 0, 0, 0, 0,
					da->allocation.width, 
					da->allocation.height);
			//gtk_widget_queue_draw_area(da,
			//	da->allocation.x,
			//	da->allocation.y,
			//	da->allocation.width,
			//	da->allocation.height); 
			//gtk_widget_show(da);
			abc=0;
		}
	
	      
		//while(gtk_events_pending())
		//gtk_main_iteration();

	} 	
	
	
	return 0;
}

