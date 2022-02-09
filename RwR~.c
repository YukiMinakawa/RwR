/*
    The Pd object of "Room within a room"

    The argument after RwR~ is the number of the windows.

    Inlets:
    1 - audio source
    2 - list of the windows in radius and angles
    3 - the energy absorbed by the wall in the reflection in decibels
    4 - X coordinate in meters of the audio source
    5 - Y coordinate in meters of the audio source
    6 - outer room size in meters

    Outlets:
    Signals
*/

#include "m_pd.h"
#include "./RwR~.h"
#include <string.h>
#include <math.h>
#define MAX_OUTLET 32
#define c 341
#define block_iter 5000

// Define a new "class"
static t_class *RwR_tilde_class;

typedef struct _RwR_tilde {
    t_object  x_obj; /// This object
    t_float f; /// dummy variable to deal with floats sent to the signal inlets

    t_float x_attenuation; // wall attenuation (dB)
    t_float x_soundsource_x; // the coordinate x of the sound source
    t_float x_soundsource_y; // the coordinate y of the sound source

    t_float x_innerroom; // length of the egde of the inner room (m)
    t_float x_outerroom; // length of the egde of the outer room (m)

    t_int x_n; // the number of the windows
    t_float *x_windows; // the coordinates of the windows, size = #window * 2

    t_inlet *x_in2; // it will get the coordinates of the windows
    t_inlet *x_in3; // gets wall attenuation
    t_inlet *x_in4; // gets the coordinate x of the sound source
    t_inlet *x_in5; // gets the coordinate y of the sound source
    t_inlet *x_in6; // gets the outer room size

    t_outlet **x_outputs; // outputs
    t_int *x_outputmapping; // output mapping array

    t_float *x_buffer; // arrays containing the delayed signals
    t_int x_bufferinmapping; // bufferin mapping array
    t_int x_bufferoutmapping; // bufferout mapping array

} t_RwR_tilde;

// This method judges the coordinates of the sound source and the windows are correct or not.
// And reset the flags that show whether delays are finished or not.
void RwR_tilde_bang(t_RwR_tilde *x){
    int i;

    post("in the method of RwR_tilde_bang");

    if(fabsf(x->x_soundsource_x - x->x_outerroom/2) > x->x_outerroom/2){
        post("soundsource_x is out of the outer room");
        x->x_soundsource_x = 0;
    }

    if(fabsf(x->x_soundsource_y - x->x_outerroom/2) > x->x_outerroom/2){
        post("soundsource_y is out of the outer room");
        x->x_soundsource_y = 0;
    }

    for(i=0; i<x->x_n; i++)
    {
        post("w%d: (%f, %f)", i, x->x_windows[i*2], x->x_windows[(i*2)+1]);
    }
    post("ss: (%.2f, %.2f)", x->x_soundsource_x, x->x_soundsource_y);
    post("a = %.2f", x->x_attenuation);
}

// This method gets the coordinates of the windows.
void RwR_tilde_winlist(t_RwR_tilde *x, t_symbol *s, int argc, t_atom *argv){
    (void)s;
    int i;

    float radius;
    float angles[32];

    // Polar coordinates
    if(argc != x->x_n+1){
        post("Unmatch with arg of RwR~ and args of 2nd inlet");
        for(i=0; i<x->x_n*2; i++){
            x->x_windows[i] = x->x_outerroom/2.0;
        }
    }
    else{
        radius = atom_getfloat(argv);
        for(i=0; i<x->x_n; i++){
            angles[i] = (-1 * atom_getfloat(argv+i+1) + 90) * (M_PI/180);
        }
        for(i=0; i<x->x_n; i++){
            x->x_windows[i*2] = radius * cos(angles[i]) + (x->x_outerroom/2);
            x->x_windows[i*2+1] = radius * sin(angles[i]) + (x->x_outerroom/2);
        }
    }
}

// This method works when the DSP is running.
t_int *RwR_tilde_perform(t_int *w){
    t_float *in, *out, *acopy;

    t_RwR_tilde *x = (t_RwR_tilde *)(w[1]);
    int n = (int)(w[2]); // length of the signal vector
    in = (t_float *)(w[3]); // input signal

    int i,j;

    Vector2 firstset[x->x_n];
    Vector2 finalset[x->x_n];
    Vector2 normal[x->x_n];
    Vector2 SoundSource = {0, 0};
    Vector2 ImageWindow = {0, 0};
    Vector2 reflect = {0, 0};
    float dist = 0.0;
    float min_dist = FLT_MAX;

    PointList hull_all = {{}, 0};
    PointList hull_direct = {{}, 0};
    PointList hull_1ref = {{}, 0};
    PointList hull_death = {{}, 0};
    PointList hull_buff = {{}, 0};
    DG hull_direct_DG[x->x_n];
    DG hull_1ref_DG[x->x_n];
    DG finalset_DG[x->x_n];

    float sr = sys_getsr();
    int smpdelay[x->x_n];
    int isSSinner = 0;
    int isDirectflag = 0;

    for(i=0; i<x->x_n; i++){
        firstset[i].x = x->x_windows[i*2];
        firstset[i].y = x->x_windows[i*2+1];
    }
    SoundSource.x = x->x_soundsource_x;
    SoundSource.y = x->x_soundsource_y;

    if(fabsf(SoundSource.x - x->x_outerroom/2) <= x->x_innerroom/2){
        if(fabsf(SoundSource.y - x->x_outerroom/2) <= x->x_innerroom/2){
            isSSinner = 1; // 1 means SoundSource is in the innerroom.
        }
    }

    buildHull(firstset, x->x_n, &hull_all);
    sortPoint(hull_all.arr, hull_all.head, x->x_outerroom/2.0);

    // Calculating normals of hull
    for(i=0; i<hull_all.head; i++){
        normal[i] = makeNormal(hull_all.arr[i], hull_all.arr[(i+1)%hull_all.head]);
    }

    // Determine if there is a direct path between the window and the sound source.
    for(i=0; i<hull_all.head; i++){
        // Direct path: isDirect returns 1.
        // Path crosses the convex hull: isDirect returns 0.
        isDirectflag = isDirect(&hull_all, normal, SoundSource, hull_all.arr[i]);

        if(isDirectflag == 1){
            hull_direct.arr[hull_direct.head] = hull_all.arr[i];
            dist = distance(SoundSource, hull_all.arr[i]);

            // delay and attenuate by the air
            getDG(&hull_direct_DG[hull_direct.head], dist, (float)c);
            hull_direct.head++;
        }
        else if(isDirectflag == 0){
            hull_buff.arr[hull_buff.head] = hull_all.arr[i];
            hull_buff.head++;
        }
        else if(isDirectflag == 2) post("isDirect returns 2");
        else if(isDirectflag == 3) post("isDirect returns 3");
        else post("isDirect error");
    }

    // Make image windows and find nearest path from four image windows.
    for(i=0; i<hull_buff.head; i++){
        min_dist = FLT_MAX;

        // Reflect on the left wall
        ImageWindow.x = hull_buff.arr[i].x * -1;
        ImageWindow.y = hull_buff.arr[i].y;
        dist = distance(SoundSource, ImageWindow);

        // The path from the sound source to the window doesn't cross the convex hull: isDirect returns 2.
        if(isDirect(&hull_all, normal, SoundSource, ImageWindow) == 2 && dist < min_dist){
            reflect.x = 0.0;
            reflect.y = (ImageWindow.y - SoundSource.y)/(ImageWindow.x - SoundSource.x)*(reflect.x - SoundSource.x)+SoundSource.y;
            if(isDirect(&hull_all, normal, reflect, hull_buff.arr[i]) == 1){
                min_dist = dist;
            }
        }

        // Reflect on the right wall
        ImageWindow.x += 2 * x->x_outerroom;
        dist = distance(SoundSource, ImageWindow);

        if(isDirect(&hull_all, normal, SoundSource, ImageWindow) == 2 && dist < min_dist){
            reflect.x = x->x_outerroom;
            reflect.y = (ImageWindow.y - SoundSource.y)/(ImageWindow.x - SoundSource.x)*(reflect.x - SoundSource.x)+SoundSource.y;
            if(isDirect(&hull_all, normal, reflect, hull_buff.arr[i]) == 1){
                min_dist = dist;
            }
        }

        // Reflect on the bottom wall
        ImageWindow.x = hull_buff.arr[i].x;
        ImageWindow.y = hull_buff.arr[i].y * -1;
        dist = distance(SoundSource, ImageWindow);

        if(isDirect(&hull_all, normal, SoundSource, ImageWindow) == 2 && dist < min_dist){
            reflect.y = 0.0;
            reflect.x = (reflect.y - SoundSource.y)*(ImageWindow.x - SoundSource.x)/(ImageWindow.y - SoundSource.y)+SoundSource.x;
            if(isDirect(&hull_all, normal, reflect, hull_buff.arr[i]) == 1){
                min_dist = dist;
            }
        }

        // Reflect on the top wall
        ImageWindow.y += 2 * x->x_outerroom;
        dist = distance(SoundSource, ImageWindow);

        if(isDirect(&hull_all, normal, SoundSource, ImageWindow) == 2 && dist < min_dist){
            reflect.y = x->x_outerroom;
            reflect.x = (reflect.y - SoundSource.y)*(ImageWindow.x - SoundSource.x)/(ImageWindow.y - SoundSource.y)+SoundSource.x;
            if(isDirect(&hull_all, normal, reflect, hull_buff.arr[i]) == 1){
                min_dist = dist;
            }
        }

        // Store the coordinate of the nearest 1st-order-reflect path.
        if(min_dist != FLT_MAX){
            hull_1ref.arr[hull_1ref.head] = hull_buff.arr[i];

            // delay and attenuate by the air
            getDG(&hull_1ref_DG[hull_1ref.head], min_dist, (float)c);

            // attenuate by the wall
            hull_1ref_DG[hull_1ref.head].attenuate_factor *= pow(10, -1*x->x_attenuation/20);

            hull_1ref.head++;
        }
        // If there is not the 1st-order-reflect path, store it in the death.
        else{
            hull_death.arr[hull_death.head] = hull_buff.arr[i];
            hull_death.head++;
        }
    }

    // Register delays and attenuate_factors to the windows.
    // At first, mute all windows.
    for(i=0;i<x->x_n;i++){
        finalset[i] = firstset[i];
        finalset_DG[i].delay = -1;
        finalset_DG[i].attenuate_factor = 0;
    }

    // If the sound source is out of the inner room, organize final set from parial set.
    for(i=0;i<hull_direct.head;i++){
        for(j=0;j<x->x_n;j++){
            if(fabsf(hull_direct.arr[i].x - firstset[j].x) <= FLT_EPSILON * fmaxf(1.f, fmaxf(fabsf(hull_direct.arr[i].x), fabsf(firstset[j].x))) && fabsf(hull_direct.arr[i].y - firstset[j].y) <= FLT_EPSILON * fmaxf(1.f, fmaxf(fabsf(hull_direct.arr[i].y), fabsf(firstset[j].y)))){
                finalset[j] = hull_direct.arr[i];
                finalset_DG[j] = hull_direct_DG[i];
            }
        }
    }

    for(i=0;i<hull_1ref.head;i++){
        for(j=0;j<x->x_n;j++){
            if(fabsf(hull_1ref.arr[i].x - firstset[j].x) <= FLT_EPSILON * fmaxf(1.f, fmaxf(fabsf(hull_1ref.arr[i].x), fabsf(firstset[j].x))) && fabsf(hull_1ref.arr[i].y - firstset[j].y) <= FLT_EPSILON * fmaxf(1.f, fmaxf(fabsf(hull_1ref.arr[i].y), fabsf(firstset[j].y)))){
                finalset[j] = hull_1ref.arr[i];
                finalset_DG[j] = hull_1ref_DG[i];
            }
        }
    }

    /*  ========================
            Output Section
        ========================
    */

    // copy input audio block
    acopy = (t_float*) getbytes(n*sizeof(t_float));
    if (!acopy) {
        post( "RwR~: cannot allocate audio copy block" );
        return 0;
    }

    memcpy(acopy, in, n*sizeof(t_float));

    for(j=0;j<n;j++){
        x->x_buffer[x->x_bufferinmapping] = *(acopy+j);
        x->x_bufferinmapping = (x->x_bufferinmapping+1)%(block_iter*n);
    }

    for (i=0; i<x->x_n; i++) {
        out = (t_float *)(w[x->x_outputmapping[i]+4]);
        smpdelay[i] = (int)(finalset_DG[i].delay * sr);

        for(j=0; j<n; j++){
            *out = fminf(x->x_buffer[(x->x_bufferinmapping+j-smpdelay[i]-n)%(block_iter*n)] * finalset_DG[i].attenuate_factor, 1);
            if(*out == 1){
                post("clip: %10.3f", x->x_buffer[(x->x_bufferinmapping+j-smpdelay[i]-n)%(block_iter*n)] * finalset_DG[i].attenuate_factor);

            }
            out++;
        }
    }
    if (acopy) freebytes(acopy, n*sizeof(t_float));

    return (w+x->x_n+4); // 4 means dsp_addvの引数 nargs, dspargs[0~2]
}

// This method regists some pointers and perform method to the DSP tree.
void RwR_tilde_dsp(t_RwR_tilde *x, t_signal **sp){
    t_int *dspargs, fi, nargs;
    dspargs = (t_int*)getbytes((x->x_n+3)*sizeof(t_int));
    dspargs[0] = (t_int)x;
    dspargs[1] = (t_int)sp[0]->s_n; // length of the signal vector
    dspargs[2] = (t_int)sp[0]->s_vec; // pointer to the input signal vector

    nargs = 3;
    for (fi=0; fi<x->x_n; fi++) {
        dspargs[3+fi] = (t_int)sp[fi+1]->s_vec; // pointer to the output signal vector
        nargs++;
    }
    dsp_addv(RwR_tilde_perform, nargs, dspargs);
    if (dspargs) {
        freebytes(dspargs, (x->x_n+3)*sizeof(t_int));
    }
}

// Constructor
void *RwR_tilde_new(t_floatarg n)
{
    int i;

    // allocate the class
    t_RwR_tilde *x = (t_RwR_tilde *)pd_new(RwR_tilde_class);

    if(!n || n > MAX_OUTLET){
        n = 2;
    }
    x->x_n = (int)n;

    post("system top level dsp block size = %d", sys_getblksize());
    post("the sample-rate of the system = %.1f", sys_getsr());
    post("#Windows = %d", x->x_n);

    x->x_in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("list"), gensym("winlist"));
    x->x_in3 = floatinlet_new(&x->x_obj, &x->x_attenuation);
    x->x_in4 = floatinlet_new(&x->x_obj, &x->x_soundsource_x);
    x->x_in5 = floatinlet_new(&x->x_obj, &x->x_soundsource_y);
    x->x_in6 = floatinlet_new(&x->x_obj, &x->x_outerroom);

    // allocate the windows
    x->x_windows = (t_float*)getbytes(x->x_n*2*sizeof(t_float));
    if(!x->x_windows){
      post("RwR~: could not allocate windows");
      return NULL;
    }

    // allocate the outputs
    x->x_outputs = (t_outlet **)getbytes(x->x_n*sizeof(t_outlet *));
    if(!x->x_outputs) {
        post( "RwR~: could not allocate outputs" );
        return NULL;
    }
    for(i=0; i<x->x_n; i++){
        x->x_outputs[i] = outlet_new(&x->x_obj, &s_signal);
    }

    // allocate the outputmapping
    x->x_outputmapping = (t_int*)getbytes(x->x_n*sizeof(t_int));
    if(!x->x_outputmapping){
        post( "RwR~: cannot allocate outputmapping array");
        return NULL; // otherwise, pd schrieks
    }
    for(i=0; i<x->x_n; i++){
        x->x_outputmapping[i] = i;
    }

    // allocate the buffer
    x->x_buffer = (t_float*) getbytes(block_iter*sys_getblksize()*sizeof(t_float)); // 5000 * 64samples / 44100Hz = 7.256seconds
    if ( !x->x_buffer ) {
        post( "RwR~: could not allocate buffer[i]" );
        return NULL;
    }

    return (void *)x;
}

// Destructor
void RwR_tilde_free(t_RwR_tilde *x)
{
    int i;

    // free any ressources associated with the given inlet
    inlet_free(x->x_in2);
    inlet_free(x->x_in3);
    inlet_free(x->x_in4);
    inlet_free(x->x_in5);
    inlet_free(x->x_in6);

    // free any ressources associated with the given outlet
    if (x->x_outputs) {
        for (i=0; i<x->x_n; i++) {
            outlet_free(x->x_outputs[i]);
        }
        freebytes(x->x_outputs, x->x_n*sizeof(t_outlet*));
    }

    // free outmapping
    if (x->x_outputmapping) {
        freebytes(x->x_outputmapping, x->x_n*sizeof(t_int));
    }

    // free windows
    if (x->x_windows) {
        freebytes(x->x_windows, x->x_n*2*sizeof(t_float));
    }

    // free buffer
    if (x->x_buffer){
        freebytes(x->x_buffer, block_iter*sys_getblksize()*sizeof(t_float));
    }
}

// Define the function-space of the class
void RwR_tilde_setup(void){
    RwR_tilde_class = class_new(gensym("RwR~"),
                                (t_newmethod)RwR_tilde_new,
                                (t_method)RwR_tilde_free,
                                sizeof(t_RwR_tilde),
                                CLASS_DEFAULT,
                                A_DEFFLOAT, 0);

    class_addbang(RwR_tilde_class, RwR_tilde_bang);

    class_addmethod(RwR_tilde_class,
                    (t_method)RwR_tilde_winlist,
                    gensym("winlist"),
                    A_GIMME, 0);

    class_addmethod(RwR_tilde_class,
                    (t_method)RwR_tilde_dsp,
                    gensym("dsp"),
                    A_DEFFLOAT, // a float argument with 0 as default
                    0);
    CLASS_MAINSIGNALIN(RwR_tilde_class, t_RwR_tilde, f); //The third argument is a (dummy-)floating point-variable of the data space
}
