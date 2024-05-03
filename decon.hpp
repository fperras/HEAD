#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <vector>
#include <gsl/gsl_multimin.h>
#include<iostream>
#include<fstream>
#include<string>
using namespace std;

struct iso_spec{
    int TD;
    int num_gauss;
    vector<double> GB;
    vector<double> intensity;
    vector<double> shift;
    vector<double> original_spectrum;
    vector<double> fitted_spectrum;
    vector<double> shift_axis;
};

void calculate_spectrum(void* params){
    struct iso_spec *spec = (struct iso_spec*) params;
    int i,j;

    for(i=0;i<spec->TD;i++){
        spec->fitted_spectrum[i]=0.;
    }

    for(j=0;j<spec->num_gauss;j++){
        for(i=0;i<spec->TD;i++){
            double offset=fabs(spec->shift[j] - spec->shift_axis[i]);
            spec->fitted_spectrum[i] += spec->intensity[j]*exp(-0.5*pow((offset)/spec->GB[j],2.));
        }
    }
}

double RMSD_G(const gsl_vector* variables, void* params){
    struct iso_spec *spec = (struct iso_spec*) params;
    int i;
    double MSD=0.;

    for(i=0;i<spec->num_gauss;i++){
        spec->GB[i]=gsl_vector_get(variables,i);
        spec->intensity[i]=gsl_vector_get(variables,i+spec->num_gauss);
        spec->shift[i]=gsl_vector_get(variables,i+2*spec->num_gauss);
    }

    calculate_spectrum(spec);

    for(i=0;i<spec->TD;i++){
        MSD += pow(spec->fitted_spectrum[i]-spec->original_spectrum[i],2.);
    }
    return sqrt(MSD);
}


void gradient_G(const gsl_vector *var, void *params, gsl_vector *df){
    //Gradient of the RMSD, used for GSL gradient optimizers

    struct iso_spec *spec = (struct iso_spec*) params;
    int var_size=3*spec->num_gauss,i;
    double cost_0 = RMSD_G(var, params), val[2];

    for(i=0;i<var_size;i++){
        val[0]=gsl_vector_get(var,i);
        val[1]=val[0]+0.0001;
        gsl_vector_set(var,i,val[1]);
        gsl_vector_set(df, i, (RMSD_G(var, params)-cost_0)/(val[1]-val[0]));
        gsl_vector_set(var,i,val[0]);
    }
}

void gradient_fdf_G(const gsl_vector *var, void *params,double *f,gsl_vector *df){
    //Simultaneously calculated the gradient and the value of the RMSD
    //For the GSL gradient optimizers

    *f = RMSD_G(var, params);
    struct iso_spec *spec = (struct iso_spec*) params;
    int var_size=3*spec->num_gauss,i;
    double cost_0 = RMSD_G(var, params), val[2];

    for(i=0;i<var_size;i++){
        val[0]=gsl_vector_get(var,i);
        val[1]=val[0]+0.001;
        gsl_vector_set(var,i,val[1]);
        gsl_vector_set(df, i, (RMSD_G(var, params)-cost_0)/(val[1]-val[0]));
        gsl_vector_set(var,i,val[0]);
    }
}

void auto_decon(){
    char buffer[256],junk;
    FILE *fp;
    struct iso_spec spec;
    double SW,X0,SF,max_intensity=0;
    int i;

    //reading the SPE fine, assuming a given format/order
    fp=fopen("isotropic_spectrum.spe","r");
    fgets(buffer,sizeof(buffer),fp);
    fscanf(fp,"%c%c%c%d\n",&junk,&junk,&junk,&spec.TD);
    fscanf(fp,"%c%c%c%lf\n",&junk,&junk,&junk,&SW);
    fscanf(fp,"%c%c%c%lf\n",&junk,&junk,&junk,&X0);
    fscanf(fp,"%c%c%c%lf\n",&junk,&junk,&junk,&SF);
    fgets(buffer,sizeof(buffer),fp);
    fgets(buffer,sizeof(buffer),fp);

    //resizing the vectors for the spectrum
    spec.original_spectrum.resize(spec.TD,0.);
    spec.fitted_spectrum.resize(spec.TD,0.);
    spec.shift_axis.resize(spec.TD,0.);

    //gathering the X and Y axes
    for(i=0;i<spec.TD;i++){
        fgets(buffer,sizeof(buffer),fp);
        sscanf(buffer,"%lf %c",&spec.original_spectrum[i],&junk);
        spec.shift_axis[i]=(X0-i*SW/spec.TD)/SF;
        if(spec.original_spectrum[i]>max_intensity)
            max_intensity=spec.original_spectrum[i];
    }
    fclose(fp);

    for(i=0;i<spec.TD;i++){
        spec.original_spectrum[i]=spec.original_spectrum[i]/max_intensity;
    }

    spec.num_gauss=0;
    do{
        max_intensity=0.;
        int max_i;
        for(i=0;i<spec.TD;i++){
            if((spec.original_spectrum[i]-spec.fitted_spectrum[i]) > max_intensity){
                max_intensity=(spec.original_spectrum[i]-spec.fitted_spectrum[i]);
                max_i=i;
            }
        }

        if(max_intensity > 0.05){
            spec.GB.push_back(0.2);
            spec.intensity.push_back(0.5);
            spec.shift.push_back(spec.shift_axis[max_i]);
            spec.num_gauss++;

        gsl_vector *variables;
        variables = gsl_vector_alloc(3*spec.num_gauss);
        for(i=0;i<spec.num_gauss;i++){
            gsl_vector_set(variables,i,spec.GB[i]);
            gsl_vector_set(variables,i+spec.num_gauss,spec.intensity[i]);
            gsl_vector_set(variables,i+2*spec.num_gauss,spec.shift[i]);
        }

        const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
        gsl_multimin_function_fdf minex_func; //declaring the minimizer
        minex_func.n = spec.num_gauss*3;  //number of variables
        minex_func.f = RMSD_G; //cost function assignment
        minex_func.df = gradient_G; //gradient calculation function
        minex_func.fdf = gradient_fdf_G; //simultaneous gradient and cost calculation function
        minex_func.params = &spec; //parameter assignment (non-optimized things)
        gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, minex_func.n); //pointer for the minimizer
        gsl_multimin_fdfminimizer_set (s, &minex_func, variables, 0.01, 1e-0);  //Assigning the initial weights, functions, etc.
        double cost,cost_min=100000000.; //initial cost is set arbitrarily high

        size_t iter = 0;
        int status;


        //Optimization loop with a maximum of 1 million steps
        do{
            iter++;
            status = gsl_multimin_fdfminimizer_iterate (s);

            if (status)
                break;

            //maximum allowable gradient convergence criterion
            status = gsl_multimin_test_gradient (s->gradient, 0.0000001);

            if (status == GSL_SUCCESS){
                break;
            }
            cost=s->f;

            //Print to the screen if there is a better solution
            if(cost<cost_min){
                printf ("%d %f\n",iter,s->f);
                cost_min=cost;
            }
        }  while (status == GSL_CONTINUE && iter < 1000000);
        }
    }while(max_intensity>0.05);

    fp=fopen("isotropic_spectrum_dcon.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nTYPE=SPE\nDATA\n",spec.TD,SW,X0,SF);
    for(i=0;i<spec.TD;i++){
        fprintf(fp,"%lf 0\n",spec.fitted_spectrum[i]);
    }
    fprintf(fp,"%s","END");
    fclose(fp);
}
