#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <vector>
#include <gsl/gsl_multimin.h>
#include "decon.hpp"
using namespace std;

struct data{
    int TD1;
    int TD2;
    double lambda;
    vector<int> index;
    vector< vector<double> > spectrum;
    vector<double> F2_sum;
    vector<double> F1_sum;
};

double RMSD(const gsl_vector* weights, void* params){
    //Cost function used to calculate the offset between the experimental F2 spectrum
    //and the predicted isotropic spectrum in F1. Also includes Tikhonov weighting.

    //gathering the relevant parameters from the data structure
    struct data *spec = (struct data*) params;
    int TD2=spec->TD2,i,j;
    int TD1=spec->TD1;
    vector<double> F1_sum(TD2,0.);
    vector<double> corr_weights(TD2,0.);
    double MSD=0., gradient=0.;

    //Derivative of the weights, used for Tikhonov
    for(i=1;i<spec->index.size();i++){
        gradient+=abs(gsl_vector_get(weights,i-1)-gsl_vector_get(weights,i));
    }
    gradient*=spec->lambda;

    //Calculation of the F1 spectrum
    //looping over the basis spectra
    for(i=0;i<spec->index.size();i++){
        int ii=spec->index[i];
        int start_j= (ii-TD1/2)*(ii>=(TD1/2));
        int end_j= (ii+TD1/2)*((ii+TD1/2)<TD2)+(TD2-1)*((ii+TD1/2)>TD2);

        //looping over the data points of each basis spectrum
        for(j=start_j;j<=end_j;j++){
            F1_sum[j]=F1_sum[j] + spec->spectrum[j][ii]*abs(gsl_vector_get(weights,i));
        }
    }

    //Calculation of the mean squared deviation between the calculated F1 spectrum and the F2 spectrum
    for(i=0;i<TD2;i++){
        MSD=MSD + pow(F1_sum[i]-spec->F2_sum[i],2.);
    }

    return sqrt(MSD)+gradient;
}

void gradient(const gsl_vector *var, void *params, gsl_vector *df){
    //Gradient of the RMSD, used for GSL gradient optimizers

    struct data *spec = (struct data*) params;
    int var_size=spec->index.size(),i;
    double cost_0 = RMSD(var, params), val[2];

    for(i=0;i<var_size;i++){
        val[0]=gsl_vector_get(var,i);
        val[1]=val[0]+0.0001;
        gsl_vector_set(var,i,val[1]);
        gsl_vector_set(df, i, (RMSD(var, params)-cost_0)/(val[1]-val[0]));
        gsl_vector_set(var,i,val[0]);
    }
}

void gradient_fdf (const gsl_vector *var, void *params,double *f,gsl_vector *df){
    //Simultaneously calculated the gradient and the value of the RMSD
    //For the GSL gradient optimizers

    *f = RMSD(var, params);
    struct data *spec = (struct data*) params;
    int var_size=spec->index.size(),i;
    double cost_0 = RMSD(var, params), val[2];

    for(i=0;i<var_size;i++){
        val[0]=gsl_vector_get(var,i);
        val[1]=val[0]+0.001;
        gsl_vector_set(var,i,val[1]);
        gsl_vector_set(df, i, (RMSD(var, params)-cost_0)/(val[1]-val[0]));
        gsl_vector_set(var,i,val[0]);
    }
}

void calc_F1sum(const gsl_vector* weights, void* params, vector<double>& F1_sum, vector<double>& iso_spec, vector<double>& wt){
    //Function used to integrate over the HEAD spectrum to produce the isotropic spectrum with appropriate intensities

    //Gathering parameters from the data structure
    struct data *spec = (struct data*) params;
    int TD2=spec->TD2,i,j,ii;
    int TD1=spec->TD1;
    vector<double> corr_weights(TD2,0.);
    std::fill(F1_sum.begin(), F1_sum.end(), 0.);

    //Looping over the basis spectra
    for(i=0;i<spec->index.size();i++){
        ii=spec->index[i];
        int start_j= (ii-TD1/2)*(ii>=(TD1/2));
        int end_j= (ii+TD1/2)*((ii+TD1/2)<TD2)+(TD2-1)*((ii+TD1/2)>=TD2);

        //looping over the data points of each basis spectrum
        //spectrum[F1][F2], i is for F2 indices
        for(j=start_j;j<=end_j;j++){
            iso_spec[ii]=iso_spec[ii] + spec->spectrum[j][ii]*abs(gsl_vector_get(weights,i));
        }

        //saving the weights
        wt[ii]=abs(gsl_vector_get(weights,i));
    }
}

int main() {
    //parameter declaration
    int TD2, TD1, i=0, j=0;
    char totxt_filename[128], buffer[256], word[24], pound;
    double left, right, width, delta,lambda;
    FILE *fp;

    //Interface to gather the totxt filename and the lambda parameter for regularization
    printf("What is the filename for the 2D CS-lineshape correlation spectrum converted using totxt?\n");
    printf("Note that the digital resolution in F1 and F2 needs to be identical.\n");
    scanf("%s",totxt_filename);

    printf("Lambda value for regularization (0 for pure least squares)\n");
    scanf("%lf",&lambda);

    //reading the header file to get TD1, TD2, and the spectral width and offset in F2.
    fp=fopen(totxt_filename,"r");
    if(fp==NULL){
        FILE *error_file;
        error_file=fopen("error.txt","a");
        fprintf(error_file, "\nERROR: totxt spectrum file '%s' not found\n", totxt_filename);
        fclose(error_file);
        exit(1);
    }

    int state=0;
    while ((fgets(buffer, sizeof(buffer), fp) != NULL)||(state<3)) {
        if(buffer[0]=='#'){
                sscanf(buffer,"%c %s",&pound,word);
                if(strcmp(word,"F2LEFT")==0){
                    sscanf(buffer,"%c %s %s %lf %s %s %s %lf",&pound,word,word,&left,word,word,word,&right);
                    sprintf(word,"void");
                    state++;
                }
                else if(strcmp(word,"NROWS")==0){
                    sscanf(buffer,"%c %s %s %d",&pound,word,word,&TD1);
                    sprintf(word,"void");
                    state++;
                }
                else if(strcmp(word,"NCOLS")==0){
                    sscanf(buffer,"%c %s %s %d",&pound,word,word,&TD2);
                    sprintf(word,"void");
                    state++;
                }
            }
    }
    fclose(fp);
    width=left-right;
    delta=width/TD2;

    //creating the data structures
    struct data spec;
    spec.spectrum.resize(TD2);
    spec.F2_sum.resize(TD2,0.);
    spec.F1_sum.resize(TD2,0.);
    spec.TD2=TD2;
    spec.TD1=TD1;
    spec.lambda=lambda;
    vector<double> F1_sum(TD2,0.);
    for(i=0;i<TD2;i++){
        spec.spectrum[i].resize(TD2,0.);
    }

    //reading the spectrum intensities and storing them as a sheared spectrum TD2xTD2
    fp=fopen(totxt_filename,"r");
    for(j=0;j<TD1;j++){
        for(i=0; i<TD2; i++){
            fgets(buffer, sizeof(buffer), fp);
            if(buffer[0]=='#'){
                i--;
            }
            else{
                int F1_index=-TD1/2+i+j;
                if((F1_index>0)&&(F1_index<TD2)){
                    sscanf(buffer,"%lf",&spec.spectrum[F1_index][i]);
                    if(spec.spectrum[F1_index][i]<0.)
                        spec.spectrum[F1_index][i]=0.;
                }
            }
        }
    }
    fclose(fp);

    //Saving the original 2D spectrum as a *.spe file
    fp=fopen("original_spectrum.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*600.0,left*600.0,600.0,TD1,width*600.0*TD1/TD2,width*300.0*TD1/TD2,600.0);
    for(i=0;i<TD1;i++){
        for(j=0;j<TD2;j++){
            int F1_index=-TD1/2+i+j;
            if((F1_index>0)&&(F1_index<TD2))
                fprintf(fp,"%lf   0.0\n",spec.spectrum[F1_index][j]);
            else
                fprintf(fp,"0.0   0.0\n");
        }
    }
    fclose(fp);

    //Saving the sheared 2D spectrum as a *.spe file
    fp=fopen("sheared_spectrum.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*600.0,left*600.0,600.0,TD2,width*600.0,left*600.0,600.0);
    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            fprintf(fp,"%lf   0.0\n",spec.spectrum[i][j]);
        }
    }
    fclose(fp);

    //normalizing the F2 spectrum, (sum of rows)
    double max_F2=0.;
    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            spec.F2_sum[j]=spec.F2_sum[j]+spec.spectrum[i][j];
            F1_sum[i]=F1_sum[i]+spec.spectrum[i][j];
            if(spec.F2_sum[j]>max_F2)
                max_F2=spec.F2_sum[j];
        }
    }

    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            spec.spectrum[i][j]=spec.spectrum[i][j]/max_F2;
        }
        spec.F2_sum[i]=spec.F2_sum[i]/max_F2;
        F1_sum[i]=F1_sum[i]/max_F2;
    }

    //deciding which datapoint intensities are to be optimized
    //Here a default minimum intensity of 1% is used, this can be changed.
    for(i=0;i<TD2;i++){
         if(spec.F2_sum[i]>0.01){
            spec.index.push_back(i);
         }
    }
    printf("%d weights to optimize\n",spec.index.size());

    //setting the initial weights for the Tikhonov minimization to 1
    gsl_vector *weights;
    weights = gsl_vector_alloc (spec.index.size());
    for(i=0;i<spec.index.size();i++){
        gsl_vector_set(weights,i,1.);
    }

    //The GSL gradient minimization is done here. This closely mirrors the example from
    //the GSL manual. BFGS, CG, and steepest descent algorithms all work well, but simplex
    //methods do not. If you chose to change the algorithm, you can uncomment the appropriate line below
    //The gsl_multimin_fdfminimizer_set parameters may need to be tweaked in that case.

    const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
   // const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
   // const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_steepest_descent;
    gsl_multimin_function_fdf minex_func; //declaring the minimizer
    minex_func.n = spec.index.size();  //number of variables
    minex_func.f = RMSD; //cost function assignment
    minex_func.df = gradient; //gradient calculation function
    minex_func.fdf = gradient_fdf; //simultaneous gradient and cost calculation function
    minex_func.params = &spec; //parameter assignment (non-optimized things)
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, minex_func.n); //pointer for the minimizer
    gsl_multimin_fdfminimizer_set (s, &minex_func, weights, 0.01, 1e-0);  //Assigning the initial weights, functions, etc.
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
            printf ("%d %f \n",iter,s->f);
            cost_min=cost;
        }
    }  while (status == GSL_CONTINUE && iter < 1000000);

    printf("Converged!\n");

    //writing out the result of the fit to a comma delimited file
    vector<double> iso_spec(TD2,0.);
    vector<double> wt(TD2,0.);
    calc_F1sum(s->x,&spec,F1_sum,iso_spec,wt);
    fp=fopen("result.csv","w");
    fprintf(fp,"index,shift,source,fit,isotropic,weights\n");
    for(i=0;i<TD2;i++){
        fprintf(fp,"%d,%lf,%lf,%lf,%lf,%lf\n",i,left-delta*i,spec.F2_sum[i],F1_sum[i],iso_spec[i],wt[i]);
    }
    fclose(fp);

    //Also saving the isotropic spectrum as a *.spe file
    fp=fopen("isotropic_spectrum.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*600.0,left*600.0,600.0);
    for(i=0;i<TD2;i++){
        fprintf(fp,"%lf 0\n",iso_spec[i]);
    }
    fprintf(fp,"%s","END");
    fclose(fp);

    //*.spe file of the 2D version of the above spectrum
    fp=fopen("isotropic_spectrum_2D.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*600.0,left*600.0,600.0,TD2,width*600.0,left*600.0,600.0);
    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            fprintf(fp,"%lf   0.0\n",spec.spectrum[i][j]*wt[j]);
        }
    }
    fclose(fp);

  //  auto_decon();

    return 0;
}
