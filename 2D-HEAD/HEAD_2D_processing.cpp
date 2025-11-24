#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <vector>
#include <gsl/gsl_multimin.h>
using namespace std;

struct data{
    int TD1;
    int TD2;
    double lambda;
    vector<int> index;
    vector< vector<double> > spectrum;
    vector< vector<double> > spectrum_scaled;
    vector<double> F2_sum;
    vector<double> F1_sum;

    vector< vector<double> > hetspectrum;
    vector<double> hetF1_sum;
    int hetTD1;
    vector<int> hetindex;

    int slice;
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
    double MSD=0., gradient=0., Euc_norm=0.;
     vector<double> slice(TD2,0.);
    for(i=0;i<TD2;i++){
        slice[i]=spec->hetspectrum[spec->hetindex[spec->slice]][i];
    }

    //Derivative of the weights, used for Tikhonov
   // for(i=1;i<spec->index.size();i++){
    //    gradient+=abs(gsl_vector_get(weights,i-1)-gsl_vector_get(weights,i));
   // }
   // gradient*=spec->lambda;

    //To avoid there being a larger lambda for the tails of the peak shapes
   // gradient*=spec->hetF1_sum[spec->hetindex[spec->slice]];

    //Calculation of the F1 spectrum
    //looping over the basis spectra
    for(i=0;i<spec->index.size();i++){
        int ii=spec->index[i];
        int start_j= (ii-TD1/2)*(ii>=(TD1/2));
        int end_j= (ii+TD1/2)*((ii+TD1/2)<TD2)+(TD2-1)*((ii+TD1/2)>=TD2);
        Euc_norm+= pow(gsl_vector_get(weights,i),2.);

        //looping over the data points of each basis spectrum
        for(j=start_j;j<=end_j;j++){
                F1_sum[j]=F1_sum[j] + spec->hetspectrum[spec->hetindex[spec->slice]][ii]*spec->spectrum_scaled[j][ii]*abs(gsl_vector_get(weights,i));
        }
    }

    //Calculation of the mean squared deviation between the calculated F1 spectrum and the F2 spectrum
    for(i=0;i<TD2;i++){
        MSD=MSD + pow(F1_sum[i]-slice[i],2.);
    }

    //return sqrt(MSD)+gradient;
    return sqrt(MSD)+spec->lambda*Euc_norm;
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
        val[1]=val[0]+0.0001;
        gsl_vector_set(var,i,val[1]);
        gsl_vector_set(df, i, (RMSD(var, params)-cost_0)/(val[1]-val[0]));
        gsl_vector_set(var,i,val[0]);
    }
}

void calc_HETCOR(const gsl_vector* weights, void* params, vector< vector<double> >& iso_spec, int F1_index){
    //Function used to integrate over the HEAD spectrum to produce the isotropic spectrum with appropriate intensities

    //Gathering parameters from the data structure
    struct data *spec = (struct data*) params;
    int TD2=spec->TD2,i,j,ii;
    int TD1=spec->TD1;
    vector<double> corr_weights(TD2,0.);
    vector<double> F2_sum(TD2,0.);

    //Looping over the basis spectra
    for(i=0;i<spec->index.size();i++){
        ii=spec->index[i];
        int start_j= (ii-TD1/2)*(ii>=(TD1/2));
        int end_j= (ii+TD1/2)*((ii+TD1/2)<TD2)+(TD2-1)*((ii+TD1/2)>=TD2);

        //looping over the data points of each basis spectrum
        for(j=start_j;j<=end_j;j++){
            F2_sum[ii]=F2_sum[ii] + spec->hetspectrum[spec->hetindex[spec->slice]][ii]*spec->spectrum_scaled[j][ii]*abs(gsl_vector_get(weights,i));
        }

        //saving the isotropic spectrum and weights
        iso_spec[F1_index][ii]=F2_sum[ii];
    }
}

int main() {
    //parameter declaration
    int TD2, TD1, TD1X, i=0, j=0,k,junk, datatype;
    char totxt_filename[128], totxt_filename2[128], buffer[256], word[24], pound;
    double left, right, width,lambda, delta;
    double left1, right1, width1, delta1;
    FILE *fp;

    //Interface to gather the totxt filename and the lambda parameter for regularization
    printf("Select the type of 2D spectrum you wish to Hahn-echo assisted deconvolute:\n");
    printf("1 1H-detected HETCOR\n");
    printf("2 SQ-SQ correlation\n");
    printf("3 SQ-DQ correlation\n");
    scanf("%d",&datatype);

    printf("\nWhat is the filename for the 2D CS-lineshape correlation spectrum converted using totxt?\n");
    printf("Note that the digital resolution in F1 and F2 needs to be identical.\n");
    scanf("%s",totxt_filename);

    printf("\nWhat is the filename for the totxtx 2D spectrum to be deconvoluted?\n");
    printf("Note that the digital resolution in F2 needs to be identical to the previous spectrum.\n");
    scanf("%s",totxt_filename2);

    printf("\nLambda value for regularization (0 for pure least squares)\n");
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
    state=0;

    fp=fopen(totxt_filename2,"r");
    if(fp==NULL){
        FILE *error_file;
        error_file=fopen("error.txt","a");
        fprintf(error_file, "\nERROR: totxt spectrum file '%s' not found\n", totxt_filename2);
        fclose(error_file);
        exit(1);
    }
    while ((fgets(buffer, sizeof(buffer), fp) != NULL)||(state<3)) {
        if(buffer[0]=='#'){
                sscanf(buffer,"%c %s",&pound,word);
                if(strcmp(word,"F1LEFT")==0){
                    sscanf(buffer,"%c %s %s %lf %s %s %s %lf",&pound,word,word,&left1,word,word,word,&right1);
                    sprintf(word,"void");
                    state++;
                }
                else if(strcmp(word,"NROWS")==0){
                    sscanf(buffer,"%c %s %s %d",&pound,word,word,&TD1X);
                    sprintf(word,"void");
                    state++;
                    if((datatype==2)&&(TD1X!=TD2)){
                        FILE *error_file;
                        error_file=fopen("error.txt","a");
                        fprintf(error_file, "\nERROR: the two TD2 values do not match\n");
                        fclose(error_file);
                        exit(1);
                    }
                    else if((datatype==3)&&(TD1X!=2*TD2)){
                        FILE *error_file;
                        error_file=fopen("error.txt","a");
                        fprintf(error_file, "\nERROR: the two TD2 values do not match\n");
                        fclose(error_file);
                        exit(1);
                    }
                }

                else if(strcmp(word,"NCOLS")==0){
                    sscanf(buffer,"%c %s %s %d",&pound,word,word,&junk);
                    sprintf(word,"void");
                    state++;
                    if(junk!=TD2){
                        FILE *error_file;
                        error_file=fopen("error.txt","a");
                        fprintf(error_file, "\nERROR: the two TD2 values do not match\n");
                        fclose(error_file);
                        exit(1);
                    }
                }
            }
    }
    fclose(fp);
    width1=left1-right1;
    delta1=width1/TD2;

    if(datatype==3)//DSSQ into a SQSQ format
        TD1X/=2;

    //creating the data structures
    struct data spec, het;
    spec.spectrum.resize(TD2);
    spec.spectrum_scaled.resize(TD2);
    spec.F2_sum.resize(TD2,0.);
    spec.F1_sum.resize(TD2,0.);
    spec.TD2=TD2;
    spec.TD1=TD1;
    spec.lambda=lambda;
    vector<double> F1_sum(TD2,0.);
    for(i=0;i<TD2;i++){
        spec.spectrum[i].resize(TD2,0.);
        spec.spectrum_scaled[i].resize(TD2,0.);
    }
    //HETCOR spectrum to be fitted
    spec.hetspectrum.resize(TD1X);
    for(i=0;i<TD1X;i++){
        spec.hetspectrum[i].resize(TD2,0.);
    }
    spec.hetTD1=TD1X;
    spec.hetF1_sum.resize(TD1X,0.);

    //reading the HEAD spectrum intensities and storing them as a sheared spectrum TD2xTD2
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

                    spec.spectrum_scaled[F1_index][i]=spec.spectrum[F1_index][i];
                }
            }
        }
    }
    fclose(fp);
    for(i=0;i<TD2;i++){
        double maximum=0;
        for(j=0;j<TD2;j++){
            if(spec.spectrum_scaled[j][i]>maximum)
                maximum=spec.spectrum_scaled[j][i];
        }
        if(maximum>0.){
        for(j=0;j<TD2;j++){
            spec.spectrum_scaled[j][i]= spec.spectrum_scaled[j][i]/maximum;
        }}
    }

    if(datatype!=3){
    //reading the HETCOR spectrum intensities
    fp=fopen(totxt_filename2,"r");
    for(j=0;j<TD1X;j++){
        for(i=0; i<TD2; i++){
            fgets(buffer, sizeof(buffer), fp);
            if(buffer[0]=='#'){
                i--;
            }
            else{
                    sscanf(buffer,"%lf",&spec.hetspectrum[j][i]);
                    if(spec.hetspectrum[j][i]<0.)
                        spec.hetspectrum[j][i]=0.;
                  /*  for(k=0;k<TD2;k++){
                        spec.spectrum_scaled[j][k][i]*=spec.hetspectrum[j][i];
                    }*/
            }
        }
    }
    fclose(fp);
    }


    else{//DQSQ reading into SQSQ
    fp=fopen(totxt_filename2,"r");
    for(j=0;j<TD1X*2;j++){
        for(i=0; i<TD2; i++){

            int F1_index=-i+j;

            fgets(buffer, sizeof(buffer), fp);
            if(buffer[0]=='#'){
                i--;
            }
            else if((F1_index>=0)&&(F1_index<TD2)){
                    sscanf(buffer,"%lf",&spec.hetspectrum[F1_index][i]);
                    if(spec.hetspectrum[F1_index][i]<0.)
                        spec.hetspectrum[F1_index][i]=0.;
                   /* for(k=0;k<TD2;k++){
                        spec.spectrum_scaled[F1_index][k][i]*=spec.hetspectrum[F1_index][i];
                    }*/
            }
        }
    }
    fclose(fp);
    }



    //Saving the original 2D spectrum as a *.spe file
    fp=fopen("original_HEAD.spe","w");
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
    fp=fopen("sheared_HEAD.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*600.0,left*600.0,600.0,TD2,width*600.0,left*600.0,600.0);
    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            fprintf(fp,"%lf   0.0\n",spec.spectrum[i][j]);
        }
    }
    fclose(fp);

    //Saving the original HETCOR spectrum as a *.spe file
    if(datatype!=3){
    fp=fopen("original_correlation.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*600.0,left*600.0,600.0,TD1X,width1*150.0,left1*150.,150.0);
    for(i=0;i<TD1X;i++){
        for(j=0;j<TD2;j++){
            fprintf(fp,"%lf   0.0\n",spec.hetspectrum[i][j]);
        }
    }
    fclose(fp);
    }

    else{
    fp=fopen("original_correlation.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,10.0*600.0,10.0*600.0,600.0,TD2*2,20.0*600.0,20.0*600.0,600.0);
    for(i=0;i<TD2*2;i++){
        for(j=0;j<TD2;j++){
            int F1_index=-j+i;
            if((F1_index>=0)&&(F1_index<TD2))
                fprintf(fp,"%lf   0.0\n",spec.hetspectrum[F1_index][j]);
            else
                fprintf(fp,"0   0.0\n");
        }}
    fclose(fp);
    }

    //normalizing the scaled HEAD spectrum
    double max_F2=0.;
    fill(spec.F2_sum.begin(), spec.F2_sum.end(), 0.);
    fill(F1_sum.begin(), F1_sum.end(), 0.);
    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            spec.F2_sum[j]=spec.F2_sum[j]+spec.spectrum_scaled[i][j];
            F1_sum[i]=F1_sum[i]+spec.spectrum_scaled[i][j];
            if(spec.F2_sum[j]>max_F2)
                max_F2=spec.F2_sum[j];
        }
    }

    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            spec.spectrum_scaled[i][j]=spec.spectrum_scaled[i][j]/max_F2;
        }
        spec.F2_sum[i]=spec.F2_sum[i]/max_F2;
        F1_sum[i]=F1_sum[i]/max_F2;
    }


        //normalizing the F2 spectrum, (sum of rows)
    fill(spec.F2_sum.begin(), spec.F2_sum.end(), 0.);
    fill(F1_sum.begin(), F1_sum.end(), 0.);
    max_F2=0.;
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


    //Normalizing the HETCOR
    max_F2=0.;
    double max_F1=0.;
    for(i=0;i<TD1X;i++){
        for(j=0;j<TD2;j++){
            if(spec.hetspectrum[i][j]>max_F2)
                max_F2=spec.hetspectrum[i][j];
    }}
    for(i=0;i<TD1X;i++){
        for(j=0;j<TD2;j++){
            spec.hetspectrum[i][j]=spec.hetspectrum[i][j]/max_F2;
            spec.hetF1_sum[i]+=spec.hetspectrum[i][j];
    }
        if(spec.hetF1_sum[i]>max_F1)
            max_F1=spec.hetF1_sum[i];
    }
    for(i=0;i<TD1X;i++){
        spec.hetF1_sum[i]/=max_F1;
    }


    //deciding which datapoint intensities are to be optimized
    //Here a default minimum intensity of 1% is used, this can be changed.
    for(i=0;i<TD2;i++){
         if(spec.F2_sum[i]>0.01){
            spec.index.push_back(i);
         }
    }
    for(i=0;i<TD1X;i++){
         if(spec.hetF1_sum[i]>0.05){
            spec.hetindex.push_back(i);
         }
    }

    vector< vector<double> > iso_spec;
    iso_spec.resize(TD1X);
    for(i=0;i<TD1X;i++){
        iso_spec[i].resize(TD2,0.);
    }

    int done=0;
    #pragma omp parallel for
    for(int J=0;J<spec.hetindex.size();J++){
    int I,K;
    struct data spec_temp=spec;
    printf("%d weights to optimize\n",spec.index.size());
    //setting the initial weights for the Tikhonov minimization to 1
    gsl_vector *weights;
    weights = gsl_vector_alloc (spec.index.size());
    spec_temp.slice=J;
    for(I=0;I<spec.index.size();I++){
        float wt=1.0;
        gsl_vector_set(weights,I,wt);
    }

    //The GSL gradient minimization is done here. This closely mirrors the example from
    //the GSL manual. BFGS, CG, and steepest descent algorithms all work well, but simplex
    //methods do not. If you chose to change the algorithm, you can uncomment the appropriate line below
    //The gsl_multimin_fdfminimizer_set parameters may need to be tweaked in that case.

    const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
    //const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
   // const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_steepest_descent;
    gsl_multimin_function_fdf minex_func; //declaring the minimizer
    minex_func.n = spec.index.size();  //number of variables
    minex_func.f = RMSD; //cost function assignment
    minex_func.df = gradient; //gradient calculation function
    minex_func.fdf = gradient_fdf; //simultaneous gradient and cost calculation function
    minex_func.params = &spec_temp; //parameter assignment (non-optimized things)
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, minex_func.n); //pointer for the minimizer
    gsl_multimin_fdfminimizer_set (s, &minex_func, weights, 0.1, 1e-0);  //Assigning the initial weights, functions, etc.
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
        status = gsl_multimin_test_gradient (s->gradient, 0.00000001);

        if (status == GSL_SUCCESS){
            break;
        }
        cost=s->f;

        //Print to the screen if there is a better solution
          if(cost<cost_min){
            printf ("%d/%d %d %f \n",done,spec.hetindex.size(),iter,s->f);
            cost_min=cost;
        }
    }  while (status == GSL_CONTINUE && iter < 10000);

    printf("Converged!\n");
    done++;
    calc_HETCOR(s->x,&spec_temp,iso_spec,spec.hetindex[J]);
    gsl_vector_free(weights);
    gsl_multimin_fdfminimizer_free(s);
    }

    //Saving the isotropic spectrum as a *.spe file
    if(datatype!=3){
    fp=fopen("isotropic_correlation_F2only.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*600.0,left*600.0,600.0,TD1X,width1*600.0,left1*600.,600.0);
    for(i=0;i<TD1X;i++){
        for(j=0;j<TD2;j++){
            fprintf(fp,"%lf   0.0\n",iso_spec[i][j]);
        }
    }
    fprintf(fp,"%s","END");
    fclose(fp);
    }

    else{
    fp=fopen("isotropic_correlation_F2only.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,10.0*600.0,10.0*600.0,600.0,TD2*2,20.0*600.0,20.0*600.0,600.0);
    for(i=0;i<TD2*2;i++){
        for(j=0;j<TD2;j++){
            int F1_index=-j+i;
            if((F1_index>=0)&&(F1_index<TD2))
                fprintf(fp,"%lf   0.0\n",iso_spec[F1_index][j]);
            else
                fprintf(fp,"0   0.0\n");
        }}
    fclose(fp);
    }

    if(datatype==1)//HETCOR, onle F2 deconvolution
        return 0;

    //transpose
    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            spec.hetspectrum[i][j]=iso_spec[j][i];
        }
    }

    done=0;
    #pragma omp parallel for
    for(int J=0;J<spec.hetindex.size();J++){
    int I,K;
    struct data spec_temp=spec;
    printf("%d weights to optimize\n",spec.index.size());
    //setting the initial weights for the Tikhonov minimization to 1
    gsl_vector *weights;
    weights = gsl_vector_alloc (spec.index.size());
    spec_temp.slice=J;
    for(I=0;I<spec.index.size();I++){
        float wt=1.0;
        gsl_vector_set(weights,I,wt);
    }

    //The GSL gradient minimization is done here. This closely mirrors the example from
    //the GSL manual. BFGS, CG, and steepest descent algorithms all work well, but simplex
    //methods do not. If you chose to change the algorithm, you can uncomment the appropriate line below
    //The gsl_multimin_fdfminimizer_set parameters may need to be tweaked in that case.

    const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
    //const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
   // const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_steepest_descent;
    gsl_multimin_function_fdf minex_func; //declaring the minimizer
    minex_func.n = spec.index.size();  //number of variables
    minex_func.f = RMSD; //cost function assignment
    minex_func.df = gradient; //gradient calculation function
    minex_func.fdf = gradient_fdf; //simultaneous gradient and cost calculation function
    minex_func.params = &spec_temp; //parameter assignment (non-optimized things)
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, minex_func.n); //pointer for the minimizer
    gsl_multimin_fdfminimizer_set (s, &minex_func, weights, 0.1, 1e-0);  //Assigning the initial weights, functions, etc.
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
        status = gsl_multimin_test_gradient (s->gradient, 0.00000001);

        if (status == GSL_SUCCESS){
            break;
        }
        cost=s->f;

        //Print to the screen if there is a better solution
          if(cost<cost_min){
            printf ("%d/%d %d %f \n",done,spec.hetindex.size(),iter,s->f);
            cost_min=cost;
        }
    }  while (status == GSL_CONTINUE && iter < 10000);

    printf("Converged!\n");
    done++;
    calc_HETCOR(s->x,&spec_temp,iso_spec,spec.hetindex[J]);
    gsl_vector_free(weights);
    gsl_multimin_fdfminimizer_free(s);
    }

    //Saving the isotropic spectrum as a *.spe file
    if(datatype!=3){
    fp=fopen("isotropic_HOMCOR_spectrum.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*600.0,left*600.0,600.0,TD2,width*600.0,left*600.,600.0);
    for(i=0;i<TD1X;i++){
        for(j=0;j<TD2;j++){
            fprintf(fp,"%lf   0.0\n",iso_spec[i][j]);
        }
    }
    fprintf(fp,"%s","END");
    fclose(fp);
    }

    else{
    fp=fopen("isotropic_HOMCOR_spectrum.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,10.0*600.0,10.0*600.0,600.0,TD2*2,20.0*600.0,20.0*600.0,600.0);
    for(i=0;i<TD2*2;i++){
        for(j=0;j<TD2;j++){
            int F1_index=-j+i;
            if((F1_index>=0)&&(F1_index<TD2))
                fprintf(fp,"%lf   0.0\n",iso_spec[F1_index][j]);
            else
                fprintf(fp,"0   0.0\n");
        }}
    fclose(fp);
    }

    return 0;
}
