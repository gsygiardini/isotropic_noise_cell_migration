#include <stdio.h>
#include <math.h>

#include <time.h>

#include <gsl/gsl_rng.h>                                 //For the program to be run it is necessary to install the Gnu Scientific Library (GSL)
#include <gsl/gsl_randist.h>                             //into the computer

//////////////////////////////////////////////////////////////////////////////

#define File_Name_Length                         100     //Maximum name lenght of the output files

#define N_Files                                  10      //Number of output files

//////////////////////////////////////////////////////////////////////////////
                                                         //The Number of cell crawlers and the K parameters are defined at each run of the program

const double dt                                 = 1e-4;  //Minimum time interval between steps

const double Gamma                              = 1.0;   //Gamma parameter of the cell crawling model
const double g                                  = 10.0;   //G parameter of the cell crawling model
const double q                                  = 0.1;   //Q parameter of the cell crawling model
const double phi                                = 0.0;  //Field parameter for theta orientation
const double theta0                             = 0.0;   //Field parameter for theta orientation
const double b                                  = 1.0;   //How much a cell wants to move forwards in the direction of its polarity instead of backwards

const double delta                              = 1e-4;  //Delta parameter for the measure of the Mean Velocity Autocorrelation Function, the Velocity Relaxation Time and the Average Velocity Histogram

const double v0_mult                            = 1.0;   //Multiplier of the initial velocity v0, where v0 is the stationary velocity

const double delay                              = 1.0;   //delay defines how slow the MVACF loop will run, if delay=0 the loop will be run proportionally to 2^n       

const int iter                                  = 1e8;    //Number of steps given by the crawling Cells
                                                         
const int n_bins                                = 1000;

const int n_bins_2d                             = 400;
                                                         //Down below, at the open_output_files(void) functions, it is necessary to define wich files will be created and
                                                         //consequently which measuring functions will be executed by the program.
                                                         //Just change the value of files[i]=1 to create and execute the corresponding file and measure

//////////////////////////////////////////////////////////////////////////////

int files[N_Files];
char file_name[File_Name_Length];
FILE *f_Pos, *f_MeanSqd, *f_Corr, *f_TruncMeanVacf, *f_OvrlpMeanVacf, *f_MeanSpd, *f_VRelax, *f_Histo2D, *f_Pdf, *f_Velocities;

int step;
int N_Cell; double k;

void open_output_files(void)
{
  for(step = 0 ; step < N_Files ; step++)
    {
      files[step] = 0;
    }

  //Set files[i]=1 if you want to activate both the file declaration and its corresponding function
  //eg: If you want that the program print the trajectories, set files[0] equal to 1: files[0]=1
  

  files[0] = 0;   //Set files[0]=1 if you want the program to generate and output to a file all the cells' positions
  if(files[0]==1)
    {
      sprintf(file_name, "Pos_%d_%g_%g_%g.dat",N_Cell,k,phi,b);
      f_Pos = fopen(file_name,"w"); 
    }

  files[1] = 1;   //Set files[1]=1 if you want the program to generate and output to a file all the cells' MSD's
  if(files[1]==1)
    {
      sprintf(file_name, "MeanSqd_%d_%g_%g_%g_%g.dat",N_Cell,k,Gamma,phi,b);
      f_MeanSqd = fopen(file_name,"w"); 
    }

  files[2] = 0;   //Set files[2]=1 if you want the program to generate and output to a file all the cells' 
  if(files[2]==1) //correlations(deprecated?)
    {
      sprintf(file_name, "Corr_%d_%g_%g_%g.dat",N_Cell,k,phi,b);
      f_Corr = fopen(file_name,"w"); 
    }

  files[3] = 0;   //Set files[3]=1 if you want the program to generate and output to a file all the cells' 
  if(files[3]==1) //Truncated Mean Velocity Autocorrelation Function 
    {
      sprintf(file_name, "TrMeanVacf_%d_%g_%g_%g.dat",N_Cell,k,phi,b);
      f_TruncMeanVacf = fopen(file_name,"w"); 
    }

  files[4] = 0;    //Set files[4]=1 if you want the program to generate and output to a file all the cells' Mean Speed
  if(files[4]==1)
    {
      sprintf(file_name, "MeanSpd_%d_%g_%g_%g.dat",N_Cell,k,phi,b);
      f_MeanSpd = fopen(file_name,"w"); 
    }

  files[5] = 0;    //Set files[5]=1 if you want the program to generate and output to a file all the cells'
  if(files[5]==1)  //<V^2> relaxation times
    {
      sprintf(file_name, "VRelax_%d_%g_%g_%g.dat",N_Cell,k,phi,b);
      f_VRelax = fopen(file_name,"w"); 
    }

  files[6] = 0;    //Set files[6]=1 if you want the program to generate and output to a file all the cells'
  if(files[6]==1)  //Overlaped Mean Velocity Autocorrelation Function
    {
      sprintf(file_name, "OvrlpMeanVacf_%d_%g_%g_%g.dat",N_Cell,k,phi,b);
      f_OvrlpMeanVacf = fopen(file_name,"w");
    }

  files[7] = 0;    //Set files[7]=1 if you want the program to generate and output to a file all the cells'
                   //Histogram functions

  //sprintf(file_name, "Velocities_%d_%g.dat",N_Cell,k);
  //f_Velocities = fopen(file_name,"w");
    
  return ;
}

//////////////////////////////////////////////////////////////////////////////

unsigned long int rseed;
gsl_rng * rand_vec;
const gsl_rng_type * T;

int i, j;
int cell;

int measure;
int step, tau;

int *cnt;

double d;
double sd;
double v0, u0;
double etapar;
double nx, ny;
double S, P, D;
double vpar_var;
double alpha, theta;
double RNaturalUnits;
double VNaturalUnits;
double V2NaturalUnits;
double MSDNaturalUnits;
double sqrt_dt, inv_dt, dt_2;

double *MSpd;
double *msd, *corr, *vacf;
double *v2, *v2par, *v2perp;

double *abs_vpar;
double *pol_angle;
double *vpar_x, *vpar_y;

double *ovacf;
double **vpar;
double **avacf;
double **rx, **ry;
double **px, **py;
double **rpar,**rperp;

//////////////////////////////////////////////////////////////////////////////
void mean_speed(void);
void histograms(void);
void set_gsl_rng(void);
void free_pointers(void);
void set_constants(void);
void open_output_files(void);
void position_print(int step);
void initial_conditions(void);
void reset_values(int vec_size);
void cell_dynamics(int i, int j);
void parallel_dynamics(int step);
void mean_square_displacement(void);
void allocate_pointers(int vec_size);
void perpendicular_dynamics(int step);
void velocity_relaxation(double delta);
void overlaped_mean_velocity_autocorrelation_function(void);
void truncated_mean_velocity_autocorrelation_function(double delta);

double wrap_angle (double theta);
double truncate(double val, double sig_algarisms);

//////////////////////////////////////////////////////////////////////////////
void set_gsl_rng(void)
{
#ifdef DEBUG
  rseed=0;
#else
  rseed=time(NULL);
#endif

  gsl_rng_env_setup();
  T    = gsl_rng_default;
  rand_vec = gsl_rng_alloc (T);
  gsl_rng_set (rand_vec, rseed);

  return;
}

//////////////////////////////////////////////////////////////////////////////

void set_constants(void)
{
  sqrt_dt      = sqrt(dt);
  inv_dt       = 1.0 / dt;
  dt_2         = dt / 2.0;  
  
  u0           = sqrt(g/(2*(Gamma + k)));

  v0           = v0_mult * u0;
 
  P            = 1.0/(Gamma + 2.0*k);
  D            = g / (2.0*(Gamma + 2.0*k)*(Gamma + k));
  S            = (2.0*q*k*(Gamma + 2.0*k)*(Gamma + k))/(g+(2.0*q*k*(Gamma + 2.0*k)*(Gamma + k)));

  MSDNaturalUnits = 2.0*D*P/(1.0-S);
  RNaturalUnits   = sqrt(MSDNaturalUnits);
  V2NaturalUnits  = 2.0*D/(P*(1.0-S));
  VNaturalUnits  = sqrt(2.0*D/(P*(1.0-S)));

  vpar_var = 0.0;
  
  return ;
}

//////////////////////////////////////////////////////////////////////////////

void cell_dynamics (int i, int j)
{
  double dx, dy;
  double dtheta;
  double xitheta;
  double xi;
  double nx_old, ny_old;
  
  etapar  = sqrt(g)*sqrt_dt*( gsl_ran_gaussian_ziggurat(rand_vec, 1.0)) + b * dt ;

  xitheta = sqrt(2.0*k)*sqrt_dt*gsl_ran_gaussian_ziggurat(rand_vec, 1.0);

  xi = sqrt(2.0*q*k)*sqrt_dt*gsl_ran_gaussian_ziggurat(rand_vec, 1.0);
  
  theta = atan2(py[i][j-1],px[i][j-1]);
  theta = (1.0 - phi * dt)*theta - phi * theta0 * dt + xitheta;
  //theta = theta; //1D model testing
  theta = wrap_angle(theta);
   
  dtheta = xitheta;

  px[i][j] = cos(theta);
  py[i][j] = sin(theta);

  nx_old = nx;
  ny_old = ny;
  
  nx =  py[i][j];
  ny = -px[i][j];
  
  vpar[i][j] = ((1.0-(Gamma)*dt)*vpar[i][j-1] + etapar)*cos(dtheta);
  
  rpar[i][j] = rpar[i][j-1] + vpar[i][j]*dt + etapar;
  rperp[i][j] = rperp[i][j-1] + sqrt(q)*xitheta;

  //dx = (vpar[i][j]*dt + xi )*px[i][j] + sqrt(q)*xitheta*nx_old;
  //dy = (vpar[i][j]*dt + xi )*py[i][j] + sqrt(q)*xitheta*ny_old;
  
  dx = (u0*dt + xi )*px[i][j] + sqrt(q)*xitheta*nx_old;
  dy = (u0*dt + xi )*py[i][j] + sqrt(q)*xitheta*ny_old;
  
  rx[i][j] = rx[i][j-1] + dx;
  ry[i][j] = ry[i][j-1] + dy; 


  //PLEASE, PLEASE, THIS IS THE SECOND TIME THAT I SPENT A LOOOT OF TIME TRYING TO FIGURE WHY THE ANALYTICAL AND NUMERICAL SOLUTIONS WERE DIFFERENT
  //ALWAYS VERIFY IF YOU ARE USING  NATURAL UNITS OR NOT, IN FACT, YOU SHOULD CRETE AN OPTION TO SELECT EITHER OF THEM...
  vpar_x[measure] = vpar[i][j]*px[i][j-1];///VNaturalUnits; 
  vpar_y[measure] = vpar[i][j]*py[i][j-1];///VNaturalUnits;
  abs_vpar[measure] = vpar[i][j];///VNaturalUnits;
  pol_angle[measure] = theta;
  
  vpar_var = vpar_var + vpar[i][j]*vpar[i][j];
  
  measure = measure + 1;
    
  return ;
}

void trajectories (void)
{
  for(step = 1 ; step < iter  ; step++)
    {      
      for(cell = 0 ; cell < N_Cell ; cell++)
	{
	  cell_dynamics(cell,step);
	}
      
      if(files[0]==1)
	{
	  position_print(step-2);
	}
    }
  printf("Finished Trajectories \n");

  printf("Starting Measures \n");

  vpar_var = vpar_var/measure;

  if(files[4]==1)
    {
      mean_speed();
    }

  if(files[1]==1)
    {
      mean_square_displacement();
    }
  
  if(files[5]==1)
    {
      velocity_relaxation(delta);
    }

  if(files[6]==1)
    {
      overlaped_mean_velocity_autocorrelation_function();
    }
  
  if(files[3]==1)
    {
      truncated_mean_velocity_autocorrelation_function(delta);
    }

  if(files[7]==1)
    {
      histograms();
    }
    
  return ;
}


//////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{

  if (argc != 3)
    {
      printf("Program needs 2 arguments:\n 1) Value of N,\n 2) Value of K.\n");
      exit(1);
    }
  N_Cell=atoi(argv[1]);
  k=atof(argv[2]);
  
  printf("\n#Setting GSL RGN...\n");
  set_gsl_rng();
  printf("\n#Done\n");

  printf("\n#Setting Constants...\n");
  set_constants();
  printf("\n#Done\n");

  printf("\n#Opening Files...\n");
  open_output_files();
  printf("\n#Done\n");

  printf("\n#Allocating Pointers...\n");
  allocate_pointers(iter);
  printf("\n#Done\n");

  printf("\n#Reseting Data...\n");
  reset_values(iter);
  printf("\n#Done\n");

  printf("\n#Setting Initial Conditions...\n");
  initial_conditions();
  printf("\n#Done\n");

  printf("\n#Starting Simulation...\n");
  trajectories();
  printf("\n#Done\n");

  printf("\n#Freeing Pointers...\n");
  //free_pointers();
  printf("\n#Done\n");

  printf("\n#Finished Program\n");
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void position_print(int step)
{
  //fprintf(f_Pos,"%.10lf \t %.10lf \t %.10lf \t", (step)*dt/P , rx[0][step]/RNaturalUnits, ry[0][step]/RNaturalUnits);

  /*
  for(cell = 1 ; cell < N_Cell ; cell++)
    {
      fprintf(f_Pos,"%.10lf \t %.10lf \t", rx[cell][step]/RNaturalUnits, ry[cell][step]/RNaturalUnits);
    }
  */

  fprintf(f_Pos,"%.10lf \t %.10lf \t %.10lf \t %.10lf", (step)*dt , rx[0][step], ry[0][step], px[0][step], py[0][step]);
  
  for(cell = 1 ; cell < N_Cell ; cell++)
    {
      fprintf(f_Pos,"%.10lf \t %.10lf \t %.10lf \t %.10lf", rx[cell][step], ry[cell][step], px[cell][step], py[cell][step]);
    }
  

  fprintf(f_Pos,"\n");

  return ;
} 

//////////////////////////////////////////////////////////////////////////////

void mean_square_displacement(void)
{

  int i, j, m, mm;

  double cnt, dr;

  double sum_corr, sum_vacf;

  
  for(i = 0 ; i < N_Cell ; i++)
    {
      for(m = 1 ; m < iter ; m = m + mm + 1)
	{

	  mm = (int)((double)(m)/delay);
	  
	  cnt = 0.0;
	  sd = 0.0;
	  sum_corr = 0.0;
	  sum_vacf = 0.0;
	  
	  for(j = 1000 ; j + m < iter ; j++)
	    {
	      sd = sd + (rx[i][j+m] - rx[i][j])*(rx[i][j+m] - rx[i][j]) + (ry[i][j+m] - ry[i][j])*(ry[i][j+m] - ry[i][j]);

	      sum_vacf = sum_vacf + vpar[i][j] * vpar[i][j+m] * (px[i][j]*px[i][j+m] + py[i][j]*py[i][j+m]);
	      
	      dr = sqrt((rx[i][j+m] - rx[i][j])*(rx[i][j+m] - rx[i][j]) + (ry[i][j+m] - ry[i][j])*(ry[i][j+m] - ry[i][j]));
	      
	      sum_corr = sum_corr + (px[i][j]*(rx[i][j+m] - rx[i][j]) + py[i][j]*(ry[i][j+m] - ry[i][j]))/dr;
	      
	      cnt++;
	      
	    }
	  msd[m]  = msd[m]  + (sd/((double)(cnt*N_Cell)));
	  corr[m] = corr[m] + (sum_corr/((double)(cnt*N_Cell)));
	  vacf[m] = vacf[m] + (sum_vacf/((double)(cnt*N_Cell)));
	  
	  if(m%1000==0)
	    {
	      printf("%d of %d\n", m, iter);
	    }
	}
    }

  for(m = 1 ; m < iter ; m = m + mm + 1)
    {
      mm = (int)((double)(m)/delay);

      //fprintf(f_MeanSqd,"%.10lf \t %.10lf \t %lf \t %lf\n", m*dt/P, msd[m]/MSDNaturalUnits, vacf[m], corr[m]);
      fprintf(f_MeanSqd,"%.10lf \t %.10lf \t %.10lf \t %lf\n", m*dt, msd[m], vacf[m], corr[m]);
    }

  return;
}
//////////////////////////////////////////////////////////////////////////////
void histogram(double *data, int data_size, char *filename)
{
  double max 	 = -1e10;
  double min 	 = 1e10;
    
  double EPSILON  = 1e-6;
  double bin_size;
  double size_factor;

  size_factor = 1.0;
  

  int hist[n_bins+100];
  
  //int data = (int)(data_size/N_Cell);
  int bin;
  int i;
	
  for(bin=0; bin<n_bins+100; ++bin)
    {
      hist[bin] = 0 ;
    }
  
  for (i=0; i < data_size; ++i)
    {
      if (data[i] > max)
	{
	  max = data[i];
	}
      if (data[i] < min)
	{
	  min = data[i];
	}
    }
  	 
  printf("passed through data\n");fflush(stdout);
  max *= (1.0 + EPSILON);
    
  max = max/size_factor;
  min = min/size_factor;
    
  bin_size = (max-min)/n_bins;
  
  printf("making histogram\n");fflush(stdout);
  for(i=0; i < data_size; ++i)
    {
      if(data[i]>min && data[i]<max)
	{
	  bin = (int)((data[i] - min)/bin_size);
	}
      hist[bin] = hist[bin] + 1;
    }
  
  printf("opening file\n");fflush(stdout);
  f_Pdf = fopen(filename, "w");
  
  for(bin=0; bin<n_bins; ++bin)
    {
      fprintf(f_Pdf,"%.15lf \t %.15lf \n", min + (bin + 0.5)*bin_size, 1.0*hist[bin]/(data_size*bin_size));
    }
  
  printf("closing file\n");fflush(stdout);
  fclose(f_Pdf);
	
  return;
}


//////////////////////////////////////////////////////////////////////////////
void histogram_2d(double *datax, double *datay, int data_size, char *filename)
{
  FILE *f_test;
  f_test = fopen("f_vx_f_vy.dat","w");
  double integral_x, integral_y;
  

  double binArea;
  double sum;
  
  double max_x 	 = -1e10;
  double min_x 	 = 1e10;
  double max_y 	 = -1e10;
  double min_y 	 = 1e10;
  
  double EPSILON  = 1e-6;
  double bin_size_x,bin_size_y;
  double size_factor;

  //double cov_xy[][];
  
  size_factor = 1.0;
  
  int binx,biny,hist[n_bins_2d+100][n_bins_2d+100];
  int i;
	
  for(binx=0; binx<n_bins_2d; ++binx)
    for(biny=0; biny<n_bins_2d; ++biny)
      hist[binx][biny] = 0;


  for (i=0; i < data_size; ++i)
    {
      if (datax[i] > max_x)
	{
	  max_x = datax[i];
	}
      if (datax[i] < min_x)
	{
	  min_x = datax[i];
	}
		
      if (datay[i] > max_y)
	{
	  max_y = datay[i];
	}
      if (datay[i] < min_y)
	{
	  min_y = datay[i];
	}
    }
  
	 
  printf("passed through data\n");fflush(stdout);
  max_x *= (1.0 + EPSILON);
  max_y *= (1.0 + EPSILON);
  
  max_x = max_x/size_factor;
  max_y = max_y/size_factor;
  min_x = min_x/size_factor;
  min_y = min_y/size_factor;

  if (max_x>max_y)
    {
      max_y=max_x;
    }
  else if(max_y>max_x)
    {
      max_x = max_y;
    }

  if (min_x>min_y)
    {
      min_x=min_y;
    }
  else if(min_y>min_x)
    {
      min_y=min_x;
    }

  //max_x = 6.0;
  //max_y = 6.0;
  //min_x = -6.0;
  //min_y = -6.0;
      
  bin_size_x = (max_x-min_x)/n_bins_2d;
  bin_size_y = (max_y-min_y)/n_bins_2d;

  binArea = bin_size_x * bin_size_y;
	
  printf("making histogram\n");fflush(stdout);
  for(i=0; i < data_size; ++i)
    {
      if(datax[i]>min_x && datax[i]<max_x)
	{
	  binx = (int)((datax[i] - min_x)/bin_size_x);
	}
      if(datay[i]>min_y && datay[i]<max_y)
	{
	  biny = (int)((datay[i] - min_y)/bin_size_y);
	}
      
      ++hist[binx][biny];
    }

  f_Pdf = fopen(filename,"w");

  sum = 0.0;

  double aux;
  
  FILE *f_Pvx;
  f_Pvx = fopen("Pvx.dat","w");
  
  printf("printing to file\n");fflush(stdout);
  for (binx=0; binx<n_bins_2d; ++binx)
    {
      for (biny=0; biny<n_bins_2d; ++biny)
	{
	  fprintf(f_Pdf,"%.10lf \t %.10lf \t %.10lf \n", min_x + (binx + 0.5)*bin_size_x, min_y + (biny + 0.5)*bin_size_y, 1.0*hist[binx][biny]/(data_size * binArea));
	  
	  sum = sum + 1.0*hist[binx][biny]/(data_size * binArea);

	}
      aux = (((double)(n_bins_2d))-1.0)/2.0;
      fprintf(f_Pvx,"%.10lf \t %.10lf \n", min_x + (binx + 0.5)*bin_size_x, 1.0*hist[binx][(int)(aux)]/(data_size * binArea));
      
      fprintf(f_Pdf,"\n");
    }

  fprintf(f_Pdf,"%lf \n", sum * binArea);
  printf("closing file\n");fflush(stdout);
  fclose(f_Pdf);

  integral_x = 0.0;
  integral_y = 0.0;
  for(i=1; i<n_bins_2d ; ++i)
    {
      for(j=1; j<n_bins_2d; ++j)
	{
	  integral_x = integral_x + 1.0*hist[i][j]/(data_size * sqrt(binArea));// integral_x = integral_x + hist[i][j] * sqrt(binArea) / (data_size * binArea)
	  integral_y = integral_y + 1.0*hist[j][i]/(data_size * sqrt(binArea));
	}
      fprintf(f_test,"%lf \t %lf \t %lf\n", min_x + (i + 0.5)*sqrt(binArea), integral_x, integral_y);
      integral_x=0.0;
      integral_y=0.0;
    }
  
  fclose(f_test);
	
  return;
}

//////////////////////////////////////////////////////////////////////////////

void histograms()
{
  char filename[200];

  int dataSize;
  dataSize = N_Cell*iter;

  double *speed, *avgVelX, *avgVelY, *avgSpeed, *avgVelII, *avgVelT;
  
  speed = (double *)malloc((dataSize) * sizeof(double));
  avgVelX = (double *)malloc((dataSize) * sizeof(double));
  avgVelY = (double *)malloc((dataSize) * sizeof(double));
  avgSpeed = (double *)malloc((dataSize) * sizeof(double));

  avgVelT = (double *)malloc((dataSize) * sizeof(double));
  avgVelII = (double *)malloc((dataSize) * sizeof(double));

  
  sprintf(filename,"hist_vpar_%d_%g.dat",N_Cell,k);
  histogram(abs_vpar,dataSize,filename);

  /*
  for(i=1 ; i<dataSize ; i++)
    {
      speed[i] = sqrt(vpar_x[i]*vpar_x[i] + vpar_y[i]*vpar_y[i]);
    }

  int dJ, l;

  l=0;
  dJ = (int)(delta/dt);

  for(i=0 ; i< N_Cell ; i++)
    {
      for(j=0 ; j<iter ; j++)
	{
	  avgVelX[l] = 0.0;
	  avgVelY[l] = 0.0;

	  avgVelII[l] = 0.0;
	  avgVelT[l] = 0.0;
	  
	  avgSpeed[l] = 0.0;
	  l=l+1;
	}
    }

  l=0;
  
  for(i=0 ; i< N_Cell ; i++)
    {
      for(j=0 ; j<iter-dJ ; j++)
	{
	  //avgVelII[l] = ((rx[i][j+dJ] - rx[i][j]) )/delta; //I still need to make the projections
	  //avgVelT[l] = ((ry[i][j+dJ] - ry[i][j]) )/delta; //I still need to make the projections
	  
	  avgVelX[l] = (rx[i][j+dJ] - rx[i][j])/delta;
	  avgVelY[l] = (ry[i][j+dJ] - ry[i][j])/delta;
	  avgSpeed[l] =  sqrt((rx[i][j+dJ] - rx[i][j])*(rx[i][j+dJ] - rx[i][j]) + (ry[i][j+dJ] - ry[i][j])*(ry[i][j+dJ] - ry[i][j]) ) /(delta*delta);
	  l=l+1;
	}
    }

  sprintf(filename,"hist_2d_%d_%g_%g.dat",N_Cell,k,delta);
  histogram_2d(avgVelX,avgVelY,dataSize,filename);

  sprintf(filename,"hist_avg_spd_%d_%g_%g.dat",N_Cell,k,delta);
  histogram(avgSpeed,dataSize,filename);
  
  */
  //sprintf(filename,"hist_spd_%d_%g.dat",N_Cell,k);
  //histogram(speed,dataSize,filename);
  
  sprintf(filename,"hist_vx_%d_%g.dat",N_Cell,k);
  histogram(vpar_x,dataSize,filename);

  //sprintf(filename,"hist_vy_%d_%g.dat",N_Cell,k);
  //histogram(vpar_y,dataSize,filename);
  
  //sprintf(filename,"pol_angle_%d_%g.dat",N_Cell,k);
  //histogram(pol_angle,dataSize,filename);

  sprintf(filename,"vpar_2d_%d_%g.dat",N_Cell,k);
  histogram_2d(vpar_x,vpar_y,dataSize,filename);
  

  /*
  for(i=1 ; i<dataSize ; i++)
    {
      vx2[i] = vpar_x[i]*vpar_x[i];
      vy2[i] = vpar_y[i]*vpar_y[i];
    }
  
  sprintf(filename,"hist_vx2_%d_%g.dat",N_Cell,k);
  histogram(vx2,N_Cell*iter,filename);

  sprintf(filename,"hist_vy2_%d_%g.dat",N_Cell,k);
  histogram(vy2,N_Cell*iter,filename);
  */
  
  return;
}
//////////////////////////////////////////////////////////////////////////////

void truncated_mean_velocity_autocorrelation_function(double delta)
{
  int i, j, m, mm;

  int n;
  
  double delta2;
  
  double cnt;

  double avx, avy;

  double sum_avacf, sum_avacf1, sum_avacf2, sum_avacf3;

  double trunc_avx1, trunc_avy1, trunc_avx2, trunc_avy2, trunc_avx3, trunc_avy3;

  delta2 = delta*delta;

  n = (int)(delta / dt);
  
  for(i = 0 ; i < N_Cell ; i++)
    {	  
      for(m = 0 ; m < iter ; m = m + mm + 1)
	{
	  mm = (int)((double)(m)/delay);
	  
	  cnt = 0.0;
	      
	  sum_avacf = 0.0;
	  sum_avacf1 = 0.0;
	  sum_avacf2 = 0.0;
	  sum_avacf3 = 0.0;
	      
	  for(j = 0 ; j + n + m < iter ; j++)
	    {
	      trunc_avx1 = ((truncate(rx[i][j+m+n],1.0) - truncate(rx[i][j+m],1.0))*(truncate(rx[i][j+n],1.0) - truncate(rx[i][j],1.0)))/delta2; 
	      trunc_avy1 = ((truncate(ry[i][j+m+n],1.0) - truncate(ry[i][j+m],1.0))*(truncate(ry[i][j+n],1.0) - truncate(ry[i][j],1.0)))/delta2;
	      
	      trunc_avx2 = ((truncate(rx[i][j+m+n],2.0) - truncate(rx[i][j+m],2.0))*(truncate(rx[i][j+n],2.0) - truncate(rx[i][j],2.0)))/delta2; 
	      trunc_avy2 = ((truncate(ry[i][j+m+n],2.0) - truncate(ry[i][j+m],2.0))*(truncate(ry[i][j+n],2.0) - truncate(ry[i][j],2.0)))/delta2;

	      trunc_avx3 = ((truncate(rx[i][j+m+n],3.0) - truncate(rx[i][j+m],3.0))*(truncate(rx[i][j+n],3.0) - truncate(rx[i][j],3.0)))/delta2; 
	      trunc_avy3 = ((truncate(ry[i][j+m+n],3.0) - truncate(ry[i][j+m],3.0))*(truncate(ry[i][j+n],3.0) - truncate(ry[i][j],3.0)))/delta2; 

	      avx = ((rx[i][j+m+n] - rx[i][j+m])*(rx[i][j+n] - rx[i][j]))/delta2;
	      avy = ((ry[i][j+m+n] - ry[i][j+m])*(ry[i][j+n] - ry[i][j]))/delta2;

	      sum_avacf  = sum_avacf  + (avx + avy);
	      sum_avacf1 = sum_avacf1 + (trunc_avx1 + trunc_avy1);
	      sum_avacf2 = sum_avacf2 + (trunc_avx2 + trunc_avy2);
	      sum_avacf3 = sum_avacf3 + (trunc_avx3 + trunc_avy3);
	      
	      cnt++;
	    }
	  avacf[0][m] = avacf[0][m] + (sum_avacf/((double)(cnt*N_Cell)));
	  avacf[1][m] = avacf[1][m] + (sum_avacf1/((double)(cnt*N_Cell)));
	  avacf[2][m] = avacf[2][m] + (sum_avacf2/((double)(cnt*N_Cell)));
	  avacf[3][m] = avacf[3][m] + (sum_avacf3/((double)(cnt*N_Cell)));
	}
    }
      
  for(m = 1 ; m < iter ; m = m + mm + 1)
    {
      mm = (int)((double)(m)/delay);
	  
      fprintf(f_TruncMeanVacf,"%lf \t %.10lf \t %.10lf \t %.10lf \t %.10lf\n", m*dt, avacf[0][m], avacf[3][m], avacf[2][m], avacf[1][m]);
    }
  
  return ;
  
}

//////////////////////////////////////////////////////////////////////////////

void overlaped_mean_velocity_autocorrelation_function(void)
{
  int i, j, m, mm, l;

  double delta2;

  double cnt;

  double avx, avy;
  
  double sum_avacf;
  

  l = (int)(delta / dt);
  delta2 = delta * delta;
  
  for(i = 0 ; i < N_Cell ; i++)
    {
      for(m = 1 ; m < iter ; m = m + mm + 1)
	{
	  mm = (int)((double)(m)/delay);
	  
	  cnt = 0.0;
	  
	  sum_avacf = 0.0;
	  
	  for(j = 0 ; j + l + m < iter ; j++)
	    {
	      avx = ((rx[i][j+l+m] - rx[i][j+m])*(rx[i][j+l] - rx[i][j]));
	      
	      avy = ((ry[i][j+l+m] - ry[i][j+m])*(ry[i][j+l] - ry[i][j]));
	      
	      sum_avacf = sum_avacf + (avx + avy) / delta2;
	      
	      cnt++;
	    }
	  ovacf[m] = ovacf[m] + (sum_avacf/((double)(cnt*N_Cell)));
	}
    }

    
  for(m = 1 ; m < iter ; m = m + mm + 1)
    {
      mm = (int)((double)(m)/delay);
      
      fprintf(f_OvrlpMeanVacf,"%lf \t %lf \n", m*dt, ovacf[m]);
     
    }
   
  return ;
}

//////////////////////////////////////////////////////////////////////////////

void mean_speed(void)
{
  int i, j, m, l;

  double delta;
  
  double sum_speed; 

  for(i = 0 ; i < N_Cell ; i++ )
    {
      delta = 0.0;
      
      for(l = 1 ; l < iter ; l=l+l)
	{
	  delta = ((double)(l)) * dt;
	  sum_speed = 0.0;
	  
	  for(j = 0 ; j + l < iter ; j++)
	    {
	      sum_speed = sum_speed + sqrt((rx[i][j+l] - rx[i][j]) * (rx[i][j+l] - rx[i][j]) + (ry[i][j+l] - ry[i][j])*(ry[i][j+l] - ry[i][j]))/(delta*((double)(iter)));
	    }
	  MSpd[l] = MSpd[l] + sum_speed/((double)(N_Cell));
	}
      printf("%d of %d\n", i, N_Cell);
    }
  
  for(m = 1 ; m < iter ; m = m + m)
    {
      fprintf(f_MeanSpd,"%lf \t %lf \n", m*dt , MSpd[m]);
    }

  return;
}

//////////////////////////////////////////////////////////////////////////////

void velocity_relaxation(double delta)
{
  int i, j, m;
  
  m = (int)(delta/dt);

  fprintf(f_VRelax, "NatTime \t NatVRelax \t NatVParRelax \t NatVPerpRelax \t Time \t VRelax \t VParRelax \t VPerpRelax \n");
  
  for(j = 0 ; j < iter ; j++)
    {
      for(i = 0 ; i < N_Cell ; i++)
	{
	  v2[j] = v2[j] + ((rx[i][j+m]-rx[i][j])*(rx[i][j+m]-rx[i][j]) + (ry[i][j+m]-ry[i][j])*(ry[i][j+m]-ry[i][j]))/(((double)(N_Cell))*delta*delta);
	  
	  v2par[j] = v2par[j] + ((rpar[i][j+m]-rpar[i][j])*(rpar[i][j+m]-rpar[i][j]))/(((double)(N_Cell))*delta*delta);

	  v2perp[j] = v2perp[j] + ((rperp[i][j+m]-rperp[i][j])*(rperp[i][j+m]-rperp[i][j]))/(((double)(N_Cell))*delta*delta);
	}
      fprintf(f_VRelax, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t%lf \t %lf\n", j*dt/P, (v2[j] - u0*u0)/V2NaturalUnits, (v2par[j]-u0*u0)/V2NaturalUnits, v2perp[j]/V2NaturalUnits, j*dt, v2[j]-u0*u0,v2par[j]-u0*u0,v2perp[j]);
    }
  
  return ;
}

//////////////////////////////////////////////////////////////////////////////

void initial_conditions(void)
{
  int i;
  double random;

  //random   = gsl_rng_uniform(rand_vec);
  
  for(i = 0 ; i < N_Cell ; i++)
    {
      random   = gsl_rng_uniform(rand_vec);
      
      px[i][0] = cos(2.0*M_PI*random);
      py[i][0] = sin(2.0*M_PI*random);
      
      px[i][1] = cos(2.0*M_PI*random);
      py[i][1] = sin(2.0*M_PI*random);
      
      rx[i][0] = 0.0;
      ry[i][0] = 0.0;

      rpar[i][1]  = cos((2.0*M_PI)*random)*v0*dt;
      rperp[i][1] = sin((2.0*M_PI)*random)*v0*dt;

      vpar[i][1] = v0;
      
      vpar[i][0] = v0 * gsl_ran_gaussian_ziggurat(rand_vec,1.0);
      vpar[i][1] = v0;
    }

  return ;
}

//////////////////////////////////////////////////////////////////////////////


void reset_values(int vec_size)
{

  for(step = 1 ; step < vec_size ; step++)
    {
      //msd[step] = 0.0;
      
      //vacf[step] = 0.0;
	  
      //corr[step] = 0.0;

      cnt[step] = 0;

      //avacf[0][step] = 0.0;
      //avacf[1][step] = 0.0;
      //avacf[2][step] = 0.0;
      
      for(i = 0 ; i < N_Cell ; i++)
	{

	  px[i][step]    = 0.0;
	  py[i][step]    = 0.0;

	  rx[i][step]    = 0.0;
	  ry[i][step]    = 0.0;

	  vpar[i][step]  = 0.0;

	  rpar[i][step]  = 0.0;
	  rperp[i][step] = 0.0;
	}
    }
  
  return ;
}

//////////////////////////////////////////////////////////////////////////////

void allocate_pointers(int vec_size)
{
  int i;

  cnt       = (int *)malloc(vec_size * sizeof(int));
  
  rx        = (double **)malloc(N_Cell * sizeof(double));
  ry        = (double **)malloc(N_Cell * sizeof(double));
  
  px        = (double **)malloc(N_Cell * sizeof(double));
  py        = (double **)malloc(N_Cell * sizeof(double));

  rpar      = (double **)malloc(N_Cell * sizeof(double));
  rperp     = (double **)malloc(N_Cell * sizeof(double));
    
  vpar      = (double **)malloc(N_Cell * sizeof(double));

  for(i = 0 ; i < N_Cell ; i++)
    {
      rx[i]        = (double *)malloc(vec_size * sizeof(double));
      ry[i]        = (double *)malloc(vec_size * sizeof(double));

      vpar[i]      = (double *)malloc(vec_size * sizeof(double));
      
      px[i]        = (double *)malloc(vec_size * sizeof(double));
      py[i]        = (double *)malloc(vec_size * sizeof(double));
      
      rpar[i]      = (double *)malloc(vec_size * sizeof(double));
      rperp[i]     = (double *)malloc(vec_size * sizeof(double));

    }

  if (files[1] == 1)
    {
      msd          = (double *)malloc(vec_size * sizeof(double));
      vacf         = (double *)malloc(vec_size * sizeof(double));
      corr         = (double *)malloc(vec_size * sizeof(double));
    }

  if(files[3]==1) //Truncated Mean Velocity Autocorrelation Function 
    {
      avacf     = (double **)malloc(4 * sizeof(double));

      for(i = 0 ; i < 4 ; i++)
	{
	  avacf[i]     = (double *)malloc(vec_size * sizeof(double));
	}
    }

  if(files[4]==1)
    {
      MSpd         = (double *)malloc(vec_size * sizeof(double));
    } 
  
  if(files[5]==1)  //<V^2> relaxation times
    {
      v2           = (double *)malloc(vec_size * sizeof(double));
      v2par        = (double *)malloc(vec_size * sizeof(double));
      v2perp       = (double *)malloc(vec_size * sizeof(double));
    }

  if(files[6]==1)  //Overlaped Mean Velocity Autocorrelation Function
    {
      ovacf     = (double *)malloc(vec_size * sizeof(double));
    }
  
  //if (files[7] == 1)
  //  {
  vpar_x       = (double *)malloc((N_Cell * vec_size ) * sizeof(double));
  vpar_y       = (double *)malloc((N_Cell * vec_size ) * sizeof(double));
  abs_vpar     = (double *)malloc((N_Cell * vec_size ) * sizeof(double));
  pol_angle    = (double *)malloc((N_Cell * vec_size ) * sizeof(double));
  //  }

  return ;
}

//////////////////////////////////////////////////////////////////////////////

double wrap_angle (double theta)
{
  return theta - 2.0*M_PI*floor(theta/(2.0*M_PI));
}

//////////////////////////////////////////////////////////////////////////////
double truncate (double val, double sig_algarisms)
{
  double power;

  power = pow(10.0,sig_algarisms);
  
  return (trunc(val * power))/power;
}

//////////////////////////////////////////////////////////////////////////////

void free_pointers (void)
{
  for(i = 0 ; i < N_Cell ; i++)
    {
      free(rx[i]);
      free(ry[i]);
      free(px[i]);
      free(py[i]);
      free(vpar[i]);
    }

  for(i = 0 ; i < 4 ; i++)
    {
      free(avacf[i]);
    }
  
  free(v2);
  free(msd); 
  free(cnt);
  free(MSpd);
  free(vacf);
  //free(ovacf);
  free(v2par);
  
  return ;
}


