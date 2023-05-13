// This code simulates the accumulation of deleterious and beneficial mutations in two clonal populations of fixed size (Npop each) that exchange migrants (nmig)
// Deleterious mutations occur at a rate U, beneficial mutations at a rate Ub
// Deleterious mutations follow an infinite site model but Beneficial mutations follow a finite site model
// the effects of deleterious mutations are taken from a Exponential distribution of mean meanSd
// the effects of beneficial mutations are taken from a Exponential distribution of mean meanSa
// Time of evolution is 400 generations
//
// to compile the code type g++ SocialMicrobiomeExpS400.cpp -o SocialMicrobiomeExpS400 -lm
// to run the code type ./SocialMicrobiomeExpS400 seed Npop Ud Ub meanSd meanSa nmig


#include<math.h>
#include<iostream>
#include<string.h>
#include<stdio.h>
#include<stdlib.h> 
#include<iomanip> 				
#include<fstream>

 
#define PI 3.1415927
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long idum; // seed for the random number generator

struct sequencia{
  double U; // rate deleterious mutations
  double Ub; //rate beneficial mutations
  double sd; // mean deleterious exp distributed
  double *F; // vector of fitnesses of each individual
  int *NMUTDEL; //number of deleterious mutations of each individual
  int *NMUTBEN; //number of beneficial mutations of each individual
  int **SITEBEN; //sites in the genome for beneficial mutations
  double Ftot; // mean population fitness
  double samedio;// mean beneficial exp distributed
    double *sa; // vector of effects of benficial mutations in each site
};

double ran2();
void fitnessinit(struct sequencia *sequence, int Npop); // inicializa o vector F
void mutationselection(struct sequencia *sequence,int Npop,int t, int Genomesite); // geracoes for mutation selection
void migration(struct sequencia *sequence,struct sequencia *sequence1, int Npop,int nmig);
double fat(int m);
int poisson(double mut);


int i1=0;

  int main(int ac, char **av)
{

  FILE *ptt1, *ptt2;
  sequencia sequence, sequence1; // two demes made of sequences, each deme has Npop individuals

 int TSample,Npop,t,tmax,config,conf,k,nmig,*auxDELM,*auxBENM,*auxDELM1,*auxBENM1,nmigmax,i,Genomesite,**auxSITEBENM,**auxSITEBENM1;
 double nmutdel,nmutben, nmutdel1,nmutben1,*auxFM,*auxFM1,MIG,r,saaux;


  char arq1[200],arq2[200];

    if(ac!=8)
        
    { std::cout << "start like this:\n" << av[0]
        << " seed Npop Ud Ua meanSd meanSa nmig \n"
        << std::endl;
        exit(-1);
            }

    tmax = 400; // maximum time of evolution
    nmigmax=500; // maximum number of migrants allowed
    config = 20;  // number of replicates of a simulation with a set of parameters
    Genomesite=200; // number of sites in the genome for beneficial mutations
    
    idum = atoi (av[1]);
    Npop = atoi (av[2]);
    sequence.U = atof (av[3]);
    sequence.Ub = atof (av[4]);
    sequence.sd = atof (av[5]);
    sequence.samedio = atof (av[6]);
    MIG = atof (av[7]);
    
    sequence.F = new double[Npop];
    sequence.NMUTDEL = new int[Npop];
    sequence.NMUTBEN = new int[Npop];
    sequence.Ftot=1.0;
    sequence.SITEBEN = new int*[Npop];
    for(k=0;k<=Npop;k++) sequence.SITEBEN[k]=new int[Genomesite];
    
    sequence1.U = sequence.U;
    sequence1.Ub = sequence.Ub;
    sequence1.sd = sequence.sd;
    sequence1.samedio = sequence.samedio;
    sequence.sa = new double[Genomesite];
    
    sequence1.F = new double[Npop];
    sequence1.NMUTDEL = new int[Npop];
    sequence1.NMUTBEN = new int[Npop];
    sequence1.Ftot=1.0;
    sequence1.SITEBEN = new int*[Npop];
    for(k=0;k<Npop;k++) sequence1.SITEBEN[k]=new int[Genomesite];
    sequence1.sa = new double[Genomesite];
    
    
    auxFM = new double[nmigmax];
    auxDELM = new int[nmigmax];
    auxBENM = new int[nmigmax];
    auxSITEBENM=new int*[nmigmax];
    for(k=0;k<nmigmax;k++) auxSITEBENM[k]=new int[Genomesite];
    
    auxFM1 = new double[nmigmax];
    auxDELM1 = new int[nmigmax];
    auxBENM1 = new int[nmigmax];
    auxSITEBENM1=new int*[nmigmax];
    for(k=0;k<nmigmax;k++) auxSITEBENM1[k]=new int[Genomesite];
    

  
  conf=0;
    TSample=400;
  //popfitaverage=0;
 
  
  while( (++conf)<=config )
    {  
        for( k=0; k<Genomesite; k++ ){sequence.sa[k]=0;sequence1.sa[k]=0;}
            
        printf(" Savector in %d config \n",conf);
        for( k=0; k<Genomesite; k++ )
            {
            r = ran2(); saaux = -log(1.-r)*sequence.samedio;
                printf(" %lf \t", saaux);
            sequence.sa[k] = saaux;
            sequence1.sa[k]= saaux;
             }
        printf(" \n");
      fitnessinit(&sequence,Npop);
      fitnessinit(&sequence1,Npop);
        for(i=0;i<Npop;i++)
            for(k=0;k<Genomesite;k++) {sequence.SITEBEN[i][k]=0;sequence1.SITEBEN[i][k]=0;}
        
        for( k=0; k<nmigmax; k++)
        {
            auxFM[k]=0;auxDELM[k]=0;auxBENM[k]=0;auxFM1[k]=0;auxDELM1[k]=0;auxBENM1[k]=0;
            for(i=0;i<Genomesite;i++){auxSITEBENM[k][i]=0;auxSITEBENM1[k][i]=0;}
        }
  
        
       
      t = 0;
      while( t<tmax )
	{
	  t++;
	  mutationselection(&sequence,Npop,t,Genomesite); // mutation and selection in one population
      mutationselection(&sequence1,Npop,t,Genomesite); // mutation and selection in the other population
        
        if((t%TSample)==0){
            // make stats for calculating frequency of parallel mutations
            int  contafreq=0;
            int contafreq1=0;
            int parallel=0;
            int fixation=0;
            int fixation1=0;
            int countben=0;
            
            sprintf(arq2,"SocialMicro400ExpSf-N%d-Ud%lf-Ua%lf-Sd%lf-Sa%lf-MIG%lf.dat",Npop,sequence.U,sequence.Ub,sequence.sd,sequence.samedio,MIG);
            ptt2 = fopen(arq2,"a");
            fprintf(ptt2," \n %d \t", t);
            for(i=0;i<Genomesite;i++)
            {contafreq=0;contafreq1=0;
                for(k=0;k<Npop;k++) {contafreq=contafreq+sequence.SITEBEN[k][i];contafreq1=contafreq1+sequence1.SITEBEN[k][i];};
                if(contafreq>0 || contafreq1>0)
                {countben=countben+1;
                    printf(" %d %d \t", contafreq, contafreq1);
                    fprintf(ptt2," %lf %lf \t", (double)contafreq/Npop, (double)contafreq1/Npop);
                   }
                if(contafreq>0 && contafreq1>0) parallel=parallel+1;
                if(contafreq==Npop) fixation=fixation+1;
                if(contafreq1==Npop) fixation1=fixation1+1;} //end cycle for genomesite counts
            
            printf(" \n ben %d parallel %d fix %d %d\n",countben, parallel, fixation,fixation1);
            fclose(ptt2);

            sequence.Ftot=0.0; nmutdel=0.0; nmutben=0.0; sequence1.Ftot=0.0; nmutdel1=0.0; nmutben1=0.0;
            for( k=0; k<Npop; k++ )
            {sequence.Ftot += sequence.F[k]/Npop; sequence1.Ftot += sequence1.F[k]/Npop;
                nmutdel += sequence.NMUTDEL[k]; nmutdel1 += sequence1.NMUTDEL[k];
                nmutben += sequence.NMUTBEN[k]; nmutben1 += sequence1.NMUTBEN[k];
            }
            
            printf("%d %d %lf %lf %lf %lf %lf %lf %lf %d %d %d\n",conf,t,sequence.U,sequence.Ub,sequence.sd,sequence.samedio,((double)sequence.Ftot), nmutdel/Npop, nmutben/Npop,parallel,fixation,countben);
            printf("%d %d %lf %lf %lf %lf %lf %lf %lf %d %d %d\n",conf,t,sequence1.U,sequence1.Ub,sequence1.sd,sequence1.samedio,((double)sequence1.Ftot), nmutdel1/Npop, nmutben1/Npop,parallel,fixation1,countben);
            if(sequence.Ftot>0)
            {
                sprintf(arq1,"SocialMicro400ExpSres-N%d-Ud%lf-Ua%lf-Sd%lf-Sa%lf.dat",Npop,sequence.U,sequence.Ub,sequence.sd,sequence.samedio);
                ptt1 = fopen(arq1,"a");
                fprintf(ptt1," %d %lf %d %lf %lf %lf %lf %lf %lf %d %d %d %d\n",conf,MIG,t,((double)sequence.Ftot), nmutdel/Npop, nmutben/Npop,((double)sequence1.Ftot), nmutdel1/Npop, nmutben1/Npop,parallel,fixation,fixation1,countben);
                fclose(ptt1);
                //             exit(-1);
            }
            
        }// close printing
        
        nmig=poisson(MIG);
        if(nmig>nmigmax)nmig=nmigmax;
        // make the first nmig individulas migrate from pop 1 to pop2 and vice versa
        // copy info of migrants from pop1
        for(k=0;k<nmig;k++){
            auxFM[k]=sequence.F[k];auxDELM[k]=sequence.NMUTDEL[k];auxBENM[k]=sequence.NMUTBEN[k];
            for(i=0;i<Genomesite;i++){auxSITEBENM[k][i]=sequence.SITEBEN[k][i];};}
       
        // copy info of migrants from pop2
        for(k=0;k<nmig;k++){auxFM1[k]=sequence1.F[k];auxDELM1[k]=sequence1.NMUTDEL[k];auxBENM1[k]=sequence1.NMUTBEN[k];for(i=1;i<=Genomesite;i++)auxSITEBENM1[k][i]=sequence1.SITEBEN[k][i];}
        
        // replace info of migrants from pop1 to pop2
        for(k=0;k<nmig;k++){sequence.F[k]=auxFM1[k];sequence.NMUTDEL[k]=auxDELM1[k];sequence.NMUTBEN[k]=auxBENM1[k];for(i=1;i<=Genomesite;i++)sequence.SITEBEN[k][i]=auxSITEBENM1[k][i];}
        // replace info of migrants from pop2 to pop1
        for(k=0;k<nmig;k++){sequence1.F[k]=auxFM[k];sequence1.NMUTDEL[k]=auxDELM[k];sequence1.NMUTBEN[k]=auxBENM[k];for(i=1;i<=Genomesite;i++)sequence1.SITEBEN[k][i]=auxSITEBENM[k][i];}
        
          
	} // close cycle for generations
  
    
    }  //close loop for runs

    delete[] auxFM;
    delete[] auxDELM;
    delete[] auxBENM;
  //  for(int k=0; k<Npop; k++) delete[] auxSITEBENM[k];
  
    
    delete[] auxFM1;
    delete[] auxDELM1;
    delete[] auxBENM1;
//    for(int k=0; k<Npop; k++) delete[] auxSITEBENM1[k];

    
    
}



void fitnessinit(sequencia *sequence,int Npop)
{
  int  k,i;
  double r, saux;
    
  for( k=0; k<Npop; k++ ){
    sequence->F[k] = 1.;
    sequence->NMUTDEL[k]=0;
    sequence->NMUTBEN[k]=0;

                           }
    
    
}



void mutationselection(sequencia *sequence, int Npop, int t, int Genomesite)
{
  int i, j,k, m, cont,mb,*auxDEL,*auxBEN,*SauxDEL,*SauxBEN,auxsite,**auxSITEBEN,**SauxSITEBEN ;
  double *auxF,*SauxF, x, soma, mut, Fmax, fitness, saux, mutb, r, sdel;
  
  auxF = new double[Npop];
  auxDEL = new int[Npop];
  auxBEN = new int[Npop];
  SauxF = new double[Npop];
  SauxDEL = new int[Npop];
  SauxBEN = new int[Npop];
  
    auxSITEBEN=new int*[Npop];
    for(k=0;k<Npop;k++) auxSITEBEN[k]=new int[Genomesite];
    SauxSITEBEN=new int*[Npop];
    for(k=0;k<Npop;k++) SauxSITEBEN[k]=new int[Genomesite];
    
    
  mut = sequence->U;
  mutb = sequence->Ub;
  
   // printf("Genomesite %d \n",Genomesite);
    
    // introduce mutation
    for( cont=0; cont<Npop; cont++ )
    {
        auxF[cont]=sequence->F[cont];
        auxDEL[cont]=sequence->NMUTDEL[cont];
        auxBEN[cont]=sequence->NMUTBEN[cont];
        for(k=0;k<Genomesite;k++) auxSITEBEN[cont][k]=sequence->SITEBEN[cont][k];
        
        // introduced deleterious mutations
        m = poisson(mut);
        auxDEL[cont]= auxDEL[cont]+m;
        fitness=auxF[cont];
        for (i=1; i<=m; i++)
        {  r = ran2(); sdel = -log(1.-r)*sequence->sd;
            fitness *= (1-sequence->sd);
        }
        auxF[cont]=fitness;
        // introduced beneficial mutations
        mb=poisson(mutb);
        fitness=auxF[cont];
        for (i=1; i<=mb; i++)
        {   auxsite= (int)(ran2()*Genomesite);
            if(auxSITEBEN[cont][auxsite]==0)
            {
                auxSITEBEN[cont][auxsite]=1;
                fitness *= (1+sequence->sa[auxsite]);auxBEN[cont]= auxBEN[cont]+1;
            }
         }
        auxF[cont]=fitness;
       
     }// end cycle of mutation for all population
    
    // calculate maximum in the population
  
	  Fmax = 0;
	  for( i=0; i<Npop; i++ )
	    {
	      if( auxF[i]>Fmax ) Fmax = auxF[i];
	    }
	  
// select individuals from the mutated population according to their fitness
    cont = 0;
    while( cont<Npop )
    {
        k = (int)(ran2()*Npop); x = Fmax*ran2();
        if( auxF[k]>x )
        {
            SauxF[cont]=auxF[k];
            SauxDEL[cont]=auxDEL[k];
            SauxBEN[cont]=auxBEN[k];
            for(i=0;i<Genomesite;i++) SauxSITEBEN[cont][i]=auxSITEBEN[k][i];
            cont++;
        }
        
    } // end cycle of selection
    
	   
	  for( i=0; i<Npop; i++ )
	    {
	      sequence->F[i] = SauxF[i]; sequence->NMUTDEL[i] = SauxDEL[i]; sequence->NMUTBEN[i] = SauxBEN[i];
            for(j=0;j<Genomesite;j++){sequence->SITEBEN[i][j]=SauxSITEBEN[i][j];
       //                         printf(" %d \t", sequence->SITEBEN[i][j]);
                                      };
     //       printf("\n");
            
	    }

  delete[] auxF;
  delete[] auxDEL;
  delete[] auxBEN;
  delete[] SauxF;
  delete[] SauxDEL;
  delete[] SauxBEN;
    for(int k=0; k<Npop; k++) delete[] auxSITEBEN[k];
    for(int k=0; k<Npop; k++) delete[] SauxSITEBEN[k];

}




int poisson(double mut)
{
    int m;
    double soma, x;
    
    soma = 0; m = 0; x = ran2();
    while( soma<x )
    {
        if( m<=60 )
            soma += pow(mut,(double)m)*exp(-mut)/fat(m);
        if( m>60 )
            soma += exp(-mut+m)*pow((mut/m),(double)m)*pow((2*m*PI),-0.5);
        if( m>2000 )
            soma = 1.0;
        m++;
    }
    m--;
    
    return(m);
    
}


double fat(int m)
{
    int i;
    double produto;
    produto = 1.0;
    for( i=1; i<=m; i++ )
        produto *= i;
    return(produto);
}



double ran2()
{

         int j;
         long k;
         static long idum2=123456789;
         static long iy=0;
         static long iv[NTAB];
         float temp;

         if (idum <= 0) {
                 if (-(idum) < 1) idum=1;
                 else idum = -(idum);
                 idum2=(idum);
                 for (j=NTAB+7; j>=0; j--) {
                         k=(idum)/IQ1;
                         idum=IA1*(idum-k*IQ1)-k*IR1;
                         if (idum < 0) idum += IM1;
                         if (j < NTAB) iv[j] = idum;

                 }
                 iy=iv[0];
         }
         k=(idum)/IQ1;
         idum=IA1*(idum-k*IQ1)-k*IR1;
         if (idum < 0) idum += IM1;
         k=idum2/IQ2;
         idum2=IA2*(idum2-k*IQ2)-k*IR2;
         if (idum2 < 0) idum2 += IM2;
         j=iy/NDIV;
         iy=iv[j]-idum2;
         iv[j]=idum;
         if (iy < 1) iy += IMM1;
         if ((temp=AM*iy) > RNMX) return RNMX;
         else return temp;
}





