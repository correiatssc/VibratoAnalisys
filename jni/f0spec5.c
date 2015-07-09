
/* 6/8/2003  => correcao de problema qdo dados temporais veem fora de ordem */
/* 8/8/2003  => retira interpolacao no tracado de F0 */
/* 8/9/2003  => Uses DFT for spectral estimation */
/* 9/9/2003  => removes DC level from F0 time series; HiPass filter */
/* 10/9/2003 => control of spectral resolution and frequency band for computation */
/* 10/9/2003 => input length of time window */
/* 14/10/203 => correction of bugs in filtering */
/* 16/10/2003 => output option: data for 3D plot */

#include "stdio.h"
#include "conio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

#define PI  3.14159265358979
#define FILT1 1
#define FILT2 2
#define NOFILTER 0
#define NO 0
#define YES 1

/* FUNCTION PROTOTYPES */
void usage(void),
     msg(void),
     initialization(int, int, int),
     windowing(int, int, int),
     CalArray(int),
     removeDC(int, int);
int dft(int);


/* FILE RELATED */
FILE *fdr = NULL;			/* data file header */
float *f = NULL;			/* processed f0 samples buffer */
					/* contents changes along processing */
float *f_orig = NULL;			/* copy of original f0 series */
float  ni;				/* nb. of samples actually read */

FILE *fpw = NULL;                      	/* header for output file */
FILE *picsFile = NULL;                 	/* header for output file (picos)*/
float *yim = NULL,			/* buffer for spectrum imag */
      *yre = NULL,                      /* ..., real */
      *ymag = NULL,			/* ..., magnitude */
      *f_filt = NULL;			/* buffer for filtering */

/* TIME INTERVALS */
int TAMX,			/* max length of time window */
    TAMY,			/* # points of dft */
    N,				/* # of pair of points in F0 series */
    window_type = 'h',  	/* default */
    n_frames = 0,
    flag = 0,
    F0filter = FILT1,
    offset, shift,
    specgram = NO;

float delta_t,
      delta_f = 0.1,
      time_window,
      Fmax,                     /* depende de sampling time */
      t1, f1,
      t2, f2,
      A, F0 = 220,
      pic_limen = 1.0,		/* noise floor for peak search */
      aux, aux1,
      file_dur,			/* length of file (sec) */
      DC_value,
      DC_count = 0.0,
      w_DC_value,		/* DC value (Mean F0) of each time window in spectrogram */
      amp_factor;		/* amplitude correction; depends on window shape */


char time[15], freq[15];

void main(argc,argv)
int argc;
char *argv[];

{
  int i, j, k;

  int input_arg = -1,
      output_arg = -1,
      pics_file_arg = -1,
      delta_f_arg = -1,
      Fmax_arg = -1,
      time_window_arg = -1,
      F0filter_arg = -1,
      window_arg = -1,
      pic_limen_arg = -1,
      calibration_arg = -1,
      delta_t_arg = -1,
      shift_arg = -1,
      specgram_arg = -1;

  /* CHECK PROGRAM CALL */
  if(argc<2)  usage();

  for(i=1; i<argc && *argv[i] == '-'; i++) {
    switch(argv[i][1]) {
      default:
	if (argc <= i+1) usage();

	switch(argv[i++][1]) {
	  case 'c':   /* calibration */
	    calibration_arg = i;
	    F0 = atof(argv[i]);
	    if(F0<0) usage();
	    break;
	  case 'i':
	    input_arg = i;
	    break;
	  case 'o':
	    output_arg = i;
	    break;
	  case 'p':
	    pics_file_arg = i;
	    break;
	  case 'r':
	    delta_f_arg = i;
	    delta_f = atof(argv[delta_f_arg]);
	    if (delta_f < 0 )  usage();
	    break;
	  case 'x':
	    time_window_arg = i;
	    time_window = atof(argv[time_window_arg]);
	    if (time_window < 0)  usage();
	    break;
	  case 'f':
	    F0filter_arg = i;
	    F0filter = atoi(argv[F0filter_arg]);
	    if(!(F0filter == FILT1 || F0filter == FILT2
				   || F0filter == NOFILTER)) usage();
	    break;
	  case 't':
	    delta_t_arg = i;
	    delta_t = atof(argv[delta_t_arg])/1000.0;
	    if ( delta_t < 0 )  usage();
	    break;
	  case 'u':
	    Fmax_arg = i;
	    Fmax = atof(argv[Fmax_arg]);
	    break;
	  case 'w':
	    window_arg = i;
	    window_type = (int) *argv[window_arg];

	    if( (window_type != 'r') && (window_type != 'h') &&
		(window_type != 't') ) usage();
	    break;
	  case 'l':
	    pic_limen_arg = i;
	    pic_limen = atof(argv[pic_limen_arg]);
	    break;

	  case 's':
	    shift_arg = i;
	    shift = atoi(argv[shift_arg]);
	    if(shift<0 || shift > 20) usage();
	    break;

	  case 'g':
	    specgram_arg = i;
	    specgram = atoi(argv[specgram_arg]);
	    if(!(specgram == 0 || specgram == 1)) usage();
	    break;

	  default:
	    usage();
	    break;
	}
	break;
    }
  }
  if ((i != argc && *argv[i] != 'i') || input_arg == -1 ||
       output_arg == -1) usage();

  /* OPEN F0 TIME SERIES */
  if((fdr = fopen(argv[input_arg], "rt")) ==  NULL) {
    printf("Error while opening file %s\n", argv[input_arg]);
    exit(0);
  }

  /* OPEN OUTPUT  FILE (spectrum) */
  if((fpw = fopen(argv[output_arg], "wt"))==NULL) {
    printf("Error while creating file (%s)\n", argv[output_arg]);
    exit(1);
  }

  /* OPEN OUTPUT  FILE (picos) */
  if(pics_file_arg != -1) {
    if((picsFile = fopen(argv[pics_file_arg], "wt"))==NULL) {
      printf("Error while creating spectral-peaks file\n");
      exit(1);
    }
  }

  /* INITIALIZE VARIABLES, BUFFERS, ETC */
  fscanf(fdr, "%s%s\n", freq, time);  /* 1st line:  0.0 dur */
  file_dur = atof(time);
  initialization(delta_t_arg, Fmax_arg, time_window_arg);

  do { /* skip initial -1 -1 */
    fscanf(fdr, "%s%s\n", time, freq);
    t1 = atof(time);  f1 = fabs(atof(freq));
  } while(t1 == -1.0);
  f[0] = f1; DC_value = f1; DC_count = 1.0;

  /* PUT MESSAGES ON THE SCREEN */
  msg();
  printf("TAMX = %i (%4.2fs), shift = %5.2fs, TAMY = %i, in = %s, out = %s,\n",
		 TAMX, TAMX*delta_t,  shift*delta_t,    TAMY, argv[input_arg], argv[output_arg]);

  printf("time_window = %5.2fs, delta_t = %5.2fs, delta_f = %5.2fHz, F0filter = %i\n",
			time_window,  delta_t,   delta_f,          F0filter);
  if(window_type == 'r') printf("rectangular window\n");
  else if(window_type == 'h') printf("hamming window\n");
  printf("file_dur = %fs\n",file_dur);

  /*---------REM0VE DC VALUE AND LOWPASS FILTER F0 SERIES ---------------*/
  j=1;
  while( (ni = (float) fscanf(fdr, "%s%s\n", time, freq) ) == 2) {
    t2 = atof(time); f2 = fabs(atof(freq));

    if(t2 == (float) -1.0) { /* end of file */
      if(j>TAMX/2) {
	if(calibration_arg != -1) CalArray(N);
	removeDC(N, F0filter);
      }
      break;
    }
    else { /* continua a leitura */
      f[j] = f2;
      f_orig[j] = f2;
      j++;

      DC_value += f2; DC_count += 1.0;

      if(j>=N) {
	printf("j=%i, break\n", j);
	if(calibration_arg != -1) CalArray(N);
	removeDC(N, F0filter);
      }
    }
  }


  /*------------------ SPECTRAL ANALSYS ------------------------------*/
  offset = 0.0;
  while(offset + TAMX < N) {
    windowing(offset, TAMX, window_type); /* windowed series in f[0,...,TAMX-1] */
    dft(TAMX);
    offset += shift;
  }

  if(offset + TAMX >= N) { /* remaining block */
    windowing(offset, TAMX, window_type); /* windowed series in f[0,...,TAMX-1] */
    dft(TAMX);
  }

  /******************** CALCULATE AVERAGE SPECTRUM ************************/
  if(n_frames>0) {
    flag = 1; /* indica que ha dados do espectro */
    /* DC value */
    ymag[0] = ymag[0] / (float) n_frames;

    for(i=1; i<TAMY; i++) {
      ymag[i] = (float) ymag[i] / (float) n_frames;
    }

    /* write to file */
    if(specgram == YES) {
      fprintf(fpw, "%5.2f ", (float) (offset+shift)*delta_t);
    }

    for(i=0; i<TAMY; i++) {
      if(specgram == YES) {
	fprintf(fpw, "%5.2f ", (float) ymag[i]);
      }
      else {
	fprintf(fpw, "%5.2f %5.4f\n", (float) i*delta_f, ymag[i] );
      }
    }
    if(specgram == YES) fprintf(fpw, "\n");
  }
  else fprintf(fpw, "%5.2f %5.2f\n", (float) -1.0, -1.0);


  /* busca picos ; yre[] armazena a magnitude */
  if(flag!=0 && pics_file_arg != -1 ) {

    for(i = (int) (2.0/delta_f); i<TAMY - 1; i++) {
      if(ymag[i] > pic_limen) {
	if(ymag[i]>ymag[i+1] && ymag[i]>ymag[i-1]) {
	  fprintf(picsFile, "%5.2f %6.4f\n", (float) i*delta_f, ymag[i]);
	}
      }
    }
  }


  free(f); free(f_filt); free(f_orig);
  free(yre); free(yim); free(ymag);
  fclose(fdr);
  fclose(fpw);
  if(pics_file_arg != -1) fclose(picsFile);
  exit(0);
}


/*---------------------------- FUNCTIONS ---------------------------*/

void usage()
{
  msg();
  printf("\n");
  printf("fm -i file1  -o file2.lts [- arg {descr (defaults)}]\n");
  printf("-i * {* = input file  (F0 time series)}\n");
  printf("-o * {* = output file (freq x magnitude or freq x mod. index)}\n");
  printf("-g * {* = output for 3D graph(specgram), 0 = no, 1 = yes (0)}\n");
  printf("-p * {* = estimated peaks file (freq x magnitude) }\n");
  printf("-r * {* = spectral beam spacing in Hz  (0.1)}\n");
  printf("-u * {* = upper frequency in Hz; * < Samp. Freq/2  (15)}\n");
  printf("-x * {* = length of time window in sec.  (2.0); 0 = File length}\n");
  printf("-s * {* = shift time window in samples, range: 1 ... 20,  (1.0)}\n");
  printf("-f * {HP filter(0: no, 1: Fc = 1Hz/25Hz, 2: Fc = 2Hz/25Hz) (1)}\n");
  printf("-w * {* = window type: 'r'ectangular, 't'riangular, or 'h'amming (h)}\n");
  printf("-l * {* = limen for spectral peak search (1.0)}\n");
  printf("-t * {* = time series sampling time (= 1/Samp. Freq), in ms (20)}\n");
  printf("-c * {calibration, FM = F0 + SUM[A(f)*fm*cos[2*pi*f*t]\n");
  printf("      * = F0 (220 Hz), (A ua,f Hz) = (5,5); (10,10); (15,15)}\n");
  exit(0);
}


void initialization(int delta_t_arg, int Fmax_arg, int time_window_arg)
{
  int i;

    printf("x4a = %f\n", time_window);
  /* set delta_t */
  if(delta_t_arg == -1) delta_t = 0.020;  		/* default */

  shift = 2; /* input: number of time intervals */
  printf("shift = %i\n", shift);

  /* set Fmax and TAMY */
  if(Fmax_arg == -1) Fmax = 15;                         /* default */
  if (Fmax > 1.0/(2.0*delta_t) )  usage();
  TAMY =  (int) ceil(Fmax/delta_f);

  /* set TAMX */
  if(time_window_arg == -1) time_window = 2;             /* default */
  if(time_window > file_dur || time_window == 0) time_window = file_dur;
  TAMX =  (int) ceil(time_window/delta_t);
  N = file_dur/delta_t;

  /* CREATE BUFFERs */
  if((f = malloc( sizeof(float) * (N+8) )) == NULL) {
    printf("out of memory in call to malloc(f).\n");
    exit(1);
  }

  if((f_orig = malloc( sizeof(float) * (N+8) )) == NULL) {
    printf("out of memory in call to malloc(f_orig).\n");
    exit(1);
  }

  if((f_filt = malloc( sizeof(float) * (N+8) )) == NULL) {
    printf("out of memory in call to malloc(f).\n");
    exit(1);
  }

  /* CREATE OUTPUT BUFFER */
  if((yre = malloc( sizeof(float) * (TAMY + 8) )) == NULL) {
    printf("out of memory in call to malloc(yre).\n");
    exit(1);
  }

  if((yim = malloc( sizeof(float) * (TAMY + 8) )) == NULL) {
    printf("out of memory in call to malloc(yim).\n");
    exit(1);
  }

  if((ymag = malloc( sizeof(float) * (TAMY + 8) )) == NULL) {
    printf("out of memory in call to malloc(ymag).\n");
    exit(1);
  }

  /* INITIALIZE BUFFERS */
  for(i=0; i < N + 8; i++) {
    f[i] = 0.0; f_filt[i] = 0.0;
  }

  for(i=0; i<TAMY + 8; i++) {
    yre[i] = yim[i] = ymag[i] = 0.0;
  }

  /*........ magnitude correction ..........*/
  if(window_type == 'r') amp_factor = 2.0;
  else if(window_type == 'h')  amp_factor = 2.0/0.54; /* teorico = 3.7037; empirico: 3.7515 */
  else amp_factor = 4.0;                             /* triang., teorico = 2.0/0.5 */

}


void msg(void)
{
  printf(" \n (C) Maurilio N. Vieira, 15 Oct 2003 \n");
  printf(" F0 time-series spectral analysis\n");
}


void windowing(int offset, int TAMX, int win_type)
/* applies time window to a block of the filtered f0 time series;
   results are stored in buffer f[]
*/
/******* global variables *****************************
 int f[], f_filt[], TAMX
 float w_DC_value
*******************************************************/

{
  int i;
  float tmp;

  w_DC_value = 0.0;
  for(i=0; i<TAMX; i++) {
    w_DC_value += f_orig[i+offset];
  }
  w_DC_value /= TAMX;

  switch(win_type) {
    case 'r':
    break;

    case 'h':
    for(i=0; i<TAMX; i++) {
      f[i] = f_filt[i+offset]*( 0.54 - 0.46*cos(2.0*PI*i/(TAMX-1.0)) );
    }
    break;

    case 't':
      for(i=0; i<TAMX/2; i++) {
	tmp = 2.0*(i+1)/TAMX;
	f[i] = f_filt[i+offset]*tmp;
	f[TAMX-1-i] = f_filt[offset+TAMX-1-i]*tmp;
      }
    break;
  }
}


/* dados no array global f[] de dimensao = length */

dft(int Length)
/******* global variables *****************************
 int f[], DC_value, n_frames, TAMY, specgram, offset
 float delta_f, delta_t, w_DC_value
*******************************************************/
{
  register int n, k;
  float M, dw;

  /*...................................................................
    delta_f foi passado c/ relacao a Fs/2 e nao Fs, como usado na DFT.
    Assim, sendo Fs = 1/delta_t,
	   delta_f = (Fs/2)/M' => M' = 1/(2*delta_f*delta_t)
	   M = 2*M' = 1/(delta_f*delta_t)
  ....................................................................*/
  M = 1.0/(delta_t*delta_f);
  dw = 2.0*PI/M;

  /*............ DFT ................*/
  /* DC value dealt with separetely */
  ymag[0] += DC_value;

  for(k=1; k<TAMY; k++) {
    yre[k] = 0.0; yim[k] = 0.0;
    for(n=0; n<Length; n++) {
      yre[k] += f[n]*cos(k*dw*n);
      yim[k] -= f[n]*sin(k*dw*n);
    }
  }

  /* specgram file format
  t1 ymag[0] ymag[1] .... ymag[TAMY-1]
  t2 ymag[0] ymag[1] .... ymag[TAMY-1]
		     ...
  tN ymag[0] ymag[1] .... ymag[TAMY-1] => last line: average values
  */


  if(specgram == YES) {
    fprintf(fpw, "%5.2f %5.2f ", (float) offset*delta_t, w_DC_value);
  }

  for(k=1; k<TAMY; k++) {
    M =  (float) amp_factor*sqrt(yre[k]*yre[k] + yim[k]*yim[k])/(float) Length;
    ymag[k] += M;

    if(specgram == YES) fprintf(fpw, "%5.2f ", (float) M);
  }

  if(specgram == YES) fprintf(fpw, "\n");

  n_frames ++;
}

void CalArray(int N)
/* replaces F0 series by simulated values for calibration purposes */
{
  /*************** global variables ***************************
     (float) f[] = array,
     (float) DC_value, (float) DC_count
  ************************************************************/
  int i;
  for(i=0; i<N; i++) {
    f[i] = (int) (F0 +  10*sin(2*PI*4.0*delta_t*i) +
			5*sin(2*PI*8.0*delta_t*i) +
			2.5*sin(2*PI*12*delta_t*i));
  }

  /* p/ funcionar com removeDC() */
  DC_count = N; DC_value = F0*DC_count;
}


void removeDC(int N, int F0filter)
{
  /*************** global variables ***************************
     (float) DC_value, (float) DC_count
     (float) f[], f_filt[]
  ************************************************************/

  int i, j;

  /*-------------- High Passfilter=> F0Hipass.m ----------------------
   % design hipass filter for f0 time series sampled at 50 Hz
   wn = 1/25; % cutoff frequency at 2 Hz
   [B,A] = cheby2(6,45,wn,'high');
   [h,w] = freqz(B, A, 256);
   mag = 20*log10(abs(h));

   for i=1:256
     w(i) = 25*i/256;
   end

   save c:\maurilio\estudantes\lorena\f0hipass.ch2 B A /ascii /double;

   plot(w,mag), grid;
   xlabel('Frequencia (Hz)');
   ylabel('Magnitude (dB)');
   title('F0Hipass.m');
  ------------------------------------------------------------------*/
  double A[6+1], B[6+1];

  /* lowpass filter, fc = 1 Hz, Fs = 50 Hz */
  double B1[6+1] = {
     7.5782251288852620e-001,
    -4.5289904018779690e+000, 1.1295665290541940e+001,
    -1.5048994709664700e+001, 1.1295665290541930e+001,
    -4.5289904018779630e+000, 7.5782251288852440e-001
  };

  double A1[6+1] = {
       1.0,
      -5.4238946867157000e+000,  1.2302437870616220e+001,
      -1.4933646178677020e+001,  1.0230251036630820e+001,
      -3.7494263866010280e+000,  5.7429496104074540e-001
  };

  /* lowpass filter, fc = 2 Hz, Fs = 50 Hz */
  double B2[6+1] = {
     5.7491756551221150e-001,
    -3.3951046323898120e+000, 8.4074507212169750e+000,
    -1.1174522744874880e+001,  8.4074507212169750e+000,
    -3.3951046323898130e+000,  5.7491756551221160e-001
  };

  double A2[6+1] = {
       1.0,
      -4.8155996555846370e+000,  9.8210339253067580e+000,
      -1.0836619424565210e+001,  6.8135759448765440e+000,
      -2.3121094254845010e+000,  3.3053020729528700e-001
  };


  /* remove DC value */
  DC_value = DC_value/DC_count;
  printf("DC = %5.2fHz\n", DC_value);

  for(i=0; i<N; i++) {
    f[i] = f[i] - DC_value;
  }

  switch(F0filter) {
    case 0: /* no filter */
    break;

    case FILT1:
      for(i=0; i<=6; i++) {
	A[i] = A1[i];
	B[i] = B1[i];
      }
    break;

    case FILT2:
      for(i=0; i<=6; i++) {
	A[i] = A2[i];
	B[i] = B2[i];
      }
    break;
  }

  if(F0filter == 0) { /* only DC removal; no LP filter */
    for(i=0; i<N; i++) {
      f_filt[i] = f[i];
    }
  }
  else if(F0filter == FILT1 || F0filter == FILT2) { /* Executa filtro */
    /* Desloca (p/ direita) amostras no buffer f[] p/ filtragem */
    for(i=N-1; i>=0; i--){
      f[i+6] = f[i];
    }

    /* zero initial conditions */
    for(i=0; i<6; i++) {
      f[i] = 0.0;
      f_filt[i] = 0.0;
    }

    /* filtragem*/
    for(i=6; i<N -1 + 6; i++) {
      f_filt[i] = 0.0;
      /*---- zeros ------*/
      for(j=0; j<=6; j++) {
	f_filt[i] = f_filt[i] + B[j]*f[i-j];
      }

      /*------poles -----*/
      for(j=1; j<=6; j++) {
	f_filt[i] = f_filt[i] - A[j]*f_filt[i-j];
      }
    }

    /* recoloca amostras filtradas com offset 0 no buffer f_filt[] */
    for(i=6; i<N+6; i++) {
      f_filt[i-6] = f_filt[i];
    }
  }
}

