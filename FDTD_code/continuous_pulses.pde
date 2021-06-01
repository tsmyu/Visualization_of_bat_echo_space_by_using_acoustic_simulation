
String name = "12000_5000";
String sim_info = "_batG_1st_";
//

int[] O_X1 = new int[]{2500, 2542, 2792, 3325};  //tate
int[] O_Y1 = new int[]{392, 825, 1758, 2600};   //yoko
float[] angle1 = new float[]{91, 95, 93, 93};//float[] angle1 = new float[]{120};
float EndTime = 21.6; // [ms]  // End time of calculation

/*  SELECT  CONDITION   */
boolean ArraymMeasurement = false;
boolean EmitGaussian = false;
boolean BeamWidth = true;
boolean AbsorbWall = false;

boolean DrawJetColormap = false;    // true: MATLAB-like jet colormap / false: Gray scale
boolean BlueColormap = false;
boolean DrawAbsoluteValue = false;      // true: abs(P), false: P
boolean DisplayHelpMessage = true;
/*----END  CONDITION----*/

float angle = 0.0;
float angle_90 = 0.0;
float rad = 0.0;
float rad_90 = 0.0;

int loop_count = 0;
int Thresh = 40;			// Threshold for Binarizing Photo Image (0 to 255)
float image_intensity = 1000;		// Brightness of Acoustic Field (>0)
float model_intensity = 10;		// Brightness of Model (>0)

int WaveformVerticalSize = 150;         // Display size of each waveform [px]
int RandomPointSourceInterval = 2000;	// Interval of Rondom Source Transmission [steps]
int MaxLineSourceNumber = 1000;		// Maximum Number of Points of Line Source

//float dx = 3.0e-4;			// Spatial Resolution [m/pixel]
//float dx = 5.0e-4;   // 40kHz  ... index: 16.59
//float dx = 4.0e-4;  //70kHz    //fine  for40kHz
float dx = 3.0e-4; // CF_FM

//float dt = 8.0e-7;			// Temporal Resolution [s/step]
//float  dt = 4.0e-7; //70kHz
//float dt = 5.0e-7; //50kHz
//float  dt = 0.5e-7; //acryl
float dt = 5.0e-7;
//float dt = 1.5e-7;

float dt_over_dx = dt/dx;		// for more efficient calculation

// setting of emmitted pulse
// pulse duration [ms]
float duration = 2.0;
// pulse frequency [kHz]
int frequency_iFM = 68;
int frequency_tFM = 50;
float frequency_CF = 110;
// sampling frequency [Hz]
int Fs = round(1/dt);
// pulse sample number
//int sample_number_pulse = round(duration/dt/1000);
float nw = 8;
int sample_number_pulse = round(nw/(frequency_CF*1000)/dt);



int intensity=0;
int amp =0;

// Acoustic properties of materials
float rho[]   = { //left=0,right=1
  1.29, 1.29e3
  //1.29, 1.29e5
  //   1.29 , 2.7e3
  //1.29 , 1.18e3   //acryl

};		// Densities [kg/m^3]
float kappa[] = {
  142.0e3, 142.0e6
  //142.0e3, 142.0e8
  //  142.0e3  , 76.0e9
  // 142.0e3  , 8.79e9   //acryl

};	        // Bulk Moduli [Pa]

float velocity = 331.8;
//float freq = 58.0e3;			// Frequency of Initial Waveform [Hz]
float freq = 46.0e3;
float D_from_O_meter = 0.01;         //  中点と両耳の距離を記入(m)
float burst = 1;                      // Wavenumber of Burst [waves]


float lambda = velocity /freq;   // To calculate distance of 2 sources
//float velocity = sqrt(kappa/rho);

//lambda = velocity /freq;

float D_of_2_sources = lambda/dx;
//float D_of_2_sources = 0.0001/dx;
//int D_from_O_2_so = round(D_of_2_sources/8);
int D_from_O_2_so = round(2.5e-3/2/dx);


int absorb_point = D_from_O_2_so -5;

boolean dipole = false;                 // Enable Dipole Transmitter (dipole point source)
int dipole_X = 0;                       // Height of Dipole Transmitter [px]
int dipole_Y = 17;                      // Width of Dipole Transmitter [px]

int NX;					// Spatial Size for X direction [pixels]
int NY;					// Spatial Size for Y direction [pixels]
int O_X;                                // set the origin point  X
int O_Y;                                // set the origin point  Y
int D_from_O;                           // set the distance from origin point

float[][] Vx;				// Particle Velocity for X direction [m/s]
float[][] Vy;				// Particle Velocity for Y direction [m/s]
float[][] P;				// Sound Pressure [Pa]
int[][] Model;				// Model
float[][] Model_intensity;

float[][] Mur_X1;			// Mur's 2nd-Order Absorption Layer
float[][] Mur_X2;			// Mur's 2nd-Order Absorption Layer
float[][] Mur_Y1;			// Mur's 2nd-Order Absorption Layer
float[][] Mur_Y2;			// Mur's 2nd-Order Absorption Layer

int PointSourceX;			// X coordinate of Point Sound Source [pixel]
int PointSourceY;			// Y coordinate of Point Sound Source [pixel]
int[] LineSourceX;			// X coordinate of Line Sound Source [pixel]
int[] LineSourceY;			// Y coordinate of Line Sound Source [pixel]
float[] SourceWaveform;                 // Source waveform

int PointReceiver0_X = 265;             // X coordinate of Point Receiver #0 [pixel]
int PointReceiver0_Y = 384;             // Y coordinate of Point Receiver #0 [pixel]
int PointReceiver1_X = 300;
int PointReceiver1_Y = 300;

int PointReceiver_center_X = 304;
int PointReceiver_center_Y = 312;

int BeamWidth0_X = 12;
int BeamWidth0_Y = 12;
int BeamWidth1_X = 12;
int BeamWidth1_Y = 12;

float[] PointReceiver0_Waveform;        // Received waveform at Point Receiver #0
float[] PointReceiver1_Waveform;
float[] PointReceiver_center_Waveform;
float[][] ArrayReceiver_Waveform;

float ReceiverAmplitude = 3000.0;       // Amplitude of displaying received waveform
boolean WaveformRoll = true;            // true: Roll / false: Sync

int n = 0;
int n_point = 0;
int n_line = 1000;
float time = 0;
float sig_point;
float sig_line;
double sig_gaussian;

int LineSourceNumber = 0;
boolean LineSourceDragging = false;
boolean LineSourceRunning = false;
float col, col_r, col_g, col_b, col_model, col_r1, col_g1, col_b1;
PImage photo;
PImage img;


float  fm_point;
int i = 0;


boolean pause = false;

PrintWriter output;

// for Array

int[] angle_Array;
float[] rad_Array = new float[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//int[] ArrayReceiver_X = new int[]{5100,5100,5100,5100,5100,5100  ,  3850,3850,3850,3850,3850,3850};
//int[] ArrayReceiver_Y = new int[]{42, 63,125,42,63,125  ,  42, 63,125,42,63,125};
int[] ArrayReceiver_X = new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
int[] ArrayReceiver_Y = new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//*/

String O_X2;
String O_Y2;
String angle2;
String output_name_sub;
String output_name;

int no;
String pulse_no;

/* Preparation *********************************************************/
void settings() {
  String photo_file  = name + ".png";      // filename to load

  /* Read Photo */
  photo = loadImage(photo_file);
  NX = photo.height;
  NY = photo.width;
  D_from_O = round(D_from_O_meter / dx );  //set the distance ()

  size(NY, NX);
  //size(NY, NX+WaveformVerticalSize*2);
}
void setup()
{
  println(D_from_O_2_so);
  n=0;
  n_point = 0;
  time = 0;

  // Array Receiver (24     //for Array

  if (ArraymMeasurement) {
    int[] angle_Array = new int[]{0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345};
    //int[] angle_Array = new int[]{180,195,210,225,240,255,270,285,300,315,330,345,0,15,30,45,60,75,90,105,120,135,150,165,180};

    for (int i = 0; i< ArrayReceiver_X.length; i++) {

      rad_Array[i] = angle_Array[i] * PI /180;
      //ArrayReceiver_X[i] = round(O_X1 + (2*D_from_O*  sin((rad_Array[i])) ) );            //sankou
      //ArrayReceiver_Y[i] = round(O_Y1 + (2*D_from_O*  cos((rad_Array[i])) ) );
    }
  }

  /* Allay Allocation */
  Vx = new float[NX+1][NY];
  Vy = new float[NX][NY+1];
  P  = new float[NX][NY];
  Model= new int[NX][NY];
  Model_intensity = new float[NX][NY];
  LineSourceX  = new int[MaxLineSourceNumber];  // plot place
  LineSourceY  = new int[MaxLineSourceNumber];  // plot place

  SourceWaveform = new float[NY];
  PointReceiver0_Waveform = new float[NY];
  PointReceiver1_Waveform = new float[NY];
  PointReceiver_center_Waveform = new float[NY];

  ArrayReceiver_Waveform = new float[NY][12];   // for Array

  Mur_X1 = new float[4][NY];
  Mur_X2 = new float[4][NY];
  Mur_Y1 = new float[NX][4];
  Mur_Y2 = new float[NX][4];



  /* Initialize Field Values */
  InitField(); // call function
  InitMur2nd(); // call function

  /* Create Model */
  for (int i=0; i<NX; i++) {
    for (int j=0; j<NY; j++) {
      color c = photo.pixels[i*NY+j];
      if ( brightness(c) > Thresh ) {
        Model[i][j] = 1;	// dense material
        //Model_intensity[i][j] = 1.0;
      } else {
        Model[i][j] = 0;	// air
        //Model_intensity[i][j] = 0.2;
      }
    }
  }

  ///////////////////////////
  ///////////////////////////

  rad = angle1[loop_count] * PI /180;
  angle = angle1[loop_count] ;//* -1;

  O_X2 = "" +O_X1[loop_count];    // '2'  is  String
  O_Y2 = "" +O_Y1[loop_count];
  angle2 = "" + angle1[loop_count];

  no = loop_count +1;
  pulse_no = "no" + no;
  //String output_name_sub = name +"_"+ O_X2 +"-" + O_Y2+ "__1";
  output_name_sub = name +sim_info +"_"+ O_Y2 +"_" + O_X2+"_deg"+angle2;
  output_name =   output_name_sub + ".csv";

  /* write data */
  //println(output_name);

  output = createWriter(output_name);   //chenge the name every time!!

  // output.println("time[s],Transmitted signal[dBSPL],Received signal1[dBSPL],Received signal2[dB]");
  //センターも録音、
  output.println("time[s],Transmitted signal[dBSPL],Received signal1[dBSPL],Received signal2[dB],center[dB]");
  println("finish setup");
  // output.println("time[s],Transmitted signal[dBSPL],Received signal1[dBSPL],Received signal2[dB],center[dB],42(1m-L), 63(1m-L),125(1m-L),42(1m-R),63(1m-R),125(1m-R)  ,  42(0.5m-L), 63(0.5m-L),(0.5m-L),42(0.5m-R),63(0.5m-R),125(0.5m-R)");
}   //end setup

///////////////////////////
///////////////////////////

/* Main Loop *********************************************************/
void draw()
{
  img = createImage( NY, NX, ALPHA );
  //img = createImage( NY, NX+WaveformVerticalSize*2, ALPHA );
  int[] ArrayReceiver_X = new int[]{5100, 5100, 5100, 5100, 5100, 5100, 3850, 3850, 3850, 3850, 3850, 3850};
  int[] ArrayReceiver_Y = new int[]{42, 63, 125, 42, 63, 125, 42, 63, 125, 42, 63, 125};

  if ( !pause ) {


    // while(loop_count <3){
    //println(output_name);

    // origin point
    O_X = O_X1[loop_count];     //set the origin pfooint
    O_Y = O_Y1[loop_count];



    //angle = angle1[loop_count] ;                // Right is "+"  ..  Left is "-"    degree   (just above is 0 degree)

    rad = angle * PI /180 ;    // convert degree to radian

    //  binaural-ish   (2 receive points)

    PointReceiver0_X = round(O_X - (D_from_O*  sin((rad)) ) );            //LEFT  mic place  (vertical)
    PointReceiver0_Y = round(O_Y - (D_from_O*  cos((rad)) ) );  // (horizontal)
    PointReceiver1_X = round(O_X + (D_from_O*  sin((rad)) ) );            //RIGHT mic place  (vertical)
    PointReceiver1_Y = round(O_Y + (D_from_O*  cos((rad)) ) );  // (horizontal)
    PointReceiver_center_X = O_X;
    PointReceiver_center_Y = O_Y;

    /* Update the Acoustic Field (This is the main part of FDTD !!)  */
    UpdateV();
    UpdateP();

    /* Mur's 2nd Order Absorption */
    Mur2nd();

    /* Initial Waveform from a Point Source (1 pulse of sinusoidal wave with Hann window) */


    /*------------------------------------------------------------------------   */
    //  Emitter Condition

    if (EmitGaussian) {
      if ( n_point < 500) {  // gaussian
        sig_gaussian = 1 * exp( -1* (sq( n_point-40.0) / 2*  sq( 0.1 )) );
        //a.^-((t).^2 / 2* (0.05.^2))
      }
    } else {
      if ( n_point <= sample_number_pulse) {
            if(n_point-round(sample_number_pulse/2)== 0.0){
            sig_point = 1.0;
          }
          else{
            sig_point = sin(2*PI*frequency_CF*1000*(n_point-round(sample_number_pulse/2))*dt)/(2*PI*frequency_CF*1000*(n_point-round(sample_number_pulse/2))*dt);
          //sig_point = sin(2*PI*(frequency_iFM*1000*n_point*dt+(((frequency_tFM-frequency_iFM)*1000/(duration/1000))*sq(n_point*dt))/2));
          //sig_point = sin(2*PI*(frequency_iFM*1000*n_point*dt+(((frequency_tFM-frequency_iFM)*1000/(duration/1000))*sq(n_point*dt))/2)) * (0.5 - 0.5*cos(2*PI*n_point/sample_number_pulse));
          //sig_point = sin(2*PI*frequency_CF*1000*n_point*dt)* (0.5 - 0.5*cos(2*PI*n_point/sample_number_pulse));
          }
      }
      else{
         sig_point = 0;
       }
    }
    println("n_point"+n_point);
    println("sig_point"+sig_point);


    if (dipole) {
      if (PointSourceX+dipole_X < NX && PointSourceY+dipole_Y < NY) {
        P[PointSourceX][PointSourceY] += sig_point;
        P[PointSourceX+dipole_X][PointSourceY+dipole_Y] += sig_point;
      }
    } else if (BeamWidth) {   //you can use here to change a sound source type ( circular piston)

      BeamWidth0_X = round(O_X - ( D_from_O_2_so*  sin((rad)) ) );
      BeamWidth0_Y = round(O_Y - ( D_from_O_2_so*  cos((rad)) ) );
      BeamWidth1_X = round(O_X + ( D_from_O_2_so*  sin((rad)) ) );
      BeamWidth1_Y = round(O_Y + ( D_from_O_2_so*  cos((rad)) ) );

      /* produce 2 points */
      if (EmitGaussian) {
        P[BeamWidth0_X][BeamWidth0_Y] += sig_gaussian;
        P[BeamWidth1_X][BeamWidth1_Y] += sig_gaussian;
      } else {
        P[BeamWidth0_X][BeamWidth0_Y] += sig_point;
        P[BeamWidth1_X][BeamWidth1_Y] += sig_point;
      }
    } else {     // SINGLE EMISSION
      PointSourceX = O_X;
      PointSourceY = O_Y;

      if (EmitGaussian) {
        P[PointSourceX][PointSourceY] += sig_gaussian;
      } else {
        P[PointSourceX][PointSourceY] += sig_point;
      }
    }


    // Record the Source Waveform
    if (EmitGaussian) {
      SourceWaveform[n%NY] = (float)sig_gaussian;    // gauss
    } else {
      SourceWaveform[n%NY] = sig_point;           // non-gauss
    }

    // End  emitter condition
    /*------------------------------------------------------------------------   */



    /* Initial Waveform from a Line Source (1 pulse of sinusoidal wave with Hann window) */
    if ( LineSourceRunning == true) {
      if ( n_line < (1.0/freq)/dt*burst ) {
        sig_line = (1.0-cos((2.0*PI*freq*n_line*dt/burst)))/2.0 * sin((2.0*PI*freq*n_line*dt)) / sqrt(LineSourceNumber); // line source
        for (int i=0; i<LineSourceNumber; i++) {
          P[LineSourceX[i]][LineSourceY[i]] += sig_line;
        }
        SourceWaveform[n%NY] += sig_line;
      } else {
        LineSourceRunning = false;
      }
    }

    /* Record the Reveived Waveform(s) */
    PointReceiver0_Waveform[n%NY] = P[PointReceiver0_X][PointReceiver0_Y];
    PointReceiver1_Waveform[n%NY] = P[PointReceiver1_X][PointReceiver1_Y];
    PointReceiver_center_Waveform[n%NY] = P[PointReceiver_center_X][PointReceiver_center_Y];

    if (ArraymMeasurement) {
      for (int i = 0; i <ArrayReceiver_X.length; i++) {

        if (n>=0) {

          // Record the ARRAY Reveived Waveform(s)
          // ArrayReceiver_X[i] = ArrayReceiver_X[i];
          if ( i <3 ) {
            ArrayReceiver_Y[i] = PointReceiver0_Y +ArrayReceiver_Y[i];
          } else if (i<6) {
            ArrayReceiver_Y[i] = PointReceiver1_Y +ArrayReceiver_Y[i];
          } else if (i<9) {
            ArrayReceiver_Y[i] = PointReceiver0_Y +ArrayReceiver_Y[i];
          } else if (i<12) {
            ArrayReceiver_Y[i] = PointReceiver1_Y +ArrayReceiver_Y[i];
          }
        }
        /*
      println(PointReceiver0_Y);
         println(PointReceiver1_Y);
         println(ArrayReceiver_X[i]);
         println(ArrayReceiver_Y[i]);
         println(i);
         */

        ArrayReceiver_Waveform[n%NY][i] = P[ ArrayReceiver_X[i] ][ ArrayReceiver_Y[i] ];
      }
    }

    //beam


    // Absorb wall

    if (AbsorbWall) {
      if ( n_point < (1.0/freq)/dt*burst +500) {
      //if (n_point < sample_number_pulse){
        for ( int i =-30; i<30; i++) {
          //println(rad);
          //println(tan(rad));
          //println(O_X);
          //println(tan(rad));
          P[ O_X + round(i* sin(rad))][O_Y-1 + round(i* cos(rad))] = 0;
          //P[ O_X + round(i* sin(rad) +1)][O_Y + round(i* cos(rad))] = 0;
          //P[ O_X + round(i* sin(rad) +2)][O_Y + round(i* cos(rad))] = 0;

          //P[ O_X +i][O_Y + round(i* tan(0) -1)] = 0;       //  when → is front(0 deg)
          //P[ O_X +i][O_Y + round(i* tan(0) -2)] = 0;
          //P[O_X +1][O_Y+i] =0;
        }
      }
    }

    /* Draw the Acoustic Field */
    DrawAcousticField();
    /* Draw the Line Source while Mouse Dragging */
    DrawLineSource();

    /* Draw the Waveforms */
    DrawWaveforms();

    /* Output the Image */
    image(img, 0, 0);

    /* Show Text Messages */
    //ShowTextMessages();

    /* Output Step Numbers to Console */
    //println(n, n_point, n_line,);

    if (ArraymMeasurement) {
      //println(time, SourceWaveform[n%NY], PointReceiver0_Waveform[n%NY], PointReceiver1_Waveform[n%NY], ArrayReceiver_Waveform[n%NY][0]);

      /* write data   (5 array)*/
      //output.println(nf(time,1,10)+","+nf(SourceWaveform[n%NY],1,10)+","+nf(PointReceiver0_Waveform[n%NY],1,10)+","+nf(PointReceiver1_Waveform[n%NY],1,10)   + ","+nf(ArrayReceiver_Waveform[n%NY][0],1,10)   +","+nf(ArrayReceiver_Waveform[n%NY][1],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][2],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][3],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][4],1,10)   );
      
      //below must active
      /// output.println(nf(time, 1, 10)+","+nf(SourceWaveform[n%NY], 1, 10)+"," +nf(PointReceiver0_Waveform[n%NY], 1, 10)+","+nf(PointReceiver1_Waveform[n%NY], 1, 10)+ ","+nf(ArrayReceiver_Waveform[n%NY][0], 1, 10)   +","+nf(ArrayReceiver_Waveform[n%NY][1], 1, 10)  +","+nf(ArrayReceiver_Waveform[n%NY][2], 1, 10)  +","+nf(ArrayReceiver_Waveform[n%NY][3], 1, 10)  +","+nf(ArrayReceiver_Waveform[n%NY][4], 1, 10)  +","+nf(ArrayReceiver_Waveform[n%NY][5], 1, 10)  +","+nf(ArrayReceiver_Waveform[n%NY][6], 1, 10)  +","+nf(ArrayReceiver_Waveform[n%NY][7], 1, 10)  +","+nf(ArrayReceiver_Waveform[n%NY][8], 1, 10)  +","+nf(ArrayReceiver_Waveform[n%NY][9], 1, 10)  +","+nf(ArrayReceiver_Waveform[n%NY][10], 1, 10)  +","+nf(ArrayReceiver_Waveform[n%NY][11], 1, 10) );
      
      /* write data   (25 array)*/
      //output.println(nf(time,1,10)+","+nf(SourceWaveform[n%NY],1,10)+"," +nf(PointReceiver0_Waveform[n%NY],1,10)+","+nf(PointReceiver1_Waveform[n%NY],1,10)+ ","+nf(ArrayReceiver_Waveform[n%NY][0],1,10)   +","+nf(ArrayReceiver_Waveform[n%NY][1],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][2],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][3],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][4],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][5],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][6],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][7],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][8],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][9],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][10],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][11],1,10) +","+nf(ArrayReceiver_Waveform[n%NY][12],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][13],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][14],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][15],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][16],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][17],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][18],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][19],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][20],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][21],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][22],1,10)  +","+nf(ArrayReceiver_Waveform[n%NY][23],1,10) +","+nf(ArrayReceiver_Waveform[n%NY][24],1,10));
    } else {
      /* Output Step Numbers to Console */
      //println(n, n_point, n_line,);
      //println(time, SourceWaveform[n%NY], PointReceiver0_Waveform[n%NY], PointReceiver1_Waveform[n%NY]);
      
      //below must active
      ///output.println(nf(time, 1, 10)+","+nf(SourceWaveform[n%NY], 1, 10)+","+nf(PointReceiver0_Waveform[n%NY], 1, 10)+","+nf(PointReceiver1_Waveform[n%NY], 1, 10     ));

      //println(time,SourceWaveform[n%NY], PointReceiver0_Waveform[n%NY], PointReceiver1_Waveform[n%NY],ArrayReceiver_Waveform[n%NY][0]);
    }


    time += dt;

    if ( !pause ) {
      n++;
      n_point++;
      n_line++;

      /* Random Source at given time intervals */
      //if ( n_point > RandomPointSourceInterval && n_line > RandomPointSourceInterval )
      //  RandomPointSource();
    }



    // if (n%(round(0.5e-3/dt)) ==1 ) {
    if (n%1000 ==1 ) {

      DrawAcousticField();
      //DrawWaveforms();
      image(img, 0, 0);
      //ShowTextMessages();
      CaptureToPNG("./cap_"+ name +"/"+ output_name +"/", output_name, i, 3, ".tif");
      i++;
    }


    if (n>=( EndTime *(1.0e-3/dt)) ) {   // break & get out of the Main Loop
      InitField();
      InitMur2nd();
      n= 0;
      n_point = 0;
      i = 0;

      loop_count += 1;

      if ( loop_count < O_X1.length  ) {
        setup();
      } else if (loop_count >= O_X1.length) {
        exit();
      }
    }

    //}//while
  }//if
}    // draw end


/* Sub Routines ************************************************************/

/* Draw the Acoustic Field */
void DrawAcousticField()
{
  float P_max = 0.001;  // 画像化にあたり最大値を設定

  for (int i=0; i<NX; i++) {
    for (int j=0; j<NY; j++) {
      col = (float)(P[i][j]*image_intensity + Model[i][j]*model_intensity);  // Value in Gray scale
      //col = (float)(P[i][j]/P_max + Model[i][j]);
      //col = (float)(P[i][j]/P_max  + Model_intensity[i][j]);


      if ( DrawJetColormap ) {  // conversion into Jet colormap
        if (DrawAbsoluteValue)
          col = min(4.0, abs(col));  // avoid black color
        else
          col = max(min(4.0, col+2.0), 0);  // avoid black color

        col_r = 255*min(max(min(col-1.5, -col+4.5), 0), 1);
        col_g = 255*min(max(min(col-0.5, -col+3.5), 0), 1);// JET color (base: green)
        col_b = 255*min(max(min(col+0.5, -col+2.5), 0), 1);
      } else if (BlueColormap) {
        col = (float)(P[i][j]*image_intensity);

        //col = max(min(4.0, col+2.0), 0);  // avoid black color
        col_r = 255*min(max(min(col-0.5, -col+3.5), 0), 1);

        col_g = 255*min(max(min(col-0.5, -col+0.5), 0), 1);  //base: blue
        col_b = 255*min(max(min(col+0.5, -col+3.5), 0), 1);

        if (Model[i][j] == 1) {

          col_r = col_g = col_b = 0;
        }
      } else {  // conversion into 255 scales
        if (DrawAbsoluteValue)
          //          col_r = 255*min(max(abs(col), 0), 1);
          col_r = 32768*col/2+P_max*32768/2;

        else
          // col_r = 255*col/2 +255/2;
          col_r =255 - 255*col;

        //col_r = 255*min(max(col+0.5, 0), 1);

        col_g = col_r;
        col_b = col_r;
      }
      //img.pixels[i*NY+j] = color(col_r, col_g, col_b);
      img.pixels[i*NY+j] = color(col_r, col_g, col_b);
    }
  }

  if ( DisplayHelpMessage ) {
    if ( 1<PointReceiver0_X && PointReceiver0_X<NX-2 && 1<PointReceiver0_Y && PointReceiver0_Y<NY-2 ) {
      for (int i=-2; i<3; i++) {
        img.pixels[(PointReceiver0_X+i)*NY+(PointReceiver0_Y)] = color(0, 0, 0);
      }
      for (int i=-2; i<3; i++) {
        img.pixels[(PointReceiver0_X)*NY+(PointReceiver0_Y+i)] = color(0, 0, 0);
      }
    }

    if ( 1<PointReceiver1_X && PointReceiver1_X<NX-2 && 1<PointReceiver1_Y && PointReceiver1_Y<NY-2 ) {
      for (int i=-2; i<3; i++) {
        img.pixels[(PointReceiver1_X+i)*NY+(PointReceiver1_Y)] = color(0, 0, 0);
      }
      for (int i=-2; i<3; i++) {
        img.pixels[(PointReceiver1_X)*NY+(PointReceiver1_Y+i)] = color(0, 0, 0);
      }
    }
  }
}   //END  DrawAcousticField()




/* Draw the Line Source while Mouse Dragging */
void DrawLineSource()
{
  if ( LineSourceDragging == true ) {
    for (int i=0; i<LineSourceNumber; i++) {
      img.pixels[LineSourceX[i]*NY+LineSourceY[i]] = color(255, 127, 127);
    }
  }
}

/* Show Text Messages */
void ShowTextMessages()
{
  fill(0);
  if ( DisplayHelpMessage ) {

    text("IMG : " + output_name, 10, 10);  //show filename
    float time_ms = n_point*dt*1.0e3;
    text(time_ms +  " [ms]", 10, 20);
    text("dx : " +dx, 10, 30);
    text("dt : " +dt, 10, 40);
    text("Freq : " +freq +" [Hz]", 10, 50);
    text("angle : " +angle + " [°]", 10, 60);
    text("Origin point : (X,Y) = (" + O_X2 + "," +O_Y2+ ")", 10, 70);
    text("(Intensity , Amp) = ("+ intensity +"," + amp +")", 10, 80);
    if (dipole)
      text("dipole", NY-40, NX-5);
    if (pause)
      text("pause ([space] to restart)", 10, 10);
    text("[ Click or Drag anywhere / 'r': reset / [Ctrl]+click: set sensor / 'd': toggle dipole-monopole ]", 10, NX+WaveformVerticalSize*2-5);
    text("transmitted wave", 10, NX+int(WaveformVerticalSize*1.0/2)-5);
    text("received wave at + (mic1)", 10, NX+WaveformVerticalSize+int(WaveformVerticalSize*0.4/2)-5);
    text("received wave at + (mic2)", 10, NX+WaveformVerticalSize+int(WaveformVerticalSize*1/2)-5);
  }
}

/* Draw the Waveforms */
void DrawWaveforms()
{
  // clear image
  noStroke();
  fill(255);
  rect(0, NX, NY, WaveformVerticalSize*2);// make suquare

  // Draw Source Waveform
  for (int i=1; i<NY-1; i++) {
    stroke(0);
    strokeWeight(2);
    if ( WaveformRoll ) {
      line( i, NX+int(WaveformVerticalSize*1.0/2)-int(60*SourceWaveform[(i+n)%NY]), 
        i+1, NX+int(WaveformVerticalSize*1.0/2)-int(60*SourceWaveform[(i+1+n)%NY]) );
    } else {
      line( i, NX+int(WaveformVerticalSize*1.0/2)-int(60*SourceWaveform[i]), 
        i+1, NX+int(WaveformVerticalSize*1.0/2)-int(60*SourceWaveform[i+1]) );
    }
  }

  // Draw Received Waveform
  for (int i=1; i<NY-1; i++) {
    stroke(0, 127, 0);
    strokeWeight(2);
    if ( WaveformRoll ) {
      line( i, NX+WaveformVerticalSize+int(WaveformVerticalSize*0.4/2)-int(ReceiverAmplitude*PointReceiver0_Waveform[(i+n)%NY]), 
        i+1, NX+WaveformVerticalSize+int(WaveformVerticalSize*0.4/2)-int(ReceiverAmplitude*PointReceiver0_Waveform[(i+1+n)%NY]) );
    } else {
      line( i, NX+WaveformVerticalSize+int(WaveformVerticalSize*0.4/2)-int(ReceiverAmplitude*PointReceiver0_Waveform[i]), 
        i+1, NX+WaveformVerticalSize+int(WaveformVerticalSize*0.4/2)-int(ReceiverAmplitude*PointReceiver0_Waveform[i+1]) );
    }
  }

  // Draw Receiving Point
  if ( DisplayHelpMessage ) {
    if ( 1<PointReceiver0_X && PointReceiver0_X<NX-2 && 1<PointReceiver0_Y && PointReceiver0_Y<NY-2 ) {
      for (int i=-2; i<3; i++) {
        img.pixels[(PointReceiver0_X+i)*NY+(PointReceiver0_Y)] = color(0, 0, 0);
      }
      for (int i=-2; i<3; i++) {
        img.pixels[(PointReceiver0_X)*NY+(PointReceiver0_Y+i)] = color(0, 0, 0);
      }
    }
  }

  //Second mic

  fill(255);
  rect(0, NX+WaveformVerticalSize*2, NY, WaveformVerticalSize);

  for (int i=1; i<NY-1; i++) {
    stroke(127, 127, 127);
    strokeWeight(2);
    if ( WaveformRoll ) {
      line( i, NX+WaveformVerticalSize+int(WaveformVerticalSize*1/2)-int(ReceiverAmplitude*PointReceiver1_Waveform[(i+n)%NY]), 
        i+1, NX+WaveformVerticalSize+int(WaveformVerticalSize*1/2)-int(ReceiverAmplitude*PointReceiver1_Waveform[(i+1+n)%NY]) );
    } else {
      line( i, NX+WaveformVerticalSize+int(WaveformVerticalSize*1/2)-int(ReceiverAmplitude*PointReceiver1_Waveform[i]), 
        i+1, NX+WaveformVerticalSize+int(WaveformVerticalSize*1/2)-int(ReceiverAmplitude*PointReceiver1_Waveform[i+1]) );
    }
  }

  if ( 1<PointReceiver1_X && PointReceiver1_X<NX-2 && 1<PointReceiver1_Y && PointReceiver1_Y<NY-2 ) {
    for (int i=-2; i<3; i++) {
      img.pixels[(PointReceiver1_X+i)*NY+(PointReceiver1_Y)] = color(0, 0, 0);
    }
    for (int i=-2; i<3; i++) {
      img.pixels[(PointReceiver1_X)*NY+(PointReceiver1_Y+i)] = color(0, 0, 0);
    }
  }

  /*

   // Draw "+" at ArrayReceive point
   //for Array
   for (int n = 0; n<25 ; n++){
   
   
   if( 1<ArrayReceiver_X[n] && ArrayReceiver_X[n]<NX-2 && 1<ArrayReceiver_Y[n] && ArrayReceiver_Y[n]<NY-2 ){
   for (int i=-2; i<3; i++) {
   img.pixels[(ArrayReceiver_X[n]+i)*NY+(ArrayReceiver_Y[n])] = color(255,0,0);
   }
   for (int i=-2; i<3; i++) {
   img.pixels[(ArrayReceiver_X[n])*NY+(ArrayReceiver_Y[n]+i)] = color(255,0,0);
   }
   }
   
   }
   */

  // draw cursor
  noStroke();
  fill(127);
  if ( WaveformRoll ) {
    //    rect(NY-3, NX+2, 4, WaveformVerticalSize*2-4);
  } else {
    rect(n%NY, NX+2, 4, WaveformVerticalSize*2-4);
  }
}

void CaptureToPNG(String fileURL, String output_name, int num, int digits, String fileType) {

  String numString = nf(num, digits);

  String pngName = fileURL + output_name + numString + fileType;
  save(pngName);
}

/* Create a Random Point Source */
void RandomPointSource()
{
  PointSourceX = int(random(NX-1));
  PointSourceY = int(random(NY-1));
  n_point = 0;
}

/* Create a Point Source at Mouse Clicking Point */
void mouseClicked()
{
  if ( mouseY < NX && mouseX < NY ) {
    if ( keyPressed == true && ( (key == 's')||(key == CODED && keyCode == CONTROL) ) ) {
      PointReceiver0_X = mouseY;
      PointReceiver0_Y = mouseX;
    } else {
      if ( n_point >= (1.0/freq)/dt ) {	// Only one point source can exist simultaneously
        if ( 0 < mouseY && mouseY < NX && 0 < mouseX && mouseX < NY ) {
          PointSourceX = mouseY;
          PointSourceY = mouseX;
          n_point = 0;
        }
      }
    }
  }
}

/* Create a Line Source on Mouse Dragging Trace */
void mouseDragged()
{
  if ( mouseY < NX && mouseX < NY ) {
    if ( n_line >= (1.0/freq)/dt ) {	// Only one line source can exist simultaneously
      // Start Dragging
      if ( LineSourceDragging == false ) {
        LineSourceDragging = true;
        LineSourceNumber = 0;
      }
      // Add a Point to the Line Source
      if ( LineSourceNumber < MaxLineSourceNumber) {
        if ( 0 < mouseY && mouseY < NX && 0 < mouseX && mouseX < NY ) {
          LineSourceX[LineSourceNumber] = mouseY;
          LineSourceY[LineSourceNumber] = mouseX;
          LineSourceNumber++;
        }
      }
    }
  }
}

/* End Dragging */
void mouseReleased()
{
  if ( LineSourceDragging == true && LineSourceRunning == false && n_line >= (1.0/freq)/dt ) {
    LineSourceDragging = false;
    LineSourceRunning = true;
    n_line = 0;
  }
}

void keyPressed()
{
  /* Toggle Pause / Run */
  if (key == 32) {
    pause = ! pause;
  }

  /* Reset Sound Field */
  if (key == 'r') {
    InitField();
    InitMur2nd();
    n_line = int((1.0/freq*burst)/dt + 1);	// Disable Point Sound Source
    n_point = int((1.0/freq*burst)/dt + 1);	// Disable Line Sound Source
  }

  /* Toggle Color map (Jet and Gray scale) */
  if (key == 'c') {
    DrawJetColormap = ! DrawJetColormap;
  }

  /* Toggle abs / original value */
  if (key == 'a') {
    DrawAbsoluteValue = ! DrawAbsoluteValue;
  }

  /* Toggle Displaying Help Message */
  if (key == 'm') {
    DisplayHelpMessage = ! DisplayHelpMessage;
  }

  /* Toggle Dipole Transmitter */
  if (key == 'd') {
    dipole = ! dipole;
  }

  /* Change the Brightness of Acoustic Field */
  if (key == 'q') {
    image_intensity *= 1.2;
    intensity = intensity +1;
  }
  if (key == 'z') {
    image_intensity /= 1.2;
    intensity = intensity-1;
  }

  /* Change the amplitude of displaying received waveform */
  if (key == 'w') {
    ReceiverAmplitude *= 1.2;
    amp = amp+1;
  }
  if (key == 'x') {
    ReceiverAmplitude /= 1.2;
    amp = amp-1;
  }

  /* Change the Wavenumber of Burst */
  if (key == '1') {
    burst = 1.0;
  }
  if (key == '2') {
    burst = 2.0;
  }
  if (key == '3') {
    burst = 3.0;
  }
  if (key == '4') {
    burst = 4.0;
  }
  if (key == '5') {
    burst = 5.0;
  }
}

void mouseWheel(MouseEvent event) {
  float e = event.getCount();
  if ( keyPressed == true && key == CODED && keyCode == CONTROL) {
    if ( e > 0)  ReceiverAmplitude /= 1.2;
    else        ReceiverAmplitude *= 1.2;
  } else {
    if ( e > 0)  image_intensity /= 1.2;
    else        image_intensity *= 1.2;
  }
}
