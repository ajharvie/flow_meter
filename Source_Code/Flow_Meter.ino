//Microcontroller firmware for flow meter. Uses a modified version of the sliding DFT.


#include <ADC.h>
#include "sliding_dft.h"
#define APin A1 //pin for A sensor
#define BPin A2 //pin for B sensor
#define discardPin 14 //connect to discard valve if using
#define ledPin 13
#define pi 3.141592654 //something to do with circles


//Variables
//Timing

int looptime = 10000; //in micros

int current_time = 0;
int old_time = 0;
int gap = 0;
bool wait = true;

//Speed calculation
float corr[200];
float bigCorr = 0;
int peakLocation;
float phase0;
float phase1;
float phaseGap;

//Stability heuristic calculation
float accVar = 0.05; //unstable with Q above this value
int sampleCheck = 150; //number of samples to wait between comparisons

float corrSum = 0; // I
float corrSumOld = 0; // Check with value from sampleCheck samples ago
int stabCheckCount = 0;
float stabValue;

//Output
float gapVolume = 0.785; //volume between detectors in ul. (~0.785 for 1 mm ID) 

float freqSpace;
float flowSpeed;
float flowRate;
float timeDelay;

//Interpolation
float peakLeft;
float peakCentre;
float peakRight;
float delta;
float preciseFreq;

float gapCentre;
float gapLeft;
float gapRight;

//Discarding waste
int stableWait = 0;
int stableThreshold = 4; //number of "stable" samples required before collecting again

//ADC object
ADC *adc = new ADC();

//Two DFT objects 
static SlidingDFT<float, 1500> dft0;
static SlidingDFT<float, 1500> dft1;

void setup() {
  adc->setResolution(12,0);
  adc->setResolution(12,1);
  //sampling speed
  adc->setSamplingSpeed(ADC_SAMPLING_SPEED::VERY_LOW_SPEED, 0);
  adc->setSamplingSpeed(ADC_SAMPLING_SPEED::VERY_LOW_SPEED, 1);
  pinMode(ledPin, OUTPUT); 
  pinMode(discardPin, OUTPUT);
  // serial port for monitoring
  Serial.begin(9600);
  digitalWrite(ledPin, HIGH);
  

  freqSpace = 100.0/1500.0;
}

void loop() {

  //wait between samples
  while(wait == true){
    current_time = micros();
    gap = current_time - old_time;
    if (gap > looptime){
      old_time = current_time;      
      wait = false;
      break;
    }
  }

  // read sensors and update sDFTS
  ADC::Sync_result result = adc->analogSyncRead(APin, BPin);

  int value0 = result.result_adc0; 
  int value1 = result.result_adc1; 

  dft0.update(value0);
  dft1.update(value1);

  //Find common peak by correlating sDFTs
  for(int i=2; i<200; ++i){ //first 200 bins
    corr[i] = abs(dft0.dft[i])*abs(dft1.dft[i]);
    corrSum = corrSum + corr[i]; //Sort of integrates DFT
  }
  
  bigCorr = 0; //peak
  for (int i=2; i<200;++i){ //first 200 bins
    if (corr[i] > bigCorr){
      bigCorr = corr[i];
      peakLocation = i;
    }
  }

  //"RMS" values of common DFT peak
  peakLeft = sqrt(corr[peakLocation-1]);
  peakCentre = sqrt(corr[peakLocation]);
  peakRight = sqrt(corr[peakLocation+1]);

  //interpolation of frequency peak

  delta = 2*((peakRight-peakLeft)/(2*peakCentre+peakLeft+peakRight)); //factor of 2 because Hann window

  preciseFreq = freqSpace*(float(peakLocation)+delta);
  
 //Linear phase interpolaton
 if(delta > 0.0){ //max to rhs of peak
    gapCentre = arg(dft1.dft[peakLocation])-arg(dft0.dft[peakLocation]);
    gapRight = arg(dft1.dft[peakLocation + 1])-arg(dft0.dft[peakLocation + 1]);

    if (gapCentre < 0.0){
      gapCentre = gapCentre + 2*pi;
    }
    if (gapRight < 0.0){
      gapRight = gapRight + 2*pi;
    }
  
    phaseGap = gapCentre + delta*(gapRight - gapCentre);

 }
 else{ //to LHS of peak
    gapCentre = arg(dft1.dft[peakLocation])-arg(dft0.dft[peakLocation]);
    gapLeft = arg(dft1.dft[peakLocation - 1])-arg(dft0.dft[peakLocation - 1]);

    if (gapCentre < 0.0){
      gapCentre = gapCentre + 2*pi;
    }
    if (gapLeft < 0.0){
      gapLeft = gapLeft + 2*pi;
    }
  
    phaseGap = gapCentre - delta*(gapLeft - gapCentre);

 }
  
  //Calculate time delays and flow rate
  timeDelay = phaseGap/(2*pi*preciseFreq);
  flowSpeed = 1/timeDelay; // mm/s (assuming distance between detectors is precisely 1 mm)
  
  flowRate = (gapVolume*flowSpeed)*60.0; //factor of 60 for units of minutes^-1

  ++stabCheckCount;
  if (stabCheckCount == sampleCheck){
    stabValue = (corrSum - corrSumOld)/corrSum; //dividing to normalise it
    corrSumOld = corrSum;
    stabCheckCount = 0;
      Serial.print(flowRate); // flow rate (depends on calibration parameter for given tubing/block)
      Serial.print(" ");
      Serial.print(flowSpeed,4); // "Speed" of droplets (in mm/s)
      Serial.print(" ");
      Serial.print(stabValue, 4); // Stability heuristic
      Serial.print(" "); 
      if (stabValue > accVar){
        Serial.println("Discard");
        digitalWrite(discardPin, LOW);
        stableWait = 0;
      }
      else if (stabValue < -1*accVar){
        Serial.println("Discard");
        digitalWrite(discardPin, LOW);
        stableWait = 0;
      }
      else{
        if(stableWait == stableThreshold){
          Serial.println("Keep");
          digitalWrite(discardPin, HIGH);
          //stableWait = 0;
        }
        else{
          Serial.println("Wait");
          ++stableWait;
        }
      }
  }
  
  corrSum = 0; 
  wait = true; 
    
}
