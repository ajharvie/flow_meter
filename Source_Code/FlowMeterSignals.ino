#include <ADC.h>

#define ledPin 13

int current_time = 0;
int old_time = 0;
int gap = 0;
int looptime = 10000;
bool wait = true;

ADC *adc = new ADC();


void setup() {
  pinMode(ledPin, OUTPUT); 
  // serial port for monitoring
  Serial.begin(9600);
  digitalWrite(ledPin, HIGH);
}

void loop() {

  while(wait == true){
    current_time = micros();
    gap = current_time - old_time;
    if (gap > looptime){
      old_time = current_time;      
      wait = false;
      break;
    }
  }
  
  ADC::Sync_result result = adc->analogSyncRead(A1, A2);

  int value0 = result.result_adc0; 
  int value1 = result.result_adc1; 

  Serial.print(value0);
  Serial.print(" ");
  Serial.println(value1);

  wait = true; 
  
}