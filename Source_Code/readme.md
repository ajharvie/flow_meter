Source code for the flow meter. Three versions are provided:
1) Flow_Meter.ino is the version used to generate the data shown in the paper, using a teensy 3.6
2) FlowMeterUPDATED.ino is updated for the current (Feb 2020) version of teensyduino (v1.51). There are a few changes, specifically removal of a few functions related to the ADC library. All settings are otherwise the same.
3) FlowMeterHighRate.ino is updated for the capability of the new Teensy 4 board - the sample rate has been increased 10x, and the number of samples considered in the calculation has been increased 2.5x. This allows for accurate measurement of much higher flow rates.
