/******************************************************
 * Written by : Usman Ali Butt                        *           
 * Dated      : 1 August 2018                         *
 * Property off:  www.microcontroller-project.com     *
 *****************************************************/

#include "oled/SSD1306Wire.h"
#include "htu/SparkFunHTU21D.h"

const int trigP = 2;  //D4 Or GPIO-2 of nodemcu
const int echoP = 0;  //D3 Or GPIO-0 of nodemcu

long duration;
float distance;
float distance_p;
float alpha = 0.1;

/************************** Pins **********************************************/
#define SDA_OLED (4)
#define SCL_OLED (15)
#define RST_OLED (16)


/******************************** OLED ****************************************/
SSD1306Wire *display;
char S1[30];
char S2[30];
char S3[30];

void display_init() {
  display = new SSD1306Wire(0x3c, SDA_OLED, SCL_OLED, RST_OLED, GEOMETRY_128_64);
  
  display->init();
  display->flipScreenVertically();
  display->setFont(ArialMT_Plain_10);
  display->drawString(0, 0, "OLED initial done!");
  display->display();
}


void drawFrame(OLEDDisplay *display, int16_t x, int16_t y) {
  display->clear();
  display->setTextAlignment(TEXT_ALIGN_LEFT);
  display->setFont(ArialMT_Plain_10);
  display->drawString(x, y, S1);
  display->drawString(x, y + 15, S2);
  display->drawString(x, y + 30, S3);
  display->display();
}



/******************************** Setup ***************************************/
HTU21D htu;

void setup() {
  /* initialise HTU32D */
  htu.begin(Wire, 4, 15);
  
  pinMode(trigP, OUTPUT);  // Sets the trigPin as an Output
  pinMode(echoP, INPUT);   // Sets the echoPin as an Input
  Serial.begin(115200);      // Open serial channel at 9600 baud rate
  display_init();
}


/******************************** Main ***************************************/
void loop() {
  // get temperature reading
  float temp = htu.readTemperature();
  float v = 331 + 0.6 * temp;

  // get ultrasonic reading
  digitalWrite(trigP, LOW);   // Makes trigPin low
  delayMicroseconds(2);       // 2 micro second delay 
  
  digitalWrite(trigP, HIGH);  // tigPin high
  delayMicroseconds(10);      // trigPin high for 10 micro seconds
  digitalWrite(trigP, LOW);   // trigPin low
  
  duration = pulseIn(echoP, HIGH);     // Read echo pin, time in microseconds
  distance = duration * (v/10000.0) / 2; // Calculating actual/real distance (in cm)
  
  float distance_est = (alpha) * distance + (1 - alpha) * distance_p;
  distance_p = distance;
  
  Serial.print("D = ");      // Output distance on arduino serial monitor 
  Serial.println(distance_est);

  sprintf(S1, "D: %0.1f cm", distance_est);
  sprintf(S2, "T: %0.1f C", temp);
  sprintf(S3, "V: %0.1f m/s", v);
  
  drawFrame(display, 0, 0);
  
  delay(100);                // Pause a little bit and start measuring distance again
}
