#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "chiper.h"

void printLicKey(unsigned char* key);
unsigned char* getLicText(char* name);
void licensetopass_(char* licName, int* nameLen, int* initDate, int* endDate, int* actDate, int* lastUse, int* maxHour, int* hourCounter);
void cactualize_(char* licName, int* nameLen, double* newhourCounter, int* today);