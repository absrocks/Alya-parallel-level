//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include "chiper.h"

#include "cwrapper.h"

void printLicKey(unsigned char* key){
	for (int i=0; i < 32; i++)
    	printf("%02X ",(unsigned)key[i]);
	printf("\n");
}

unsigned char* getLicText(char* name){

	//Open license file
	//printf("-name: -%s-\n",name );
	//FILE *licFile = fopen("/home/testsuite/Alya/Executables/unix/Alya-License.lic", "rb");
	FILE *licFile = fopen(name, "rb");


	//printf("filePointer: %p",licFile);

	//Initialize vars
	int licNum;										//license number necessary to get the key and the IV
	int len = 0;									//Lenght of encripted data.

	fread(&licNum, sizeof(int), 1, licFile); 
	unsigned char* licKey = generateKeyFromLicense(licNum); 	//License Key
	unsigned char* licIV = generateIV(licNum);					//License IV
	unsigned char* licText; 									//Cripted data buffer.


	//Read Data from the file
	fread(&len, sizeof(int), 1, licFile);
	//printf("len: %d\n", len);
	licText = (unsigned char*) malloc(sizeof(unsigned char)*len*10);
	fread(licText, sizeof(unsigned char), len, licFile);

	//printf("licNum: %d\n" , licNum);
	//printLicKey(licKey);
	//printf("licIV: %s\n", licIV);
	//printf("Text: %s\n", licText);

	unsigned char* plaintext;	//data buffer
	plaintext = (char*) malloc(sizeof(char)*len*10);

	decrypt(licText, 
			len, 
			licKey,
			licIV, 
			plaintext);


	//Set the end of the line character.
	plaintext[len] = '\0';

	//finallize
	fclose(licFile);
	free(licText);
	free(licKey);
	free(licIV);

	return plaintext;

}

void licensetopass_(char* licName, int* nameLen, int* initDate, int* endDate, int* actDate, int* lastUse, int* maxHour, int* hourCounter){


	//Convert Fortran string to C string
	char* name;
	name = (char*)malloc(sizeof(char)*(*nameLen));
	for(int i=0; i<*nameLen; i++)name[i] = licName[i];
	name[*nameLen] = '\0';

	//printf("licName: -%s-\n", licName);
	//printf("Name: -%s-\n", name);
	//printf("len: %d\n", *nameLen);

	//Get the text after decripted
	unsigned char* plaintext = getLicText(name);
	
	//printf("decrypted Text: %s\n", plaintext);

	//Process the text and "send" it to fortran.
	*initDate = atoi(plaintext+2);
	*endDate = atoi(plaintext+13);
	if(plaintext[23]=='A'){
		*actDate = atoi(plaintext+24);
		*lastUse = -1;
	}
	else if(plaintext[23] == 'U'){
		*actDate = -1;
		*lastUse = atoi(plaintext+24);
	}
	for(int i=0; i<strlen(plaintext); i++){
		if(plaintext[i] == 'H' && plaintext[i+1] == 'M')
			*maxHour = atoi(plaintext+i+2);
		if(plaintext[i] == 'H' && plaintext[i+1] == 'C')
			*hourCounter = atoi(plaintext+i+2);
	}


	//printf("initDate: %d, endDate: %d actDate: %d lastUse: %d maxHour: %d hourCounter: %d\n", *initDate, *endDate, *actDate, *lastUse, *maxHour, *hourCounter);

	//finallize
	free(name);
	free(plaintext);

}

void cactualize_(char* licName, int* nameLen, double* newhourCounter, int* today){

	//Convert Fortran string to C string
	char* name;
	name = (char*)malloc(sizeof(char)*(*nameLen));
	for(int i = 0; i<*nameLen; i++) name[i] = licName[i];
	name[*nameLen] = '\0';

	//Get the text after decripted	
	unsigned char* plaintext = getLicText(name);
	
	//Read license Number get IV and KEY
	int licNum; 									// license number
	FILE *licFile = fopen(name, "rb");
	fread(&licNum, sizeof(int), 1, licFile);
	fclose(licFile);
	unsigned char* licKey = generateKeyFromLicense(licNum);
	unsigned char* licIV = generateIV(licNum);

	//Get the number of consumed hours and the max hours avaiable
	int maxHour;
	int hourCounter;
	for(int i=0; i<strlen(plaintext); i++){
		if(plaintext[i] == 'H' && plaintext[i+1] == 'M')
			maxHour = atoi(plaintext+i+2);
		if(plaintext[i] == 'H' && plaintext[i+1] == 'C')
			hourCounter = atoi(plaintext+i+2);
	}


	unsigned char licText[500];
	unsigned char licText2[500];

	//Generate the new license Actualized
	memcpy(licText, plaintext, 23);
	sprintf(licText2, "%sU%d:HM%d:HC%d", licText, *today, maxHour, (int)(hourCounter + *newhourCounter));

	//Last Testing stuff
	printf("hourCounter: %d, newhourCounter: %f\n", hourCounter, *newhourCounter);
	printf("hourCounter: + newhourCounter: %d\n", hourCounter + *newhourCounter);
	printf("licText: %s\n", licText2);

	unsigned char licTextCipher[500];

	//Encrypt the license.
	int len = encrypt(	licText2,
						strlen(licText2),
						licKey,
  						licIV,
  						licTextCipher);

	//Write the license with the modifications
	licFile = fopen(name, "wb");
	fwrite(&licNum, sizeof(int), 1, licFile);			//License number
	fwrite(&len, sizeof(int), 1, licFile);				//Lenght of the chiper data
	fwrite(licTextCipher, sizeof(unsigned char), len, licFile);	//Chiper data

	//finallize
	fclose(licFile);
	free(plaintext);
	free(name);
	free(licKey);
	free(licIV);


}
