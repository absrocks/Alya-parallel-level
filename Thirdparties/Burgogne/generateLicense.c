#include <stdio.h>
#include <stdlib.h>
#include "chiper.h"


void printLicKey(unsigned char* key){
	for (int i=0; i < 32; i++)
    	printf("%02X ",(unsigned)key[i]);
	printf("\n");
}

void help(){
	printf("\n----------------License Generator for ALYA:----------------\n\n");
	printf("Usage:  \n\n./generateLicense initDate finalDate activationDate MaxHour\n");
	printf("Date format: YYYYMMDD\n\n");
	abort();
}

int main(int argc, char *argv[]){


	//Check Parameteers
	if(argc < 5)help();

	for(int i=1; i<4; i++)
		if(strlen(argv[i]) != 8)
			help();

	//Initialize with a random seed.
	srand(time(NULL));
	//srand(9);

	//Initializing vars.
	int licNum = generateLicenseNumber();
	unsigned char* licKey = generateKeyFromLicense(licNum);
	printLicKey(licKey);
	char* licIV = generateIV(licNum);
	char licName[50];
	sprintf(licName, "Alya-License.lic"	);

	//char licText[200] = "From:20170230 To:20170422";
	char licText[200];
	sprintf(licText, "DI%s:DF%s:DA%s:HM%s:HC00", argv[1], argv[2], argv[3], argv[4]);
	printf("licText: %s\n", licText);
	char licTextCipher[200];

	//Encrypt the license.
	int len = encrypt(	licText,
						strlen(licText),
						licKey,
  						licIV,
  						licTextCipher);

	
	printf("lenText: %d\n",len);

	//Write the license file
	FILE *licFile = fopen(licName, "wb");
	fwrite(&licNum, sizeof(int), 1, licFile);			//License number
	fwrite(&len, sizeof(int), 1, licFile);				//Lenght of the chiper data
	fwrite(licTextCipher, sizeof(char), len, licFile);	//Chiper data

	//Finalize
	fclose(licFile);
	//free(licKey);
	//free(licIV);

	printf("Generated Successful\n" );

	return 1;
}
