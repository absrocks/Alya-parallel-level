#include "chiper.h"

int generateLicenseNumber(){
  int r = rand();
  return r;
}

unsigned char* generateKeyFromLicense(int n){

  int v[10];
  int auxn = n;
  int nmod10;
  unsigned char* randkey;
  randkey = (unsigned char*)malloc(sizeof(unsigned char*)*32);

  nmod10 = n%10;
  for(int i=0; i<10; i++){
    v[(nmod10+i)%10] = auxn%10 + 1;
    auxn /=10;
  }

  srand(v[0]);

  for (int i = 0; i < 32; ++i)
    randkey[i] = (rand()*v[i%10])%255;

  for (int i = 0; i < 32; ++i)
    randkey[i] = (randkey[i]*v[3])%255;

  for (int i = 0; i < 32; ++i)
    randkey[i] = (randkey[i]/v[5])%255;

  for (int i = 0; i < 32; ++i)
    randkey[i] = (randkey[i]*randkey[i]*v[6])%255;

  return randkey;

}

char* generateIV(int licNum){
  int v[10];
  int auxn = licNum;
  int nmod10;
  unsigned char* randIV;
  randIV = (unsigned char*)malloc(sizeof(unsigned char*)*32);

  nmod10 = (licNum/10)%10;
  for(int i=0; i<10; i++){
    v[(nmod10+i)%10] = auxn%10 + 1;
    auxn /=10;
  }

  srand(v[3]*v[2]);

  for (int i = 0; i < 32; ++i)
    randIV[i] = (rand()*v[i%10])%255;

  return randIV;
}

void handleErrors(void)
{
  ERR_print_errors_fp(stderr);
  abort();
}

int decrypt(unsigned char *ciphertext, int ciphertext_len, unsigned char *key,
  unsigned char *iv, unsigned char *plaintext)
{
  EVP_CIPHER_CTX *ctx;

  int len;

  int plaintext_len;

  /* Create and initialise the context */
  if(!(ctx = EVP_CIPHER_CTX_new())) handleErrors();

  /* Initialise the decryption operation. IMPORTANT - ensure you use a key
   * and IV size appropriate for your cipher
   * In this example we are using 256 bit AES (i.e. a 256 bit key). The
   * IV size for *most* modes is the same as the block size. For AES this
   * is 128 bits */
  if(1 != EVP_DecryptInit_ex(ctx, EVP_aes_256_cbc(), NULL, key, iv))
    handleErrors();

  /* Provide the message to be decrypted, and obtain the plaintext output.
   * EVP_DecryptUpdate can be called multiple times if necessary
   */
  if(1 != EVP_DecryptUpdate(ctx, plaintext, &len, ciphertext, ciphertext_len))
    handleErrors();
  plaintext_len = len;

  /* Finalise the decryption. Further plaintext bytes may be written at
   * this stage.
   */
  if(1 != EVP_DecryptFinal_ex(ctx, plaintext + len, &len)) handleErrors();
  plaintext_len += len;

  /* Clean up */
  EVP_CIPHER_CTX_free(ctx);

  return plaintext_len;
}

int encrypt(unsigned char *plaintext, int plaintext_len, unsigned char *key,
  unsigned char *iv, unsigned char *ciphertext)
{
  EVP_CIPHER_CTX *ctx;

  int len;

  int ciphertext_len;

  /* Create and initialise the context */
  if(!(ctx = EVP_CIPHER_CTX_new())) handleErrors();

  /* Initialise the encryption operation. IMPORTANT - ensure you use a key
   * and IV size appropriate for your cipher
   * In this example we are using 256 bit AES (i.e. a 256 bit key). The
   * IV size for *most* modes is the same as the block size. For AES this
   * is 128 bits */
  if(1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_cbc(), NULL, key, iv))
    handleErrors();

  /* Provide the message to be encrypted, and obtain the encrypted output.
   * EVP_EncryptUpdate can be called multiple times if necessary
   */
  if(1 != EVP_EncryptUpdate(ctx, ciphertext, &len, plaintext, plaintext_len))
    handleErrors();
  ciphertext_len = len;

  /* Finalise the encryption. Further ciphertext bytes may be written at
   * this stage.
   */
  if(1 != EVP_EncryptFinal_ex(ctx, ciphertext + len, &len)) handleErrors();
  ciphertext_len += len;

  /* Clean up */
  EVP_CIPHER_CTX_free(ctx);

  return ciphertext_len;
}
