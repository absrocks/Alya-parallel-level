/*-------------------------------------------------------------------------
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from
the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not
   be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution.
-------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include "cgnslib.h"
#include "cgns_header.h"
#include "ADF.h"

char cgns_error_mess[200] = "no CGNS error reported";

void cgi_error(char *format, ...) {
    va_list arg;
    va_start(arg, format);
    vsprintf(cgns_error_mess,format, arg);
    va_end(arg);
}

void cgi_warning(char *format, ...) {
    va_list arg;
    fprintf(stdout,"*** Warning:");
    va_start(arg, format);
    vfprintf(stdout,format,arg);
    va_end(arg);
    fprintf(stdout," ***\n");
}

char const *cg_get_error() {
    return cgns_error_mess;
}

void cg_error_exit() {
    fprintf(stderr,"%s\n",cgns_error_mess);
    exit(1);
}

void cg_error_print() {
    fprintf(stderr,"%s\n",cgns_error_mess);
}

void adf_error(char *routine_name, int ier) {
    char adf_err[80];

    ADF_Error_Message(ier, adf_err);
    cgi_error("Error in routine '%s':\n '%s'",routine_name,adf_err);
}

