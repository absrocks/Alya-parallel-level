/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 ******************************************************************************/
/// \file

#include "xwin_scoped_error_handler.h"

#include <X11/Xutil.h>
#include <X11/Xproto.h>
#include <X11/Xlibint.h>
#undef min
#undef max

#include <cstdio>
#include <map>

#include "common/forwarding_logger.h"

namespace MI {
namespace GPU {


// A map of error ids to their occurence
std::map<int, int>* error_map = 0;


int glx_error_handler(Display *dp, XErrorEvent *evt)
{
    // This is heavily based on the source code of Xmu library's
    // XmuSimpleErrorHandler function. Original copyright notice follows:


    // Copyright 1988, 1998  The Open Group
    //
    // Permission to use, copy, modify, distribute, and sell this software and its
    // documentation for any purpose is hereby granted without fee, provided that
    // the above copyright notice appear in all copies and that both that
    // copyright notice and this permission notice appear in supporting
    // documentation.
    //
    // The above copyright notice and this permission notice shall be included in
    // all copies or substantial portions of the Software.
    //
    // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
    // OPEN GROUP BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
    // AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    // CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    //
    // Except as contained in this notice, the name of The Open Group shall not be
    // used in advertising or otherwise to promote the sale, use or other dealings
    // in this Software without prior written authorization from The Open Group.

    // ignore errors for XQueryTree, XGetWindowAttributes,
    // and XGetGeometry, print a message for everything else.
    switch (evt->request_code) {
        case X_QueryTree:
        case X_GetWindowAttributes:
            if (evt->error_code == BadWindow) return 0;
            break;
        case X_GetGeometry:
            if (evt->error_code == BadDrawable) return 0;
            break;
    }

    unsigned long error_id = (evt->error_code   << 2*8*sizeof(evt->error_code)) |
                             (evt->request_code <<   8*sizeof(evt->error_code)) |
                              evt->minor_code;

    // We only need the output buffer. This is necessary because there is no way to
    // override _XExtension::error_values with XESetPrintErrorValues, this is specific
    // to each extension so it is the only way to capture the full output.
    FILE* fp = fopen("/dev/null", "w");
    static char stream_buffer[20*BUFSIZ];
    setbuf(fp, stream_buffer);

    static char buffer[BUFSIZ];
    static char mesg[BUFSIZ];
    static char number[32];
    static const char *mtype = "XlibMessage";
    _XExtension *ext  = NULL;
    _XExtension *bext = NULL;
    XGetErrorText(dp, evt->error_code, buffer, BUFSIZ);
    XGetErrorDatabaseText(dp, mtype, "XError", "X Error", mesg, BUFSIZ);
    fputs("\n", fp);
    fprintf(fp, "%s:  %s\n", mesg, buffer);
    XGetErrorDatabaseText(dp, mtype, "MajorCode", "Request Major code %d", mesg, BUFSIZ);
    fprintf(fp, mesg, evt->request_code);
    if (evt->request_code < 128) {
        snprintf(number, sizeof(number), "%d", evt->request_code);
        XGetErrorDatabaseText(dp, "XRequest", number, "", buffer, BUFSIZ);
    }
    else {
        // XXX this is non-portable
        for (ext = dp->ext_procs; ext && (ext->codes.major_opcode != evt->request_code); ) {
            ext = ext->next;
        }
        if (ext) {
            snprintf(buffer, sizeof(buffer), "%s", ext->name);
        }
        else {
            buffer[0] = '\0';
        }
    }
    fprintf(fp, " (%s)", buffer);
    fputs("\n", fp);
    if (evt->request_code >= 128) {
        XGetErrorDatabaseText(dp, mtype, "MinorCode", "Request Minor code %d",
            mesg, BUFSIZ);
        fprintf(fp, mesg, evt->minor_code);
        if (ext) {
            snprintf(mesg, sizeof(mesg),
                "%s.%d", ext->name, evt->minor_code);
            XGetErrorDatabaseText(dp, "XRequest", mesg, "", buffer, BUFSIZ);
            fprintf(fp, " (%s)", buffer);
        }
        fputs("\n", fp);
    }
    if (evt->error_code >= 128) {
        // kludge, try to find the extension that caused it
        buffer[0] = '\0';
        for (ext = dp->ext_procs; ext; ext = ext->next) {
            if (ext->error_string) 
                (*ext->error_string)(dp, evt->error_code, &ext->codes, buffer, BUFSIZ);
            if (buffer[0]) {
                bext = ext;
                break;
            }
            if (ext->codes.first_error &&
                ext->codes.first_error < evt->error_code &&
                (!bext || ext->codes.first_error > bext->codes.first_error))
            {
                bext = ext;
            }
        }
        if (bext) {
            snprintf(buffer, sizeof(buffer), "%s.%d",
                bext->name, evt->error_code - bext->codes.first_error);
        }
        else {
            strcpy(buffer, "Value");
        }
        XGetErrorDatabaseText(dp, mtype, buffer, "", mesg, BUFSIZ);
        if (mesg[0]) {
            fputs("  ", fp);
            fprintf(fp, mesg, evt->resourceid);
            fputs("\n", fp);
        }
        // let extensions try to print the values
        for (ext = dp->ext_procs; ext; ext = ext->next) {
            if (ext->error_values) {
                (*ext->error_values)(dp, evt, fp);
            }
        }
    }
    else {
        if ((evt->error_code == BadWindow  ) ||
            (evt->error_code == BadPixmap  ) ||
            (evt->error_code == BadCursor  ) ||
            (evt->error_code == BadFont    ) ||
            (evt->error_code == BadDrawable) ||
            (evt->error_code == BadColor   ) ||
            (evt->error_code == BadGC      ) ||
            (evt->error_code == BadIDChoice) ||
            (evt->error_code == BadValue   ) ||
            (evt->error_code == BadAtom    ))
        {
            if (evt->error_code == BadValue) {
                XGetErrorDatabaseText(dp, mtype, "Value", "Value 0x%x", mesg, BUFSIZ);
            }
            else {
                if (evt->error_code == BadAtom) {
                    XGetErrorDatabaseText(dp, mtype, "AtomID", "AtomID 0x%x", mesg, BUFSIZ);
                }
                else {
                    XGetErrorDatabaseText(dp, mtype, "ResourceID", "ResourceID 0x%x", mesg, BUFSIZ);
                }
            }
            fprintf(fp, mesg, evt->resourceid);
            fputs("\n", fp);
        }
    }

    XGetErrorDatabaseText(dp, mtype, "ErrorSerial", "Error Serial #%d", mesg, BUFSIZ);
    fprintf(fp, mesg, evt->serial);
    fputs("\n", fp);
    XGetErrorDatabaseText(dp, mtype, "CurrentSerial", "Current Serial #%d",mesg, BUFSIZ);
    fprintf(fp, mesg, NextRequest(dp)-1);
    fputs("\n", fp);

    // close the dummy file
    fclose(fp);

    if ((*error_map)[error_id]++ < 3) {
        ERROR_LOG << stream_buffer;
    }

    if ((*error_map)[error_id] == 3) {
        ERROR_LOG << "Logging of subsequent errors of this kind will be suppressed.";
    }

    if (evt->error_code == BadImplementation)
        return 0;
    return 1;
}


Scoped_error_handler::Scoped_error_handler()
{
    if (!error_map)
    {
        error_map = new std::map<int, int>;
    }
    m_old_error_handler = XSetErrorHandler(glx_error_handler);
}


Scoped_error_handler::~Scoped_error_handler()
{
    if (error_map) 
    {
        delete error_map;
        error_map = 0;
    }
    //XSetErrorHandler(m_old_error_handler);
}

}} // GPU, MI
