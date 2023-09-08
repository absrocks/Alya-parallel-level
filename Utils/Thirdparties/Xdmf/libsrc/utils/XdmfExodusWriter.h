/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfExodusWriter.h,v 1.2 2010-03-24 20:03:48 kwleiter Exp $  */
/*  Date : $Date: 2010-03-24 20:03:48 $ */
/*  Version : $Revision: 1.2 $ */
/*                                                                 */
/*  Author:                                                        */
/*     Kenneth Leiter                                              */
/*     kenneth.leiter@arl.army.mil                                 */
/*     US Army Research Laboratory                                 */
/*     Aberdeen Proving Ground, MD                                 */
/*                                                                 */
/*     Copyright @ 2009 US Army Research Laboratory                */
/*     All Rights Reserved                                         */
/*     See Copyright.txt or http://www.arl.hpc.mil/ice for details */
/*                                                                 */
/*     This software is distributed WITHOUT ANY WARRANTY; without  */
/*     even the implied warranty of MERCHANTABILITY or FITNESS     */
/*     FOR A PARTICULAR PURPOSE.  See the above copyright notice   */
/*     for more information.                                       */
/*                                                                 */
/*******************************************************************/

#include <Xdmf.h>

class XdmfExodusWriterNameHandler;

/*!
 * @brief XdmfExodusWriter class encapsulates the operation of writing an
 *        ExodusII file from an Xdmf file.
 */

class XdmfExodusWriter
{
  public:
    /*!
     * Constructor.
     */
    XdmfExodusWriter();

    /*!
     * Destructor.
     */
    ~XdmfExodusWriter();

    /*!
     * Write XdmfGrid to exodus file
     *
     * @param fileName a const char * containing the path of the exodus ii file to write
     * @param gridToWrite Pointer to XdmfGrid to write.
     *
     */
    void write(const char * fileName, XdmfGrid * gridToWrite);

  private:
    /*!
     * Convert from xdmf to exodus ii cell types
     *
     * @param XdmfInt32 containing the xdmf topology to convert
     *
     * @return a std::string of the exodus ii topology type equivalent to the xdmf topology type.
     */
    std::string DetermineExodusCellType(XdmfInt32 xdmfElementType);
    XdmfExodusWriterNameHandler * nameHandler;
};
