<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text"/>  
  <xsl:template match="/log">
    --------ALYA SVN COMMIT INFO:
    ----Message:"<xsl:value-of select="logentry/msg"/>"
    ----Revision:"<xsl:value-of select="logentry/@revision"/>"
    ----Author:"<xsl:value-of select="logentry/author"/>"
    ----Date:"<xsl:value-of select="logentry/date"/>"
    ----Changes:
    <xsl:for-each select="logentry/paths/path">
    --File:<xsl:value-of select="."/> Action:<xsl:value-of select="@action"/>  
    </xsl:for-each>
  </xsl:template>  
</xsl:stylesheet>