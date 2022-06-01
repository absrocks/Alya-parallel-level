<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text"/>  
  <xsl:template match="/log">
    <xsl:for-each select="logentry">
        ----LOGENTRY
        ----REVISION:<xsl:value-of select="@revision"/>
        ----AUTHOR:<xsl:value-of select="author"/>
        ----CHANGES:<xsl:for-each select="paths/path"><xsl:value-of select="."/>#</xsl:for-each>
        ----MESSAGE:<xsl:value-of select="msg"/>
    </xsl:for-each>
  </xsl:template>  
</xsl:stylesheet>