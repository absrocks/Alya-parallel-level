<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="html"/>  
  <xsl:template match="/testResult">   
    <html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">   
      <head>        
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"></meta>     
        <title>ALYA TESTSUITE</title>
        <meta name="description" content="alya testSuite results"></meta>
        <link href="../../libs/application.css" media="all" rel="stylesheet" type="text/css"></link>          
        <script src="../../libs/prototype.js" type="text/javascript">
          <xsl:text> </xsl:text>
        </script>
        <script src="../../libs/effects.js" type="text/javascript">
          <xsl:text> </xsl:text>
        </script>
        <script src="../../libs/dragdrop.js" type="text/javascript">
          <xsl:text> </xsl:text>
        </script>
        <script src="../../libs/controls.js" type="text/javascript">
          <xsl:text> </xsl:text>
        </script>
        <script src="../../libs/application.js" type="text/javascript">
          <xsl:text> </xsl:text>
        </script>
          
        <!--[if IE 6]>
        <style type="text/css">
        * html body{ width: expression( document.documentElement.clientWidth < 900 ? '900px' : '100%' ); }
        body {behavior: url(/redmine/stylesheets/csshover.htc?1331454347);}
        </style>
        <![endif]-->         
    </head>    
      <br></br>
      <body class="controller-projects action-settings">
        <div id="wrapper">
          <div id="wrapper2">            
            <div id="header">              
              <h1>ALYA TestSuite</h1>              
            </div>            
            <div class="nosidebar" id="main">
              <div id="sidebar">
                <xsl:text> </xsl:text>
              </div>              
              <div id="content">
                <h2><xsl:value-of select="module"/> Module <xsl:value-of select="name"/> Test</h2>                
                <xsl:for-each select="results/result">
                  <div class="tabs">
                    <ul>
                      <li><xsl:value-of select="processors"/> processors</li>      
                    </ul>  
                  </div>                               
                  <div class="tab-content" id="table1" style="">
                    <table class="list versions" border="1px">
                      <thead><tr>
                        <th>File</th>
                        <th>Differences</th>
                      </tr></thead>
                      <tbody>
                        <xsl:for-each select="file">
                          <xsl:if test="passed = 'false'">
                            <xsl:choose>
                              <xsl:when test="error != ''">
                                <tr class="version odd open ">
                                  <td class="name" align="center"><xsl:value-of select="fileName"/></td>
                                  <td class="description" align="center" colspan="2">ERROR: <xsl:value-of select="error"/></td>
                                </tr>
                              </xsl:when>
                              <xsl:when test="diff != ''">                                
                                <tr class="version odd open ">
                                  <td class="name" align="center"><xsl:value-of select="fileName"/></td>
                                  <td class="diff" align="center" colspan="2"><!--code><xsl:call-template name="insertBreaks">
                                    <xsl:with-param name="pText" select="diff"/>
                                  </xsl:call-template></code--><xsl:value-of select="diff"/></td>
                                </tr>
                              </xsl:when>
                              <xsl:when test="section != ''">
                                <tr>
                                  <td class="name" align="center"><xsl:value-of select="fileName"/></td>
                                  <td class="description" align="center">
                                    <div class="tab-content" id="table11" style="">
                                      <table class="list versions" border="1px">
                                      <thead>
                                        <tr>
                                          <th colspan="3"><xsl:value-of select="section"/></th>
                                        </tr>
                                      </thead>
                                      <tbody>
                                        <tr class="version odd open ">
                                          <td class="name" align="center">Column</td>
                                          <td class="description" align="center">Original</td>
                                          <td class="description" align="center">Obtained</td>
                                        </tr>
                                        <xsl:for-each select="diffValue">
                                          <tr class="version odd open ">
                                            <td class="name" align="center"><xsl:value-of select="column"/></td>
                                            <td class="description" align="center"><xsl:value-of select="original"/></td>
                                            <td class="description" align="center"><xsl:value-of select="obtained"/></td>
                                          </tr>                                          
                                        </xsl:for-each>
                                      </tbody>
                                    </table>
                                    </div>  
                                  </td>
                                </tr>
                              </xsl:when>
                              <xsl:otherwise>
                                <tr class="version odd open ">
                                  <td class="name" align="center"><xsl:value-of select="fileName"/></td>
                                  <td class="description" align="center"><xsl:value-of select="original"/></td>
                                  <td class="description" align="center"><xsl:value-of select="obtained"/></td>
                                </tr>
                              </xsl:otherwise>
                            </xsl:choose>
                          </xsl:if>                                                   
                        </xsl:for-each>                                            
                      </tbody>
                    </table>                                    
                  </div>                  
                </xsl:for-each>
              </div>
            </div>            
            <div id="footer">
              <div class="bgl"><div class="bgr">
                <a href="http://www.bsc.es/">BSC-CNS</a>
              </div></div>
            </div>
          </div>
        </div>        
      </body>
    </html>    
  </xsl:template>
  <xsl:template name="string-replace-all">
    <xsl:param name="text" />
    <xsl:param name="replace" />
    <xsl:param name="by" />
    <xsl:choose>
      <xsl:when test="contains($text, '&#xA;')">
        <xsl:value-of select="substring-before($text,'&#xA;')" />
        <br/>
        <xsl:call-template name="string-replace-all">
          <xsl:with-param name="text"
            select="substring-after($text,'&#xA;')" />
          <xsl:with-param name="replace" select="'&#xA;'" />
          <xsl:with-param name="by" select="$by" />
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$text" />
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>
  
  <xsl:template name="insertBreaks">
    <xsl:param name="pText" select="."/>    
    <xsl:choose>
      <xsl:when test="not(contains($pText, '&#xA;'))">
        <xsl:copy-of select="$pText"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="substring-before($pText, '&#xA;')"/>
        <br> </br>
        <xsl:call-template name="insertBreaks">
          <xsl:with-param name="pText" select=
            "substring-after($pText, '&#xA;')"/>
        </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>
</xsl:stylesheet>