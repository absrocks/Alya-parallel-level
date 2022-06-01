<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="html"/>  
  <xsl:template match="/testSuite">   
    <html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">   
      <head>        
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"></meta>     
        <title>ALYA TESTSUITE</title>
        <meta name="description" content="alya testSuite results"></meta>
        <link href="./libs/application.css" media="all" rel="stylesheet" type="text/css"></link>          
        <script src="./libs/prototype.js" type="text/javascript">
          <xsl:text> </xsl:text>
        </script>
        <script src="./libs/effects.js" type="text/javascript">
          <xsl:text> </xsl:text>
        </script>
        <script src="./libs/dragdrop.js" type="text/javascript">
          <xsl:text> </xsl:text>
        </script>
        <script src="./libs/controls.js" type="text/javascript">
          <xsl:text> </xsl:text>
        </script>
        <script src="./libs/application.js" type="text/javascript">
          <xsl:text> </xsl:text>
        </script>
          
        <!--[if IE 6]>
        <style type="text/css">
        * html body{ width: expression( document.documentElement.clientWidth < 900 ? '900px' : '100%' ); }
        body {behavior: url(/redmine/stylesheets/csshover.htc?1331454347);}
        </style>
        <![endif]-->         
    </head>      
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
                <h2>Summary</h2>                
                <div class="tabs">
                  <ul>
                    <li>Information</li>      
                  </ul>  
                </div>                
                <div class="tab-content" id="tab-content-versions" style="">
                  <table class="list versions" border="1">
                    <thead><tr>
                      <th>Platform</th>
                      <th>SVN version</th>
                      <th>Compiler</th>    
                      <th style="width:15%">Date</th>
                    </tr></thead>
                    <tbody>                      
                      <tr class="version odd open ">
                        <td class="name" align="center"><xsl:value-of select="information/platform"/></td>
                        <td class="description" align="center"><xsl:value-of select="information/svnRevision"/></td>
                        <td align="center"><xsl:value-of select="information/compiler"/></td>
                        <td class="date" align="center">      
                          <xsl:value-of select="information/date"/>
                        </td>
                      </tr>
                    </tbody>
                  </table>
                </div>    
                <!-- COMPILATION RESULTS-->
                <div class="tabs">
                  <ul>
                    <li>COMPILATION</li>      
                  </ul>  
                </div>                
                <div class="tab-content" id="table1" style="">
                  <table class="list versions" border="1px">
                    <thead><tr>
                      <th>Module</th>
                      <th>Result</th>    
                      <th style="width:15%"><xsl:text> </xsl:text></th>
                    </tr></thead>
                    <tbody>
                      <xsl:for-each select="testSuiteJobResult/moduleCompilation/moduleCompilationResult">
                        <tr class="version odd open ">
                          <td class="name" align="center"><xsl:value-of select="moduleName"/></td>
                          <td  align="center">
                            <xsl:choose>
                              <xsl:when test="passed = 'true'"><img src="./libs/check.png" /></xsl:when>
                              <xsl:otherwise><img src="./libs/x.png" /></xsl:otherwise>
                            </xsl:choose>
                          </td>
                          <td class="buttons"  align="center">
                            <xsl:if test="passed = 'false'"><a class="icon icon-edit">
                              <xsl:attribute name="href">
                                <xsl:value-of select="file" />
                              </xsl:attribute>
                              Open</a></xsl:if>                                  
                          </td>
                        </tr>
                      </xsl:for-each>                                            
                    </tbody>
                  </table>                                    
                </div>
                <div class="tab-content" id="table1" style="">
                  <table class="list versions" border="1px">
                    <thead><tr>
                      <th>Service</th>
                      <th>Result</th>    
                      <th style="width:15%"><xsl:text> </xsl:text></th>
                    </tr></thead>
                    <tbody>
                      <xsl:for-each select="testSuiteJobResult/serviceCompilation/serviceCompilationResult">
                        <tr class="version odd open ">
                          <td class="name" align="center"><xsl:value-of select="serviceName"/></td>
                          <td  align="center">
                            <xsl:choose>
                              <xsl:when test="passed = 'true'"><img src="./libs/check.png" /></xsl:when>
                              <xsl:otherwise><img src="./libs/x.png" /></xsl:otherwise>
                            </xsl:choose>
                          </td>
                          <td class="buttons"  align="center">
                            <xsl:if test="passed = 'false'"><a class="icon icon-edit">
                              <xsl:attribute name="href">
                                <xsl:value-of select="file" />
                              </xsl:attribute>
                              Open</a></xsl:if>                                  
                          </td>
                        </tr>
                      </xsl:for-each>                                            
                    </tbody>
                  </table>                                    
                </div>
                <!-- TESTS RESULTS-->
                <div class="tabs">
                  <ul>
                    <li>TESTS</li>      
                  </ul>  
                </div>                
                <div class="tab-content" id="table1" style="">
                  <table class="list versions" border="1px">
                    <thead><tr>
                      <th>Module</th>
                      <th>Description</th>
                      <th>Result</th>    
                      <th style="width:15%"><xsl:text> </xsl:text></th>
                    </tr></thead>
                    <tbody>
                      <xsl:for-each select="testSuiteResult/moduleResults/moduleResult">
                        <tr class="version odd open ">
                          <td class="name" align="center"><xsl:value-of select="moduleName"/></td>
                          <td class="description"  align="center"><xsl:value-of select="description"/></td>
                          <td  align="center">
                            <xsl:choose>
                              <xsl:when test="passed = 'true'"><img src="./libs/check.png" /></xsl:when>
                              <xsl:otherwise><img src="./libs/x.png" /></xsl:otherwise>
                            </xsl:choose>
                          </td>
                          <td class="buttons"  align="center">
                            <xsl:if test="passed = 'false'"><a class="icon icon-edit">
                              <xsl:attribute name="href">
                                <xsl:value-of select="moduleName" />/<xsl:value-of select="moduleName" />.html
                              </xsl:attribute>
                              Open</a></xsl:if>                                  
                          </td>
                        </tr>
                      </xsl:for-each>                                            
                    </tbody>
                  </table>                                    
                </div>
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
</xsl:stylesheet>