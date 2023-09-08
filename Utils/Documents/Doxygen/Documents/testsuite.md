# TestSuite {#testsuite}

## Basic information

The **TestSuite** is the tool that we use for testing, This tool compile Alya many times with different configurations and modules in order to check the indepence of the different Alya Modules, the external libraries such as openmp, ompss, libple and agmg and some special compilation flags like integer 8 or NO\_MPI.
Then for checking the results we do regresions tests, every regresion test is an @ref alyaCase and we check the numbers found in the outputfiles and compare it with some previous verified base files. That's called regresion testing

**Regression testing** is a type of software testing which verifies that software which was previously developed and tested still performs the same way after it was changed or interfaced with other software. Changes may include software enhancements, patches, configuration changes, etc. During regression testing, new software bugs or regressions may be uncovered. Alya TestSuite impact analysis is performed to determine what areas could be affected by the proposed changes. These areas may include functional and non-functional areas of the system.

<span style="color:orange; font-size:18px;"> The purpose of regression testing is to ensure that changes such as those mentioned above have not introduced new faults. One of the main reasons for regression testing is to determine whether a change in one part of the software affects other parts of the software.
</span>

### How to execute my TestSuite

The TestSuite 2.0 has been thought to be easy and simple for Alya programmers and also to be very usefull and configurable for your own tests. So let's start!

<ol>
	<li>First checkout the TestSuite 2.0:</li>
	<br/>From MN4:
	\code
		svn co file:///gpfs/projects/bsc21/svnroot/AlyaTS/Trunk
	\endcode
	Local Machine:
	\code
		svn co svn+ssh://bsc21***@mn3.bsc.es/gpfs/projects/bsc21/svnroot/AlyaTS/Trunk
	\endcode
	<li>Second configure the testSuite:</li>
	\code
		./configure
	\endcode
	<li>Finally execute the testSuite:</li>
	\code
		./testSuite
	\endcode
</ol>

### TestSuite 2.0 finished and now what?

**Where can I find the report?**
You will find a new directory named generalReport and inside an index.html file, open it with your favourite navigator.</p>

Once you have oppened the documentation you will see something like this:
	
<img src="../img/TSReport.png" width="70%" style="border: 10px solid; padding: 30px;">

Where:
<ol>
	<li>General information</li>
	<li>Details on each subexecution</li>
	<li>Name of the subexectuion</li>
	<li>Loaded modules on Marenostrum</li>
	<li>Compilation summary</li>
	<li>Test summary</li>
	<li>Compilation details</li>
	<li>Tests details</li>
	<li>Link to the used config.in for this subexecution</li>

</ol>

### Configure Options

By the moment the configure is quite simple, in the future the options will increase.

<ul>
 <li><b>-m</b> To select the modules that will be compiled.</li>
 <li><b>-e</b> To change the user email.</li>
 <li><b>-s</b> To select the path in marenostrum where the testsuite will be executed.</li>
 <li><b>-b</b> To select the branch to be tested.</li>
 <li><b>-t</b> To select the tests to be executed.</li>
 <li><b>-all</b> To do the general configure, is the only way to reconfigure your MN user or MN login.</li>
</ul>

### Launch Options

<ul>
 <li><b>./testSuite </b> -> Normal execution </li>
 <li><b>./testSuite launch </b>-> Execute the TS and finish after the initialization.</li>
 <li><b>./testSuite check </b>-> After "launch" check if the TS has finished</li>
 <li><b>./testSuite report </b>-> After launch check if the TS has finished and generate the report </li>
</ul>

You can do `./testSuite launch`, switch off your PC and then `./testSuite check/report`

### How to add a test in the TestSuite

http://tannat.bsc.es/TestSuite/addTestForm.php

## How does the TestSuite work internally?

The aim of this secction is to explain what does the TestSuite do to test ALYA.
Before starting you must know that the TestSuite 2.0 is able to launch mutiple SubTestsSuites with different configuratioons, for exemple one for intel compilers, another one for gfortran compilers and also for example a TestSuite for integer 8. **What we will explain in the next section is what each of this text do.**

### First: Checkout
The first thing that the TestSuite 2.0 does is to checkout the following:
<ul>
	<li>Alya</li>
	<li>TestSuite files needed to execute</li>
	<li>TestSuite tests to execute</li>
	<li>Also it copy all the necessary files to execute the job, configurations etc...</li>
</ul>

### Second: Compilation
<ul>
	<li>First the TestSuite compile the externat libraries such as metis or libple.</li>
	<li>The second step is to execute single module compilation that compile each module one by one, this test is made to detect compilation errors in each module and also to test the indepence of this module from the others.</li>
	<li>Then is time to do some extra compilations test called "special compilation tests" these are to compile alya with special flags such as MPI_OFF, openmp or integer 8. <b>Important: This compilation tests only compile the modules that has been passed the previous compilation.</b></li>
	<li>Finally we compile Alya without any special flags and all the modules that has been passed the single compilation.</li>
</ul>


### Third: Test Execution
This step only execute all the tests with the different configuration MPI/OMP, every excution is done inside the correspondig folder following this rule: MPI\_OMP. At the end we have all the diffente executions in a different folder so we can check later the result


### Fourth: Comparison
First we check if the execution finish correctfully, the rule to determine this is if the output print "ALYA  CALCULATIONS CORRECT" then we supose that everything worked well.
Then we compare the output files specified in the test configuration:


#### What do we exactly compare?

We extract from the output files all the numbers that match with the following regular expresion.xº

    -?\d\.\d\d\d\d\d\d\d\dE[+-]\d\d\d

Then we can specify which **rows** and **columns** do we want to compare in the configuration file for each test.

#### Method "albert"

**Note: the method described below is obsolete and must be replaced by the methods rounding, absolute or relative**

There are two kind of tolerancy in the TestSuite test configuration:
<ul>
	<li><b>Tolerancy:</b> Tolerancy means how tolerant is a number from errors, it use a integer value that means ignore x decimals beggining from the left</li>
	<p><b>For example:</b> for this given value: <b><span style="color: green;">0.11102<span style="color: red;">260</span>E+001</span></b>, with tolerancy 3 we will compare only the green values but not the red ones</p>
	<li><b>Precision:</b> This is thought to compare the FINITE ELEMENT ERRORS but also can be used for other comparision, in this case we also have an integer value but it is ralted to the exponent. It means how different can be the exponenent between the base value and de obtained</li>
	<p><b>For example: </b><b>0.11102260E+001 (base) , 0.11102260E+003 (test)</b> with a tolerancy value of 2 will pass the test.</p>
</ul>

#### Method "absolute"

This method enables an absolute comparison between two numbers.
It evaluates if the absolute difference between the test result and the reference is inferior than a threshold we call _tolerancy_.
Tolerancy is expressed as a float number.
This method is ideal when all the numbers have the same order of magnitude.

#### Method "relative"

This method enables a relative comparison between two numbers.
It evaluates if the absolute difference between the test result and the reference values divided by the reference value is inferior to a ratio _diff_.
Diff is expressed as a float number.
This method is ideal when the numbers have different orders of magnitude.

#### Method "rounding"

This method enables to specify the float precision, according to the value of _tolerancy_, expressed as an integer. Tolerancy corresponds to the number of digit after the comma that we take into account.
For instance, for a number 0.11102260E+001 with tolerancy 5, we will consider only the differences greater than 0.00001E+001.
This method is ideal whatever the order of magnitude amplitude.

### Fifth: Postprocess
<ul>
	<li>Copy files from Marenostrum.</li>
	<li>Generate the report.</li>
	<li>Send Emails.</li>
	<li>Insert to database (NOT IMPLEMNTED)</li>
</ul>

## Frequently asked questions (FAQ)
<ol>
	<li><b>Why my test fail allways even when there are no errors?</b></li><br/>
	Well this is probably because your test is not well thought, but also can be because the machine precision due to the floating point numbers it's not exact! In the second situation you should ajust the tolereancy of your test.<br/><br/>

	<li><b>Then obvious question, how should be my test?</b></li><br/>
	The objective of our Test Suite is to assure code stability: the code keeps producing the same results although it is constantly evolving. The TS is a differential tool, so its power relies on both an extensive code coverage and a high frequency of testing. Think of it as a computing a partial derivative in a discrete way! 
	Considering that a change in the code results is not necessarily wrong, the idea of the TS is to inform of any change. Once it happens, an IMMEDIATE decision must be taken: either accept it or reject it. In the first case, the TS must be modified to account for the new output. In the second case, the source of error can be easily tracked back following the code changes and eliminated.
	<br/>
	<br/>

	With all these in mind, the TS should include tests following these guidelines (ideally):

	<br/>

	- Tests must be short and small, targeting EXCLUSIVELY to find differences between versions.
	- By short we mean taking a few seconds, with no more than 2 or 3 time steps.
	- By small we mean anything that can be solved in a few seconds in one MN node (say with
	fewer than 1000 elements) using a few MPI tasks. 
	- Sequential tests are also fine, depending on what you want to check.
	- It is much better to use 10 different tests to cover 10 different potentially changing features than 1 test to cover all of them. So don't be neither shy nor modest: prepare as many tests as you want. 
	- They are not benchmarks (which is a different kind of test, not covered in the current TS version). Benchmarks are excellent, but they are a different issue.

	<br/>
	<li><b>Where can I find the TestSuite tests?</b></li><br/>
	\code
		svn co file:///gpfs/projects/bsc21/svnroot/AlyaTS/tests
	\endcode
	<p></p>

</ol>

