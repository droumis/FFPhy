<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" dir="ltr" lang="en"><head>


  <meta http-equiv="Content-Type" content="text/html; charset=utf-8"><title>data_analysis [Frank Lab Wiki]</title>
  

  <meta name="generator" content="DokuWiki Release 2006-11-06">
<meta name="robots" content="noindex,nofollow">
<meta name="date" content="2008-08-01T09:32:14-0700">
<meta name="keywords" content="data_analysis">
<link rel="start" href="http://hippo/wiki/">
<link rel="contents" href="http://hippo/wiki/doku.php?id=data_analysis&amp;do=index" title="Index">
<link rel="alternate" type="application/rss+xml" title="Recent Changes" href="http://hippo/wiki/feed.php">
<link rel="alternate" type="application/rss+xml" title="Current Namespace" href="http://hippo/wiki/feed.php?mode=list&amp;ns=">
<link rel="alternate" type="text/html" title="Plain HTML" href="http://hippo/wiki/doku.php?do=export_xhtml&amp;id=data_analysis">
<link rel="alternate" type="text/plain" title="Wiki Markup" href="http://hippo/wiki/doku.php?do=export_raw&amp;id=data_analysis">
<link rel="stylesheet" media="screen" type="text/css" href="wikipreprocessing_files/css_002.php">
<link rel="stylesheet" media="print" type="text/css" href="wikipreprocessing_files/css.php">
<script type="text/javascript" charset="utf-8" src="wikipreprocessing_files/js.php"></script>

  <link rel="shortcut icon" href="http://hippo/wiki/lib/tpl/default/images/favicon.ico"></head><body>
<div class="dokuwiki">
  
  <div class="stylehead">

    <div class="header">
      <div class="pagename">
        [[<a href="http://hippo/wiki/doku.php?id=data_analysis&amp;do=backlink">data_analysis</a>]]
      </div>
      <div class="logo">
        <a href="http://hippo/wiki/doku.php?id=" name="dokuwiki__top" id="dokuwiki__top" accesskey="h" title="[ALT+H]">Frank Lab Wiki</a>      </div>

      <div class="clearer"></div>
    </div>

    
    <div class="bar" id="bar__top">
      <div class="bar-left" id="bar__topleft">
        <form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="rev" value="" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit this page" class="button" accesskey="e" title="Edit this page [ALT+E]" type="submit"></div></form>        <form class="button" method="get" action="/wiki/doku.php"><div class="no"><input name="do" value="revisions" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Old revisions" class="button" accesskey="o" title="Old revisions [ALT+O]" type="submit"></div></form>      </div>

      <div class="bar-right" id="bar__topright">
        <form class="button" method="get" action="/wiki/doku.php"><div class="no"><input name="do" value="recent" type="hidden"><input name="id" value="" type="hidden"><input value="Recent changes" class="button" accesskey="r" title="Recent changes [ALT+R]" type="submit"></div></form>        <form action="/wiki/doku.php?id=" accept-charset="utf-8" class="search" id="dw__search"><div class="no"><input name="do" value="search" type="hidden"><input id="qsearch__in" accesskey="f" name="id" class="edit" title="[ALT+F]" type="text"><input value="Search" class="button" title="Search" type="submit"><div id="qsearch__out" class="ajax_qsearch JSpopup"></div></div></form>&nbsp;
      </div>

      <div class="clearer"></div>
    </div>

        <div class="breadcrumbs">
      Trace: <span class="bcsep">»</span> <a href="http://hippo/wiki/doku.php?id=data_filtering" class="breadcrumbs" title="data_filtering">data_filtering</a> <span class="bcsep">»</span> <a href="http://hippo/wiki/doku.php?id=start" class="breadcrumbs" title="start">start</a> <span class="bcsep">»</span> <span class="curid"><a href="http://hippo/wiki/doku.php?id=data_analysis" class="breadcrumbs" title="data_analysis">data_analysis</a></span>          </div>
    
    
  </div>
  
  
  <div class="page">
    <!-- wikipage start -->
    <div class="toc">
<div class="tocheader toctoggle" id="toc__header"><img alt="-" src="wikipreprocessing_files/arrow_up.gif" id="toc__hide"><img style="display: none;" alt="+" src="wikipreprocessing_files/arrow_down.gif" id="toc__show">Table of Contents</div>
<div id="toc__inside">

<ul class="toc">
<li class="level1"><div class="li"><span class="li"><a href="#frank_lab_analysis_instructions" class="toc">FRANK LAB ANALYSIS INSTRUCTIONS</a></span></div>
<ul class="toc">
<li class="level2"><div class="li"><span class="li"><a href="#moving_and_extracting_data_from_the_rig" class="toc">1. MOVING AND EXTRACTING DATA FROM THE RIG</a></span></div>
<ul class="toc">
<li class="level3"><div class="li"><span class="li"><a href="#a._transfer_the_data" class="toc">A. Transfer the data</a></span></div></li>
<li class="level3"><div class="li"><span class="li"><a href="#b._extract_the_data" class="toc">B.  Extract the data</a></span></div></li>
</ul>
</li>
<li class="level2"><div class="li"><span class="li"><a href="#reconstructing_the_animal_s_position" class="toc">2. RECONSTRUCTING THE ANIMAL'S POSITION</a></span></div></li>
<li class="level2"><div class="li"><span class="li"><a href="#clustering_spike_data_using_matclust" class="toc">3. CLUSTERING SPIKE DATA USING MATCLUST</a></span></div>
<ul class="toc">
<li class="level3"><div class="li"><span class="li"><a href="#a._setting_up_matclust_on_a_network_computer" class="toc">A. Setting up matclust on a network computer</a></span></div></li>
<li class="level3"><div class="li"><span class="li"><a href="#b._making_the_parameter_files" class="toc">B. Making the parameter files</a></span></div></li>
</ul>
</li>
<li class="level2"><div class="li"><span class="li"><a href="#processing_the_matclust_files" class="toc">4. PROCESSING THE MATCLUST FILES</a></span></div></li>
<li class="level2"><div class="li"><span class="li"><a href="#linearizing_the_data" class="toc">5. LINEARIZING THE DATA</a></span></div></li>
<li class="level2"><div class="li"><span class="li"><a href="#filtering_the_eeg_data" class="toc">6. FILTERING THE EEG DATA</a></span></div></li>
<li class="level2"><div class="li"><span class="li"><a href="#extracting_ripples_and_high_amplitude_theta_times" class="toc">7. EXTRACTING RIPPLES AND HIGH AMPLITUDE THETA TIMES</a></span></div></li></ul>
</li></ul>
</div>
</div>
<h1><a name="frank_lab_analysis_instructions" id="frank_lab_analysis_instructions">FRANK LAB ANALYSIS INSTRUCTIONS</a></h1>
<div class="level1">

<p>
Instructions on how to perform the initial data processing after data has been collected.  
</p>

</div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="1-143" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="FRANK LAB ANALYSIS INSTRUCTIONS" type="submit"></div></form></div>
<h2><a name="moving_and_extracting_data_from_the_rig" id="moving_and_extracting_data_from_the_rig">1. MOVING AND EXTRACTING DATA FROM THE RIG</a></h2>
<div class="level2">

</div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="144-200" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="1. MOVING AND EXTRACTING DATA FROM THE RIG" type="submit"></div></form></div>
<h3><a name="a._transfer_the_data" id="a._transfer_the_data">A. Transfer the data</a></h3>
<div class="level3">

<p>
On your computer, go to the local data drive (use the ‘df’ command to
list the available directories and choose the data directory whose
entry starts with “dev”). On that hardrive, create a directory for the
animal (e.g. (/datadrive/username/charley) Within this directory create
a subdirectory for each day of recording (e.g.
/datadrive/username/charley/charley01,
/datadrive/username/charley/charley02, etc.). sftp the data from the
data acquisition machines to the appropriate subdirectory
</p>

</div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="201-746" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="A. Transfer the data" type="submit"></div></form></div>
<h3><a name="b._extract_the_data" id="b._extract_the_data">B.  Extract the data</a></h3>
<div class="level3">

<p>
cd to one of the day directories for the animal (e.g. /datadrive/username/charley/charley01). 
</p>
<pre class="code"> &gt;&gt; nspike_extract -all charley*******.dat</pre>

<p>
 This will run the extraction program on one of the .dat files. (Here, <strong></strong>*
means whatever you called your .dat file). If you’re using a 64-bit
machine make sure that /usr/local/bin64 is included in your path and
that it comes before /usr/local/bin. Run this program for each .dat
file in each day’s directory. </p>

<p>
*note there is a 2G limit on the files that can be created using
nspike_extract. if you get a “file size limit exceded”, extract the
offending file separately in two parts. For instance use nspike_extract
-pos file.ike.dat in one directory, and nspike_extract -pos
file.ike.dat.2 in another directory, then reconstruct the position and
put the .p files in the same directory.
</p>

<p>
Run this for each day.  This will create a number of files from each .dat file:  </p><ol><li class="level1"><p><strong>a posimage file</strong>
</p>

<p>
containing the position frames (eg. *.mpeg) </p></li><li class="level1"><p><strong>Tetrode subdirectories and tt files</strong>: 
</p>

<p>Data from each tetrode will be in it’s own directory (e.g. 01-110)
where 01 is the number of the tetrode and 110 is the depth. Within that
directory the data are stored in a ttfile (e.g. 01-110.tt) containing
timestamps and spike waveforms. </p></li><li class="level1"><p><strong>EEG files</strong>: 
</p>

<p>
Each eeg channel will have a file (eg. 01-110.eeg for tetrode 1) </p></li><li class="level1"><p><strong>.config files</strong>: 
</p>

<p>
These files contain the configuration information for the software. </p></li><li class="level1"><p><strong>.event files</strong>: 
</p>

<p>
This file contains the digital input and output events from the
recording. Note that where there are multiple data files from the same
data acquisition computer, the software will automatically append the
data onto the same output file (as it should). If you need to reextract
data by hand you should delete the output files created by the
extraction and then reextract the data. If, for some reason you reset
the clock between data files and you wish to append one file with a
time offset onto another you can use the -toffset argument (here to add
one hour and thirty minutes to each timestamp): </p>
<pre class="code">   &gt;&gt; nspike_extract -all -toffset 1:30:00 datafilename</pre></li></ol></div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="747-2920" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="B.  Extract the data" type="submit"></div></form></div>
<h2><a name="reconstructing_the_animal_s_position" id="reconstructing_the_animal_s_position">2. RECONSTRUCTING THE ANIMAL'S POSITION</a></h2>
<div class="level2">

<p> Go to the first day’s data (e.g.
/datadrive/username/charley/charley01) There will be a file called
“charley*****.mpeg“(*** indicating the date information). For each
epoch of data collection, run nspike_fixpos as follows: </p><ol><li class="level1"><p>Sleep epoch (usually epochs 1, 3 and 5, below is for epoch 1)
</p>
<pre class="code"> &gt;&gt; nspike_fixpos -p charley****.mpeg -t  charley****.postimestamp -o charley01_1.p -tstart epochstarttime -tend epochendtime -skip 10 -f charley***.mpegoffset</pre></li><li class="level1"><p>Run epoch 
</p>
<pre class="code"> For the run epochs, the -skip option should be set to 3</pre></li></ol><p> There are four arguments that change for each epoch:</p><ul><li class="level1"><p><code>-o filename</code> – where filename is day (charley01) %”_”% epoch number (1, 2, etc)</p></li><li class="level1"><p><code>-tstart epochstarttime</code>  – where epochstarttime is the starting time you wrote down in your notes for the epoch</p></li><li class="level1"><p><code>-tend epochendtime</code>  –  where epochendtime is the end time you wrote down in your notes for the epoch</p></li><li class="level1"><p><code>-skip #</code>
– where # is the number of frames that the program is allowed to skip
before stopping for user input when it can’t figure out where the
diodes are. This is usually 3 for run sessions and 10 for sleep
sessions.</p></li></ul><p> You will have started the data collection
before putting the animal on the track. As such, for run sessions it is
very important that you mark down the time the animal actually started
running on the track. </p>

<p>
Once you have started spike_fixpos, set the size of the exclude ring
for the diode positions (use a value of about 10 for runs and about 15
for sleeps)
</p>

<p>Using the left and right mouse buttons, enter the position of the
front and back diode for frames where the playback stops. Press return
to tell the computer to accept the positions. </p>

<p>
The following is additional spike_fixpos information from the top of spike_fixpos.c:
</p>

<p>
<strong>  Defining Diode Positions</strong>
</p>

<p>
When the position markers are hollow, the user can define where the
front and back diode positions are. The left mouse button is for the
front diode, and the right button is for the back diode. When return is
pressed, the position is saved and playback is resumed.
</p>

<p>
<strong>  Main Menu</strong>
</p>

<p>
The main menu pops up when the middle mouse button is pressed. Keep the button pressed down and release to make a selection.
</p>

<p>
<strong>  Fast Keys </strong>
</p>

<p>
You can press keys at any time to make quick changes. These keys are: 
</p>
<pre class="code">spacebar - pause playback
enter - execute command
' - revert to the last user-picked diode positions
t - change threshold
r - change exclude ring size
b - step back one frame
n - step foreward one frame
c - clear excluded pixels
x - add pixel to exclude list
S - foreward to next full image
B - back to first frame
Q - quit</pre>

<p>
 <strong>  Command Line </strong>
</p>

<p>
Type <code>spike_fixpos -a</code> for command line options.
</p>

</div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="2921-5696" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="2. RECONSTRUCTING THE ANIMAL'S POSITION" type="submit"></div></form></div>
<h2><a name="clustering_spike_data_using_matclust" id="clustering_spike_data_using_matclust">3. CLUSTERING SPIKE DATA USING MATCLUST</a></h2>
<div class="level2">

</div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="5697-5752" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="3. CLUSTERING SPIKE DATA USING MATCLUST" type="submit"></div></form></div>
<h3><a name="a._setting_up_matclust_on_a_network_computer" id="a._setting_up_matclust_on_a_network_computer">A. Setting up matclust on a network computer</a></h3>
<div class="level3">

<p>
  (each user needs to do this the first time he/she uses matclust)
</p>

<p>If you have not analyzed data in matlab yet, then you’ll need to set
up it’s path so all the needed programs can be found. This is easily
done by creating ~/matlab/startup.m with the following text: </p>
<pre class="code">p=path; path(userpath,p);
p=path; path('~/Src/Matlab/',p);
p=path; path('/usr/local/matclust',p);
p=path; path('/usr/local/matclust/datacompile',p);
p=path; path('/usr/local/linearization',p);
p=path; path('/usr/local/mpgread',p);</pre>

<p>
 You may need other directories in there later, but it’s a good start for now.
</p>

<p>
All of the required matclust files are located at /usr/local/matclust
on hippo. /usr/local should be mounted on your computer already. Before
you set up matclust on a computer, you need the matlab path to include
/usr/local/matclust (above). </p>

<p>
Start matlab. At the prompt, run <code>setuplocal</code>. This
should bring up a window asking for two directories to create different
files. The first box asks where the matclust folder should be created.
You should put the folder somewhere on your home directory, like
/home/yourname. Here, a folder named ‘matclust’ will be created that
contains your personal startup defaults and your personal helper
programs for matclust. The second box asks for a directory on your
local harddrive for quick writing while matclust is running. Anywhere
on the local hardrive is fine, like /data??/yourname. You will probably
never use the files put here– matclust just needs to know where the
local drive is to write temporary files. Then press ‘ok’. Finally, it
is important that the new matclust folder created on your home
directory is on the matlab path. <strong>ADD ANOTHER LINE IN STARTUP.M FOR THE NEW DIRECTORY</strong>.
You should now be able to run matclust from the matlab prompt. Just
type ‘matclust’. If the program opens with no errors, you are good to
go. Now you just need some parameter files to open.
</p>

</div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="5753-7747" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="A. Setting up matclust on a network computer" type="submit"></div></form></div>
<h3><a name="b._making_the_parameter_files" id="b._making_the_parameter_files">B. Making the parameter files</a></h3>
<div class="level3">

<p>The parameter generator currently assumes that you have extracted
the data and that it is organized in the following way (with ‘charley’
as the example animal name): </p>
<pre class="code">animal's raw data directory &gt;&gt;
                     charley01 &gt;&gt;
                           01-120 &gt;&gt;
                               01-120.tt
	</pre>

<p>
and that you have already analyzed the position data to produce a
separate .p file for each epoch of recording. These .p files should be
located in each day’s folder. For an animal named ‘charley’ on
recording day 1, the files should be called charley01_1.p, charley01_2,
..., charley01_6.p for six epochs (you can have more or less). If this
is all true, then in matlab, change your directory to the folder
containing the folders for all the recording days (’raw data directory’
above). Run makedayparms and input the name of the folder that contains
the data you want to cluster: </p>
<pre class="code">&gt;&gt; makedayparms(dayname) 
&gt;&gt; makedayparms(dayname,thresh)
&gt;&gt; makedayparms(dayname,thresh, maxallowedamp)
&gt;&gt; makedayparms(dayname,thresh, maxallowedamp, system)</pre>

<p> The program will read the .tt files and create a paramter .mat file
in the same subdirectory as the .tt file. It will also create the .mat
files containing the waveform info.
</p>
<pre class="code"> THRESH - the threshold that at least one spike must excede in order to be included (default 0 microvolts)
 
 MAXALLOWEDAMP - if the amplitude on any channel excedes this level (in microvolts), exclude the spike. Default 2500.
 
 SYSTEM = 1 for the old rig with spike
 SYSTEM = 2 for the new rigs with nspike (default)</pre>

<p> The parameter files are titled xxxxxx_params.mat. In matclust, go
to file→open→raw parameter file and open up a parameter file. You can
now cluster away. Use the help menu to get some specifics on how to use
the program. NOTE: in order for dayprocess.m (below) to work correctly,
you must name each epoch in the following way: sleep1, sleep2,
sleep3,... for the sleep epochs, and run1, run2, run3, ... for the runs
and export the time filters.
</p>

<p>
<strong>IMPORTANT: THE LAB USES STRICT RULES ON CLUSTERING PROCEDURES.  IF YOU HAVE NOT HAD 
A COMPLETE TUTORIAL FROM SOMEONE, ASK!</strong>
</p>

</div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="7748-9950" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="B. Making the parameter files" type="submit"></div></form></div>
<h2><a name="processing_the_matclust_files" id="processing_the_matclust_files">4. PROCESSING THE MATCLUST FILES</a></h2>
<div class="level2">

<p>The next step is to get the clustered data from the matclust files
into a more manageable state. Before you compile the matclust files,
you need to create a folder where the processed data will go
(’ANIMDIRECT’ below). This folder should be in your local harddrive
(not on hippo), and separate from the animal’s raw data folder.
</p>

<p>Next, run dayprocess.m in matlab to extract compile the matclust
files. Run it from the animal’s raw data directory (example:
/datadrive/username/charley/) Note: if you want to go through all the
days with one function, you can use <code>processdays.m</code> instead.  Type <code>help processdays</code> for info. Here is some info on how to run dayprocess:
</p>
<pre class="code">    &gt;&gt; dayprocess(DAYDIRECT, ANIMDIRECT, FILEPREFIX, DAYNUM, OPTIONS)</pre><ul><li class="level1"><p>DAYDIRECT – folder name where the day’s data is stored (example: ‘charley01’)</p></li><li class="level1"><p>ANIMDIRECT – the path to where the animal’s processed data will be stored – example ‘/datadrive/username/Cha’</p></li><li class="level1"><p>FILEPREFIX
– also the first three letters of the animal’s name (example ‘cha’),
which will attach to the beginning of the .mat files containing the
variables. I recommend using this. If not, use the empty string (<code>''</code>).</p></li><li class="level1"><p>DAYNUM – the day number for the experiment (starting with 1)</p></li><li class="level1"><p>OPTIONS – optional input in the form of: ‘option’, value, ‘option’, value optional arguments are:</p><ul><li class="level2"><p>DIODEPOS
– 0 uses back diode for pos, 1 uses front diode for pos, and values in
between use the proportionate distance between the two diodes. (default
0.5, ie, average of diodes)</p></li><li class="level2"><p>CMPERPIX – size of each pixel in centimeters (default 1)</p></li><li class="level2"><p>POSFILT – filter to use for pos smoothing when computing velocity (default gaussian(30*0.5, 60))</p></li><li class="level2"><p>PROCESSEEG – 1 to process EEG, 0 to skip EEG processing (default 1)</p></li><li class="level2"><p>SYSTEM – 1 for old rig, 2 for nspike (default 2)</p></li><li class="level2"><p>VARPREFIX
– the first three letters of the animals name, which is attached to all
variable names. I recommend using the empty string (<code>''</code>) for this, which will not include any name in the variables. (Default <code>''</code>)</p></li></ul></li></ul><p> Example: 
</p>
<pre class="code">&gt;&gt; dayprocess('charley01', 'Cha', 'cha', 1, 'cmperpix', .9)</pre>

<p> To Calculate CMPERPIX: Put /usr/local/mpgread in your path (if you
used the path list above, this is already done). In matlab, go to a
folder with one of your mpegs, ie /data99/student/charley/Charley01/ </p>
<pre class="code">&gt;&gt;M = mpgread('day1.mpeg', frame, 'truecolor') </pre>

<p>
 Frame is the number of a frame during your run.  To calculate use 30frames/sec and time to first run.  
</p>
<pre class="code">&gt;&gt;a=frame2im(M) to turn frame into image
&gt;&gt;image(a) to show selected frame
&gt;&gt;ginput(2) to create cross hairs on image</pre>

<p> Click on endpoints that you want to measure. Now you can calculate
the number of pixels that equate to a known length (like your track) </p>

<p>
In the processed data directory, you will now have files containing raw
position data, interpolated position data, spike data, and EEG data.
Generally, there is one file for each day of recording, like
chapos01.mat will contain the position data for the first day. If you
load one of these file into matlab, you’ll see that each variable is
organized as a nested cell matrix. For position data, pos{day}{epoch}
will give you the cell containing data for that data and epoch. Spike
information has two more nested cells -
spikes{day}{epoch}{tetrode}{cell}.
</p>

</div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="9951-13333" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="4. PROCESSING THE MATCLUST FILES" type="submit"></div></form></div>
<h2><a name="linearizing_the_data" id="linearizing_the_data">5. LINEARIZING THE DATA</a></h2>
<div class="level2">

<p> If you ran your animal on a track, you will probably want to
linearize the position data. The programs that do this require some
specific information about the track geometry. These programs are
located in /usr/local/linearization. </p><ol><li class="level1"><p>You must run CREATETASKSTRUCT in order to get the track inforamtion into a standardized format.
</p>
<pre class="code">&gt;&gt; createtaskstruct(animdirect,fileprefix,index,coordprogram, options) &lt;
</pre>

<p>
Here, INDEX is a matrix where each row describes an epoch - [day epoch]. 
</p>

<p>
COORDPROGRAM - a string with the name of the function that will provide
the track segment location info. This program should take three inputs
- 1) directoryname 2) pos.data for an epoch, and 3) the epoch index
[day epoch] The output of the program is a cell array, where each cell
describes one trajectory on the track (do not count foreward and
backward movement as 2 trajectories). Inside each cell is a three
dimentional M-by-2-by-N matrix. The matrix gives the x and y coodinates
for the M segment endpoints, across N timeframes. N should equal the
length of pos.data. The segment end points should start with the
position that the user wants to be 0 on the linear scale and progress
to the last segment. All endpoints that are shared by two segments
should only exist once per trajectory. The function is called like
this: coords = coordprogram(direcoryname,pos{day}{epoch}.data,[day
epoch]). If you are using a w track, this program has already been
written: ‘getcoord_wtrack’. Type ‘help getcoord_wtrack’ for info on
what to click with the mouse. </p></li><li class="level1"><p>Now, run LINEARDAYPROCESS for the days you want to linearize. 
</p>
<pre class="code">&gt;&gt; lineardayprocess(animdirect,fileprefix,days, options)</pre>

<p> Type ‘help linearizeposition’ for a detailed description of the
options and the outputs. This program will save a file called ‘linpos’
for each day you specified. For each epoch, a matrix called
‘statematrix’ describes different aspects of the linear behavior for
each time step. &lt; </p></li><li class="level1"><p>To
calculate which trajectory the animal was on for each time step and to
find invalid time points, run GETBEHAVESTATE using the data in linpos.
Type <code>help getbehavestate</code> for more info on this.  </p></li><li class="level1"><p>To
calculate linear occupancy normalized rate maps, run CALCLINFIELDS.
Type ‘help calclinfields’ for more info on this. The output of this
function can be used to plot a linear rate map for each trajectory.</p></li></ol>

</div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="13334-15754" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="5. LINEARIZING THE DATA" type="submit"></div></form></div>
<h2><a name="filtering_the_eeg_data" id="filtering_the_eeg_data">6. FILTERING THE EEG DATA</a></h2>
<div class="level2">

<p>
 There are three eeg filtering functions in /usr/local/filtering:
</p>

<p>
thetadayprocess.m gammadayprocess.m rippledayprocess.m
</p>

<p>
These functions will create new files in your EEG/ directory for theta,
gamma and ripple band activity. All of them have a format similar to
dayprocess, but there are some options that are unique to each. The
following is the recommended command line that can be used for all of
them.
</p>

<p>
xxxxxdayprocess(’/dataxx/user/animaldatafolder/’, ‘animalprefix’, [1:ndays]);
</p>

<p>
where xxxxx is theta, gamma or ripple.
</p>

<p>Note that this will take an hour or more of processing time per day
of data. One suggestion would be to open two matlab windows and run the
functions on half of the day in one window and on the other half in the
other window.
</p>

<p>
The final result of this will be files such as 
</p>

<p>
chatheta-01-2-03.mat chagamma-01-2-03.mat charipple-01-2-03.mat
</p>

<p>
Where 01 is the day, 2 is the epoch and 03 is the tetrode number.
</p>

</div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="15755-16732" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="6. FILTERING THE EEG DATA" type="submit"></div></form></div>
<h2><a name="extracting_ripples_and_high_amplitude_theta_times" id="extracting_ripples_and_high_amplitude_theta_times">7. EXTRACTING RIPPLES AND HIGH AMPLITUDE THETA TIMES</a></h2>
<div class="level2">

<p> The next step is to extract all of the ripple and high amplitude
theta events. We may also want to extract high amplitude gamma at some
point, but that code is not yet written.
</p>

<p>To extract ripples you should run the extractripples() function on
each day of data. You can write a short script that does this: </p>
<pre class="code">directoryname = '/dataxx/name/Cha';
fileprefix = 'cha';
days = 1:10;
tetrodes = -1; % specify all tetrodes
mindur = 0.015; % 15 ms minimum duration
nstd = 2; % 2 std
for d = days
    extractripples(directoryname, fileprefix, d, tetrodes, mindur, nstd);
end</pre>

<p>
 Use the mindur and nstd of 0.015 and 2 above unless you have good reason to do something else.
</p>

<p>
This will create a new data structure in the animal data directory called
</p>

<p>
animalpreripples01.mat for day 1, and so on.
</p>

<p>
 Exactly the same sort of script should be used to extract the high theta amplitude events: 
</p>
<pre class="code">directoryname = '/dataxx/name/Cha';
fileprefix = 'cha';
days = 1:10;
tetrodes = -1; % specify all tetrodes
mindur = 1; % 1 s minimum duration
nstd = 1; % 1 std
for d = days
    extracthightheta(directoryname, fileprefix, d, tetrodes, mindur, nstd); 
end</pre>

<p>
 Once again, use the mindur and nstd of 1 second unless you have reason to do otherwise.
</p>

<p>
 This will create a new data structure in the animal data directory called
</p>

<p>
chahightheta01.mat for day 1, and so on.
</p>

</div>
<div class="secedit"><form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="lines" value="16733-" type="hidden"><input name="rev" value="1217608334" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit" class="button" title="7. EXTRACTING RIPPLES AND HIGH AMPLITUDE THETA TIMES" type="submit"></div></form></div>
    <!-- wikipage stop -->
  </div>

  <div class="clearer">&nbsp;</div>

  
  <div class="stylefoot">

    <div class="meta">
      <div class="user">
              </div>
      <div class="doc">
        data_analysis.txt · Last modified: 2008/08/01 09:32 by 10.1.1.10      </div>
    </div>

   
    <div class="bar" id="bar__bottom">
      <div class="bar-left" id="bar__bottomleft">
        <form class="button" method="post" action="/wiki/doku.php"><div class="no"><input name="do" value="edit" type="hidden"><input name="rev" value="" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Edit this page" class="button" accesskey="e" title="Edit this page [ALT+E]" type="submit"></div></form>        <form class="button" method="get" action="/wiki/doku.php"><div class="no"><input name="do" value="revisions" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Old revisions" class="button" accesskey="o" title="Old revisions [ALT+O]" type="submit"></div></form>      </div>
      <div class="bar-right" id="bar__bottomright">
                                <form class="button" method="get" action="/wiki/doku.php"><div class="no"><input name="do" value="login" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Login" class="button" title="Login" type="submit"></div></form>        <form class="button" method="get" action="/wiki/doku.php"><div class="no"><input name="do" value="index" type="hidden"><input name="id" value="data_analysis" type="hidden"><input value="Index" class="button" accesskey="x" title="Index [ALT+X]" type="submit"></div></form>        <a class="nolink" href="#dokuwiki__top"><input class="button" value="Back to top" onclick="window.scrollTo(0, 0)" title="Back to top" type="button"></a>&nbsp;
      </div>
      <div class="clearer"></div>
    </div>

  </div>

</div>

<div class="footerinc">
  <a href="http://hippo/wiki/feed.php" title="Recent changes RSS feed"><img src="wikipreprocessing_files/button-rss.png" alt="Recent changes RSS feed" height="15" width="80"></a>

  <a href="http://creativecommons.org/licenses/by-nc-sa/2.0/" rel="license" title="Creative Commons License"><img src="wikipreprocessing_files/button-cc.gif" alt="Creative Commons License" height="15" width="80"></a>

  <a href="https://www.paypal.com/xclick/business=andi%40splitbrain.org&amp;item_name=DokuWiki+Donation&amp;no_shipping=1&amp;no_note=1&amp;tax=0&amp;currency_code=EUR&amp;lc=US" title="Donate"><img src="wikipreprocessing_files/button-donate.gif" alt="Donate" height="15" width="80"></a>

  <a href="http://www.php.net/" title="Powered by PHP"><img src="wikipreprocessing_files/button-php.gif" alt="Powered by PHP" height="15" width="80"></a>

  <a href="http://validator.w3.org/check/referer" title="Valid XHTML 1.0"><img src="wikipreprocessing_files/button-xhtml.png" alt="Valid XHTML 1.0" height="15" width="80"></a>

  <a href="http://jigsaw.w3.org/css-validator/check/referer" title="Valid CSS"><img src="wikipreprocessing_files/button-css.png" alt="Valid CSS" height="15" width="80"></a>

  <a href="http://wiki.splitbrain.org/wiki:dokuwiki" title="Driven by DokuWiki"><img src="wikipreprocessing_files/button-dw.png" alt="Driven by DokuWiki" height="15" width="80"></a>



<!--

<rdf:RDF xmlns="http://web.resource.org/cc/"
    xmlns:dc="http://purl.org/dc/elements/1.1/"
    xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
<Work rdf:about="">
   <dc:type rdf:resource="http://purl.org/dc/dcmitype/Text" />
   <license rdf:resource="http://creativecommons.org/licenses/by-nc-sa/2.0/" />
</Work>

<License rdf:about="http://creativecommons.org/licenses/by-nc-sa/2.0/">
   <permits rdf:resource="http://web.resource.org/cc/Reproduction" />
   <permits rdf:resource="http://web.resource.org/cc/Distribution" />
   <requires rdf:resource="http://web.resource.org/cc/Notice" />
   <requires rdf:resource="http://web.resource.org/cc/Attribution" />
   <prohibits rdf:resource="http://web.resource.org/cc/CommercialUse" />
   <permits rdf:resource="http://web.resource.org/cc/DerivativeWorks" />
   <requires rdf:resource="http://web.resource.org/cc/ShareAlike" />
</License>

</rdf:RDF>

-->


</div>

<div class="no"><img src="wikipreprocessing_files/indexer.php" alt="" height="1" width="1"></div>
</body></html>