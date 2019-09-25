Getting Started Using Matclust
******************************


1.Setting up matclust on a network computer

All of the required files are located at /usr/local/matclust on hippo. /usr/local should be mounted on your computer already.  Before you set up matclust on a computer, you need the matlab path to include /usr/local/matclust. Most people in the lab have a startup.m file located on their home directories in the matlab folder.  Add: 
p=path; path('/usr/local/matclust',p);
after the other path definitions.  If no path definitions exist, you also need the line
path(pathdef); before the previous line.  Save the altered startup.m and start matlab.  At the prompt, run setuplocal.  This should bring up a window asking for two directories to  create different files.  The first box asks where the matclust folder should be created.  You should put the folder somewhere on your home directory, like /home/yourname.  Here, a folder named 'matclust' will be created that contains your personal startup defaults and your personal helper programs for matclust.  The second box asks for a directory on your local harddrive for quick writing while matclust is running.  Anywhere on the local hardrive is fine, like /data??/yourname.  You will probably never use the files put here-- matclust just needs to know where the local drive is to write temporary files.  Then press 'ok'.  Finally, it is important that the new matclust folder created on your home directory is on the matlab path.  Add another line in startup.m for the new directory.  You should now be able to run matclust from the matlab prompt.  Just type 'matclust'.  If the program opens with no errors, you are good to go.  Now you just need some parameter files to open.

2.Making the parameter files

The parameter generator currently assumes that you have extracted the data and that it is organized in the following way:

animal directory >>
	day01 >>
		01-120 >>
			01-120.tt
		
and that you have already analyzed the position data to produce a separate .p file for each epoch of recording.  These .p files should be located in each day's folder. For an animal named 'dude' on recording day 1, the files should be called dude01_1.p, dude01_2, ..., dude01_6.p for six epochs.  If this is all true, then in matlab, change your directory to the folder containing the folders for all the recording days ('animal directory' above). Then run makeallparms.  [If you do not want to cluster all your data at once, run makedayparms and input the name of the folder tha contains the data you want to cluster, ie if the folder is called day01 and you used the new rig at the prompt type in makedayparms('day01',2).  If you used the old rig type in makedayparms('day01',1).]  This will create parameter and waveform .mat files in each tetrode's folder.  The parameter files are titled xxxxxx_params.mat.  In matclust, go to file->open->raw parameter file and open up a parameter file. You can now cluster away.  

		





    

