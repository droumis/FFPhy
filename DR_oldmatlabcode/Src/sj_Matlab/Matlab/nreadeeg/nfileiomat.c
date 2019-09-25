#include "mexheader.h"
#include "mex.h"
#include "fileiomat.h"


int readposrec(FILE *file, PosRec *posrec)
   /* returns 0 on error */
{
   u32 time;
   short s[4];

   if (fread(&time, sizeof(u32), 1, file) != 1) 
       return 0;
   if (fread(s, sizeof(short), 4, file) != 4) 
       return 0;
   
   posrec->time = time;
   posrec->x1 = s[0];
   posrec->y1 = s[1];
   posrec->x2 = s[2];
   posrec->y2 = s[3];
   return 1;
}

int readttrec(FILE *file, TTRec *ttrec)
   /* returns 0 on error */
{
	u32 time;

   	if (fread(&(ttrec->time), sizeof(u32), 1, file) != 1) 
       	return 0;
   	if (fread(ttrec->x, sizeof(unsigned short), TTRECSIZE, file) != TTRECSIZE) 
       	return 0;
  	 
   return 1;
}

int readcontrec(FILE *file, ContRec *contrec)
{
	u32 time;

   	if (fread(&(contrec->timestamp), sizeof(u32), 1, file) != 1) 
       	return 0;
   	if (fread(&(contrec->numsamples), sizeof(int), 1, file) != 1) 
       	return 0;
   	if (fread(&(contrec->sampfreq), sizeof(double), 1, file) != 1) 
       	return 0;

	/* allocate space for the data */
	contrec->data = (short *) mxCalloc(sizeof(short), contrec->numsamples);
	if (fread(contrec->data, sizeof(short), contrec->numsamples, file) != 
		contrec->numsamples) 
		return 0;
  	 
   return 1;
}

int readnewcontrec(FILE *file, NewContRec *contrec, int nsamp)
{
	u32 time;

   	if (fread(&(contrec->timestamp), sizeof(u32), 1, file) != 1) 
       	return 0;

	/* allocate space for the data if necessary*/
	if (contrec->data == NULL) {
	    contrec->data = (short *) mxCalloc(sizeof(short), nsamp);
	}
	if (fread(contrec->data, sizeof(short), nsamp, file) != 
		nsamp) 
		return 0;
  	 
   return 1;
}

u32 ParseTimestamp(char *s)
{
char	*ptr;
char	*ptr2;
u32	time;
u32 hour;
u32 min;
u32 sec;
float fracsec;
int	ncolons;
char	*fracptr;
char	timestr[100];	

    if(s == NULL){
	return(0);
    }
    /*
    ** copy the passed argument to the timestring for
    ** manipulation
    */
    strcpy(timestr,s);
    /*
    ** check for hr:min:sec.fracsec format vs min:sec
    */
    ncolons = strcount(timestr,':');
    fracsec = 0;
    if((fracptr = strchr(timestr,'.')) != NULL){
	sscanf(fracptr,"%f",&fracsec);
	*fracptr = '\0';
    };
    switch(ncolons){
    case 0:
	if(fracptr){
	    sscanf(timestr,"%d",&sec);
	    time = (u32) (sec*1e4 + (fracsec*1e4 + 0.5)); 
	} else {
	    /*
	    ** straight timestamp
	    */
	    sscanf(timestr,"%d",&time);
	}
	break;
    case 1:
	/*
	** find the colon
	*/
	ptr = strchr(timestr,':');
	/*
	** separate the minutes and the seconds into two strings
	*/
	*ptr = '\0';
	/*
	** read the minutes before the colon
	*/
	sscanf(timestr,"%d",&min);
	/*
	** read the seconds after the colon
	*/
	sscanf(ptr+1,"%d",&sec);
	/*
	** compute the timestamp
	*/
	time = (u32) (min*6e5 + sec*1e4 + (fracsec*1e4 + 0.5)); 
	break;
    case 2:
	/*
	** find the first colon
	*/
	ptr = strchr(timestr,':');
	/*
	** find the second colon
	*/
	ptr2 = strchr(ptr+1,':');
	/*
	** separate the hours, minutes and the seconds into strings
	*/
	*ptr = '\0';
	*ptr2 = '\0';
	/*
	** read the hours before the first colon
	*/
	sscanf(timestr,"%d",&hour);
	/*
	** read the minutes before the second colon
	*/
	sscanf(ptr+1,"%d",&min);
	/*
	** read the seconds after the colon
	*/
	sscanf(ptr2+1,"%d",&sec);
	/*
	** compute the timestamp
	*/
	time = (u32) (hour*36e6 + min*6e5 + sec*1e4 + (fracsec*1e4 + 0.5)); 
	break;
    default:
	fprintf(stderr,"unable to parse timestamp '%s'\n",timestr);
	return(0);
    }
    return(time);
}

int IsStringEmpty(char *str)
{
    if(str == NULL) return(1);
    /*
    ** scan the string to see if there are any non-white space
    ** characters in it
    */
    while(str && (*str != '\0')){
	if((*str != ' ') && (*str != '\t') && (*str != '\n')){
	    /*
	    ** found a non-white space character
	    */
	    return(0);
	}
	str++;
    }
    /*
    ** all white space
    */
    return(1);
}

char *TimestampToString(u32 timestamp)
{
int	hour;
int	min;
int	sec;
static char string[20];
double	fracsec;

    hour = (int)((timestamp/1e4)/3600);
    min = (int)((timestamp/1e4)/60 - 60*hour);
    sec = (int)(timestamp/1e4 - 60*min - 3600*hour);
    fracsec = ((double)timestamp/1.0e4 - (int)(timestamp/1e4))*1e4 + 0.5;
    if(hour > 0){
	sprintf(string,"%d:%02d:%02d.%04d",hour,min,sec,(int)fracsec);
    } else {
	sprintf(string,"%02d:%02d.%04d",min,sec,(int)fracsec);
    }
    return(string);
}

void FormatTime(u32 timestamp,int *min, int *sec)
{
    *min = (int)((timestamp/1e4)/60);
    *sec = (int)(timestamp/1e4 - 60*(int)((timestamp/1e4)/60));
}

int strcount(char *s,char c)
{
int	count;

    if(s == NULL) return(0);
    count = 0;
    while(*s != '\0'){
		if(*s == c) count++;
		s++;
    }
    return(count);
}

