u32 ParseTimestamp(char *s);
int IsStringEmpty(char *s);
char *TimestampToString();
int strcount(char *s, char c);
int readcontrec(FILE *infile, ContRec *contrec);
int readnewcontrec(FILE *infile, NewContRec *contrec, int nsamp);
int readposrec(FILE *infile, PosRec *posrec);
int readttrec(FILE *ttfile, TTRec *ttrec);

