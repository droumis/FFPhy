u32 ParseTimestamp(char *s);
int IsStringEmpty(char *s);
char *TimestampToString();
int strcount(char *s, char c);
int readcontrec(FILE *infile, ContRec *contrec);
int readposrec(FILE *infile, PosRec *posrec);
int readttrec(FILE *ttfile, TTRec *ttrec);

