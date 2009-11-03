/*  Last edited: Apr  3 13:19 1995 (esr) */

/******************************************************************************** 
*
*   EFETCH.C
*
*   Fetches a sequence from a database using index files in the EMBL CD-ROM format
*
*   Erik Sonnhammer, 930210
*
*********************************************************************************


  Date        Modification
--------   ----------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * ----------------------------------------------------
93-11-17  Added database: prefix of queries
93-12-05  Added Mini output format and Mixed db format
93-12-21  Added -A and ProDom support
*/

/* $Id: efetch.c,v 1.1 2009-11-03 18:28:23 edgrif Exp $ */

#include "dbidx.h"

extern char dbsource[32];
char  *env, dbdir[128]="", *seqName;
int    wormpep2swiss = 0;

void swissprot(void)
{
    if (env = getenv("SWDIR")) strcpy(dbdir, env);
    else strcpy(dbdir, DEFAULT_SWDIR);
    strcpy(dbsource, "S"); 
}
void pir(void)
{
    if (env = getenv("PIRDIR")) strcpy(dbdir, env);
    else strcpy(dbdir, DEFAULT_PIRDIR);
    strcpy(dbsource, "P"); 
}
void wormpep(void)
{
    if (env = getenv("WORMDIR")) strcpy(dbdir, env);
    else strcpy(dbdir, DEFAULT_WORMDIR);
    strcpy(dbsource, "W");
    fasta = 1;
}
void embl(void)
{
    if (env = getenv("EMBLDIR")) strcpy(dbdir, env);
    else strcpy(dbdir, DEFAULT_EMBLDIR);
    strcpy(dbsource, "E");
}
void genbank(void)
{
    if (env = getenv("GBDIR")) strcpy(dbdir, env);
    else strcpy(dbdir, DEFAULT_GBDIR);
    strcpy(dbsource, "G");
}
void prodom(void)
{
    if (env = getenv("PRODOMDIR")) strcpy(dbdir, env);
    else strcpy(dbdir, DEFAULT_PRODOMDIR);
    strcpy(dbsource, "D");
}
void prosite(void)
{
    if (env = getenv("PROSITEDIR")) strcpy(dbdir, env);
    else strcpy(dbdir, DEFAULT_PROSITEDIR);
    strcpy(dbsource, "R");
    fasta = 0;
}

void main(int argc, char *argv[])
{
  FILE        *idx, *db, *actrg, *achit;
  static char  idxfile[128]="", dbfile[128]="", divfile[128]="", 
               actrgfile[128]="", achitfile[128]="", dbprefix[32]="";
  static char *file, customName[128]="";
  char         query[32], *tmpstr, *cp, fetchstr[256]; 
  int          optc, pos, i;
  int          Startseq = -1, Endseq = -1;
  extern int   optind;
  extern char *optarg;
  char        *optstring="haAB:fiqs:S:e:E:D:d:xl:HopbrmMn:w"; 

  Header      *idxhead, *actrghead, *achithead;
  Entrynam    *rec;
  Acnum       *ACrec;

  static char *cc_date = 
#if defined(__DATE__)
      __DATE__
#else
      ""
#endif
      ;

  char *usage;
  static char usageText[] = "\n\
  EFETCH - retrieve entries from sequence databases.\n\
\n\
  Usage: efetch -options [database:]<query> \n\
\n\
    Databases:  SWissprot/SP, PIR, WOrmpep/WP, EMbl, GEnbank/GB, ProDom, ProSite\n\
\n\
    Options:\n\
    -a            Search with Accession number\n\
    -f            Fasta format output\n\
    -q            Sequence only output (one line)\n\
    -s <#>        Start at position #\n\
    -e <#>        Stop at position #\n\
    -o            More options and info...\n\
\n\
  by Erik Sonnhammer (esr@sanger.ac.uk)\n\
  Version 2.1, ";

  if (!(usage = malloc(strlen(usageText) + 102))) {
      fprintf(stderr, "Malloc error\n");
      exit(1);
  }

  strcpy(usage, usageText);
  strncat(usage, cc_date, 100);
  strcat(usage, "\n");

  /* parse command line */
  while ((optc = getopt(argc, argv, optstring)) != -1)
  switch (optc) {
   
  case 'h':
      fprintf(stderr, "\n%s\n", usage); exit(1);
  case 'B':
      strcpy(dbsource, optarg); break;
  case 'a':
      ACsearch = 1; break;
  case 'A':
      AC2entryname = 1; ACsearch = 1; break;
  case 'b': 
      Reversebytes = 0; break;
  case 'f':
      fasta = 1; break;
  case 'i':
      strcpy(idxfile, optarg); break;
  case 'q':
      fasta = 1;
      seqonly = 1; break;
  case 'D':
      strcpy(dbdir, optarg); break;
  case 'd':
      strcpy(dbfile, optarg); break;
  case 'H':
      display_head=1; break;
  case 'o':
      printf("\n\
         -D <dir>      Specify database directory\n\
         -H            Display index header data\n\
         -p            Display entrynames in search path\n\
         -r            Print sequence in 'raw' format\n\
         -m            Fetch from mixed mini database\n\
         -M            Mini format output\n\
         -b            Do NOT reverse the order of bytes\n\
                              (SunOS, IRIX do reverse, Alpha not)\n\
         -d <dbfile>   Specify database file (avoid this)\n\
	 -i <idxfile>  Specify index file (avoid this)\n\
         -l <divfile>  Specify division lookup table (avoid this)\n\
         -B <database> Specify database (archaic)\n\
         -A            Only return entryname for accession number\n\
         -n <name>     Give the sequence this name\n\
         -x            Don't require query to match entry's name exactly (avoid)\n\
         -w            For Wormpep: also fetch cross-referenced SwissProt entry\n\
\n\
\n\
  Environment variables:\n\
  SWDIR      = SwissProt  directory - database and EMBL index files\n\
  PIRDIR     = PIR        -- \" --\n\
  WORMDIR    = Wormpep    -- \" --\n\
  EMBLDIR    = EMBL       -- \" --\n\
  GBDIR      = Genbank    -- \" --\n\
  PRODOMDIR  = ProDom     -- \" --\n\
  PROSITEDIR = ProSite    -- \" --\n\
  DBDIR      = User's own -- \" -- (fasta format)\n\
\n\
  SEQDB    database file (default SwissProt)\n\
  SEQDBIDX index file\n\
  DIVTABL  division lookup table\n\
\n\
  Ex. setenv DBDIR /pubseq/seqlibs/embl/\n\
\n\
  Note that Prodom family consensus seqs can be fetched by PD:_#\n\n");
      exit(0);
  case 'p':
      search = 1; break;
  case 'l':
      strcpy(divfile, optarg); break;
  case 's':
  case 'S':
      Startseq = atoi(optarg);
      fasta = 1;               break;
  case 'e':
  case 'E':
      Endseq = atoi(optarg); 
      fasta = 1;               break;
  case 'r':
      raw_seq = 1;             break;
  case 'x':
      Exact_match = 0;         break;
  case 'm':
      Mixdb = 1;               break;  /* Mixed source database */
  case 'M':
      Mini = 1;  /* Mini output format */
      fasta = 1;               break;
  case 'n':
      strcpy(customName, optarg); break;
  case 'w':
      wormpep2swiss = 1;       break;
  default:
      fprintf(stderr, "\n%s\n", usage); exit(1);      
  }

  if (argc - optind != 1)
  { 
    fprintf(stderr, "%s\n", usage); exit(1);      
  }
  else strcpy(query, argv[argc - 1]);

  for (i=0; i<strlen(query); i++) query[i] = toupper(query[i]);
  if (tmpstr = (char *)strchr(query, '$')) *tmpstr = '_';
  if (tmpstr = (char *)strchr(query, '*')) {
      *tmpstr = 0;
      Exact_match = 0;
  }

  for (i=0; i<strlen(dbsource); i++) dbsource[i] = toupper(dbsource[i]);

#ifdef NOREVBYTES
  Reversebytes = 0;
#endif

#ifdef ALPHA
  Reversebytes = 0;
#endif

/***********************************  Database directory *******************************/

/* Parse the database source */

#ifdef MIXDB
  Mixdb = 1;
#endif

  if (Mixdb)
  {
      if (env = getenv("DBDIR")) strcpy(dbdir, env);
      dbsource[0] = 'M';
  }
  else if (tmpstr = (char *)strchr(query, ':'))
  {
      strcpy(dbprefix, query);
      dbprefix[tmpstr-query+1] = 0;

      if (!strncmp(query, "SP", 2)) swissprot();
      else if (!strncmp(query, "SW", 2)) swissprot();
      else if (!strncmp(query, "PIR", 3)) pir();
      else if (!strncmp(query, "WP", 2)) wormpep();
      else if (!strncmp(query, "WO", 2)) wormpep();
      else if (!strncmp(query, "EM", 2)) embl();
      else if (!strncmp(query, "GB", 2)) genbank();
      else if (!strncmp(query, "GE", 2)) genbank();
      else if (!strncmp(query, "PROD", 4)) prodom();
      else if (!strncmp(query, "PD", 2)) prodom();
      else if (!strncmp(query, "PROS", 4)) prosite();
      else if (!strncmp(query, "PS", 2)) prosite();
      else fprintf(stderr, "Sorry, I don't support the database prefix %s:\n", dbprefix);

      strcpy(query, tmpstr+1);
      if (!Mixdb && !Mini) *dbprefix = 0;
  }
  else /* Use -B specifier (Old way - keep for compatibility) */
  if (strcmp(dbsource, ""))
  {
      if (!strncmp(dbsource, "S", 1)) swissprot();
      else if (!strncmp(dbsource, "P", 1)) pir();
      else if (!strncmp(dbsource, "W", 1)) wormpep();
      else if (!strncmp(dbsource, "E", 1)) embl();
      else if (!strncmp(dbsource, "G", 1)) genbank();
      else if (!strncmp(dbsource, "D", 1)) prodom();
      else {
	  fprintf(stderr, "Sorry, I don't support the database %s\n", dbsource);
	  exit(-11);
      }
  }
  else /* Use -D database dir */
      if (strcmp(dbdir, ""));
  else /* Use environment DBDIR database dir */
      if (env = getenv("DBDIR")) {
	  strcpy(dbdir, env);
	  strcpy(dbsource, "W");
  }
  else /* Try to make an intelligent guess which database the user wants */
  {
      if (ACsearch) swissprot();
      else if (strchr(query, '_') || strlen(query) < 5 ) swissprot();
      else if (strchr(query, '.')) wormpep();
      else pir();
  }


/*********************************************************************************************/

  /* Open files
     (1) from command line seqdb
     (2) from environment SEQDB
     (3) defaults from dbidx.h
  */

  /* Entryname index file */
  if (!strcmp(idxfile, "")) 
  {
    env = getenv("SEQDBIDX");
    if (env) strcpy(idxfile, env);
      else sprintf(idxfile, "%s%s", dbdir, DEFAULT_IDXFILE);
  }
  if ((idx = fopen(idxfile, "r")) == NULL)
  {
          fprintf(stderr,"Could not open the index file named:  %s\n", idxfile);
	  exit(1);
  }

  /* Division lookup table file */
  if (!strcmp(divfile, ""))
  {
    env = getenv("DIVTABL");
    if (env != NULL) strcpy(divfile, env);
      else sprintf(divfile, "%s%s", dbdir, DEFAULT_DIVFILE);
  }
    

  idxhead = readheader(idx);
  if (!idxhead) {
    fprintf(stderr,"\nError reading index header - possibly corrupted\n");
    exit(1);
  }
  if (display_head) printheader(idxhead);


  
  if (ACsearch)
  {
    /* Accession TARGET file (contains AC numbers */
    sprintf(actrgfile, "%s%s", dbdir, DEFAULT_ACTRGFILE);
    if ((actrg = fopen(actrgfile, "r")) == NULL) 
	{
	    fprintf(stderr,"Couldn't open Accession Target file %s\n", actrgfile);
	    exit(1);
	}
    actrghead = readheader(actrg);
    if (!actrghead) {
	fprintf(stderr,"\nError reading index header - possibly corrupted\n");
	exit(1);
    }
    if (display_head) printheader(actrghead);

    /* Accession HIT file (contains pointers to Entrynam.idx) */
    sprintf(achitfile, "%s%s", dbdir, DEFAULT_ACHITFILE);
    if ((achit = fopen(achitfile, "r")) == NULL) 
	{
	    fprintf(stderr,"Couldn't open Accession Hit file %s\n", achitfile);
	    exit(1);
	}
    achithead = readheader(achit);
    if (!achithead) {
	fprintf(stderr,"\nError reading index header - possibly corrupted\n");
	exit(1);
    }
    if (display_head) printheader(achithead);


/********************** Find Accession Number ************************/

    if ((pos = bin_search(query, actrg, actrghead, getAcname)) == -1) 
    {
	fprintf(stderr,"\n%s not found in index %s\n", query, actrgfile);
	exit(1);
    }
    if (!(ACrec = readActrg(actrg, actrghead, pos)))
    {
	fprintf(stderr,"\nError: no record at position %d in index %s\n", 
		pos, actrgfile);
	exit(1);
    }
    if (display_head) printActrg(ACrec);

    if (AC2entryname || ACrec->nhits > 1) 
    {
	if (!AC2entryname) /* Multiple hits - test with P12713 */
	    printf("Warning: Non-primary accession number, shared by these entries:\n");

	for (i=0; i < ACrec->nhits; i++)
	{
	    pos = getAcentry_rec_num(achit, achithead, (ACrec->hit_rec_num+i-1) );
	    
	    if (!(rec = readEntrynam(idx, idxhead, pos)))
	    { 
		fprintf(stderr,"\nError: no record at position %d in index %s\n", 
			pos, idxfile);
		exit(1);
	    }
	    printf("%s\n", rec->entry_name);
	}
	exit(-1);
    }
    else /* Primary AC# */
    {
	pos = getAcentry_rec_num(achit, achithead, (int)(ACrec->hit_rec_num-1) );
    }
    fclose(actrg);
    free(actrghead);
    fclose(achit);
    free(achithead);
  }
  else /* Normal, primary index search */
  {
    /* Find the record closest to query in idxhead by binary search */
    pos = bin_search(query, idx, idxhead, getname);
    if (pos == -1)
    {
	fprintf(stderr,"\n%s not found in index %s\n", query, idxfile);
	exit(1);
    } 
  }

  if (!(rec = readEntrynam(idx, idxhead, pos)))
  { 
	fprintf(stderr,"\nError: no record at position %d in index %s\n", 
		pos, idxfile);
	exit(1);
  }
  if (display_head) printEntrynam(rec, idxhead);

  /* Now we have the division code and can open the database file
     Use division code only if dbfile is not specified on command line 
     or environment.  [ Cmdline - env - div.lkp - default ]
  */
  if (!strcmp(dbfile, "")) 
  {
    env = getenv("SEQDB");
    if (env != NULL) strcpy(dbfile, env);
    else
    {
      file = getdbfile(rec->div_code, divfile);
      if (file != NULL) sprintf(dbfile, "%s%s", dbdir, file);
      else
      {
        sprintf(dbfile, "%s%s", dbdir, DEFAULT_DBFILE);
        fprintf(stderr, "Warning: couldn't open division lookup table %s \n"
		"Using Default Database file %s\n", divfile, dbfile);
      }
    }
  }
  if (display_head) printf("Database file: %s \n", dbfile);

  /* Now the database itself */
  if ((db = fopen(dbfile, "r")) == NULL) {
          fprintf(stderr,"Could not find or open the database file named:  %s\n", dbfile);
	  exit(1);
  }

/* If database is chosen with -B -f, make sure the correct format is engaged.
     
     -B PIR  -f -> printPIRseq
     -B GB   -f -> printGBseq
     -B SW   -f -> printEMBLseq
     -B EMBL -f -> printEMBLseq
*/

  if (*customName) seqName = customName;
  else seqName = rec->entry_name;

  if (fasta)
  {
      if(!strncmp(dbsource, "P", 1))
      {
        if (!seqonly) printf(">%s%s ", dbprefix, seqName);
        printPIRseq(db, rec, Startseq, Endseq);
      }
      else if(*dbsource == 'W' || *dbsource == 'M' || *dbsource == 'D')
      {
	  if (!seqonly) printf(">%s%s ", dbprefix, seqName);
	  tmpstr = printFastaseq(db, rec, Startseq, Endseq);

	  /* Fetch cross-referenced swissprot entry */
	  if (wormpep2swiss && *dbsource == 'W' && (cp = (char *)strstr(tmpstr, "SW:"))) {
	      if (tmpstr = (char *)strchr(cp, ' ')) *tmpstr = 0;
	      strcpy(fetchstr, "efetch -a ");
	      strcat(fetchstr, cp);
	      printf("\nFetching cross-referenced swissprot entry: %s\n", fetchstr);
	      system(fetchstr);
	  }
      }
      else if(!strncmp(dbsource, "G", 1))
      {
        if (!seqonly) printf(">%s%s ", dbprefix, seqName);
        printGBseq(db, rec, Startseq, Endseq);
      }
      else  /* Assume anything else is EMBL format (same for SwissProt) */
      {
	if (!seqonly) printf(">%s%s ", dbprefix, seqName);
        printEMBLseq(db, rec, Startseq, Endseq);
      }
  }
  else printAnn(db, rec);

  /* Clean up */
  fclose(idx);
  fclose(db);
  free(idxhead);
}
