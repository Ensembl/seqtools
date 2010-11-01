/*
 *  blxparser.h
 *  blixem
 *
 *  Created by Gemma Barson on 02/09/2010.
 *  Copyright 2010 Sanger Institute. All rights reserved.
 *
 *  Description: A parser for old-style exblx or feature-series
 *		 type files. Kept for backwards compatability but
 *		 GFF3 files should be used going forward.
 */


void           parseFS(MSP **MSPlist, FILE *file, char *opts, GList* featureLists[], GList **seqList, GSList *supportedTypes, GSList *styles,
		       char **seq1, char *seq1name, char **seq2, char *seq2name, const int qOffset) ;
