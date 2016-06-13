/*  File: blxcontext.hpp
 *  Author: Gemma Barson, 2016-04-06
 *  Copyright (c) 2016 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * SeqTools is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
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
 * ---------------------------------------------------------------------------
 * This file is part of the SeqTools sequence analysis package, 
 * written by
 *      Gemma Barson      (Sanger Institute, UK)  <gb10@sanger.ac.uk>
 * 
 * based on original code by
 *      Erik Sonnhammer   (SBC, Sweden)           <Erik.Sonnhammer@sbc.su.se>
 * 
 * and utilizing code taken from the AceDB and ZMap packages, written by
 *      Richard Durbin    (Sanger Institute, UK)  <rd@sanger.ac.uk>
 *      Jean Thierry-Mieg (CRBM du CNRS, France)  <mieg@kaa.crbm.cnrs-mop.fr>
 *      Ed Griffiths      (Sanger Institute, UK)  <edgrif@sanger.ac.uk>
 *      Roy Storey        (Sanger Institute, UK)  <rds@sanger.ac.uk>
 *      Malcolm Hinsley   (Sanger Institute, UK)  <mh17@sanger.ac.uk>
 *
 * Description: A Blixem context class, containing all status information 
 *              required for a blixem instance.
 *----------------------------------------------------------------------------
 */

#ifndef DEF_BLXCONTEXT_H
#define DEF_BLXCONTEXT_H

#include <gtk/gtk.h>
#include <blixemApp/blixem_.hpp>



class BlxContext
{
public:
  // Constructors
  BlxContext(CommandLineOptions *options,
             const IntRange* const refSeqRange_in,
             const IntRange* const fullDisplayRange_in,
             const char *paddingSeq_in,
             GArray* featureLists_in[],
             GList *seqList_in,
             GSList *supportedTypes_in,
             GtkWidget *widget_in,
             GtkWidget *statusBar_in,
             const gboolean External_in,
             GSList *styles_in);

  ~BlxContext();

  // Access
  double charWidth() const;
  double charHeight() const;
  BlxStrand activeStrand() const;

  // Query
  bool isSeqSelected(const BlxSequence *seq) const;
  SequenceGroup *getSequenceGroup(const BlxSequence *seqToFind) const;
  GList *getSelectedSeqsByType(const BlxSequenceType type) const;
  BlxSequence* getSelectedTranscript(int *num_transcripts) const;

  void highlightBoxCalcBorders(GdkRectangle *drawingRect, GdkRectangle *highlightRect, 
                               const IntRange *fullRange, const IntRange *highlightRange, 
                               const int yPadding);

  // Modify
  void saveSettingsFlags(GKeyFile *key_file);
  void killAllSpawned();

  void calculateDepth(const int numUnalignedBases);
  int calculateTotalDepth(const IntRange *range, const BlxStrand strand);
  int getDepth(const int coord, const char *base_char = NULL, const BlxStrand strand = BLXSTRAND_NONE);
  int getDepthForCounter(const int coord, const DepthCounter counter);

  SequenceGroup* getQuickGroup(const bool isFilter);
  void destroySequenceGroup(SequenceGroup **seqGroup);
  void deleteAllSequenceGroups();
  void deleteAllQuickGroups();


  GtkWidget *statusBar;                   /* The Blixem window's status bar */
    
  char *refSeq;                           /* The reference sequence (always forward strand, always DNA sequence) */
  const char *refSeqName;                 /* The name of the reference sequence */
  IntRange refSeqRange;                   /* The range of the reference sequence */
  IntRange fullDisplayRange;              /* The range of the displayed sequence */
  int refSeqOffset;                       /* how much the coordinate system has been offset from the input coords */

  BlxBlastMode blastMode;                 /* The type of blast matching that was used */
  BlxSeqType seqType;                     /* The type of the match sequences, e.g. DNA or peptide */
  char **geneticCode;                     /* The genetic code used to translate DNA <-> peptide */
  int numFrames;                          /* The number of reading frames */

  GArray *bulkFetchDefault;               /* The default method(s) of bulk fetching sequences (can be overridden by an MSPs data-type properties) */
  GArray *userFetchDefault;               /* The default method(s) for interactively fetching individual sequences */
  GArray *optionalFetchDefault;           /* The default method(s) for bulk fetching optional sequence data */
  GHashTable *fetchMethods;               /* List of fetch methods, keyed on name as a GQuark */

  char* dataset;                          /* the name of a dataset, e.g. 'human' */
  gboolean optionalColumns;               /* load data for optional columns on startup to populate additional info like tissue-type */
  gboolean saveTempFiles;                 /* whether to save temporary files created to store results of file conversions */

  MSP *mspList;                           /* List of all MSPs. Obsolete - use featureLists array instead */
  GArray* featureLists[BLXMSP_NUM_TYPES]; /* Array indexed by the BlxMspType enum. Each array entry contains a zero-terminated array of all the MSPs of that type. */
    
  GList *matchSeqs;                       /* List of all match sequences (as BlxSequences). */
  GSList *supportedTypes;                 /* List of supported GFF types */
  const char *paddingSeq;                 /* A sequence of padding characters, used if the real sequence could not be found. All padded MSPs
                                           * use this same padding sequence - it is constructed to be long enough for the longest required seq. */
    
  gboolean displayRev;                    /* True if the display is reversed (i.e. coords decrease as you read from left to right, rather than increase). */
  gboolean external;                      /* True if Blixem was run externally or false if it was run internally from another program */
    
  GList *selectedSeqs;                    /* A list of sequences that are selected (as BlxSequences) */
  GList *sequenceGroups;                  /* A list of SequenceGroups */
    
  DotterRefType dotterRefType;            /* Whether to dotter a ref seq range or a transcript */
  DotterMatchType dotterMatchType;        /* Saved type of match to call dotter on */
  char *dotterAdhocSeq;                   /* Saves the sequence text the user pastes into the dotter dialog */
  gboolean dotterHsps;                    /* Whether the dotter "HSPs only" option is on by default */
  gboolean dotterSleep;                   /* Whether the sleep-dotter option is on by default */
  int dotterStart;                        /* Start coord to call dotter on, or UNSET_INT to calculate automatically */
  int dotterEnd;                          /* End coord to call dotter on, or UNSET_INT to calculate automatically */
  int dotterZoom;                         /* Zoom param to call dotter with, if using manual params */
    
  GArray *defaultColors;                  /* Default colors used by Blixem */
  gboolean usePrintColors;                /* Whether to use print colors (i.e. black and white) */
  char *windowColor;                      /* If not null, background color for the window */

  GList *columnList;                      /* A list of details about all the columns in the detail view (might have been better to use an array here but it's a short list so not important) */
  GSList *styles;
    
  gboolean flags[BLXFLAG_NUM_FLAGS];              /* Array of all the flags the user can toggle. Indexed by the BlxFlags enum. */
  GtkWidget *dialogList[BLXDIALOG_NUM_DIALOGS];   /* Array of all the persistent dialogs in the application */
  GSList *spawnedProcesses;                       /* List of processes spawned by Blixem */
  BlxModelId modelId;                             /* which tree model to use (i.e. normal or squashed) */
  
  /* This array holds the depth (num alignments) at each coord of the ref seq for
   * different types of counter. For each counter there is an array the same lenth as the
   * reference sequence. This could cause memory problems with a large region. Generally
   * though the length of the reference sequence is not that great and we are more concerned with
   * fast access as the user interacts with the display.*/
  int *depthArray[DEPTHCOUNTER_NUM_ITEMS];
  int minDepth;                           /* minimum value in the depthArray */
  int maxDepth;                           /* maximum value in the depthArray */

  long ipresolve;                         /* specify whether curl should use ipv4/ipv6 */
  const char *cainfo;                     /* specify location of curl cainfo file */
  bool fetch_debug;                       /* enable verbose debug output in fetch methods */

private:
  void createColors(GtkWidget *widget);
  void initialiseFlags(CommandLineOptions *options);
  void loadSettings();

  double m_charWidth;
  double m_charHeight;
};


#endif /* DEF_BLXCONTEXT_H */
