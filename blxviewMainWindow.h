/*
 *  blxviewMainWindow.h
 *  acedb
 *
 *  Created by Gemma Barson on 24/11/2009.
 *
 */

#ifndef _blxview_main_window_included_
#define _blxview_main_window_included_

#include <gtk/gtk.h>
#include <SeqTools/blixem_.h>

typedef struct _MainWindowProperties
  {
    GtkWidget *bigPicture;
    GtkWidget *detailView;
    
    MSP *mspList;		/* Linked list of match sequences */
    const char const *refSeq;	/* The reference sequence */
    
    gboolean strandsToggled;	/* If true, the reverse strand becomes the 'main' or 'top' strand */
    BlxBlastMode blastMode;	/* The type of blast matching that was used */
    BlxSeqType seqType;		/* The type of sequence, e.g. DNA or peptide */
  } MainWindowProperties;


/* Public function declarations */
void blxHelp(void);

MainWindowProperties*	  mainWindowGetProperties(GtkWidget *widget);
gboolean		  mainWindowGetStrandsToggled(GtkWidget *mainWindow);
GtkWidget*		  mainWindowGetBigPicture(GtkWidget *mainWindow);
GtkWidget*		  mainWindowGetDetailView(GtkWidget *mainWindow);
BlxBlastMode		  mainWindowGetBlastMode(GtkWidget *mainWindow);

GtkWidget*		  createMainWindow(char *refSeq, 
					   MSP *mspList, 
					   BlxBlastMode blastMode,
					   BlxSeqType seqType, 
					   int numReadingFrames);


#endif /* _blxview_main_window_included_ */