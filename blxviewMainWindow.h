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
#include <SeqTools/blxview.h>

typedef struct _MainWindowProperties
  {
    GtkWidget *bigPicture;
    GtkWidget *detailView;
    
    gboolean strandsToggled; /* If true, views look at the rev strand if their default is fwd or v.v. */
  } MainWindowProperties;


/* Public function declarations */
MainWindowProperties*	  mainWindowGetProperties(GtkWidget *widget);
GtkWidget*		  createMainWindow(char *refSeq, MSP *mspList, int numReadingFrames);


#endif /* _blxview_main_window_included_ */