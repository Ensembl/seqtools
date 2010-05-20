/*
 *  blxdotter.h
 *  acedb
 *
 *  Created by Gemma Barson on 03/02/2010.
 *
 *  Description:
 *    Provides functions for a blixem window to be able to call dotter. Includes
 *    utility functions to filter matches, find a suitable range to initiate dotter
 *    with etc.
 */

#ifndef _blx_dotter_h_included_
#define _blx_dotter_h_included_


#include <SeqTools/blxview.h>
#include <gtk/gtk.h>


void			showDotterDialog(GtkWidget *blxWindow);
gboolean		callDotter(GtkWidget *blxWindow, const gboolean hspsOnly, GError **error);


#endif /* _blx_dotter_h_included_ */
