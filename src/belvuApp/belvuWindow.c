/*  File: belvuWindow.c
 *  Author: Gemma Barson, 2011-04-11
 *  Copyright (c) 2009 - 2010 Genome Research Ltd
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
 * Description: see belvuWindow.h
 *----------------------------------------------------------------------------
 */

#include <belvuApp/belvuWindow.h>
#include <gtk/gtk.h>


#define DEFAULT_WINDOW_BORDER_WIDTH      1    /* used to change the default border width around the blixem window */
#define DEFAULT_WINDOW_WIDTH_FRACTION	 0.6  /* what fraction of the screen size the blixem window width defaults to */
#define DEFAULT_WINDOW_HEIGHT_FRACTION	 0.3  /* what fraction of the screen size the blixem window height defaults to */
#define DEFAULT_FONT_SIZE_ADJUSTMENT	 -2   /* used to start with a smaller font than the default widget font */


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

/* Set various properties for the belvu window */
static void setStyleProperties(GtkWidget *widget)
{
  /* Set the initial window size based on some fraction of the screen size */
  GdkScreen *screen = gtk_widget_get_screen(widget);
  const int width = gdk_screen_get_width(screen) * DEFAULT_WINDOW_WIDTH_FRACTION;
  const int height = gdk_screen_get_height(screen) * DEFAULT_WINDOW_HEIGHT_FRACTION;
  
  gtk_window_set_default_size(GTK_WINDOW(widget), width, height);
  
  gtk_container_set_border_width (GTK_CONTAINER(widget), DEFAULT_WINDOW_BORDER_WIDTH); 
  gtk_window_set_mnemonic_modifier(GTK_WINDOW(widget), GDK_MOD1_MASK); /* MOD1 is ALT on most systems */
  
  /* Set the default font size to be a bit smaller than usual */
  int origSize = pango_font_description_get_size(widget->style->font_desc) / PANGO_SCALE;
  const char *origFamily = pango_font_description_get_family(widget->style->font_desc);

  char parseString[500];
  sprintf(parseString, "gtk-font-name = \"%s %d\"", origFamily, origSize + DEFAULT_FONT_SIZE_ADJUSTMENT);
  gtk_rc_parse_string(parseString);
}


gboolean createBelvuWindow(BelvuContext *bc)
{
  gboolean ok = TRUE;
  
  GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setStyleProperties(window);
  

//  graphRegister(PICK, boxPick);
//  graphRegister(MIDDLE_DOWN, middleDown);
//  graphRegister(RESIZE, belvuRedraw);
//  graphRegister(KEYBOARD, keyboard);
//  graphRegister(DESTROY, belvuDestroy) ;
//  

  
  gtk_widget_show_all(window);
  

  return ok;
}
