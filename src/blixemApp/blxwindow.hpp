/*  File: blxWindow.h
 *  Author: Gemma Barson, 2009-11-24
 *  Copyright [2018-2023] EMBL-European Bioinformatics Institute
 *  Copyright (c) 2006-2017 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
 * Description: Creates the main Blixem window. Also creates a "context",
 *              which contains all of the variables associated with a Blixem
 *              session.
 *
 *              The context could live somewhere else, but was just included
 *              here because there is one context for each Blixem window.
 *              Ideally the context would be created earlier (perhaps in the
 *              main function) but so far that has not been practical.
 *----------------------------------------------------------------------------
 */

#ifndef _blxwindow_included_
#define _blxwindow_included_

#include <gtk/gtk.h>
#include <blixemApp/blixem_.hpp>
#include <seqtoolsUtils/utilities.hpp>


/* Public function declarations */
BlxContext*		  blxWindowGetContext(GtkWidget *widget);
GList*                    blxWindowGetColumnList(GtkWidget *blxWindow);
gboolean		  blxWindowGetDisplayRev(GtkWidget *blxWindow);
GtkWidget*		  blxWindowGetBigPicture(GtkWidget *blxWindow);
GtkWidget*		  blxWindowGetDetailView(GtkWidget *blxWindow);
GtkWidget*                blxWindowGetBigPictureCoverageView(GtkWidget *blxWindow);
GtkWidget*                blxWindowGetDetailViewCoverageView(GtkWidget *blxWindow);
GtkWidget*		  blxWindowGetMainMenu(GtkWidget *blxWindow);
GtkWidget*		  blxWindowGetSeqHeaderMenu(GtkWidget *blxWindow);
BlxBlastMode		  blxWindowGetBlastMode(GtkWidget *blxWindow);
IntRange*		  blxWindowGetFullRange(GtkWidget *blxWindow);
IntRange*		  blxWindowGetRefSeqRange(GtkWidget *blxWindow);
const char*		  blxWindowGetRefSeqName(GtkWidget *blxWindow);
BlxSeqType		  blxWindowGetSeqType(GtkWidget *blxWindow);
char**			  blxWindowGetGeneticCode(GtkWidget *blxWindow);
char*			  blxWindowGetRefSeq(GtkWidget *blxWindow);
int			  blxWindowGetNumFrames(GtkWidget *blxWindow);
int			  blxWindowGetDotterStart(GtkWidget *blxWindow);
int			  blxWindowGetDotterEnd(GtkWidget *blxWindow);
int			  blxWindowGetDotterZoom(GtkWidget *blxWindow);
MSP*			  blxWindowGetMspList(GtkWidget *blxWindow);
GList*			  blxWindowGetAllMatchSeqs(GtkWidget *blxWindow);
GList*			  blxWindowGetSequenceGroups(GtkWidget *blxWindow);
SequenceGroup*		  blxWindowGetSequenceGroup(GtkWidget *blxWindow, const BlxSequence *seqToFind);
const char*		  blxWindowGetPaddingSeq(GtkWidget *blxWindow);
int			  blxWindowGetOffset(GtkWidget *blxWindow);
BlxStrand		  blxWindowGetActiveStrand(GtkWidget *blxWindow);
gboolean                  blxWindowGetNegateCoords(GtkWidget *blxWindow);

GList*                    blxWindowGetSelectedSeqs(GtkWidget *blxWindow);
GList*                    blxWindowGetSelectedSeqsByType(GtkWidget *blxWindow, const BlxSequenceType type);
BlxSequence*              blxWindowGetSelectedTranscript(GtkWidget *blxWindow);
void                      blxWindowSelectSeq(GtkWidget *blxWindow, BlxSequence *seq);
void                      blxWindowSetSelectedSeqList(GtkWidget *blxWindow, GList *seqList);
void                      blxWindowDeselectSeq(GtkWidget *blxWindow, BlxSequence *seq);
void                      blxWindowDeselectAllSeqs(GtkWidget *blxWindow);
gboolean                  blxWindowIsSeqSelected(GtkWidget *blxWindow, const BlxSequence *seq);
void                      blxWindowSetSeqSelected(GtkWidget *blxWindow, BlxSequence *seq, const gboolean selected);
void                      blxWindowSelectionChanged(GtkWidget *blxWindow);
BlxSequence*              blxWindowGetLastSelectedSeq(GtkWidget *blxWindow);

int                       sequenceGetGroupOrder(GtkWidget *blxWindow, const BlxSequence *seq);
void                      copySelectionToClipboard(GtkWidget *blxWindow);
void                      findSeqsFromClipboard(GtkClipboard *clipboard, const char *clipboardText, gpointer data);
void                      findAndSelectSeqsFromClipboard(GtkClipboard *clipboard, const char *clipboardText, gpointer data);

void                      refreshDialog(const BlxDialogId dialogId, GtkWidget *blxWindow);
void                      showHelpDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void                      showSettingsDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void                      showSortDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void                      showViewPanesDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void                      showGroupsDialog(GtkWidget *blxWindow, const gboolean editGroups, const gboolean bringToFront);
void                      showFindDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void                      showAboutDialog(GtkWidget *blxWindow);
void                      showInfoDialog(GtkWidget *blxWindow);

void                      blxWindowRedrawAll(GtkWidget *blxWindow);

GtkWidget*                createBlxWindow(CommandLineOptions *options,
                                          const char *paddingSeq,
                                          GArray* featureLists[],
                                          GList *seqList,
                                          GSList *supportedTypes,
                                          const gboolean External,
                                          GSList *styles);


#endif /* _blxwindow_included_ */
